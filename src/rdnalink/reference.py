"""A reusable reference: windows cut once, identical windows merged, one MMseqs2 index per marker.

`rdnalink ref` builds a folder that `search` and `report` take in place of a window prefix:

  REF/ref.{ssu,its,lsu}.fasta   the reference windows (headers keep the labels)
  REF/ref.log.tsv               how each read was cut (when `ref` cut them)
  REF/<marker>/target.fasta     one sequence per distinct window, named by the first read ID (in sort order) that
                                carries it
  REF/<marker>/members.tsv      every read ID that carries each distinct window
  REF/<marker>/tdb*             MMseqs2 database and index of target.fasta
  REF/manifest.json             inputs and their SHA-256, tool versions, parameters and counts; written last, so a
                                folder without it is an unfinished build

A search against a built reference gives the same hits as a search against its windows: the distinct windows,
their names and the index are the ones `search` would build itself. It only skips building them for every search.
"""
import datetime
import json
import os
import shutil
import subprocess
import sys

from . import __version__
from .seqio import read_fasta, read_id, sha256, write_fasta

MANIFEST = 'manifest.json'
PREFIX = 'ref'
MARKERS = ('ssu', 'its', 'lsu')


def resolve(ref):
    """A folder built by `rdnalink ref` -> (its window prefix, its manifest); anything else is a window prefix."""
    path = os.path.join(ref, MANIFEST)
    if os.path.isdir(ref) and os.path.exists(path):
        with open(path) as fh:
            return os.path.join(ref, PREFIX), json.load(fh)
    return ref, None


def target_db(ref, marker):
    return os.path.join(ref, marker, 'tdb')


def mmseqs_version(mmseqs='mmseqs'):
    try:
        p = subprocess.run([mmseqs, 'version'], capture_output=True, text=True)
        return (p.stdout or p.stderr).strip().splitlines()[0]
    except (OSError, IndexError):
        return 'not found'


def _window_file(prefix, marker):
    for path in (f'{prefix}.{marker}.fasta', f'{prefix}.{marker}.fasta.gz'):
        if os.path.exists(path):
            return path
    sys.exit(f'no {prefix}.{marker}.fasta(.gz)')


def build(out, reads=None, windows_prefix=None, hmm=None, threads=4, chunk=4000, min_length=300, kmer=0,
          itsxrust='itsxrust', barrnap='barrnap', mmseqs='mmseqs'):
    """Build a reference folder from labelled long reads (cut here) or from windows already cut. Returns the
    manifest."""
    from .search import build_target, distinct_windows
    if (reads is None) == (windows_prefix is None):
        sys.exit('give either the reads or --windows, not both')
    if os.path.exists(os.path.join(out, MANIFEST)):
        sys.exit(f'{out} already holds a built reference; remove it or choose another folder')
    if shutil.which(mmseqs) is None and not os.path.exists(mmseqs):
        sys.exit(f'{mmseqs} not found on PATH')
    os.makedirs(out, exist_ok=True)
    prefix = os.path.join(out, PREFIX)
    man = {'rdnalink': __version__}

    if reads is not None:
        from .windows import BARRNAP_ARGS, ITSXRUST_ARGS, cut, tool_versions
        if not hmm:
            sys.exit('--hmm is needed to cut windows')
        status = cut(reads, prefix, hmm, threads, chunk, min_length, itsxrust, barrnap)
        work = prefix + '.work'
        if os.path.exists(os.path.join(work, 'tools.log')):
            shutil.copy(os.path.join(work, 'tools.log'), os.path.join(out, 'windows_tools.log'))
        shutil.rmtree(work, ignore_errors=True)
        man['source'] = {'reads': os.path.abspath(reads), 'reads_sha256': sha256(reads),
                         'reads_n': sum(status.values()), 'hmm': os.path.abspath(hmm), 'hmm_sha256': sha256(hmm)}
        man['windows'] = {'tools': tool_versions(itsxrust, barrnap), 'chunk': chunk, 'min_length': min_length,
                          'itsxrust_args': ITSXRUST_ARGS, 'barrnap_args': BARRNAP_ARGS,
                          'status': {f'{st}:{src}': n for (st, src), n in sorted(status.items())}}
    else:
        src = {m: _window_file(windows_prefix, m) for m in MARKERS}
        for m, path in src.items():
            with open(f'{prefix}.{m}.fasta', 'w') as fh:
                for h, s in read_fasta(path):
                    write_fasta(fh, h, s)
        man['source'] = {'windows': {m: os.path.abspath(p) for m, p in src.items()},
                         'windows_sha256': {m: sha256(p) for m, p in src.items()}}

    log = open(os.path.join(out, 'mmseqs.log'), 'a')
    man['markers'] = {}
    for m in MARKERS:
        t = {read_id(h): s.upper() for h, s in read_fasta(f'{prefix}.{m}.fasta')}
        members = distinct_windows(t)
        build_target(os.path.join(out, m), t, members, mmseqs, threads, kmer, log)
        with open(os.path.join(out, m, 'members.tsv'), 'w') as fh:
            fh.write('representative\tn\treads\n')
            for rep in sorted(members):
                fh.write(f'{rep}\t{len(members[rep])}\t{",".join(members[rep])}\n')
        man['markers'][m] = {'windows': len(t), 'distinct': len(members),
                             'windows_sha256': sha256(f'{prefix}.{m}.fasta')}
    log.close()
    man['mmseqs'] = {'version': mmseqs_version(mmseqs), 'kmer': kmer,
                     'index': 'createindex --search-type 3 --strand 1' + (f' -k {kmer}' if kmer else '')}
    man['created'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    with open(os.path.join(out, MANIFEST + '.tmp'), 'w') as fh:
        json.dump(man, fh, indent=2)
        fh.write('\n')
    os.replace(os.path.join(out, MANIFEST + '.tmp'), os.path.join(out, MANIFEST))
    return man
