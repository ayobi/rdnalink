import gzip
import json
import os
import random
import shutil

import pytest

from rdnalink import reference, search
from rdnalink.seqio import write_fasta


def test_distinct_windows_names_by_first_id():
    t = {'b': 'ACGT', 'a': 'ACGT', 'c': 'TTTT'}
    assert search.distinct_windows(t) == {'a': ['a', 'b'], 'c': ['c']}


def test_resolve_folder_and_prefix(tmp_path):
    assert reference.resolve(str(tmp_path / 'win')) == (str(tmp_path / 'win'), None)
    (tmp_path / 'manifest.json').write_text(json.dumps({'mmseqs': {'kmer': 13}}))
    prefix, man = reference.resolve(str(tmp_path))
    assert prefix == str(tmp_path / 'ref') and man['mmseqs']['kmer'] == 13


def _mutate(rng, s, n):
    s = list(s)
    for i in rng.sample(range(len(s)), n):
        s[i] = rng.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)


def _windows(tmp_path):
    """30 labelled reference reads (two pairs with identical windows) and 10 queries, each a mutated copy."""
    rng = random.Random(1)
    lens = {'ssu': 500, 'its': 300, 'lsu': 600}
    base = {m: [''.join(rng.choice('ACGT') for _ in range(n)) for _ in range(6)] for m, n in lens.items()}
    ref, qry = {m: {} for m in lens}, {m: {} for m in lens}
    for i in range(30):
        g = i % 6
        for m in lens:
            ref[m][f'R{i:02d};g__G{g};s__e{i % 3}'] = _mutate(rng, base[m][g], 1 + i % 7)
    for m in lens:                                      # identical windows, to be merged
        ref[m]['R30;g__G0;s__e0'] = ref[m]['R00;g__G0;s__e0']
        ref[m]['R31;g__G1;s__e1'] = ref[m]['R01;g__G1;s__e1']
    for i in range(10):
        for m in lens:
            qry[m][f'Q{i:02d}'] = _mutate(rng, base[m][i % 6], 10)
    for name, d in (('ref', ref), ('qry', qry)):
        for m, seqs in d.items():
            with open(tmp_path / f'{name}.{m}.fasta', 'w') as fh:
                for h, s in seqs.items():
                    write_fasta(fh, h, s)
    return str(tmp_path / 'ref'), str(tmp_path / 'qry')


def _ranked(path):
    """Hit rows without bit score and E-value: MMseqs2's bit scores can differ by a few bits between runs of the
    same search, and ranking never uses them."""
    with gzip.open(path, 'rt') as fh:
        return [line.split('\t')[:10] + line.split('\t')[12:] for line in fh]


@pytest.mark.skipif(shutil.which('mmseqs') is None, reason='needs mmseqs on PATH')
def test_built_reference_gives_the_same_hits(tmp_path):
    ref_prefix, qry_prefix = _windows(tmp_path)
    man = reference.build(str(tmp_path / 'REF'), windows_prefix=ref_prefix, threads=1, kmer=13)
    assert man['markers']['ssu'] == {'windows': 32, 'distinct': 30,
                                     'windows_sha256': man['markers']['ssu']['windows_sha256']}
    search.search(qry_prefix, str(tmp_path / 'REF'), str(tmp_path / 'hits_ref'), threads=1)
    search.search(qry_prefix, ref_prefix, str(tmp_path / 'hits_prefix'), threads=1, kmer=13)
    for m in ('ssu', 'its', 'lsu'):
        a, b = (_ranked(tmp_path / d / f'{m}.hits.tsv.gz') for d in ('hits_ref', 'hits_prefix'))
        assert a == b and len(a) > 10
    assert not os.path.exists(tmp_path / 'hits_ref' / 'work' / 'ssu' / 'tdb')      # the built index was used
    with pytest.raises(SystemExit):
        search.search(qry_prefix, str(tmp_path / 'REF'), str(tmp_path / 'hits_k'), threads=1, kmer=15)
