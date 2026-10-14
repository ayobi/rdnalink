"""rdnalink: species calls and checks from the SSU, ITS and LSU of one long rDNA read (prototype).

  rdnalink windows READS.fasta -o PREFIX --hmm F.hmm      cut SSU, ITS and LSU windows
  rdnalink search -q PREFIX -r REF_PREFIX -o HITS_DIR     one MMseqs2 search per window
  rdnalink report -q PREFIX -r REF_PREFIX --hits HITS_DIR -o report.tsv
                                                          per-read calls, agreement and flags

The reference is a set of labelled long reads cut with the same `windows` step (headers with k__ ... g__ s__ ranks,
as in EUKARYOME or UNITE). Run `rdnalink <command> -h` for options.
"""
import argparse
import sys

from . import __version__


def _min_identity(text):
    out = {}
    for part in (text or '').split(','):
        if part.strip():
            m, v = part.split('=')
            out[m.strip()] = float(v)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog='rdnalink', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=f'rdnalink {__version__}')
    sub = ap.add_subparsers(dest='cmd', required=True)

    w = sub.add_parser('windows', help='cut SSU, ITS and LSU windows (ITSxRust, then barrnap)')
    w.add_argument('reads')
    w.add_argument('-o', '--out', required=True, help='output prefix')
    w.add_argument('--hmm', required=True, help="ITSxRust's F.hmm")
    w.add_argument('-t', '--threads', type=int, default=4)
    w.add_argument('--chunk', type=int, default=4000, help='reads per tool run (E-values scale with it)')
    w.add_argument('--min-length', type=int, default=300, help='shortest window kept (nt)')
    w.add_argument('--itsxrust', default='itsxrust')
    w.add_argument('--barrnap', default='barrnap')

    s = sub.add_parser('search', help='one MMseqs2 search per window')
    s.add_argument('-q', '--query', required=True, help='query window prefix')
    s.add_argument('-r', '--ref', required=True, help='reference window prefix')
    s.add_argument('-o', '--out', required=True, help='output folder')
    s.add_argument('-t', '--threads', type=int, default=4)
    s.add_argument('--keep', type=int, default=300, help='hits stored per query and marker')
    s.add_argument('--min-cov', type=float, default=0.8, help='coverage of the shorter window')
    s.add_argument('--mmseqs', default='mmseqs')
    s.add_argument('--kmer', type=int, default=0, help="k-mer size (default: MMseqs2's, 15 for nucleotides, which "
                   'needs about 10 GB of RAM for a large reference; 13 fits a laptop)')

    r = sub.add_parser('report', help='per-read calls, agreement and flags')
    r.add_argument('-q', '--query', required=True)
    r.add_argument('-r', '--ref', required=True)
    r.add_argument('--hits', required=True)
    r.add_argument('-o', '--out', required=True)
    r.add_argument('--top', type=int, default=300, help='ranks read per marker')
    r.add_argument('--lead', type=float, default=1.0, help='points a confident genus must lead by')
    r.add_argument('--min-identity', default='', help='e.g. ssu=99,its=97,lsu=98 (default: none)')
    r.add_argument('--keep-self', action='store_true', help='let a read be its own hit')
    a = ap.parse_args(argv)

    if a.cmd == 'windows':
        from .windows import cut, tool_versions
        st = cut(a.reads, a.out, a.hmm, a.threads, a.chunk, a.min_length, a.itsxrust, a.barrnap)
        print('; '.join(f'{k}: {v}' for k, v in tool_versions(a.itsxrust, a.barrnap).items()), file=sys.stderr)
        for (status, src), n in sorted(st.items(), key=lambda x: -x[1]):
            print(f'{status}\t{src}\t{n}')
    elif a.cmd == 'search':
        from .search import search
        for m, (nq, uq, nt, ut) in search(a.query, a.ref, a.out, a.mmseqs, a.threads, a.keep, a.min_cov,
                                          kmer=a.kmer).items():
            print(f'{m}: {nq} queries ({uq} distinct) against {nt} windows ({ut} distinct)')
    elif a.cmd == 'report':
        from .report import report
        c = report(a.hits, a.query, a.ref, a.out, a.top, a.lead, _min_identity(a.min_identity), not a.keep_self)
        for k, v in sorted(c.items()):
            print(f'{k}\t{v}')


if __name__ == '__main__':
    main()
