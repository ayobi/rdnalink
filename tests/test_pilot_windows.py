"""Integration test: `rdnalink windows` reproduces Amendment 1's windows on the pilot.

Runs the real ITSxRust and barrnap on the 3,994 pilot reads (`pilot/pilot.fasta`) and compares every read's
status, source, strand and ITS bounds with `pilot/cut_hybridfastjm.log.tsv`, the pilot log of the rule Amendment 1
fixed (ITSxRust first, then barrnap --fast with the join and 18S-end rules). Takes about 10 minutes on 8 threads.
Windows only: no search, no report, nothing scored. Opt in with

    RDNALINK_HMM=/path/to/F.hmm pytest tests/test_pilot_windows.py

(RDNALINK_THREADS sets threads, default 8). On a mismatch the reads that differ are written to
pilot_windows_diff.tsv in the current folder.
"""
import csv
import os

import pytest

from rdnalink import windows

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PILOT = os.path.join(ROOT, 'pilot', 'pilot.fasta')
LOG = os.path.join(ROOT, 'pilot', 'cut_hybridfastjm.log.tsv')
HMM = os.environ.get('RDNALINK_HMM')
KEYS = ('status', 'source', 'strand', 'its_start', 'its_end')


def _rows(path):
    with open(path) as fh:
        return {r['id']: r for r in csv.DictReader(fh, delimiter='\t')}


@pytest.mark.skipif(not HMM, reason='set RDNALINK_HMM to F.hmm to run the pilot comparison')
def test_pilot_windows_match_amendment1(tmp_path):
    windows.cut(PILOT, str(tmp_path / 'pilot'), HMM, threads=int(os.environ.get('RDNALINK_THREADS', '8')))
    got, want = _rows(str(tmp_path / 'pilot.log.tsv')), _rows(LOG)
    assert set(got) == set(want), 'read IDs differ from the pilot log'
    diff = []
    for k, w in want.items():
        g = dict(got[k])
        if g['status'] == 'short_window':          # the pilot log kept short windows as ok
            g['status'] = 'ok'
        if any(g[c] != w[c] for c in KEYS):
            diff.append([k] + [w[c] for c in KEYS] + [g[c] for c in KEYS])
    if diff:
        with open('pilot_windows_diff.tsv', 'w') as fh:
            fh.write('id\t' + '\t'.join(f'pilot_{c}' for c in KEYS) + '\t' + '\t'.join(f'now_{c}' for c in KEYS)
                     + '\n')
            for d in diff:
                fh.write('\t'.join(d) + '\n')
    assert not diff, f'{len(diff)} of {len(want)} reads differ; see pilot_windows_diff.tsv'
