"""Sequence motifs that flag known paralog types in the ITS window.

short58: the short-5.8S paralog of Rhizophagus and [Rhizoglomus]. The short type differs from the main type in one
hairpin near the 3' end of the 5.8S (main ...CACTCTCTGG(T/C)AACCCGGGGAGTATG..., short ...CACTCTCTT(T/C)GGAGAGTATG...).
The motif YTYGGRGARTATG marks all 41 short copies of these genera among the Stefani et al. (2025) culture copies
and none of their 113 other copies; in UNITE (February 2025) it splits Rhizophagus SHs almost perfectly. It is also
the normal 5.8S of Oehlia diaphana and occurs in many unplaced Glomeraceae, so it identifies the paralog only in
reads placed in Rhizophagus. A hit counts when it sits where the hairpin is: 30-45 nt after the fITS7 primer site
(GTGARTCATCGAATCTTTG, up to two mismatches); in the culture copies the offset is 37 nt.
"""
import re

SHORT58 = re.compile('[CT]T[CT]GG[AG]GA[AG]TATG')
FITS7 = 'GTGARTCATCGAATCTTTG'
SHORT58_OFFSET = (30, 45)
SHORT58_GENERA = {'Rhizophagus', 'Rhizoglomus'}


def primer_sites(seq, primer, maxmm=2):
    """0-based starts of primer in seq with at most maxmm mismatches (R matches A or G)."""
    out, n = [], len(primer)
    for i in range(len(seq) - n + 1):
        mm = 0
        for j in range(n):
            p, c = primer[j], seq[i + j]
            if not (c == p or (p == 'R' and c in 'AG') or (p == 'Y' and c in 'CT')):
                mm += 1
                if mm > maxmm:
                    break
        if mm <= maxmm:
            out.append(i)
    return out


def short58(its_seq):
    """True if the ITS window carries the short-5.8S hairpin at its place in the 5.8S."""
    s = its_seq.upper()
    hits = [m.start() for m in SHORT58.finditer(s)]
    if not hits:
        return False
    sites = primer_sites(s, FITS7)
    lo, hi = SHORT58_OFFSET
    return any(lo <= h - i <= hi for h in hits for i in sites)


def paralog_flag(its_seq, genus):
    """'short5.8S' when the motif is present and the read is placed in Rhizophagus or Rhizoglomus; 'short5.8S_motif'
    when the motif is present but the read is placed elsewhere or nowhere; '' otherwise."""
    if not short58(its_seq):
        return ''
    return 'short5.8S' if genus in SHORT58_GENERA else 'short5.8S_motif'
