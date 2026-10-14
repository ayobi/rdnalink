"""Read labels from FASTA headers, and which of them count as names.

Headers carry ranks as k__ p__ c__ o__ f__ g__ s__ (EUKARYOME, UNITE), separated by ';', '|' or ','. A named
genus has no placeholder marks (no 'unclassified', 'sp', 'incertae', digits or dots); a named species is a named
genus plus a lower-case epithet, written Genus_epithet. The rules are those of the rdnalink test protocol
(section 2), so names mean the same thing in the tool and in the test.
"""
import re

TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')
EPITHET = re.compile(r'^[a-z][a-z-]+$')
BINOMIAL = re.compile(r'^[A-Z][a-z]+[_ ][a-z][a-z-]+$')
PLACEHOLDER = re.compile(r'(?i)(^|[_ .])(sp|spp|cf|aff|nov|unidentified|uncultured|unclassified|incertae|'
                         r'sedis|environmental|metagenome)([_ .]|$)|\d|\.')


def ranks(header):
    """{rank letter: value} for the k..s fields of a header (first value of each rank)."""
    out = {}
    for r, v in TAX.findall(header):
        out.setdefault(r, v.strip())
    return out


def names(header):
    """(named species or None, named genus or None, order or None) for one header."""
    t = ranks(header)
    g, s, o = t.get('g', ''), t.get('s', ''), t.get('o', '')
    named_g = bool(g) and g[0].isupper() and not PLACEHOLDER.search(g)
    sp = None
    if s and named_g:
        if EPITHET.match(s):
            s = f'{g}_{s}'
        if BINOMIAL.match(s) and not PLACEHOLDER.search(s):
            sp = s.replace(' ', '_')
    named_o = bool(o) and o[0].isupper() and not PLACEHOLDER.search(o)
    return sp, (g if named_g else None), (o if named_o else None)
