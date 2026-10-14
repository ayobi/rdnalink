import math

from rdnalink import motifs, report, taxonomy, windows
from rdnalink.seqio import rc

# ------------------------------------------------------------------ taxonomy
def test_names_eukaryome_header():
    h = 'EUK1;k__Fungi;p__Glomeromycota;o__Glomerales;f__Glomeraceae;g__Rhizophagus;s__irregularis'
    assert taxonomy.names(h) == ('Rhizophagus_irregularis', 'Rhizophagus', 'Glomerales')


def test_names_placeholders():
    assert taxonomy.names('X;g__Densospora;s__unclassified')[0] is None
    assert taxonomy.names('X;g__Glomeraceae_gen_Incertae_sedis;s__x')[1] is None
    assert taxonomy.names('X;o__Densosporales;g__Densosporales.fam01.gen01;s__y') == (None, None, 'Densosporales')


def test_names_unite_header():
    h = 'MT123|k__Fungi;p__Glomeromycota;g__Rhizophagus;s__Rhizophagus_irregularis|SH1.10FU'
    assert taxonomy.names(h)[:2] == ('Rhizophagus_irregularis', 'Rhizophagus')


# ------------------------------------------------------------------ windows
def test_bounds_itsxrust_both_strands():
    ssu, its, lsu = 'GATTACAGG', 'CATCGATGAAG', 'TTGACCTCAAA'
    read = ssu + its + lsu
    b = windows.bounds_itsxrust({'r': read, 's': rc(read)}, {'r': its, 's': its})
    assert b['r'] == ('ok', 'itsxrust', '+', len(ssu), len(ssu) + len(its))
    assert b['s'][0:3] == ('ok', 'itsxrust', '-')
    assert windows.split_read(rc(read), b['s']) == (ssu, its, lsu)


def test_bounds_barrnap_join_and_multi_operon():
    read = 'A' * 1000
    feats = {'r': [('18S_rRNA', 1, 300, '+', 1e-50), ('18S_rRNA', 320, 400, '+', 1e-5),
                   ('28S_rRNA', 600, 1000, '+', 1e-60)],
             'm': [('18S_rRNA', 1, 300, '+', 1e-50), ('28S_rRNA', 400, 600, '+', 1e-60),
                   ('18S_rRNA', 700, 900, '+', 1e-30)]}
    b = windows.bounds_barrnap({'r': read, 'm': read}, feats, join=True, fix=False)
    assert b['r'][3:] == (400, 599)              # ITS starts after the joined second 18S hit
    assert b['m'][0] == 'multi_operon'


def test_fix_ssu_end_moves_start():
    ssu = 'C' * 200 + windows.SSU_END_MOTIF + 'G' * 13      # motif 27 nt before the true ITS start
    insert = 'A' * 60
    seq = 'C' * 100 + 'T' * 40 + insert + ssu + 'T' * 300    # barrnap stopped before the insert
    i = 200                                                   # barrnap's ITS start, before insert + motif
    j = len(seq) - 100
    n = windows.fix_ssu_end(seq, i, j)
    assert n == seq.index(windows.SSU_END_MOTIF) + windows.MOTIF_TO_ITS


# ------------------------------------------------------------------ motifs
SHORT = ('AAACTTTCAACAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAGTGTGAATTGCAGAACTTTG'
         'TGAATCATCGAATCTTTGAACGCAAATTGCACTCTCTTCGGAGAGTATGCCTGTTTGAGGGTCAGTA')
MAIN = SHORT.replace('CACTCTCTTCGGAGAGTATG', 'CACTCTCTGGTAACCCGGGGAGTATG')


def test_short58():
    assert motifs.short58(SHORT)
    assert not motifs.short58(MAIN)
    assert motifs.paralog_flag(SHORT, 'Rhizophagus') == 'short5.8S'
    assert motifs.paralog_flag(SHORT, 'Oehlia') == 'short5.8S_motif'
    assert motifs.paralog_flag(MAIN, 'Rhizophagus') == ''


def test_short58_needs_its_place():
    far = 'GTGAATCATCGAATCTTTG' + 'A' * 200 + 'CTTCGGAGAGTATG'
    assert not motifs.short58(far)


# ------------------------------------------------------------------ report helpers
def test_unique_best_and_lead():
    sp = {'a': 'X_x', 'b': 'X_x', 'c': 'Y_y'}
    ge = {'a': 'X', 'b': 'X', 'c': 'Y'}
    hl = [(99.9, 'a'), (99.8, 'b'), (99.5, 'c')]
    assert report.unique_best(hl, sp) == ('X_x', 99.9)
    assert report.unique_best([(99.9, 'a'), (99.75, 'c')], sp) == (None, 99.9)
    g, best, lead = report.genus_lead(hl, ge, full=False)
    assert (g, best) == ('X', 99.9) and abs(lead - 0.4) < 1e-9
    assert report.genus_lead([(99.0, 'a'), (98.5, 'b')], ge, full=False)[2] == math.inf
    assert abs(report.genus_lead([(99.0, 'a'), (98.5, 'b')], ge, full=True)[2] - 0.5) < 1e-9
