"""FASTA input and output, read IDs and checksums."""
import gzip
import hashlib
import re

_COMP = str.maketrans('ACGTUNacgtunRYKMSWBDHVrykmswbdhv', 'TGCAANtgcaanYRMKSWVHDByrmkswvhdb')


def opener(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def read_fasta(path):
    """Yield (header, sequence) pairs; the header without '>' and the sequence on one line."""
    head, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(seq)
                head, seq = line[1:].rstrip('\n'), []
            else:
                seq.append(line.strip())
    if head is not None:
        yield head, ''.join(seq)


def write_fasta(fh, header, seq, width=80):
    fh.write(f'>{header}\n')
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + '\n')


def read_id(header):
    """The read ID: the header up to the first space, ';' or '|' (EUKARYOME and UNITE style headers)."""
    return re.split(r'[\s;|]', header, maxsplit=1)[0]


def rc(seq):
    return seq.translate(_COMP)[::-1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()
