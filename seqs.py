import sys

class Fasta:
    def __init__(self, name, seq=None):
        self.name = name
        self.seq = seq
    def __repr__(self):
        return ">%s\n%s".format(self.name, self.seq)

def readFastaByChrom(filename):
    """Read a file in fasta format. Return an iterator of Fasta objects where
    each object is a fasta sequence in the file. Uses a generator for
    memory-efficient reading in of the sequences"""
    f = open(filename)
    file_pos = f.tell()
    line = f.readline()
    while line.startswith('#'):
        file_pos = f.tell()
        line = f.readline()
    f.seek(file_pos)
    while True:
        yield parseFastaMatch(f)
    f.close()

def parseFastaMatch(filehandle):
    "Read in lines representing a single fasta sequence"
    readName = filehandle.readline().rstrip('\n')[1:]
    if readName == '':
        raise StopIteration
    fastaObj = Fasta(name=readName)
    filePos = filehandle.tell()
    line = filehandle.readline().rstrip('\n')
    seq = [line]
    while True:
        filePos = filehandle.tell()
        line = filehandle.readline().rstrip('\n')
        if line == '': break
        elif line.startswith('>'):
            filehandle.seek(filePos)
            break
        else:
            seq.append(line)
    fastaObj.seq = ''.join(seq)
    return fastaObj

def pretty_print(seq, colwidth=60, filehandle=sys.stdout):
    """Pretty prints a sequence -- i.e. prints the sequence with a newline
    every colwidth characters. The sequence is printed to the filehandle
    specified in the function call (stdout by default)."""
    l = len(seq)
    for i in range(0, l, colwidth):
        filehandle.write("%s\n" % seq[i:i+colwidth])
