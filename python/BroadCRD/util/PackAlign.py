
class Align:
    """Align - implement a class to read and represent an "align" struct from PackAlign.h

    a=Align("filename.align")

    yields a.pos1, a.pos2, a.nblocks, a.gaps, a.lens

    a.nblocks = len(a.gaps) = len(a.lens)
    """

    def __init__(self, filename=""):
        if filename != "":
            self.read_from(filename)
        else:
            pos1=0
            pos2=0
            nblocks=0
            gaps=list()
            lens=list()

    def read_from(self,filename):
        import struct
        align=open(filename,'r').read()
        offset=len("BINWRITE")
        if align[:offset] != "BINWRITE":
            raise Exception("file " + filename + " seems not from the BinaryWriter")
        align=align[offset:]
        fmt="iii"
        fsize=struct.calcsize(fmt)
        assert fsize <= len(align)
        self.pos1,self.pos2,self.nblocks=struct.unpack_from( fmt, align, 0 )

        afmt="i"*self.nblocks
        afsize=struct.calcsize(afmt)
        assert ( fsize + afsize ) <= len(align)
        self.gaps=list(struct.unpack_from( afmt, align, fsize ))

        assert ( fsize + afsize + afsize ) <= len(align)
        self.lens=list(struct.unpack_from( afmt, align, fsize+afsize ))

