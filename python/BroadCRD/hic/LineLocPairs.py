class LineLocPairs:
    import struct
    def __init__(self, filename):
        f=open(filename,'r')
        if f.read(len('BINWRITE')) != 'BINWRITE':
            raise Exception('not a CompRD binary file')
        length=self.struct.unpack('Q',f.read(self.struct.calcsize('Q')))
        self.data=list()
        for i in range(length[0]):
            self.data.append( self.struct.unpack('iiiiii', \
                    f.read(self.struct.calcsize('iiiiii') ) ) )

    def frame(self):
        import pandas as pd
        f=dict()
        columns=['line1','left1','right1','line2','left2','right2']
        for name in columns:
            f[name]=list()
        for rec in self.data:
            for i in range(len(columns)):
                f[columns[i]].append(rec[i])
        return pd.DataFrame(f,columns=columns)
