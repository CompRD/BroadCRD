# class VCF
#
# okay, so create a VCF object like this:
#
# v = VCF('/path/to/file.vcf')
#
# v gets the following:
#
# metadata - list of initial metadata lines in the VCF prior to the #CHROM line
# sample_names - in-order list of sample names from the #CHROM line
#
#
class VCFException(Exception):
    def __init__(self, msg):
        Exception.__init__(self,msg)

class VCFLine:
    """Takes a data line from a VCF, parses it, and encapsulates it"""
    def __init__(self, line):
        self.line=line.strip()
        fields=line.split()
        (chr,pos,id,ref,alt,qual,filter,info,format)=fields[:9]
        self.chr = chr
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.samples = fields[9:]

    @property
    def alt_list(self):
        """Returns the ALT bases as a list"""
        return self.alt.split(',')

    def get_sample_dict(self, sample):
        """Given a zero-based sample index, returns a python dictionary with the sample
        information.  For instance:
            v=VCF("vcf_file.vcf")
            for line in v.lines():
                print "Genotype={}".format( line.get_sample_dict(0)["GT"] )

        Obviously, it's up to you to ensure that the dictionary contains
        what you want.
        """

        if sample >= len(self.samples):
            raise VCFException(\
            "requested sampled index (zero-based) {}, but we only have {} samples".format( sample, len(self.samples) ) )

        names = self.format.split(":")
        vals = self.samples[sample].split(":")

        return dict(zip(names, vals))

    def __repr__(self):
       temp=[ self.chr, str(self.pos), self.id, self.ref, \
            self.alt, self.qual, self.filter, self.info, self.format \
            ]
       temp.extend( self.samples )
       return "\t".join(temp)



class VCF:
    """Simple VCF class which parses the header, saves the sample names
    as a list, saves the metadata lines, and parses the file
    line-by-line through a generator returned by the lines() method.
    VCF lines are returned as VCFLine objects"""

    def __init__(self, filename, make_index = False ):
        self.filename = filename
        self._fd = open( filename, 'r' )
        self.metadata = []
        self._index = None

        line = self._fd.readline().strip()
        while line[:len('#CHROM')] != '#CHROM':
            if line[0] != '#':
                raise VCFException("haven't seen the #CHROM line, so I expected a line starting with #, but I saw ".format(line))
            self.metadata.append( line )
            line = self._fd.readline().strip()

        self._tell = self._fd.tell()      # record the position of the first data line

        self.metadata.append( line )      # keep the #CHROM line here, too

        fields=line.split()
        self.sample_names = fields[9:]

        if make_index:
            last_chr = None
            last_tell = self._fd.tell()
            for line in self._fd:
                chr=line.split()[0]
                if chr != last_chr:
                    last_chr = chr
                    if not self._index: self._index = dict()
                    self._index[chr] = last_tell
                    last_tell = self._fd.tell()
            self._fd.seek( self._tell )

    def seek(self, contig):
        """returns found seek offset or -1 if not found
        seeks to the beginning if there is no index.  Returns the
        position to which we seek"""
        if self._index:
            if self._index.has_key( contig ):
                self._fd.seek( self._index[contig] )
                return self._index[contig]
            else:
                return -1
        else:
            self._fd.seek(self._tell)
            return 0


    def lines(self, reset = False, filter = None):
        """returns a generator that provides data lines from the VCF
        file, e.g.:

        for line in VCF( filename ).lines():
            process(line)

        if reset is True, then the file scanning starts over from beginning
        if filter is provided, then each (raw text) line is passed to
        filter first and only lines returning True from filter are
        returned as VCFLines
        """
        if reset: self._fd.seek( self._tell )

        if filter != None:
            return ( VCFLine( line ) for line in self._fd if filter(line) )
        else:
            return ( VCFLine( line ) for line in self._fd )


