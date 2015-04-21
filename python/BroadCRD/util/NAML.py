import re

class NAML(dict):
    """NAML -- Not Another Markup Language
    somewhat of a misnomer, but a response to YAML.  The point here was
    a dead-simple, self-contained class for reading key-value pairs from
    a file.  There are no types -- everything read is a string.  A few
    caveats:

    1. An input file comprises lines of the form:

    keyname: value

    2. A line starting with a # symbol in the first column is a comment.
    Comments read in from a file will be written back out, all bunched
    up at the beginning.
    3. If multiple lines specify the same key, the value is a list of
    the values aggregated from all such lines.
    4. Whitespace between the colon and the first non-whitespace
    character in the "value" is ignored.  Whitespace later in the line
    is kept intact.

    neilw - july 2014
    """
    def __init__(self, filename=None):
        super(NAML, self).__init__()
        self.comments=list()
        if filename: self.read_from_file(filename)

    class ParseFailure(Exception):
        def __init__(self, filename, line):
            msg="failed to parse line\n {}\n in file {}".format(line.rstrip(), filename)
            super(NAML.ParseFailure,self).__init__(msg)

    def read_from_file(self, filename):
        with open(filename, 'r') as input:
            for line in input:
                line=line.rstrip('\r\n')
                if line[0] == '#':
                    self.comments.append(line)
                    continue
                m=re.match('(?P<key>[^:]+):\s+(?P<val>.*)', line)
                if not m:
                    raise NAML.ParseFailure(filename, line)
                else:
                    k=m.group('key')
                    v=m.group('val')
                    if k in self:
                        try:
                            self[k].append(v)
                        except AttributeError:
                            self[k]=[self[k]]
                            self[k].append(v)
                    else:
                        self[k]=v

    def write_to_file(self, filename):
        with open( filename, 'w' ) as output:
            for c in self.comments: output.write(c+'\n')
            for k,v in self.iteritems():
                if type(v) != list: v = [v]
                for val in v: output.write( "{}: {}\n".format(k,val) )
