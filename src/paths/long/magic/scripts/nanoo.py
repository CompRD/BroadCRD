#!/usr/bin/env python
# nanoo
#
# A noo version of nanoc.  Really a complementary program for bulk converting a lot of reads, hopefully a
# bit faster.
#
# note: my hopes of faster-ness have been somewhat dashed by the fact that the to_csv function
# is Pandas is slow as molasses. But we're still better off doing multiple files at once and
# we have new file formats
#
# neilw - Nov 14, 2014

import argparse
import sys
import h5py as h
import os
import BroadCRD.util.Nanopore as Nanopore

def write_fastq(filename, fq_lines):
    open(filename, 'w').write('\n'.join(fq_lines) )
    return filename

def write_data(filename, data, attrs=None):
    with open(filename, 'w') as fd:
        try:
            if attrs is not None:
                fd.write("# ATTRS\n")
                for key, val in attrs.iteritems():
                    fd.write("{}: {}\n".format(key,val))
            fd.write("# HEAD\n")
            fd.write(",".join(data.columns.values)+"\n")
            fd.write("# DATA\n")
            data.to_csv(path_or_buf=fd, header=False, index=False)
        except:
            os.unlink(filename)
            raise
    return filename


class Fast5Files(list):
    def __init__(self,directory):
        import os
        super(Fast5Files, self).__init__()
        self.extend([os.path.join(directory,filename) \
                     for filename in os.listdir(directory) \
                     if filename[-6:] == '.fast5'])

if __name__ == '__main__':
    PREF_COMPLEMENT="comp_"
    PREF_TEMPLATE="temp_"
    PREF_CONSENSUS="cons_"


    parser = argparse.ArgumentParser(description='bulk conversion of Nanoporetech reads')
    parser.add_argument('dir', metavar='fast5-dir', help='directory of .fast5 files' )
    parser.add_argument('odir', metavar='output-dir', help='directory for output files')
    args = parser.parse_args(sys.argv[1:])


    count = 0
    for fast5 in Fast5Files(args.dir):
#        input=h.File(fast5,'r')

        (dirname,filename)=os.path.split(fast5)
        (basename,_)=os.path.splitext(filename)

        try:
            read=Nanopore.Read(fast5)
            if read.cfastq == None or read.tfastq == None or read.fastq == None:
                print "skipping {}, no full 2D data".format(fast5)
                continue
        except:
            print >>sys.stderr, "error reading fast5 file: "+fast5
            raise

        # write the consensus, complement, and template
        try:
            names=[]
            for fq_lines, prefix in [ \
                                ( read.fastq, PREF_CONSENSUS ), \
                                ( read.cfastq, PREF_COMPLEMENT ), \
                                ( read.tfastq, PREF_TEMPLATE ) \
                                ]:
                names.append( write_fastq( os.path.join(args.odir,prefix + basename+".fastq")
                                , fq_lines))

            names.append( write_data( os.path.join(args.odir,PREF_COMPLEMENT+basename+".events" ),
                                read.cevents))
            names.append( write_data( os.path.join(args.odir, PREF_COMPLEMENT+basename+".model" ),
                                read.cmodel, read.cmodel_attrs))
            names.append( write_data( os.path.join(args.odir,PREF_TEMPLATE+basename+".events" ),
                                read.tevents))
            names.append( write_data( os.path.join(args.odir, PREF_TEMPLATE+basename+".model" ),
                                read.tmodel, read.tmodel_attrs))
            names.append( write_data( os.path.join(args.odir, PREF_CONSENSUS+basename+".align"),
                                read.align))
            names.append( write_data( os.path.join(args.odir, basename+".hairpin"), read.hairpin))
            print "{}: completed FASTQ conversion of {}".format(count,fast5)
        except:
            print >>sys.stderr, "failed writing one of the files for "+fast5
            print >>sys.stderr, "removing any intermediate results"
            for name in names: os.unlink(name)
            raise

        count+=1
