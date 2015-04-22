
class Read:
    import h5py as h5py
    import pandas as pd

    def __init__(self, filename):
        self._file = self.h5py.File(filename, 'r')
        self._filename = filename

        # it's super-ugly that we need to know below whether a member variable is a DataFrame
        # (and is .empty if not present) or if it's another datatype that might be None.
        #
        # This is inelegant, but servicable
        #

        # do we have a template
        self.tevents = self.df_if_exists('/Analyses/Basecall_2D_000/BaseCalled_template/Events')
        if not self.tevents.empty:
            self.tfastq = self._file['/Analyses/Basecall_2D_000/BaseCalled_template/Fastq'].value.split('\n')
            self.tmodel = self.pd.DataFrame(self._file['/Analyses/Basecall_2D_000/BaseCalled_template/Model'].value)
            self.tmodel_attrs = self._file['/Analyses/Basecall_2D_000/BaseCalled_template/Model'].attrs
        else:
            self.tfastq = self.tmodel_attrs = None
            self.tmodel = self.pd.DataFrame()

        # do we have a complement
        self.cevents = self.df_if_exists('/Analyses/Basecall_2D_000/BaseCalled_complement/Events')
        if not self.cevents.empty:
            self.cfastq = self._file['/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq'].value.split('\n')
            self.cmodel = self.pd.DataFrame(self._file['/Analyses/Basecall_2D_000/BaseCalled_complement/Model'].value)
            self.cmodel_attrs = self._file['/Analyses/Basecall_2D_000/BaseCalled_complement/Model'].attrs
        else:
            self.cfastq = self.cmodel_attrs = None
            self.cmodel = self.pd.DataFrame()

        self.align = self.df_if_exists('/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment')
        self.hairpin = self.df_if_exists('/Analyses/Basecall_2D_000/HairpinAlign/Alignment')
        self.fastq = self.value_if_exists('/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq')
        if self.fastq: self.fastq=self.fastq.split('\n')

        read0name = self._file['/Analyses/EventDetection_000/Reads'].values()[0].name
        self.events = self.pd.DataFrame(self._file[read0name + '/Events'].value)

    def df_if_exists(self, key):
        try: return self.pd.DataFrame( self._file[key].value )
        except KeyError: return self.pd.DataFrame()     # check with .empty

    def value_if_exists(self, key):
        try: return self._file[key].value
        except KeyError: return None

    def attr_if_exists(self, key):
        try: return self._file[key].attrs
        except KeyError: return None

    def __repr__(self):
        return "Nanopore read from file: " + self._filename

    def classify(self, value, model):
        import math
        k = self.pd.DataFrame(model.kmer)
        k['p'] = 0.
        for i in range(len(k)):
            nval = (value - model.level_mean[i]) / model.level_stdv[i]
            k[i:i + 1].p = 1. / math.sqrt(2.*math.pi) * model.level_stdv * math.exp(-nval * nval / 2.)
        return k


    def tclassify(self, value): return self.classify(value, self.tmodel)
    def cclassify(self, value): return self.classify(value, self.tmodel)
