from __future__ import print_function, division
import os
import numpy as np
import pandas as pd
from pyutils import memo

class Annotation(object):
    def __init__(self, stem_chr, signed=True):
        self.stem_chr = stem_chr
        self.signed = signed

    def filestem(self, chrnum=''):
        return '{}{}'.format(self.stem_chr, chrnum)
    def annot_filename(self, chrnum):
        return self.filestem(chrnum) + '.annot.gz'
    def sannot_filename(self, chrnum):
        return self.filestem(chrnum) + '.sannot.gz'
    def sqnorm_filename(self, chrnum):
        return self.filestem(chrnum) + '.sqnorm'
    def size_filename(self, chrnum):
        return self.filestem(chrnum) + '.M'
    def ldscores_filename(self, chrnum):
        return self.filestem(chrnum) + '.l2.ldscore.gz'
    def conv_filename(self, chrnum, full=False):
        if not full:
            return self.filestem(chrnum) + '.conv.gz'
        else:
            return self.filestem(chrnum) + '.fullconv.gz'
    @classmethod
    def isfullconv(cls, filename):
        return filename.endswith('.fullconv.gz')

    @memo.memoized
    def annot_df(self, chrnum):
        return pd.read_csv(self.annot_filename(chrnum),
                compression='gzip', header=0, sep='\t')
    def sannot_df(self, chrnum):
        return pd.read_csv(self.sannot_filename(chrnum),
                compression='gzip', header=0, sep='\t')
    @memo.memoized
    def sqnorms(self, chrnum):
        return pd.read_csv(self.sqnorm_filename(chrnum), names=self.names(chrnum), sep='\t')
    @memo.memoized
    def sizes(self, chrnum):
        return pd.read_csv(self.size_filename(chrnum), names=self.names(chrnum), sep='\t')

    @memo.memoized
    def names(self, chrnum):
        if os.path.exists(self.annot_filename(chrnum)):
            return self.annot_df(chrnum).columns.values[4:]
        else: # assuming sannot exists
            return self.sannot_df(chrnum).columns.values[6:]
    @classmethod
    def names_observed(cls, names):
        return [n + '.O' for n in names]
    @classmethod
    def names_conv(cls, names, observed=True):
        O = ('.O' if observed else '')
        return [n + O + '.conv' for n in names]

    @memo.memoized
    def num_snps(self, chrnum):
        if os.path.exists(self.annot_filename(chrnum)):
            return len(self.annot_df(chrnum))
        else: # assuming sannot exists
            return len(self.sannot_df(chrnum))

    def total_sqnorms(self, chromosomes):
        return sum([self.sqnorms(c) for c in chromosomes])
    def total_sizes(self, chromosomes):
        return sum([self.sizes(c) for c in chromosomes])


if __name__ == '__main__':
    a = Annotation('/groups/price/yakir/data/annot/1000G3.wim5u/mock/')
    print(a.names(22))
    print(a.sannot_df(22))
    print(a.sqnorms(22))
    print(a.sizes(22))
    print(a.total_sqnorms([1,22]))
    print(a.total_sizes([1,22]))
