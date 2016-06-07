from __future__ import print_function, division
import os
import pandas as pd
from pysnptools.snpreader import Bed
from pybedtools import BedTool
from pyutils import memo

class Dataset(object):
    def __init__(self, bfile_chr, assembly='hg19'):
        self.bfile_chr=bfile_chr
        self.assembly=assembly

    @property
    def path(self):
        return os.path.dirname(self.bfile_chr) + '/'

    @memo.memoized
    def data(self, chrnum):
        return Bed(self.bfile(chrnum))
    def stdX(self, chrnum, r):
        return self.stdX_it(chrnum, range(r[0],r[1]))
    def stdX_it(self, chrnum, it):
        genotypes = self.data(chrnum)[:, it].read()
        genotypes.standardize(); genotypes.standardize()
        return genotypes.val

    def bfile(self, chrnum):
        return self.bfile_chr + str(chrnum)
    def bimfile(self, chrnum):
        return self.bfile(chrnum)+'.bim'
    def ucscbedfile(self, chrnum):
        return self.bfile(chrnum)+'.ucscbed'
    def frq_file(self, chrnum):
        return self.bfile(chrnum) + '.frq'
    @memo.memoized
    def frq_df(self, chrnum):
        return pd.read_csv(self.frq_file(chrnum), delim_whitespace=True, header=0)

    @memo.memoized
    def bim_df(self, chrnum):
        return pd.read_csv(self.bfile(chrnum)+'.bim',
                names=['CHR','SNP','CM','BP','A1','A2'],
                sep='\t')

    def M(self, chrnum):
        return len(self.bim_df(chrnum))
    def totalM(self, chroms=None):
        if chroms is None:
            chroms = [c for c in range(1,23) if os.path.exists(self.bimfile(c))]
        return sum(self.M(c) for c in chroms)
    def N(self):
        chrom = min([c for c in range(1,23) if os.path.exists(self.bimfile(c))])
        return self.data(chrom).iid_count

    @memo.memoized
    def ucscbed(self, chrnum):
        if os.path.exists(self.ucscbedfile(chrnum)):
            bt = BedTool(self.ucscbedfile(chrnum))
            if len(bt) == self.M(chrnum):
                return BedTool(self.ucscbedfile(chrnum))
            else:
                print('NOTE: ucscbed file for', self.bfile(chrnum), 'is wrong length.')
        print('NOTE: ucscbed file for',self.bfile(chrnum),'not found/incomplete. creating...')
        ucscbed = []
        for row in self.bim_df(chrnum).iterrows():
            ucscbed.append('chr{}\t{}\t{}'.format(
                row[1]['CHR'],
                row[1]['BP'],
                int(row[1]['BP'])+1))
        print('DONE')
        return BedTool('\n'.join(ucscbed), from_string=True).saveas(
                self.ucscbedfile(chrnum))


if __name__ == '__main__':
    d = Dataset('/groups/price/yakir/data/datasets/1000G3.wim5u/1000G3.wim5u.')
    print(d.path)
    print(d.bfile(22))
    print(d.frq_file(22))
    print(d.frq_df(22).columns)
    print(d.bim_df(22).columns)
    print(d.M(22))
    print(d.totalM())
    print(d.totalM(chroms=[21,22]))
    print(d.N())
