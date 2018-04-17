from __future__ import print_function, division
import os
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from pybedtools import BedTool
from ypy import memo
import ypy.iter as pyit
import time


# Wrapper class for plink-formatted genotype data sets
class Dataset(object):
    def __init__(self, bfile_chr, assembly='hg19'):
        self.bfile_chr=bfile_chr
        self.assembly=assembly

    @property
    def path(self):
        return os.path.dirname(self.bfile_chr) + '/'

    #TODO: investigate effect of changing count_A1 to True to conform with plink
    #   standard (shouldn't change LD estimates)
    @memo.memoized
    def data(self, chrnum):
        return Bed(self.bfile(chrnum), count_A1=False)
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

    # return ldblock info, standardized X, snp metadata
    def block_data(self, ldblocks, c, meta=None, chunksize=15, genos=True, verbose=2):
        # restrict to ld blocks in this chr and process them in chunks
        chr_blocks = ldblocks[ldblocks.chr=='chr'+str(c)]
        # get metadata about snps
        snps = self.bim_df(c)

        t0 = time.time()
        for block_nums in pyit.grouper(chunksize, range(len(chr_blocks))):
            # get ld blocks in this chunk, and indices of the snps that start and end them
            chunk_blocks = chr_blocks.iloc[block_nums]
            blockstarts_ind = np.searchsorted(snps.BP.values, chunk_blocks.start.values)
            blockends_ind = np.searchsorted(snps.BP.values, chunk_blocks.end.values)
            if verbose >= 1:
                print('{} : chr {} snps {} - {}'.format(
                    time.time()-t0, c, blockstarts_ind[0], blockends_ind[-1]))

            # read in refpanel for this chunk, and find the relevant annotated snps
            if genos:
                Xchunk = self.stdX(c, (blockstarts_ind[0], blockends_ind[-1]))
                print('read in chunk')
            else:
                Xchunk = None

            if meta is not None:
                metachunk = meta.iloc[blockstarts_ind[0]:blockends_ind[-1]]

            # calibrate ld block starts and ends with respect to the start of this chunk
            blockends_ind -= blockstarts_ind[0]
            blockstarts_ind -= blockstarts_ind[0]
            for i, start_ind, end_ind in zip(
                    chunk_blocks.index, blockstarts_ind, blockends_ind):
                if verbose >= 2:
                    print(time.time()-t0, ': processing ld block',
                            i, ',', end_ind-start_ind, 'snps')
                if genos:
                    X = Xchunk[:, start_ind:end_ind]
                else:
                    X = None
                if meta is not None:
                    metablock = metachunk.iloc[start_ind:end_ind]
                else:
                    metablock = None
                yield (chr_blocks.loc[i], X, metablock,
                    metachunk.iloc[start_ind:end_ind].index)


# basic tests
if __name__ == '__main__':
    d = Dataset('/groups/price/yakir/data/datasets/1000G3.wim5unm/1000G3.wim5unm.')
    print(d.path)
    print(d.bfile(22))
    print(d.frq_file(22))
    print(d.frq_df(22).columns)
    print(d.bim_df(22).columns)
    print(d.M(22))
    print(d.totalM())
    print(d.totalM(chroms=[21,22]))
    print(d.N())
