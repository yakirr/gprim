from __future__ import print_function, division
import os
import numpy as np
import pandas as pd
import itertools
from pyutils import memo, fs


# complementary bases
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# bases
BASES = COMPLEMENT.keys()
# true iff strand ambiguous
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
        for x in itertools.product(BASES, BASES)
        if x[0] != x[1]}
# SNPS we want to keep (pairs of alleles)
VALID_SNPS = {x for x in map(lambda y: ''.join(y), itertools.product(BASES, BASES))
        if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {
        x for x in map(lambda y: ''.join(y), itertools.product(VALID_SNPS, VALID_SNPS))
        # strand and ref match
        if ((x[0] == x[2]) and (x[1] == x[3])) or
        # ref match, strand flip
        ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
        # ref flip, strand match
        ((x[0] == x[3]) and (x[1] == x[2])) or
        # ref flip, strand flip
        ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}
# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {''.join(x):
        # ref flip, strand match
        ((x[0] == x[3]) and (x[1] == x[2])) or
        # ref flip, strand flip
        ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
        for x in MATCH_ALLELES}

# Checks if SNP columns are equal. If so, saves time by using concat instead of merge.
def smart_merge(x, y, how='inner'):
    if len(x) == len(y) and (x.SNP == y.SNP).all():
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True).drop('SNP', axis=1)
        out = pd.concat([x, y], axis=1)
    else:
        out = pd.merge(x, y, how=how, on='SNP')
    return out

# keeps only snps in ref, flips alleles in df to match those in ref, and sets value at
# non-matching snps to a missing_val (0 by default). 
def reconciled_to(ref, df, colnames, signed=True, missing_val=0):
    result = smart_merge(ref[['CHR','BP','SNP','CM','A1','A2']],
        df[['SNP','A1','A2']+list(colnames)].rename(columns={'A1':'A1_df','A2':'A2_df'}),
        how='left')

    # snps in ref but not in df
    missing = result.A1_df.isnull()
    result.loc[missing,'A1_df'] = result.loc[missing,'A1']
    result.loc[missing,'A2_df'] = result.loc[missing,'A2']
    result.loc[missing,colnames] = missing_val

    # assign to zero SNPs in ref whose alleles don't match df
    a1234 = (result.A1+result.A2+result.A1_df+result.A2_df).apply(lambda y: y.upper())
    match = a1234.apply(lambda y: y in MATCH_ALLELES)
    n_mismatch = (~match).sum()
    print('of', len(match), 'snps,', n_mismatch, 'snps do not have same alleles,')
    result.loc[~match,colnames] = missing_val
    result.loc[~match,'A1_df'] = result.loc[~match,'A1']
    result.loc[~match,'A2_df'] = result.loc[~match,'A2']

    if signed:
        a1234 = (result.A1+result.A2+result.A1_df+result.A2_df).apply(lambda y: y.upper())
        a1234[~match] = list(MATCH_ALLELES)[0] # in case ref has strand-ambiguous alleles
        flip = a1234.apply(lambda y: FLIP_ALLELES[y])
        n_flip = flip.sum()
        print('of', len(match), 'snps,', n_flip, 'snps need flipping')
        result.loc[flip,colnames] *= -1

    return result.drop(['A1_df', 'A2_df'], axis=1)


class Annotation(object):
    def __init__(self, stem_chr, signed=True):
        self.stem_chr = stem_chr
        self.signed = signed

    def filestem(self, chrnum='', mkdir=False):
        fname = '{}{}'.format(self.stem_chr, chrnum)
        if mkdir:
            fs.makedir_for_file(fname)
        return fname
    def annot_filename(self, chrnum):
        return self.filestem(chrnum) + '.annot.gz'
    def sannot_filename(self, chrnum, mkdir=False):
        return self.filestem(chrnum, mkdir) + '.sannot.gz'
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
    a = Annotation('/groups/price/yakir/data/simannot/1000G3.wim5unm/mock/')
    print(a.names(22))
    print(a.sannot_df(22))
    print(a.sqnorms(22))
    print(a.sizes(22))
    print(a.total_sqnorms([1,22]))
    print(a.total_sizes([1,22]))

    ref = pd.read_csv('/groups/price/yakir/data/simannot/1000G3.wim5unm/mock/22.sannot.gz',
            delimiter='\t').drop('ANNOT', axis=1).iloc[:20]

    ref.ix[0,'A1'] = a.sannot_df(22).ix[0,'A2']
    ref.ix[0,'A2'] = a.sannot_df(22).ix[0,'A1']
    ref.ix[1,'A1'] = 'A'; ref.ix[1,'A2'] = 'C'
    print(ref)

    print(reconciled_to(ref, a.sannot_df(22), ['ANNOT']))
