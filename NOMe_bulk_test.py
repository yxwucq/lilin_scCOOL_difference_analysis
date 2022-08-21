# -*- coding: utf-8 -*-

import fractions
import gzip
import os
import sys
from turtle import position

import numpy as np
import pandas as pd
import scipy.stats

help_text = """
Usage: python "input_bedgraph.gz" "output.bedgraph"
Example: python NOMe_bulk_test.py /datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/singleC/chr10.ACG.TCG.bed.gz 
Method: Using shift window to slide on chr, and make calls on each region
"""

def show_help():
    print(help_text, file=sys.stderr)

class NDR(object):
    """
    Peaks calling within a chromosome.
    """
    def __init__(self, chrom, methyl_fraction, p_value_cutoff=1e-20, bin_len=100, step_len=20, dist_len=140):
        self.chrom = chrom # Specifying the chromosome working on
        self.bin_len = bin_len # Peaks calling bin length
        self.step_len = step_len # Sliding step length
        self.dist_len = dist_len # Peaks width lower limit
        self.p_value_cutoff = p_value_cutoff # P-value cutoff for pearson's chi-squared test
        
        self.np_RATIO = np.array([1-methyl_fraction, methyl_fraction])
        
        self.NDR_info = {} # Contains the total information of a NDR
        
        self.l_umt_met = {} # Contains the methylation ratio of a region on Chr

    def parse_line(self, position_beginning, unmet_cnt, met_cnt):
        self.__get_bins(position_beginning)
        self.__NDR_append(unmet_cnt, met_cnt, position_beginning)
        self.__NDR_bin_detect()
        
    def last_line(self):
        if len(self.NDR_info) != 0:
            if self.__NDR_len() > self.dist_len:
                self.__NDR_report()

    def __NDR_bin_detect(self):
        """
        Check if there are existing significant peaks
        """
        for idx in sorted(self.l_umt_met):
            self.__bin_info(idx)
            
            """对这种可能显著并且不会继续落 reads 的 bin"""
            if idx < self.bin_idx_min:
                """确定NDR和现有bin 不overlap，看看是否足够长，打印出来，否则可以扔掉"""
                if len(self.NDR_info) != 0 and self.bin_begin > self.NDR_info['end']:
                    if self.__NDR_len() >= self.dist_len-1:
                        self.__NDR_report()
                            
                    self.__NDR_remove()
                        
                    del self.l_umt_met[idx]
                    continue
                                    
                pval = self.__chisq_bin(idx)
                """如果bin显著，无NDR则该bin作为NDR的基点，有NDR则和现有bin必有overlap（否则上一步被截住），对NDR elongate即可"""
                if pval > self.p_value_cutoff:
    #                    self.__NDR_remove()
                    del self.l_umt_met[idx]
                    continue
                else:
                    if 'beg' not in self.NDR_info:
                        self.NDR_info = {
                            'beg': self.bin_begin, 'end': self.bin_endin, 'idx':[idx], 'log_pval':[np.log10(pval)], 
                            'umt':[self.l_umt_met[idx][0]], 'met':[self.l_umt_met[idx][1]]
                        }
                    else:
                        if idx not in self.NDR_info['idx']:
                            self.NDR_info['end'] = self.bin_endin
                            self.NDR_info['idx'].append(idx)
                            self.NDR_info['log_pval'].append(np.log10(pval))
                            self.NDR_info['umt'].append(self.l_umt_met[idx][0])
                            self.NDR_info['met'].append(self.l_umt_met[idx][1])

    def __NDR_remove(self):
        if 'idx' in self.NDR_info:
            for idx_rm in self.NDR_info['idx']:
                if idx_rm in self.l_umt_met:
                    del self.l_umt_met[idx_rm]
        
        self.NDR_info = {}

    def __NDR_append(self, unmet_cnt, met_cnt, position_beginning):
        """
        Load count in the bins
        """
        for i in range(self.bin_idx_min, self.bin_idx_max+1):
            if i not in self.l_umt_met:
                self.l_umt_met[i] = np.zeros(2)
            self.l_umt_met[i][0] += unmet_cnt
            self.l_umt_met[i][1] += met_cnt

    def __NDR_len(self):
        return self.NDR_info['end'] - self.NDR_info['beg']

    def __NDR_report(self):
#        str_p = ",".join( ["%1.2e" % (pval) for pval in self.NDR_info['pval']] )

#        l_val = np.unique(np.array(sorted(["%d" % (idx) for idx in self.NDR_info['idx']]),dtype="int"))
#        np_idx = np.array([ self.NDR_info['idx'].index(i) for i in l_val ],dtype="int")

        np_idx = np.array(self.NDR_info['idx'], dtype="int")

        np_pval= np.array(self.NDR_info['log_pval'])
        np_umt = np.array(self.NDR_info['umt'])
        np_met = np.array(self.NDR_info['met'])
        np_length = self.NDR_info['end'] - self.NDR_info['beg'] +1

        str_p = ",".join( ["%1.2f" % (pval) for pval in np_pval] )
        str_u = ",".join( ["%d" % (umt)  for umt  in np_umt ] )
        str_m = ",".join( ["%d" % (met)  for met  in np_met ] )
        print(self.chrom, self.NDR_info['beg'], self.NDR_info['end'], np_length ,str_p, str_u, str_m, sep='\t')

    def __chisq_bin(self, idx):
        np_obs    = self.l_umt_met[idx]
        np_exp    = self.np_RATIO * np.sum(np_obs)
        chisquare = np.sum((np_obs - np_exp)**2/np_exp)
        pval = 1
        '''
            只考虑 umet的reads 足够少的情况，这种情况下 reads 全被甲基化，即未得到组蛋白足够的保护，NDR
        '''
        if np_obs[0] < np_exp[0]:  
            pval = scipy.stats.chi2.sf(chisquare, 1)
#        print self.bin_begin, self.bin_endin, np_obs, np_exp, pval, idx
        return pval

    def __get_bins(self, pos):
        """
        Get bin index between the site
        """
        self.bin_idx_min = int(np.ceil((pos - self.bin_len) / self.step_len))
        self.bin_idx_max = int(pos/self.step_len)
        if self.bin_idx_min < 0: self.bin_idx_min = 0
    
    def __bin_info(self, in_idx):
        """bin to pos"""
        self.bin_begin = (in_idx) * self.step_len + 1
        self.bin_endin = (in_idx+1) * self.step_len + self.bin_len


def calculate_met_fraction(input_file_path):
    df = pd.read_csv(input_file_path, sep=' ', header=None)
    df_sum = df.groupby(0)[[3,4]].sum()
    df_sum['fraction'] = df_sum[4] / (df_sum[3]+df_sum[4])
    return df_sum

def main():
    try:
        DEPTH = None
        input_file_path = sys.argv[1] # Input bedgraph file for calling
        if len(sys.argv) == 3:
            DEPTH = int(sys.argv[2]) # Specify the depth
    except IndexError:
        show_help()
        sys.exit(1)
    # log_p_value_cutoff = -10
    met_fraction_df = calculate_met_fraction(input_file_path)
    working_chr = None
    m_NDR = None
    input_file = open(input_file_path, "r")
    for line in input_file:
        line = line.strip("\n").split()
        if (not working_chr) or line[0] != working_chr:
            working_chr = line[0]
            methyl_fraction = met_fraction_df.loc[working_chr,'fraction']
            if m_NDR:
                m_NDR.last_line() 
            m_NDR = NDR(working_chr, methyl_fraction, p_value_cutoff=1e-20, bin_len=100, step_len=20, dist_len=140)
        
        position_beginning = int(line[1])
        unmet_cnt, met_cnt = float(line[3]), float(line[4])
        
        if DEPTH and (met_cnt+unmet_cnt) < DEPTH: continue
        m_NDR.parse_line(position_beginning, unmet_cnt, met_cnt)
    
    input_file.close()
    m_NDR.last_line()
    
if __name__ == '__main__':
    main()
