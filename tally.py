#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

INDEXES1 = frozenset(['GGCTGTTG', 'AACACCAT', 'CGTGTCCA', 'GGGGAGGT',
   'CTGCACAC', 'CGAAACAG', 'GCGAGCGG', 'GGCAAGTC', 'GAATTCAA',
   'AGTGGATC', 'TCGACAAT', 'GTGCGAGA', 'GCGCGTCT', 'ATTTTTAG',
   'ACTTAGGG', 'TAGGGCTA', 'ACACGACA', 'CCATCACT', 'TATGCCCC',
   'ACTCCGTT', 'TCGTTTTC', 'GGTAAAGA', 'AAATTTGT', 'TCACCCAC',
   'ATGAACGA', 'GCAGCTAA', 'TTCTCACA', 'AGGTGCCC', 'ACCTCTTA',
   'CAATGCTC', 'CGCATTTA', 'CCCCAGAG', 'ATCTGTGC', 'AGTACCTA',
   'ACCTGCAG', 'CGTACTCG', 'CCAAAATA', 'TGTTATAA', 'TTAATCGG',
   'TCTAGCCT', 'GGTGCAAT', 'AGGCCGCG', 'CAAGACGT', 'TGGGGGAG',
   'AAAAGGGG', 'ATAAAACT', 'TAATGAAG', 'CATGCGGG', 'TTCCTTCT',
   'ACGCTATG', 'TATTTGGA', 'TGGACTGC', 'CAGAATCT', 'GTGGCACG',
   'GCGTTGAT', 'CTGAGTTG', 'AAGTCAGC', 'TTGTATGT', 'TTCGATTA',
   'CTCGTAGG', 'ACGGGTAC', 'CTATCGTA', 'CCGAGACC', 'GAGGATAG',
   'ATGGGGTT', 'AAGCAAAA', 'TATTGTTT', 'GCTACATC', 'GGGATCTT',
   'CACTATAC', 'GCTTCCGA', 'AGCTTACT', 'TCACAGGA', 'ATTCCAAC',
   'GACAGGCT', 'GATTTACG', 'GTTGTCTC', 'GTCCTGAC', 'CGCCGGCA',
   'TCCAACTG', 'TTCGCGGT', 'GTAGAGCC', 'TGTCCAGG', 'TGAGTCAT',
   'GCCATTGC', 'CCTCTGCC', 'CGCCCCGT', 'ATACATTG', 'TAACGTCC',
   'GAAGTGTG', 'AAGATGCA', 'GGACGATT', 'TACAGATC', 'GATCGCTG',
   'CTTGGCAT', 'CGGTTCGG',])

INDEXES2 = frozenset(['ATTCCAGT', 'TACCACAT', 'ACAAATAC', 'CGATAACA',
   'GCTTGTGC', 'GCAGCTCA', 'CAGGACAG', 'AATGGAAT', 'CATTCGGC',
   'GTCCTTAA', 'AGTTGCTT', 'GTGTCGCT', 'CTAATGCA', 'CCATAGTC',
   'GTCCGCGG', 'GACGATTA', 'ATTGATGA', 'TTCTTACG', 'TGGAAAGG',
   'TGGGAGCT', 'TCTTTGTT', 'AGAGGTCG', 'TCATGAAC', 'TAAACCTG',
   'GACTCTCG', 'AACTTGAC', 'CATACCCT', 'TAAGGGTA', 'TTTAGCAT',
   'GGAACTTC', 'CAGTCATA', 'CTGCTTGG', 'ACCGTCCT', 'AGGAATCA',
   'GTAGTGAG', 'GCACAACT', 'TGCGGTAA', 'CTGAAGTT', 'GTTAAACC',
   'CCAACGGT', 'ATGGTGGT', 'AGATTAGC', 'TTGCGATT', 'CACGGCTT',
   'CGGTGGCG', 'ACCAGGAA', 'TCTGCAGA', 'GTACGGTC', 'TTTATGGG',
   'GCGTCCAC', 'TCACTCCC', 'TTGCCCAG', 'CCCTTTCA', 'GATTTAAG',
   'GTGGTATC', 'GATGCCTC', 'GGCGCGAT', 'AACCAACC', 'AATTATTG',
   'ACGCGGGG', 'CCTCCTAT', 'AAACGAGA', 'TGAATTGA', 'CTTGGCCC',
   'GGTTGGAA', 'TCGTACTA', 'AGTCTTCC', 'ACTCACTC', 'TGTCATGT',
   'CTAGCTAC', 'AAATACGT', 'CCGGGGAT', 'ATAATTTT', 'GTTACCGA',
   'GAGAAGGC', 'TATATATC', 'CCCAAAGA', 'CAACAATG', 'AGGTCCGG',
   'TGCACGCG', 'AGCTAAAT', 'CCGAGACC', 'GCCAGTTG', 'AGGCCTTT',
   'CGTGTGTA', 'CACCCGAG', 'GAAAGCAA', 'ATGACAAC', 'TCCGATCG',
   'TATTTCCA', 'TAGTGTCT', 'AGCCGCAC', 'TTTGAGTC', 'GGGCTACG',
   'CTATCAGG', 'CGACGTGC',])


def parse_stc(f):
   canonical = dict()
   for line in f:
      brcd,count,variants = line.split()
      for variant in variants.split(','):
         canonical[variant] = brcd
   return canonical


def main(canonical, f):
   dict_of_barcodes = \
      defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
   for line in f:
      brcd,index1,index2,umi = line.split()
      if brcd not in canonical:
         continue
      if index1 not in INDEXES1 or index2 not in INDEXES2:
         continue
      brcd = canonical[brcd]
      dict_of_barcodes[brcd][(index1, index2)][umi] += 1
   for brcd,dict_of_indexes in dict_of_barcodes.items():
      for (index1, index2), C in dict_of_indexes.items():
         S = [a for a in C.values() if a > 2]
         if len(S) > 0:
            print brcd, index1, index2, len(S)


if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      canonical = parse_stc(f)
   with open(sys.argv[2]) as f:
      main(canonical, f)
