#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

# Reverse complement of the primer sequence.
RT_INDEXES = frozenset(['GGCTGTTG', 'AACACCAT', 'CGTGTCCA',
   'GGGGAGGT', 'CTGCACAC', 'CGAAACAG', 'GCGAGCGG', 'GGCAAGTC',
   'GAATTCAA', 'AGTGGATC', 'TCGACAAT', 'GTGCGAGA', 'GCGCGTCT',
   'ATTTTTAG', 'ACTTAGGG', 'TAGGGCTA', 'ACACGACA', 'CCATCACT',
   'TATGCCCC', 'ACTCCGTT', 'TCGTTTTC', 'GGTAAAGA', 'AAATTTGT',
   'TCACCCAC', 'ATGAACGA', 'GCAGCTAA', 'TTCTCACA', 'AGGTGCCC',
   'ACCTCTTA', 'CAATGCTC', 'CGCATTTA', 'CCCCAGAG', 'ATCTGTGC',
   'AGTACCTA', 'ACCTGCAG', 'CGTACTCG', 'CCAAAATA', 'TGTTATAA',
   'TTAATCGG', 'TCTAGCCT', 'GGTGCAAT', 'AGGCCGCG', 'CAAGACGT',
   'TGGGGGAG', 'AAAAGGGG', 'ATAAAACT', 'TAATGAAG', 'CATGCGGG',
   'TTCCTTCT', 'ACGCTATG', 'TATTTGGA', 'TGGACTGC', 'CAGAATCT',
   'GTGGCACG', 'GCGTTGAT', 'CTGAGTTG', 'AAGTCAGC', 'TTGTATGT',
   'TTCGATTA', 'CTCGTAGG', 'ACGGGTAC', 'CTATCGTA', 'CCGAGACC',
   'GAGGATAG', 'ATGGGGTT', 'AAGCAAAA', 'TATTGTTT', 'GCTACATC',
   'GGGATCTT', 'CACTATAC', 'GCTTCCGA', 'AGCTTACT', 'TCACAGGA',
   'ATTCCAAC', 'GACAGGCT', 'GATTTACG', 'GTTGTCTC', 'GTCCTGAC',
   'CGCCGGCA', 'TCCAACTG', 'TTCGCGGT', 'GTAGAGCC', 'TGTCCAGG',
   'TGAGTCAT', 'GCCATTGC', 'CCTCTGCC', 'CGCCCCGT', 'ATACATTG',
   'TAACGTCC', 'GAAGTGTG', 'AAGATGCA', 'GGACGATT', 'TACAGATC',
   'GATCGCTG', 'CTTGGCAT', 'CGGTTCGG',])

# Reverse complement of the primer sequence.
PCR_INDEXES = frozenset(['ATTCCAGT', 'TACCACAT', 'ACAAATAC',
   'CGATAACA', 'GCTTGTGC', 'GCAGCTCA', 'CAGGACAG', 'AATGGAAT',
   'CATTCGGC', 'GTCCTTAA', 'AGTTGCTT', 'GTGTCGCT', 'CTAATGCA',
   'CCATAGTC', 'GTCCGCGG', 'GACGATTA', 'ATTGATGA', 'TTCTTACG',
   'TGGAAAGG', 'TGGGAGCT', 'TCTTTGTT', 'AGAGGTCG', 'TCATGAAC',
   'TAAACCTG', 'GACTCTCG', 'AACTTGAC', 'CATACCCT', 'TAAGGGTA',
   'TTTAGCAT', 'GGAACTTC', 'CAGTCATA', 'CTGCTTGG', 'ACCGTCCT',
   'AGGAATCA', 'GTAGTGAG', 'GCACAACT', 'TGCGGTAA', 'CTGAAGTT',
   'GTTAAACC', 'CCAACGGT', 'ATGGTGGT', 'AGATTAGC', 'TTGCGATT',
   'CACGGCTT', 'CGGTGGCG', 'ACCAGGAA', 'TCTGCAGA', 'GTACGGTC',
   'TTTATGGG', 'GCGTCCAC', 'TCACTCCC', 'TTGCCCAG', 'CCCTTTCA',
   'GATTTAAG', 'GTGGTATC', 'GATGCCTC', 'GGCGCGAT', 'AACCAACC',
   'AATTATTG', 'ACGCGGGG', 'CCTCCTAT', 'AAACGAGA', 'TGAATTGA',
   'CTTGGCCC', 'GGTTGGAA', 'TCGTACTA', 'AGTCTTCC', 'ACTCACTC',
   'TGTCATGT', 'CTAGCTAC', 'AAATACGT', 'CCGGGGAT', 'ATAATTTT',
   'GTTACCGA', 'GAGAAGGC', 'TATATATC', 'CCCAAAGA', 'CAACAATG',
   'AGGTCCGG', 'TGCACGCG', 'AGCTAAAT', 'CCGAGACC', 'GCCAGTTG',
   'AGGCCTTT', 'CGTGTGTA', 'CACCCGAG', 'GAAAGCAA', 'ATGACAAC',
   'TCCGATCG', 'TATTTCCA', 'TAGTGTCT', 'AGCCGCAC', 'TTTGAGTC',
   'GGGCTACG', 'CTATCAGG', 'CGACGTGC',])


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
      brcd,idx1,idx2,umi = line.split()
      if brcd not in canonical:
         continue
      if idx1 not in RT_INDEXES or idx2 not in PCR_INDEXES:
         continue
      brcd = canonical[brcd]
      dict_of_barcodes[brcd][(idx1, idx2)][umi] += 1
   for brcd, dict_of_indexes in dict_of_barcodes.items():
      for (idx1, idx2), C in dict_of_indexes.items():
         S = [(a,b) for (a,b) in C.items() if b > 2]
         if len(S) > 0:
            UMIs = ' '.join(sorted([a for (a,b) in S]))
            sys.stdout.write('%s %s %s %d %s\n' % \
                  (brcd, idx1, idx2, len(S), UMIs))


if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      canonical = parse_stc(f)
   with open(sys.argv[2]) as f:
      main(canonical, f)
