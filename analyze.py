#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

import seeq

from gzopen import gzopen

def main(f):
   constant = seeq.compile('TATAGTGAGTCGTATTAAAAGCGAAAGGGAAACCAGAGGAGC', 5)
   for lineno,line in enumerate(f):
      if lineno % 4 == 0:
         index2 = re.sub(r'.*\+', '', line.rstrip())
      elif lineno % 4 == 1:
         m = constant.match(line.rstrip())
         if m is None:
            continue
         try:
            barcode, ignore, tail = m.tokenize()
         except ValueError:
            continue
         if len(tail) < 8 or len(barcode) < 14:
            continue
         UMI = tail[:4]
         index1 = tail[4:12]
         print barcode, index1, index2, UMI

if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      main(f)
