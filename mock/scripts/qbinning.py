#!/usr/bin/env python

'''
Script para aplicar binning a los valores de calidad de un fastq.
Toma como entrada un fastq (opcionalmente comprimido con gzip) y devuelve un fastq con menos valores de calidad.
Usage: qbinning.py file.fq > output.fq
Lucía Balestrazzi - 2021
'''

import sys
import os
import gzip

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

try:
    file = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")

if not os.path.exists(file):
    raise SystemError("Error: File does not exist\n")

if is_gz_file(file):
    fh = gzip.open(file, 'rt')
else:
    fh = open(file, 'r')

lines = []
for line in fh:
    lines.append(line.rstrip())
    if len(lines) == 4:
        record = process(lines)
        qual = record['quality']
        newqual=""
        for base in qual:
          qsc=ord(base)-33
          if qsc <= 7:
            newqual = newqual + chr(5+33)
          elif 8 <= qsc <= 13:
            newqual = newqual + chr(12+33)
          elif 14 <= qsc <= 19:
            newqual = newqual + chr(18+33)
          elif qsc >= 20:
            newqual = newqual + chr(24+33)
        record['quality'] = newqual
        for key in record.keys():
            print(record[key])
        lines = []
