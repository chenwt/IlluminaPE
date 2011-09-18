import os
import sys
import csv
from DigitalFingerprint.DF import DFReader, DFWriter

DF_file = 'DS19175toDS19354.composite.subsample500000.DF'
name_dict = dict((r['Sample'],r['Name']) for r in csv.DictReader(open('sample.annotation.txt'), delimiter='\t'))
note_dict = dict((r['Sample'],r['Notes']) for r in csv.DictReader(open('sample.annotation.txt'), delimiter='\t'))

pyrooverlap = [line.strip() for line in open('sample.overlapWithPyroNames.txt')]

df_to_use = []
for df in DFReader(open(DF_file)):
	name = name_dict[df.name]
	if note_dict[df.name] in pyrooverlap:
		df.name = note_dict[df.name] + '-' + name
		df_to_use.append(df)

w = DFWriter(sys.stdout)
w.writes(df_to_use)
