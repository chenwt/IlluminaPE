import os,sys,random
from miscBowTie import FastqReader

input = sys.argv[1]
sum = 0
count = 0
for r in FastqReader(input):
	if random.random() >= 0.0001:
		continue
	for q in r['qual']:
		sum += ord(q)-33
		count += 1
print sum*1./count, sum, count
