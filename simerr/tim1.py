r1={'seq': 'ACTCCTACGGGAGNCAGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTANTTNTGTATGTAAAGCTCN', 'qual': '<@@DDDDDFFFCF#######################################################################################', 'offset': 0, 'ref': 'Unc06znl', 'ID': 'HWI-ST700693:106:D05FJACXX:1:1101:21360:2964 1:N:0:ATCACG/1', 'strand': '+'}
r2={'seq': 'AGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGAAAGAAAATGACGGTACCTGACTAAGAAGCCNCGGCTAACTACGTGCCAGCAGCCGCGGTAATAN', 'qual': '####################################################################################################', 'offset': 75, 'ref': 'Unc06znl', 'ID': 'HWI-ST700693:106:D05FJACXX:1:1101:21360:2964 2:N:0:ATCACG/2', 'strand': '-'}


import time
import composite
import hello
from cPickle import *
base_freq = {'A':.25,'T':.25,'C':.25,'G':.25}

s1 = time.time()
for i in xrange(10**3):
	composite.compose2(r1, r2, base_freq)
print time.time() - s1

s1 = time.time()
for i in xrange(10**3):
	hello.compose2(r1, r2, base_freq)
print time.time() - s1

