import time

seq = 'A'*10**5 + 'T'*10**5

d = {'A':0, 'T':0}

s1 = time.time()
for i in xrange(10**5): d[seq[i]] += 1
print time.time() - s1

s1 = time.time()
d['A'] += seq.count('A')
d['T'] += seq.count('T')
print time.time() - s1
		
