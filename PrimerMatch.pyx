"""
Primer matching using suffix arrays allowing mismatches and indels
"""
import os, sys
import pdb
import c_tools_karkkainen_sanders as tks
from collections import namedtuple

MatchResult = namedtuple('MatchResult', ['primer_start', 'match_len', 'miss'])

def mismatch_exceeded(s1, s2, max_mm_allowed):
    mm = 0
    for i in xrange(min(len(s1),len(s2))):
        if s1[i]!=s2[i]:
            mm += 1
            if mm > max_mm_allowed:
                return mm
    return mm

def match_primer_len(seq, primer, max_mm_allowed, min_overlap, is_reverse):
    l = len(primer)
    if not is_reverse: # forward primer
        for i in xrange(l-min_overlap+1):
            mm = mismatch_exceeded(primer[i:], seq, max_mm_allowed)
            if mm <= max_mm_allowed:# found a match!
                return l-i, mm
    else:
        for i in xrange(l-min_overlap+1):
            mm = mismatch_exceeded(primer[:l-i], seq[-(l-i):], max_mm_allowed)
            if mm <= max_mm_allowed:
                return l-i, mm
    return 0, -1


def uperm(s, symbols, acc=''):
    """
    Unique permutation
    
    ex: if we want unique permutation of 'AABBC'
    call with uperm([1,2,1], 'ABC')
    """
    if sum(s) == 0: yield acc
    else:
        for i, x in enumerate(s):
            if x > 0:
                s[i] -= 1
                for p in uperm(s, symbols, acc + symbols[i]): yield p
                s[i] += 1

class PrimerMatch:
    def __init__(self, primer):
        """
        primer -- should be in DNA form (ATCG)
        """
        self.primer = primer
        self.N = len(self.primer) # length of primer
        self.isprimer = lambda x: x < self.N
        
        self.s = None # primer$input$
        self.sa = None # suffix array
        self.lcp = None # lcp array
        self.match_result = None # MatchResult
        
    def clear(self):
        self.s = None
        self.sa = None
        self.lcp = None
        self.match_result = None
        
    def make_suffix(self, seq): 
        """
        seq -- should be a DNA sequence
        """
        seq = seq[:self.N]
        M = self.N * 2 + 2 # length of S
        s = self.primer + '$' + seq + '$'
        sa = tks.simple_kark_sort(s)
        lcp = [[0] * (M + 1) for i in xrange(M)]
        new_lcp = [[0] * (M + 1) for i in xrange(M)]
        # LCP[i][j] = longest common prefix for suffix S_i and S_j
        for i, x in enumerate(tks.LCP(s, sa)):
            lcp[i][i + 1] = x
            new_lcp[sa[i]][sa[i + 1]] = x
            new_lcp[sa[i + 1]][sa[i]] = x
        for i in xrange(2, M - 1): # can ignore first two cuz must be $ and $seq
            for j in xrange(i + 2, M):
                lcp[i][j] = min(lcp[i][j - 1], lcp[j - 1][j])
                new_lcp[sa[i]][sa[j]] = lcp[i][j]
                new_lcp[sa[j]][sa[i]] = lcp[i][j]
                if lcp[i][j] == 0: # no need to do j+1... since it'll all be zero
                    break

        self.s = s
        self.sa = sa
        self.lcp = new_lcp
        #return s, sa, new_lcp, lcp
    
    def match_by_pattern(self, pattern, min_match_len):
        """
        pattern --- list of {'M', 'I', 'D'} #mismatch, insertion, deletion
        min_match_len --- this much of primer should at least be in the seq
        
        Should be called by match()
        
        Returns MatchResult if has match, otherwise None
        """
        lcp = self.lcp
        n = self.N
        m = self.N * 2 + 2
        p_i = 0 # index on pattern
        for start in xrange(n - min_match_len + 1):
            p_i = 0
            i = start
            k = self.N + 1 # the index at which seq starts in S
            miss = [] # (type[M|I|D], position)
            while i < n and k < m:
                #pdb.set_trace()
                if lcp[i][k] > 0: # primer and seq matches for this long
                    i, k = i + lcp[i][k], k + lcp[i][k]
                    if i >= n - 1: # it's a match! return it!
                        # sometimes the match includes the $s, remember to subtract
                        match_len = k - self.N - 1 - (k==m)
                        return MatchResult(primer_start=start, match_len=match_len,\
                                           miss=miss)
                # if we're still here, no finished matching
                try:
                    p = pattern[p_i] 
                    p_i += 1
                except IndexError: # this means pattern is empty
                    break
                if p == 'M': # we're allowing a mismatch, advance i & k
                    miss.append(('M', i))
                    i += 1
                    k += 1
                elif p == 'D': # we're allowing an deletion w.r.t to primer[i]
                    miss.append(('D', k))
                    i += 1
                else: # we're allowing an insertion at seq[k]
                    miss.append(('I', i))
                    k += 1
            #print start, i, miss
        return None
            
    def match(self, min_match_len, max_mm, max_in, max_de):
        """
        min_match_len -- minimum match length to primer (probably 10)
        max_mm -- max mismatch allowed (for primer match, should be 0-2)
        max_in -- max insertion w.r.t primer allowed (probably 0)
        max_de -- max deletion w.r.t primer allowed (probably 0-1)
        
        Given S, SA, LCP which should already be computed via make_suffix
        Find the smallest i s.t. primer[i:] matches seq[:N-i] 
        satisfying mismatch/indel pattern
        
        Returns: MatchResult if has match, otherwise None
        
        NOTE: will find the longest match which sometimes means
             includes extra mismatches/indels which shouldn't matter
             if we're just trying to remove primers from reads
        """
        # todo: check that S, SA, LCP has been computed?
        for pattern in uperm([max_mm, max_de, max_in], 'MDI'):
            match_result = self.match_by_pattern(pattern, min_match_len)
            if match_result is not None: # yeah! a match!
                self.match_result = match_result
                return match_result
            else:
                self.match_result = None
                     
    def __print__(self):
        """
        Prints the current S, SA, LCP
        """
        for i in xrange(len(self.s) - 1):
            print i, self.sa[i], self.lcp[self.sa[i]][self.sa[i + 1]], self.s[self.sa[i]:], self.isprimer(self.sa[i])
                                                                                                    
    def print_match(self):
        """
        match_result should exist after calling match()
        """
        if self.match_result is None:
            print "Match is empty"
        else:
            miss = self.match_result.miss
            miss.sort(key=lambda x:x[1])
            i = self.match_result.primer_start
            j = self.N + 1
            s1 = ''
            s2 = ''
            for type, pos in miss:
                if pos > self.N:
                    if type == 'D':
                        s2 += self.s[j:pos] + '-'
                        j = pos
                elif type == 'I':
                    s1 += self.s[i:pos] + '-'
                    i = pos
            s1 += self.s[i:self.N]
            s2 += self.s[j:]
        print "Primer:", s1
        print "Input :", s2
            
                    
                
            

