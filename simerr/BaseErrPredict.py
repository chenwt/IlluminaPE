import os,sys
from csv import DictReader
from Bio import SeqIO
import pdb
import rpy2
import rpy2.robjects as robj

class BaseErrPredict():
    def __init__(self, Rdata):
        """
        Rdata --- the .Rdata file that contains the ksvm 'model'
        """
        self.Rdata = Rdata
        robj.packages.importr("kernlab")
        robj.r.load(self.Rdata)
        self.model = robj.globalenv["model"]
        
    def predict(self, df, ignorefirst):
        """
        df --- a DataFrame object for the model
              which must contain the keys 
              'Phred', 'Cycle', 'B2', 'B1', 'B0'
        ignorefirst --- ignore the first X results
              because I added them to account for all A/T/C/G
        Returns: an iterator of +(correct)/-(error)
        """
        outcome = robj.r.predict(self.model, df)
        # outcome.levels '+' means correct base '-' means error
        # don't forget that the index in R is 1-based, so we need to -1
        for o in outcome[ignorefirst:]:
            yield outcome.levels[o-1]
            
    def compare_seq(self, seq1, seq2, qual1, qual2, cycle1, cycle2, primer, max_allowed_remain_mismatch, justmm):
        """
        Remember to add dummy values in the beginning?
        
        Return: <True|False>, <list of mismatch positions>, <update list of (pos, base, qual, cycle)>
        
        Will return True (i.e. yes to merge seq2 into seq1) if
        (0) no mismatches (after accounting for 'N's), or
        (1) justmm is True and # of mismatches <= max allowed, or
        (2) after correction, # of mismatches <= max allowed
        """
        correct = []
        poses = []
        for pos in xrange(min(len(seq1), len(seq2))):
            if seq1[pos] == 'N':
                correct.append((pos, seq2[pos], qual2[pos], cycle2[pos]))
            elif seq2[pos] == 'N':
                correct.append((pos, seq1[pos], qual1[pos], cycle1[pos]))
            elif seq1[pos]!=seq2[pos]:
                poses.append(pos)
            #if len(poses) >= 10:
            #    print >> sys.stderr, "HARDCODE FOR NOW: MISMATCH EXCEEDS 10. ABORT"
            #    return False, poses, []

        if len(poses) == 0: return True, poses, correct
        
        # justmm is True if we're comparing (super)seqs that are quite abundance
        # in this case we don't predict base errors, just merge if mm is small
        if justmm and len(poses) <= max_allowed_remain_mismatch:
            #raw_input("hey it's just mm!")
            return True, poses, correct 
        
        phreds = [0, 0, 0, 0]
        cycles = [0, 0, 0, 0]
        B2s = ['A','T',"C",'G']
        B1s = ['A','T','C','G']
        B0s = ['A','T','C','G']
        for pos in poses:
            phreds.append(qual1[pos])
            phreds.append(qual2[pos])
            cycles.append(cycle1[pos])
            cycles.append(cycle2[pos])
            B2s.append(seq1[pos-2] if pos>=2 else primer[-2])
            if B2s[-1] == 'N':
                B2s[-1] = seq2[pos-2] if seq2[pos-2]!='N' else 'A'
            B2s.append(seq2[pos-2] if pos>=2 else primer[-2])
            if B2s[-1] == 'N':
                B2s[-1] = seq1[pos-2] if seq1[pos-2]!='N' else 'A'            
            B1s.append(seq1[pos-1] if pos>=1 else primer[-1])
            if B1s[-1] == 'N':
                B1s[-1] = seq2[pos-1] if seq2[pos-1]!='N' else 'A'
            B1s.append(seq1[pos-1] if pos>=1 else primer[-1])
            if B1s[-1] == 'N':
                B1s[-1] = seq1[pos-1] if seq1[pos-1]!='N' else 'A'           
            B0s.append(seq1[pos])
            B0s.append(seq2[pos])
        
        df = robj.DataFrame({\
                        "Phred": robj.IntVector(phreds),\
                        "Cycle": robj.IntVector(cycles),\
                        "B2": robj.StrVector(B2s),\
                        "B1": robj.StrVector(B1s),\
                        "B0": robj.StrVector(B0s)})
        
        #print poses
        #print df
        
        outcome = [x for x in self.predict(df, 4)]
        #print outcome
        remain_mismatch = 0
        for i, pos in enumerate(poses):
            ans = outcome[2*i] + outcome[2*i+1]
            if ans == '+-': # don't really need to do anything since we're just using seq1 as ref
                pass#correct.append((pos, seq1[pos], qual1[pos], cycle1[pos]))
            elif ans == '-+' and not justmm: # will only accept the less abundant one if the size ratio is not big
                correct.append((pos, seq2[pos], qual2[pos], cycle2[pos]))
            elif ans == '--':
                correct.append((pos, 'N', 2, cycle1[pos]))
            else: # '++'
                remain_mismatch += 1
                #return False, correct
        if remain_mismatch <=0: #<= max_allowed_remain_mismatch:
            return True, poses, correct
        else:
            return False, poses, correct
            

    

