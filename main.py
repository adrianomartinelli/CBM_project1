#%%
def create_sam_header(filename, RNAME, RLEN, version=1.4, sorting_order="unsorted"):

    # creating a new SAM file
    sam_file = open(filename, 'w')

    # adding headers
    sam_file.write("@HD\tVN:%.1f\tSO:%s\n" % (version, sorting_order))
    sam_file.write("@SQ\tSN:%s\tLN:%d\n" % (RNAME, RLEN))
    sam_file.close()


def append_sam_alignment(filename, reverse=False, secondary=False, QNAME="*", RNAME="*", POS=0, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*"):

    # calculating the flag
    FLAG = 0
    if reverse: FLAG += 16
    if secondary: FLAG += 256

    # open existing SAM file append one linear alignment as a line to the alignment body
    sam_file = open(filename, "a")
    sam_file.write("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL))
    sam_file.close()

class Genome:

    def __init__(self, file):
        self.read_genome(file)
        self.size = len(self.genome)
        self.alphabet = ['$','A','C','G','T']
        self.inverse_alphabet = ['$','T','G','C','A']
        self.bwt_index()
        self.get_C()
        self.get_Occ()

    def read_genome(self, file):
        with open(file, 'r') as f:
            line = f.readline()
            while line:
                if(line.startswith('>')):
                    self.genome = f.read().replace('\n', '')
                    line = f.readline()
                else:
                    line = f.readline()

    def rot_left(self,s):
        return(s[1:] + s[0])

    def bwt_index(self):
        s = self.genome + '$'
        m = [s]
        for i in range(len(s)-1):
            m.append(self.rot_left(m[i]))

        self.A = list(np.argsort(m))
        m = [m[i] for i in self.A]
        self.L = [i[-1] for i in m]
        self.F = [i[0] for i in m]

        #update alphabet size which is given by character in L
        #TODO this section seems not to be necessary
        #self.__alphabet = list(set(G.L))
        #idx = np.argsort(self.__alphabet)
        #self.__alphabet = [self.__alphabet[i] for i in idx]

    def get_C(self):
        self.C={}
        for c in self.alphabet:
            occ=0
            for l in self.L:
                if l<c:
                    occ+=1
            self.C[c] = occ
    
    def get_Occ(self):
        self.Occ = np.full((len(self.C), len(self.L)),-1)
        for i, c in enumerate(self.alphabet):
            occ = 0
            for j, l in enumerate(self.L):
                if l == c:
                    occ+=1
                self.Occ[(i, j)] = occ
    
    def get_FM(self, P):
        #check if pattern is feasible
        #TODO delete in final version
        for i in range(len(P)):
            if P[i] not in self.alphabet:
                return (-1, -1)

        keys = list(self.C.keys())
        i = len(P)
        c = P[i-1]
        sp = self.C[c] + 1

        idx = keys.index(c)
        if idx+1 == len(keys):
            ep = len(self.L)
        else:
            ep = self.C[keys[idx+1]]

        while sp <= ep and i >= 2:
            c = P[i-2]
            idx = keys.index(c)
            sp = self.C[c] + self.Occ[idx,sp-2] + 1
            ep = self.C[c] + self.Occ[idx, ep-1]
            i -= 1
        return sp-1, ep-1

    def get_genomeIdx(self, sp, ep):
        #if ep > sp range will be empty
        #this corresponds to no hits in the genome
        return([self.A[i] for i in range(sp,ep+1)])

    def map_kmers(self, kmers):
        del_kmer = []
        for kmer in kmers.keys():
            sp, ep = self.get_FM(kmer)
            if(sp > ep):
                del_kmer.append(kmer)
            else:
                kmers[kmer][1].extend(self.get_genomeIdx(sp,ep))
        
        for kmer in del_kmer:
            del kmers[kmer]
        
        return kmers

    def align(self, substring, string, cost, reverse_complement = False, constrain = 2):
        alphabet = self.alphabet[1:]

        C = np.full((len(substring)+1, len(string)+1), 10000)
        C[0,:] = 0
        C[:, 0] = np.arange(0, len(substring)+1)

        for row in range(1,C.shape[0]):
            for col in range(max(row-constrain, 1),min(C.shape[1], row+constrain+1)):
                idx1 = alphabet.index(substring[row-1])
                idx2 = alphabet.index(string[col-1])

                match = C[row-1,col-1] + cost[idx1, idx2]
                mismatch = C[row-1,col-1] + cost[idx1, idx2]
                gap1 = C[row-1, col] + 2
                gap2 = C[row, col-1] + 2

                C[row,col] = min(match, mismatch, gap1, gap2)
        return min(C[len(substring),:])
        #return C

def read_line(f):
    line = f.readline()
    
    if(line.startswith('@')):
        #Process id line
        line = line.strip('@').strip('\n')
        identifier, origin = line.split('/')
        origin = int(origin)

        #Process seq and skip next line which is a +
        seq = f.readline().strip('\n')
        f.readline()

        #Process phread
        phread = f.readline().strip('\n')

    else:
        seq = ''
        phread = ''
        identifier = ''
    #return line
    return seq, phread, identifier


class Func:
    def __init__(self):
        self.alphabet = ['A','C','G','T']
        self.inverse_alphabet = ['T','G','C','A']

    def complement(self, string):
        idx = [self.inverse_alphabet[self.alphabet.index(i)] for i in string]
        return ''.join(idx)

    def generate_kmers(self, seq, k, shift = 1, reverse_complement = False):
        #generates a tuple of ([list1], [list2]) where the list1 contains
        #all the position of this kmer in the read and list2 is empty
        #list2 will ultimatively filled by the function map_kmers() and will
        #contain all the positions in the genome where this kmer is mapped to
        if reverse_complement:
            s = self.complement(seq[::-1])
        else:
            s = seq

        kmer = {}
        for i in range(0,len(s)-k+1, shift):
            if s[i:i+k] in kmer:
                kmer[s[i:i+k]][0].append(i)
            else:
                kmer[s[i:i+k]] = ([i],[])
        
        return kmer

    def get_start_index(self, kmers):
        #Given the position of kmers in the read and where they mapped in the genome, this function computes
        #what would be the starting position of the read in the genome.
        idx = {}
        for kmer in kmers.keys():
            #All combinations of starting positions
            comb = list(itertools.product(kmers[kmer][0], kmers[kmer][1]))
            for i in comb:
                start = i[1]-i[0]
                if str(start) in idx:
                    idx[str(start)] += 1
                else:
                    idx[str(start)] = 1

        #Generate sorted tulple list out of dict
        idx = sorted(idx.items(), key=lambda kv: kv[1], reverse=True)
        idx = [int(i[0]) for i in idx]
        return idx

    def get_feasible_comb(self, start1, start2, max_distance):
        comb = []
        #look at all combinations within a range
        for i in range(len(start1)):
            for j in range(len(start2)):
                if (abs(start1[i] - start2[j]) < max_distance): #TODO should we also check that start2 > start1??
                    comb.append((i,j))



        return comb
#%%
import os
import numpy as np
import itertools
import collections
import time
#%%
f = Func()
#G = Genome('genome.txt')

start = time.time()
G = Genome('./data_small/genome.chr22.5K.fa')
end = time.time()
print('Genome read in (duration:)', end-start)

cost = np.full((4,4),1)
cost[np.diag_indices(len(cost))] = 0

#f1 = open('read.txt', 'r')
f1 = open('./data_small/output_tiny_30xCov1.fq', 'r')
#f2 = open('read2.txt', 'r')
f2 = open('./data_small/output_tiny_30xCov2.fq', 'r')

seq1, phread1, QNAME1 = read_line(f1)
seq2, phread2, QNAME2 = read_line(f2)

filename = "test.sam"
ref_name = "22_5K"
create_sam_header(filename, RNAME=ref_name, RLEN=100)

######################################
#Optimise those parameters for run time optimisation

max_dist = 500 #the maximum distance between two kmers such that their are considered to come from the same DNA fragment
k = 16 #length of a kmer
shift = int(k/2) #shift between kmers
constrain = 5 #maximum number of allowed gaps during alignment

######################################

start = time.time()
while seq1 and seq2:
    # print(seq1, phread1, QNAME1)
    # print(seq2, phread2, QNAME2)

    ######################################
    #In this part we invert seq2
    ######################################

    #Generate kmers for both reads
    kmers1 = f.generate_kmers(seq1, k = k, shift = shift)
    kmers2 = f.generate_kmers(seq2, k = k, shift = shift, reverse_complement = True)

    #Map kmers to the genome (exact with FM index)
    kmers1 = G.map_kmers(kmers1)
    kmers2 = G.map_kmers(kmers2)

    # print(kmers1)
    # print(kmers2)

    #Find the starting positions of each 
    start1 = f.get_start_index(kmers1)
    start2 = f.get_start_index(kmers2)

    # print(start1)
    # print(start2)

    comb = f.get_feasible_comb(start1, start2, max_dist)
    #print(I, 'Length comb:' ,len(comb))

    if len(comb) > 0:
        #Find unique start sites such that a read is not aligned to the same position multiple times
        comb1_unique = list(set([i[0] for i in comb]))
        comb2_unique = list(set([i[1] for i in comb]))
        #print(I, 'Length unique:', len(comb1_unique), len(comb2_unique))

        for i in comb1_unique:
            genome_string = G.genome[start1[i]:(start1[i]+len(seq1))]
            score = G.align(seq1, genome_string, cost, constrain = constrain)
            start1[i] = (start1[i], score)

        for i in comb2_unique:
            genome_string = G.genome[start2[i]:(start2[i]+len(seq2))]
            score = G.align(f.complement(seq2[::-1]), genome_string, cost, constrain = constrain)
            start2[i] = (start2[i], score)

        # print(start1)
        # print(start2)

        #Calculate alignment score for all start sites in start1_unique and start2_unique
        alignments = []
        for i in comb:
            alignments.append(start1[i[0]][1] + start2[i[1]][1]) #add the scores of alignment of read1 and read2 together to obtain overall score of this read mapping

        # print(alignments)

    ######################################
    #In this part we invert seq1
    ######################################

    #Generate kmers for both reads
    kmers1 = f.generate_kmers(seq1, k = k, shift = shift, reverse_complement = True)
    kmers2 = f.generate_kmers(seq2, k = k, shift = shift)

    #Map kmers to the genome (exact with FM index)
    kmers1 = G.map_kmers(kmers1)
    kmers2 = G.map_kmers(kmers2)

    #Find the starting positions of each 
    start1_ = f.get_start_index(kmers1)
    start2_ = f.get_start_index(kmers2)

    comb_ = f.get_feasible_comb(start1_, start2_, max_dist)

    #print(I, 'Length comb:' ,len(comb_))
    if len(comb_) > 0:
        #Find unique start sites such that a read is not aligned to the same position multiple times
        comb1_unique = list(set([i[0] for i in comb_]))
        comb2_unique = list(set([i[1] for i in comb_]))
        #print(I, 'Length unique:', len(comb1_unique), len(comb2_unique))

        for i in comb1_unique:
            genome_string = G.genome[start1_[i]:(start1_[i]+len(seq1))]
            score = G.align(f.complement(seq1[::-1]), genome_string, cost, constrain = constrain)
            start1_[i] = (start1_[i], score)

        for i in comb2_unique:
            genome_string = G.genome[start2_[i]:(start2_[i]+len(seq2))]
            score = G.align(seq2, genome_string, cost, constrain = constrain)
            start2_[i] = (start2_[i], score)

        #Calculate alignment score for all start sites in start1_unique and start2_unique
        alignments_ = []
        for i in comb_:
            alignments_.append(start1_[i[0]][1] + start2_[i[1]][1]) #add the scores of alignment of read1 and read2 together to obtain overall score of this read mapping


    #Decide which inversion option is best
    #Inversion of read1 did not yield any feasible starting positions or resulted in vorse alignments
    if (len(comb) and len(comb_) == 0) or (len(comb) and len(comb_) and min(alignments) < min(alignments_)):
        #Get best alignment
        min_idx = alignments.index(min(alignments))
        POS1 = start1[comb[min_idx][0]][0]+1
        POS2 = start2[comb[min_idx][1]][0]+1
        append_sam_alignment(filename, reverse=False, secondary=False, QNAME=QNAME1, RNAME=ref_name, POS=POS1, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*")
        append_sam_alignment(filename, reverse=True, secondary=False, QNAME=QNAME2, RNAME=ref_name, POS=POS2, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*")
    #Inversion of read2 did not yield any feasible starting positions or resulted in worse alignments
    elif (len(comb) == 0 and len(comb_)) or (len(comb) and len(comb_) and min(alignments) > min(alignments_)):
        min_idx = alignments_.index(min(alignments_))
        POS1 = start1_[comb_[min_idx][0]][0]+1
        POS2 = start2_[comb_[min_idx][1]][0]+1
        append_sam_alignment(filename, reverse=True, secondary=False, QNAME=QNAME1, RNAME=ref_name, POS=POS1, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*")
        append_sam_alignment(filename, reverse=False, secondary=False, QNAME=QNAME2, RNAME=ref_name, POS=POS2, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*")
    #Not able to map
    elif (len(comb) == 0 and len(comb_) == 0):
        POS1 = 0
        POS2 = 0
        append_sam_alignment(filename, reverse=False, secondary=False, QNAME=QNAME1, RNAME=ref_name, POS=POS1, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*")
        append_sam_alignment(filename, reverse=False, secondary=False, QNAME=QNAME2, RNAME=ref_name, POS=POS2, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*")
    #TODO what do we do if min(alignment_) == min(alignments_)?

    seq1, phread1, QNAME1 = read_line(f1)
    seq2, phread2, QNAME2 = read_line(f2)
    

end = time.time()
print('Processing finished (duration:)', end-start)
f1.close()
f2.close()

# %%
