#%%
class Reads:
    pair1 = []
    pair2 = []

    def __init__(self, file1=None, file2 = None):
            self.pair1 = []
            self.pair2 = []
            self.read_file(file1)
            self.read_file(file2)

    def read_file(self, file):
        with open(file, 'r') as f:
            line = f.readline()
            while line:
                if(line.startswith('@')):
                    #Process id line
                    line = line.strip('@').strip('\n')
                    idx, origin = line.split('/')
                    origin = int(origin)

                    #Process seq and skip next line which is a +
                    seq = f.readline().strip('\n')
                    f.readline()

                    #Process phread
                    phread = f.readline().strip('\n')

                    if origin == 1:
                        self.pair1.append(Read(idx, seq, phread, origin))
                    else:
                        self.pair2.append(Read(idx, seq, phread, origin))
                    
                    #Read next line
                    line = f.readline()
                else:
                    line = f.readline()

    def print_reads(self):
        for i in self.pair1:
            print('Pair1 reads:')
            print(i.get_seq, '\n', i.get_phread, sep = '')

        for i in self.pair2:
            print('Pair2 reads:')
            print(i.get_seq, '\n', i.get_phread, sep = '')

class Read:
    def __init__(self, id, seq, phread, origin, k = 0):
        self.__id = id
        self.length = len(seq)
        self.__seq = seq
        self.__phread = phread
        self.__origin = origin
        self.__alphabet = ['A','C','G','T']
        self.__inverse_alphabet = ['T','G','C','A']
        self.mapped = []

    def complement(self, string):
        idx = [self.__inverse_alphabet[self.__alphabet.index(i)] for i in string]
        return ''.join(idx)

    def generate_kmers(self, k, shift = 1, reverse_complement = False):
        #generates a tuple of ([list1], [list2]) where the list1 contains
        #all the position of this kmer in the read and list2 is empty
        #list2 will ultimatively filled by the function map_kmers() and will
        #contain all the positions in the genome where this kmer is mapped to
        if reverse_complement:
            s = self.complement(self.__seq[::-1])
        else:
            s = self.__seq

        self.kmer = {}
        for i in range(0,len(s)-k+1, shift):
            if s[i:i+k] in self.kmer:
                self.kmer[s[i:i+k]][0].append(i)
            else:
                self.kmer[s[i:i+k]] = ([i],[])

    def get_start_index(self):
        self.idx = {}
        for kmer in self.kmer.keys():
            #All combinations of starting positions
            comb = list(itertools.product(self.kmer[kmer][0], self.kmer[kmer][1]))
            for i in comb:
                start = i[1]-i[0]
                if str(start) in self.idx:
                    self.idx[str(start)] += 1
                else:
                    self.idx[str(start)] = 1
    
        #Generate sorted tulple list out of dict
        self.idx = sorted(self.idx.items(), key=lambda kv: kv[1], reverse=True)
        self.idx = [(int(i[0]), i[1]) for i in self.idx]

    def get_id(self):
        return self.__id

    def get_seq(self):
        return self.__seq

    def get_phread(self):
        return self.__phread

    def get_origin(self):
        return self.__origin

class Genome:

    def __init__(self, file):
        self.read_genome(file)
        self.__size = len(self.__genome)
        self.__alphabet = ['$','A','C','G','T']
        self.__inverse_alphabet = ['$','T','G','C','A']
        self.bwt_index()
        self.get_C()
        self.get_Occ()

    def read_genome(self, file):
        with open(file, 'r') as f:
            line = f.readline()
            while line:
                if(line.startswith('>')):
                    self.__genome = f.read().replace('\n', '')
                    line = f.readline()
                else:
                    line = f.readline()

    def get_genome(self):
        return(self.__genome)

    def get_alphabet(self):
        return(self.__alphabet)

    def rot_left(self,s):
        return(s[1:] + s[0])

    def bwt_index(self):
        s = self.__genome + '$'
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
        for c in self.__alphabet:
            occ=0
            for l in self.L:
                if l<c:
                    occ+=1
            self.C[c] = occ
    
    def get_Occ(self):
        self.Occ = np.full((len(self.C), len(self.L)),-1)
        for i, c in enumerate(self.__alphabet):
            occ = 0
            for j, l in enumerate(self.L):
                if l == c:
                    occ+=1
                self.Occ[(i, j)] = occ
    
    def get_FM(self, P):
        #check if pattern is feasible
        #TODO delete in final version
        for i in range(len(P)):
            if P[i] not in self.__alphabet:
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

    def map_kmers(self, read):
        del_kmer = []
        for kmer in read.kmer.keys():
            sp, ep = self.get_FM(kmer)
            if(sp > ep):
                del_kmer.append(kmer)
            else:
                read.kmer[kmer][1].extend(self.get_genomeIdx(sp,ep))
        
        for kmer in del_kmer:
            del read.kmer[kmer]

    def align(self, substring, string, cost):
        alphabet = self.__alphabet[1:]
        
        C = np.full((len(substring)+1, len(string)+1), -1)
        C[0,:] = 0
        C[:, 0] = np.arange(0, len(substring)+1)

        for row in range(1,C.shape[0]):
            for col in range(1,C.shape[1]):
                idx1 = alphabet.index(substring[row-1])
                idx2 = alphabet.index(string[col-1])

                match = C[row-1,col-1] + cost[idx1, idx2]
                mismatch = C[row-1,col-1] + cost[idx1, idx2]
                gap1 = C[row-1, col] + 2
                gap2 = C[row, col-1] + 2

                C[row,col] = min(match, mismatch, gap1, gap2)
        return min(C[len(substring),:])

    def extract_genome_region(self, genome_idx, read_idx, read_length, margin):
        idx1 = genome_idx - read_idx - margin
        idx2 = genome_idx + read_length - read_idx + margin
        print(idx1, idx2, self.__genome)
        #return string and start position of genome
        return self.__genome[idx1:idx2], genome_idx - read_idx

    def map_reads(self, R, k, max_distance, margin, cost):

        #For every pair read1 and read2 compute the optimal alignment
        #and inversion status.
        for I in range(len(R.pair1)):
            #Get current reads
            read1 = R.pair1[I]
            read1.generate_kmers(k)
            self.map_kmers(read1)            
            read1.get_start_index()
            idx1 = read1.idx
            #print(read1.kmer)
            #print(read1.idx)

            #Invert read 2
            read2 = R.pair2[I]
            read2.generate_kmers(k, reverse_complement = True) 
            self.map_kmers(read2)
            read2.get_start_index()
            idx2_inv = read2.idx
            #print(read2.kmer)
            #print(read2.idx)

            comb1 = []
            #look at all combinations within a range
            for n in range(len(read1.idx)):
                for m in range(len(read2.idx)):
                    if (abs(read1.idx[n][0] - read2.idx[m][0]) < max_distance):
                        comb1.append((n,m))

            #print('Comb1: ',comb1)
            alignment1 = []
            if len(comb1) > 0:
                for c in comb1:
                    #Map first read1
                    #print(read1.length)
                    genome_string1 = self.__genome[read1.idx[c[0]][0]:(read1.idx[c[0]][0]+read1.length)]
                    score1 = self.align(read1.get_seq(), genome_string1, cost)
                    #print(genome_string1, read1.get_seq(), score1)

                    #Map second read2
                    genome_string2 = self.__genome[read2.idx[c[1]][0]:(read2.idx[c[1]][0]+read2.length)]
                    s = read2.complement(read2.get_seq()[::-1])
                    score2 = self.align(s, genome_string2, cost)
                    #print(genome_string2, read2.get_seq(), score2)

                    alignment1.append(score1+score2)
            #if len(alignment1) > 0:
            #    print('Alignment: ',alignment1, alignment1.index(min(alignment1)))

            #Get current reads
            #read1 = R.pair1[i]
            read1.generate_kmers(k, reverse_complement = True) 
            self.map_kmers(read1)
            read1.get_start_index()
            idx1_inv = read1.idx
            #print(read1.kmer)
            #print(read1.idx)

            #Invert read 2
            #read2 = R.pair2[i]
            read2.generate_kmers(k, reverse_complement = False) 
            self.map_kmers(read2)
            read2.get_start_index()
            idx2 = read2.idx
            #print(read2.kmer)
            #print(read2.idx)
            
            comb2 = []
            #look at all combinations within a range
            for n in range(len(read1.idx)):
                for m in range(len(read2.idx)):
                    if (abs(read1.idx[n][0] - read2.idx[m][0]) < max_distance):
                        comb2.append((n,m))
            #print('Comb2: ',comb2)
            print(len(comb1), len(comb2))
            alignment2 = []
            if len(comb2) > 0:
                for c in comb2:
                    #Map first read1
                    #print(read1.length)
                    genome_string1 = self.__genome[read1.idx[c[0]][0]:(read1.idx[c[0]][0]+read1.length)]
                    s = read1.complement(read1.get_seq()[::-1])
                    score1 = self.align(s, genome_string1, cost)
                    #print(genome_string1, read1.get_seq(), score1)

                    #Map second read2
                    genome_string2 = self.__genome[read2.idx[c[1]][0]:(read2.idx[c[1]][0]+read2.length)]
                    score2 = self.align(read2.get_seq(), genome_string2, cost)
                    #print(genome_string2, read2.get_seq(), score2)

                    alignment2.append(score1+score2)
            #if len(alignment2) > 0:
            #    print('Alignment2: ',alignment2)
            #    print(alignment2.index(min(alignment2)))

            #Read1 & read2 inverted 
            if((len(alignment1) and len(comb2) == 0) or (len(alignment1) and len(comb2) and min(alignment1) < min(alignment2))):
                s = sorted(range(len(alignment1)), key=lambda k: alignment1[k])
                comb1_sorted = [comb1[i] for i in s]
                #print(idx1)
                #print(idx2_inv)
                for idx, i in enumerate(comb1_sorted):
                    read1.mapped.append((idx1[i[0]][0]+1, alignment1[s[idx]], 0)) #genome index, score, inversion falag
                    read2.mapped.append((idx2_inv[i[1]][0]+1, alignment1[s[idx]], 1)) #genome index, score, inversion falag
            
            elif((len(alignment1) == 0 and len(comb2)) or (len(alignment1) and len(comb2) and min(alignment1) > min(alignment2))):
                s = sorted(range(len(alignment2)), key=lambda k: alignment2[k])
                comb2_sorted = [comb2[i] for i in s]
                for i in comb2_sorted:
                    for idx, i in enumerate(comb2_sorted):
                        read1.mapped.append((idx1_inv[i[0]][0]+1, alignment2[s[idx]], 1)) #genome index, score, inversion falag
                        read2.mapped.append((idx2[i[1]][0]+1, alignment2[s[idx]], 0)) #genome index, score, inversion falagndex, score, inversion falag
            else:
                print(I,'Not mapped')
                read1.mapped = 'Not mapped'
                read2.mapped = 'NOt mapped'
            print(I, '/', len(R.pair1))

#%%
import numpy as np
import itertools
import collections
#%%

# G = Genome('genome.txt')
# R = Reads('read.txt', 'read2.txt')

G = Genome('./data_small/genome.chr22.5K.fa')
R = Reads('./data_small/output_tiny_30xCov1.fq', './data_small/output_tiny_30xCov2.fq')

cost = np.full((4,4),1)
cost[np.diag_indices(len(cost))] = 0

G.map_reads(R, k = 16, max_distance = 500, margin = 0, cost = cost)

#%%
P = 'TCT'
sp, ep = G.get_FM(P)
if(sp != -1 and ep != -1):
    idx = G.get_genomeIdx(sp, ep)
    print(idx)
    print([G._Genome__genome[i:i+len(P)] for i in idx])
    print(G.get_genome())

#%%
for i in R.pair1:
    print(i.get_seq())
    i.generate_kmers(k = 2)
    print(i.kmer)
    G.map_kmers(i)
    print(i.kmer)

#%%
for i in R.pair1:
    for kmer in i.kmer.keys():
        print(kmer)
        for idx in i.kmer[kmer][1]:
            print(G._Genome__genome[idx:idx+len(kmer)], end = ' ')
        print()

#%%
print(R.pair1[0].kmer)

# %%
s = G.extract_genome_region(R.pair1[0].kmer['GA'][1], R.pair1[0].kmer['GA'][0], len(R.pair1[0].get_seq()), 0)

# %%
G.align(s, R.pair1[0].get_seq(), cost)

# %%
R = Reads('read.txt', 'read2.txt')

#%%
G = Genome('genome.txt')



# %%
comb = itertools.product(R.pair1[0].kmer['GT'][0], R.pair1[0].kmer['GT'][1])

# %%
