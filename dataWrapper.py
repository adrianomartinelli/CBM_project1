#%%
class Reads:
    pair1 = []
    pair2 = []

    def __init__(self, file=None):
            self.pair1 = []
            self.pair2 = []
            self.read_file(file)
            # if file:
            #     self.read_file(file+'1.fq')
            #     self.read_file(file+'2.fq')

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
        self.__seq = seq
        self.__phread = phread
        self.__origin = origin
        if k != 0:
            self.generate_kmers(k)

    def generate_kmers(self, k):
        #generates a tuple of ([list1], [list2]) where the list1 contains
        #all the position of this kmer in the read and list2 is empty
        #list2 will ultimatively filled by the function map_kmers() and will
        #contain all the positions in the genome where this kmer is mapped to
        self.kmer = {}
        for i in range(len(self.__seq)-k+1):
            if self.__seq[i:i+k] in self.kmer:
                self.kmer[self.__seq[i:i+k]][0].append(i)
            else:
                self.kmer[self.__seq[i:i+k]] = ([i],[])

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
        self.__alphabet = list(set(G.L))
        idx = np.argsort(self.__alphabet)
        self.__alphabet = [self.__alphabet[i] for i in idx]

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
        for kmer in read.kmer.keys():
            sp, ep = self.get_FM(kmer)
            read.kmer[kmer][1].extend(self.get_genomeIdx(sp,ep))
        

#%%
import numpy as np

G = Genome('genome.txt')
R = Reads('read.txt')

#%%
G.bwt_index()
G.get_C()
G.get_Occ()
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
