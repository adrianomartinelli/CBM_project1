import numpy as np
import matplotlib.pyplot as plt
import os
import random

alphabet     = ['A','C','G','T']
inv_alphabet = ['T','G','C','A']


class Read():
  def __init__(self, idx, origin, string, quality, num_seeds=0, len_seed=0, inv=False):
    self.idx = idx
    self.origin = origin
    self.string = string
    self.quality = quality
    self.inv = inv
    self.num_seeds = num_seeds
    self.len_seed = len_seed
    self.num_match = 0
    self.num_inv_match = 0
    self.seed_pos = random.sample(range(len(string)-len_seed), num_seeds)
    self.subs_seed_pos = [x for x in range(len(string)-len_seed) if x not in self.seed_pos]
    self.generate_seeds()
    if inv:
      self.inv_complement()

  def generate_seeds(self):
    self.seeds = {}
    for i in self.seed_pos:
      while self.string[i:i+self.len_seed] in self.seeds:
        i = random.sample(self.subs_seed_pos, 1)[0]
      self.seeds[self.string[i:i+self.len_seed]] = {'ReadPos':[i],'GenomePos':[],'NumMatch':0}

  def generate_inv_seeds(self):
    self.inv_seeds = {}
    for i in self.seed_pos:
      while self.inv_string[i:i+self.len_seed] in self.inv_seeds:
        i = random.sample(self.subs_seed_pos, 1)[0]
      self.inv_seeds[self.inv_string[i:i+self.len_seed]] = {'ReadPos':[i],'GenomePos':[], 'NumMatch':0}

  def complement(self, inv):
    inv_read = ''
    for base in inv:
        inv_base = inv_alphabet[alphabet.index(base)]
        inv_read += inv_base
    return inv_read

  def inv_complement(self):
    #generates inverse complement of the read
    self.inv = not self.inv
    inv = self.string[::-1]
    self.inv_string = self.complement(inv)
    self.inv_quality = self.quality[::-1]
    self.generate_inv_seeds()

class Reads():
  def __init__(self, file=None, num_seeds=0, len_seed=0):
    self.pair1 = []
    self.pair2 = []
    self.num_seeds = num_seeds
    self.len_seed = len_seed
    if file:
      self.get_read_string(file+'1.fq')
      self.get_read_string(file+'2.fq')

  def get_read_string(self, file):
    with open(file, 'r') as f:
      line = f.readline()
      while line:
        if(line.startswith('@')):
          #Process ID line
          line = line.lstrip('@').rstrip('\n')
          #Get id of the genome and pair origin
          idx, origin = line.split('/')
          #Get read string
          string = f.readline().rstrip('\n')
          #Skip line with a '+'
          f.readline()
          #Get quality of the read pair
          quality = f.readline().rstrip('\n')

          #If pair origin is 1, assign the read to pair1, else assign the read to pair2
          if(origin == '1'):
              self.pair1.append(Read(idx, origin, string, quality, self.num_seeds, self.len_seed))
          else:
              self.pair2.append(Read(idx, origin, string, quality, self.num_seeds, self.len_seed))

          #Get next ID line
          line = f.readline()
        else:
          line = f.readline()

  def print_reads(self):
    print('Pair1 reads:')
    for r in self.pair1:
      print(r.idx)
      print(r.string, '+', r.quality, sep = '\n')

    print('Pair2 reads:')
    for r in self.pair2:
      print(r.idx)
      print(r.string, '+', r.quality, sep = '\n')


class Genome():
  def __init__(self, file=None):
    self.genome_string = ''
    self.genome_length = 0
    self.alphabet = ['$','A','C','G','T']
    if file:
      self.read_genome(file)
    self.bwt_index()
    self.occurence_array()
    self.occurence_matrix()

  def read_genome(self, file):
    with open(file, 'r') as f:
      self.genome_header = f.readline()
      for line in f:
        line = line.rstrip('\n')
        self.genome_string = self.genome_string + line
    self.genome_length = len(self.genome_string)

  def bwt_index(self):
    s = self.genome_string + '$'
    m = [s]
    for _ in range(len(s)-1):
      s = s[1:] + s[0]
      m.append(s)

    idx = np.argsort(m)

    self.s = [m[i] for i in idx]
    self.i = list(idx)
    self.L = [i[-1] for i in self.s]
    self.F = [i[0]  for i in self.s]

  def occurence_array(self):
    self.C={}
    for c in self.alphabet:
      occ=0
      for l in self.L:
        if l<c:
          occ+=1
      self.C[c] = occ

  def occurence_matrix(self):
    self.Occ={}
    for c in self.alphabet:
      occ=0
      for idx, l in enumerate(self.L):
        if l==c:
          occ+=1
        self.Occ[(c, idx)] = occ

  def next_char(self, char):
    idx = self.alphabet.index(char)
    if idx == (len(self.alphabet) - 1):
      return False
    else:
      return self.alphabet[idx+1]

  def find_pos(self, P):
    #Calculate initial values for FM index algorithm
    i = len(P)
    c = P[i-1]
    sp = self.C[c] + 1

    #Calculate 'ep' by finding the next character of 'c'
    if self.next_char(c):
      ep = self.C[self.next_char(c)]
    else:
      ep = len(self.L)

    while sp <= ep and i >= 2:
      c = P[i-2]
      sp = self.C[c] + self.Occ[(c, sp - 2)] + 1
      ep = self.C[c] + self.Occ[(c, ep - 1)]
      i=i-1

    if sp > ep:
      return []
    else:
      match_idx = [self.i[i] for i in range(sp-1,ep)]
      match_idx.sort()
      return match_idx

  def seed_positions(self, reads):
    for pair in [reads.pair1, reads.pair2]:
      for read in pair:
        for seed in read.seeds:
          pos = self.find_pos(seed)
          if pos:
            read.seeds[seed]['GenomePos'] = pos
            read.seeds[seed]['NumMatch'] = len(pos)
            read.num_match += len(pos)
        if read.inv:
          for inv_seed in read.inv_seeds:
            inv_pos = self.find_pos(inv_seed)
            if inv_pos:
              read.inv_seeds[inv_seed]['GenomePos'] = inv_pos
              read.inv_seeds[inv_seed]['NumMatch'] = len(inv_pos)
              read.num_inv_match += len(inv_pos)

  def edit_distance(self, x, y):
    score = [[0, 4, 2, 4, 8],
             [4, 0, 4, 2, 8],
             [2, 4, 0, 4, 8],
             [4, 2, 4, 0, 8],
             [8, 8, 8, 8, 8]]
    D = []
    for i in range(len(x) + 1):
      D.append([0] * (len(y)+1))

    for i in range(len(x) + 1):
      D[i][0] = i + score[alphabet.index(x[i-1])][-1]
    for i in range(len(y) + 1):
      D[0][i] = i + score[-1][alphabet.index(y[i-1])]

    for i in range (1, len(x)+1):
      for j in range (1, len(y)+1):
        distHor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
        distVer = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
        if x[i-1] == y[j-1]:
          distDiag = D[i-1][j-1]
        else:
          distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

        D[i][j] = min(DistHor, distVer, distDiag)

    return D[-1][-1]


##Read an example genome
G = Genome('ExampleGenome.txt')
#Print the example genome string
print(G.genome_string)

##Pattern matching example
P = 'AGT'
P_match_idx = G.find_pos(P)
if P_match_idx:
  print("Matching locations in the genome:", np.array(P_match_idx)+1, sep='\n')
  print("Patterns at matching locations:",[G.genome_string[i:i+len(P)] for i in P_match_idx],  sep='\n')
else:
  print('No match :(')

##Load the example reads and generate seeds
example_read = 'ExampleRead'
R = Reads(example_read, num_seeds=5, len_seed=3)

print("Reads and generated seeds with their positions in the reads:\n")
print(R.pair1[0].string)
print(R.pair1[0].seeds)
#Calculate inverse complement of the read and generate seeds
R.pair1[0].inv_complement()
print(R.pair1[0].inv_string)
print(R.pair1[0].inv_seeds)

##Find the locations of read seeds in genome
G.seed_positions(R)

print("Reads and generated seeds with their positions in the read and genome as well as number of matches:\n")
print(R.pair1[0].string)
print(R.pair1[0].seeds)
print(R.pair1[0].inv_string)
print(R.pair1[0].inv_seeds)