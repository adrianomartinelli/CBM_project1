#%%
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import itertools
import math
import argparse
import time

alphabet     = ['A','C','G','T']
inv_alphabet = ['T','G','C','A']

def create_sam_header(filename, RNAME, RLEN, version=1.4, sorting_order="unsorted"):
    #Creating a new SAM file
    sam_file = open(filename, 'w')

    #Adding headers
    sam_file.write("@HD\tVN:%.1f\tSO:%s\n" % (version, sorting_order))
    sam_file.write("@SQ\tSN:%s\tLN:%i\n" % (RNAME, RLEN))
    sam_file.close()

def append_sam_alignment(filename, reverse=False, secondary=False, QNAME="*", RNAME="*", POS=0, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*"):
    #Calculating the flag
    FLAG = 0
    if reverse: FLAG += 16
    if secondary: FLAG += 256

    #Open existing SAM file append one linear alignment as a line to the alignment body
    sam_file = open(filename, "a")
    sam_file.write("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL))
    sam_file.close()


#Genome class
class Genome():
  def __init__(self, k, file=None):
    """
    Create Genome class with given file
    """
    self.genome_string = ''
    self.genome_length = 0
    self.k = k
    self.alphabet = ['$','A','C','G','T']
    self.cost_matrix = [[0, 4, 2, 4, 8],
                        [4, 0, 4, 2, 8],
                        [2, 4, 0, 4, 8],
                        [4, 2, 4, 0, 8],
                        [8, 8, 8, 8, 8]]
    if file:
      print('Read in Genome file')
      self.read_genome(file)
      print('Read in Genome file: Done')

    print('Start indexing Genome with k=', self.k)
    self.index()
    print('Finished indexing Genome')

  def read_genome(self, file):
    """
    Read genome file
    """
    with open(file, 'r') as f:
      line = f.readline()
      while line:
        if(line.startswith('>')):
          print('\tRead line')
          self.genome_name = line.rstrip('\n')[1:]
          self.genome_string = f.read().replace('N', '').replace('\n', '')
          line = f.readline()
        else:
          line = f.readline()
      self.genome_length = len(self.genome_string)
      #print('Genome length', self.genome_length)

    # with open(file, 'r') as f:
    #   self.genome_name = f.readline().rstrip('\n')[1:]
    #   for line in f:
    #     line = line.rstrip('\n')
    #     if 'N' not in line:
    #       self.genome_string = self.genome_string + line
    # self.genome_length = len(self.genome_string)

  def occurence_array(self, L):
    """
    Calculate occurence array
    """
    self.C={}
    for c in self.alphabet:
      occ=0
      for l in L:
        if l<c:
          occ+=1  
      self.C[c] = occ

  def occurence_matrix(self, L):
    """
    Calculate occurence matrix
    """
    self.Occ={}
    for c in self.alphabet:
      occ=0
      for idx, l in enumerate(L):
        if l==c:
          occ+=1
        self.Occ[(c, idx)] = occ
  
  def index(self):
    """
    Calculate necessary vectors for FM-Index algorithm
    """

    if self.genome_length < self.k:
      k = self.genome_length+1
    else:
      k = self.k

    s = self.genome_string + '$'
    ss = s + s[0:k]
    #print(s)
    m = []
    l = []
    #print(self.genome_length+1)
    for i in range(self.genome_length+1):
      if not i % 100000:
        print('\t',i, '/', self.genome_length, '(',i/self.genome_length*100, '%)')
      m.append(ss[i:(i+k)])
      l.append(s[-1 + i])

    # self.m = m
    # self.l = l

    # s = self.genome_string + '$'
    # mm = [s[0:k]]
    # ll = [s[-1]]
    # cur_rot = s
    # print(len(s)-1)
    # for i in range(len(s)-1):
    #   if not (i % 10000):
    #     print('Rotation', i,'/',len(s))
    #   cur_rot = cur_rot[1:] + cur_rot[0]
    #   #s = s[1:] + s[0]
    #   mm.append(cur_rot[0:k])
    #   ll.append(cur_rot[-1])
    
    # print(mm == m)
    # print(ll == l)

    # # self.M = M
    # self.ll = ll
    #print(any(np.where(M != m)))
    #print(any(np.where(ll != l)))

    print('\tSort index')
    idx = np.argsort(m)

    self.i = list(idx)
    print('\tConstruct s')
    s = [m[i] for i in idx]
    #print(s)
    #L = [i[-1] for i in s]
    print('\tConstruct L')
    L = [l[i] for i in idx]
    #print(L)

    print('\tConstruct Occ Array')
    self.occurence_array(L)
    print('\tConstruct Occ Matrix')
    self.occurence_matrix(L)

  def next_char(self, char):
    """
    Calculate next char in the alphabet
    """
    idx = self.alphabet.index(char)
    if idx == (len(self.alphabet) - 1):
      return False
    else:
      return self.alphabet[idx+1]

  def find_pos(self, P):
    """
    Calculate the position of pattern 'P' in genome using FM-index algorithm
    """
    #Calculate initial values for FM index algorithm
    i = len(P)
    c = P[i-1]
    sp = self.C[c] + 1

    #Calculate 'ep' by finding the next character of 'c'
    if self.next_char(c):
      ep = self.C[self.next_char(c)]
    else:
      ep = self.genome_length + 1

    #FM-Index algorithm
    while sp <= ep and i >= 2:
      c = P[i-2]
      sp = self.C[c] + self.Occ[(c, sp - 2)] + 1
      ep = self.C[c] + self.Occ[(c, ep - 1)]
      i=i-1

    #If sp>ep return an empyt list, otherwise return the corresponding positions in the genome
    if sp > ep:
      return []
    else:
      match_idx = [self.i[i] for i in range(sp-1,ep)]
      match_idx.sort()
      return match_idx
    
  def seed_positions(self, read1, read2):
    """
    Calculate the genome positions of seeds of given reads
    """
    for read in [read1, read2]:
      inv_flag = False
      for seed in list(read.seeds.keys()):
        pos = self.find_pos(seed)
        if pos:
          read.seeds[seed]['GenomePos'] = pos
        else:
          del read.seeds[seed]
      
      if len(read.seeds) < read.num_seeds/5:
        return not inv_flag
    
    return inv_flag
  
  def edit_distance(self, read, start, constraint=12, margin=4):
    """
    Calculate the edit distance matrix for alignment and keep track of the best path
    """
    DELETION, INSERTION, MISMATCH, MATCH = range(4)
    cost = self.cost_matrix
    x, y = read.string, self.genome_string[max(start-margin, 0): min(start+read.len_string+margin, self.genome_length)]
    
    D = []
    for i in range(len(x) + 1):
      D.append([math.inf] * (len(y)+1))
    
    for i in range(len(y) + 1):
      D[0][i] = 0
    
    B = []
    for i in range(len(x) + 1):
      B.append([0] * (len(y)+1))
    
    for i in range(1, len(x) + 1):
      D[i][0] = D[i-1][0] + cost[alphabet.index(x[i-1])][-1]
      B[i][0] = DELETION

    for i in range (1, len(x)+1):
      for j in range(max(i-constraint, 1), min(len(y)+1,i+constraint+1)):
        distHor = (D[i][j-1] + cost[-1][alphabet.index(y[j-1])], INSERTION)
        distVer = (D[i-1][j] + cost[alphabet.index(x[i-1])][-1], DELETION)
        
        if x[i-1] == y[j-1]:
          distDiag = (D[i-1][j-1], MATCH)
        else:
          distDiag = (D[i-1][j-1] + cost[alphabet.index(x[i-1])][alphabet.index(y[j-1])], MISMATCH)
        
        D[i][j], B[i][j] = min(distDiag, distHor, distVer)
    
    best_pos = np.argmin(D[-1][:])
    cost = min(D[-1][:])
    
    alignment = []
    i, j = len(x), best_pos
    
    #Backtrack from the best position
    if B[i][j] == DELETION:
      i -= 1
      alignment.append('D')
    elif B[i][j] == INSERTION:
      j -= 1
      alignment.append('I')
    elif B[i][j] == MISMATCH:
      i -= 1
      j -= 1
      alignment.append('X')
    elif B[i][j] == MATCH:
      i -= 1
      j -= 1
      alignment.append('M')
    
    alignment.append(1)
      
    while i > 0:
      assert i >= 0
      if B[i][j] == DELETION:
        i -= 1
        if alignment[-2] == 'D':
          alignment[-1] += 1
        else:
          alignment.extend(['D', 1])
      elif B[i][j] == INSERTION:
        j -= 1
        if alignment[-2] == 'I':
          alignment[-1] += 1
        else:
          alignment.extend(['I', 1])
      elif B[i][j] == MISMATCH:
        i -= 1
        j -= 1
        if alignment[-2] == 'X':
          alignment[-1] += 1
        else:
          alignment.extend(['X', 1])
      elif B[i][j] == MATCH:
        i -= 1
        j -= 1
        if alignment[-2] == 'M':
          alignment[-1] += 1
        else:
          alignment.extend(['M', 1])
      else:
        assert(False)
    
    genome_pos = j + start - margin
    seq_alignment = ''.join(str(i) for i in list(reversed(alignment)))

    #Fill in the necessary fields in read class
    read.genome_pos[start]['Pos'] = genome_pos
    read.genome_pos[start]['Cost'] = cost
    read.genome_pos[start]['Alignment'] = seq_alignment
    
    return cost
    

class Read():
  def __init__(self, idx, origin, string, quality, len_seed=0, inv=False):
    """
    Construct read class with necessary variables
    """
    self.idx = idx
    self.origin = origin
    self.string = string
    self.quality = quality
    self.len_string = len(string)
    self.len_seed = len_seed
    self.num_seeds =  len(range(0, len(self.string) - self.len_seed, int(self.len_seed/2)))
    if inv:
      self.inv = not inv
      self.inverse()
    else:
      self.inv = inv
      self.generate_seeds()

  def generate_seeds(self):
    """
    Generate seeds for a given read
    """
    self.seeds = {}
    for i in range(0, len(self.string) - self.len_seed, int(self.len_seed/2)):
      kmer = self.string[i:i+self.len_seed]
      if kmer in self.seeds:
         self.seeds[kmer]['ReadPos'].append(i)
      else:
         self.seeds[kmer] = {'ReadPos':[i], 'GenomePos':[], 'AlignmentPos':[], 'Alignment':[]}
  
  def complement(self, inv):
    """
    Change the each base of read with its complement
    """
    inv_read = ''
    for base in inv:
        inv_base = inv_alphabet[alphabet.index(base)]
        inv_read += inv_base
    return inv_read

  def inverse(self):
    """
    Generates inverse complement of the read and its seeds
    """
    self.inv = not self.inv
    self.string = self.complement(self.string[::-1])
    self.quality = self.quality[::-1]
    self.generate_seeds()
  
  def start_index(self, G):
    """
    Calculate potential positions of reads in the genome
    """
    idx = {}
    for seed in self.seeds:
      combinations = list(itertools.product(self.seeds[seed]['ReadPos'], self.seeds[seed]['GenomePos']))
      for pair in combinations:
        #All combinations of starting positions
        start = pair[1] - pair[0]
        if start>=0 and start<=(G.genome_length - self.len_string):
          if str(start) in idx:
            idx[str(start)] += 1
          else:
            idx[str(start)] = 1
          
    #Find the max number of matches among the matching positions
    max_match = max(idx.values())
    #Generate sorted tuple list out of dict
    idx = sorted(idx.items(), key=lambda item: item[1], reverse=True)
    start_idx = [int(i[0]) for i in idx if i[1]>max_match/5]
    self.genome_pos = {i:{'Pos:':int, 'Cost':int, 'Alignment':''} for i in start_idx}
  

def read_line(f):
  line = f.readline()
  if(line.startswith('@')):
    #Process ID line
    line = line.lstrip('@').rstrip('\n')
        
    #Get id of the genome and pair origin
    idx, origin = line.split('/')
        
    #Get read string
    string = f.readline().rstrip('\n')
        
    #Skip line with '+'
    f.readline()
          
    #Get quality of the read pair
    quality = f.readline().rstrip('\n')
    
    return idx, origin, string, quality
  else:
    return '', '', '', ''

def feasible_combinations(r1, r2, dist_threshold=500):
  """
  Calculate feasible combinations of seed positions with given threshold distance between them
  """
  combinations = []
  if r1.inv == True:
    for i1 in list(r1.genome_pos.keys()):
      for i2 in list(r2.genome_pos.keys()):
        if 0 <= (i1 - i2) <= dist_threshold:
          combinations.append((i1, i2))
          
  else:
    for i1 in list(r1.genome_pos.keys()):
      for i2 in list(r2.genome_pos.keys()):
        if 0 <= (i2 - i1) <= dist_threshold:
          combinations.append((i1, i2))
          
  return combinations

#%%

if __name__ == '__main__':
  #Take the input and output paths with parser
  parser = argparse.ArgumentParser(description='Map pair-end read to reference genome')
  parser.add_argument("--genome",
                      action='store',
                      help='The path to genome file.',\
                      default='./data/genome.chr22.fa')
  parser.add_argument("--read1",
                      action='store',
                      help='The path to first paired end read file.',\
                      default='./data/output_5xCov1.fq')
  parser.add_argument("--read2",
                      action='store',
                      help='The path to second paired end read file.',\
                      default='./data_small/output_5xCov2.fq')
  parser.add_argument("--k",
                      action='store',
                      type = int,
                      help='Length of k-mers.',\
                      default='16')
  parser.add_argument("--outpath",
                      action='store',
                      help='The path to where the SAM file should be stored.',
                      default='./output.mod.sam')

  args = parser.parse_args()

  #############################################
  # Set params
  #############################################
  #Assign necessary constants to variables
  len_seed = args.k
  dist_threshold = 500
  constraint=2
  margin=0
  #############################################

  #Read the input and output files
  print('Read genome:')
  tstart = time.time()
  G = Genome(file = args.genome, k = len_seed)
  tend = time.time()
  print('Time to read genome:', tend - tstart)

  f1 = open(args.read1, 'r')
  f2 = open(args.read2, 'r')

  tstart = time.time()
  load2 = read_line(f2)
  tend = time.time()
  print('Time to read read2:', tend - tstart)

  tstart = time.time()
  load1 = read_line(f1)
  tend = time.time()
  print('Time to read read1:', tend - tstart)

  #Create the header of output file
  create_sam_header(args.outpath, G.genome_name, G.genome_length)

  #Load and process reads until the end of the read files
  tstart = time.time()
  iteration = 0
  while load1[0] and load2[0]:
    #Assign reads by using Read class
    r1 = Read(*load1, len_seed, inv=False)
    r2 = Read(*load2, len_seed, inv=True)
    
    #Calculate seed positions and take the inverse-complement of reads if necessary
    inv_flag = G.seed_positions(r1, r2)
    if inv_flag:
      r1.inverse()
      r2.inverse()
      G.seed_positions(r1, r2)
    
    #Calculate potential start positions of reads in genome
    r1.start_index(G)
    r2.start_index(G)

    #Calculate feasible combinations
    combs = feasible_combinations(r1, r2, dist_threshold=dist_threshold)
    best_comb = None
    min_cost = math.inf

    #If there is a feasible combination
    if len(combs) > 0:
      #Eliminate unnecessary start points of reads
      comb1 = list(set([i[0] for i in combs]))
      for start1 in list(r1.genome_pos.keys()):
        if start1 not in comb1:
          del r1.genome_pos[start1]
      comb2 = list(set([i[1] for i in combs]))
      for start2 in list(r2.genome_pos.keys()):
        if start2 not in comb2:
          del r2.genome_pos[start2]

      #Calculate the best combination of start indexes
      for start1 in r1.genome_pos:
        G.edit_distance(r1, start1, constraint=constraint, margin=margin)
      for start2 in r2.genome_pos:
        G.edit_distance(r2, start2, constraint=constraint, margin=margin) 
      for comb in combs:
        start1 = comb[0]
        start2 = comb[1]
        cost = r1.genome_pos[start1]['Cost'] + r2.genome_pos[start2]['Cost']
        if cost < min_cost:
          best_comb = comb
      
      #Take necessary variables for output file
      start1, start2 = best_comb[0], best_comb[1]
      pos1, pos2 = r1.genome_pos[start1]['Pos'], r2.genome_pos[start2]['Pos']
      alignment1, alignment2 = r1.genome_pos[start1]['Alignment'], r2.genome_pos[start2]['Alignment']
      inv1, inv2 = r1.inv, r2.inv

      #Write to the output file
      append_sam_alignment(args.outpath, reverse=r1.inv, secondary=False, 
                          QNAME=r1.idx, RNAME=G.genome_name, POS=pos1+1, 
                          MAPQ=0, CIGAR=alignment1, RNEXT="*", PNEXT=0, 
                          TLEN=0, SEQ="*", QUAL="*")
      append_sam_alignment(args.outpath, reverse=r2.inv, secondary=False, 
                          QNAME=r2.idx, RNAME=G.genome_name, POS=pos2+1, 
                          MAPQ=0, CIGAR=alignment2, RNEXT="*", PNEXT=0, 
                          TLEN=0, SEQ="*", QUAL="*")
    #If there is no feasible combinations
    else:
      append_sam_alignment(args.outpath, QNAME=r1.idx, RNAME=G.genome_name)
      append_sam_alignment(args.outpath, QNAME=r2.idx, RNAME=G.genome_name)

    iteration += 1
    if not (iteration % 100):
      print('Read', iteration,' Run-time:', time.time()-tstart)

    #Read the next 2 reads
    load1 = read_line(f1)
    load2 = read_line(f2)

#%%
#G = Genome(file = 'genome.txt', k = 16) 

# %%
