import numpy as np
import matplotlib.pyplot as plt
import os

class Read():
  def __init__(self, idx, origin, string, quality):
    self.idx = idx
    self.origin = origin
    self.string = string
    self.quality = quality


class Reads():
  def __init__(self, file=None):
    self.pair1 = []
    self.pair2 = []
    if file:
      self.get_read_string(file+'1.fq')
      self.get_read_string(file+'2.fq')

  def get_read_string(self, file):
    with open(file, 'r') as f:
      line = f.readline()
      while line:
        if(line.startswith('@')):
          line = line.lstrip('@').rstrip('\n')
          idx, origin = line.split('/')

          string = f.readline().rstrip('\n')
          f.readline()

          quality = f.readline().rstrip('\n')

          if(origin == '1'):
              self.pair1.append(Read(idx, origin, string, quality))
          else:
              self.pair2.append(Read(idx, origin, string, quality))

          #Read next line
          line = f.readline()
        else:
          line = f.readline()


class Genome():
  def __init__(self, file=None):
    self.genome_string = ''
    self.genome_length = 0
    self.chars = ['$','A','C','G','T']
    if file:
      self.read_genome(file)

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
    for c in self.chars:
      occ=0
      for l in self.L:
        if l<c:
          occ+=1
      self.C[c] = occ

  def occurence_matrix(self):
    self.Occ={}
    for c in self.chars:
      occ=0
      for idx, l in enumerate(self.L):
        if l==c:
          occ+=1
        self.Occ[(c, idx)] = occ

  def next_char(self, char):
    idx = self.chars.index(char)
    if idx == len(self.chars) - 1:
      return -1
    else:
      return self.chars[idx+1]

  def fm_index(self, P):
    i = len(P)
    c = P[i-1]
    sp = self.C[c] + 1
    if self.next_char(c)==-1:
      ep = len(self.L)
    else:
      ep = self.C[self.next_char(c)]

    while sp <= ep and i >= 2:
      c = P[i-2]
      sp = self.C[c] + self.Occ[(c, sp - 2)] + 1
      ep = self.C[c] + self.Occ[(c, ep-1)]
      i=i-1

    range1 = self.i[ep-1]
    range2 = self.i[sp-1]

    return range1, range2

  def calculate_best_position(self, read, left, right):
    best_position = left
    best_score    = -1000
    read_length   = len(read)
    for position in range (right-left+1):
      current_position = left + position
      score = 0
      for idx, c in enumerate(read):
        if c == self.genome_string[current_position+idx]:
          score = score + 1
        else:
          score = score - 1
      if score > best_score:
        best_score = score
        best_position = current_position
    return best_position


example_genome = 'ExampleGenome.txt'
G = Genome(example_genome)
G.bwt_index()
G.occurence_array()
G.occurence_matrix()

print(G.genome_string)

P = 'AGT'
range1, range2 = G.fm_index(P)

# Print low range - top range
if range1>range2:
  range1, range2 = range2, range1
print('Search Range: ', range1+1, range2+1)