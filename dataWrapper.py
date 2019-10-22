class Reads:
    pair1 = []
    pair2 = []

    def read_file(self, file, file_type):
        def __init__(self, file=None):
            self.pair1 = []
            self.pair2 = []
            if file:
                self.get_read_string(file+'1.fq')
                self.get_read_string(file+'2.fq')

        with open(file, 'r') as f:
            line = f.readline()
            while line:
                if(line.startswith('@')):
                    #Process id line
                    line = line.strip('@').strip('\n')
                    idx, origin = line.split('/')

                    #Process seq and skip next line which is a +
                    seq = f.readline().strip('\n')
                    f.readline()

                    #Process phread
                    phread = f.readline().strip('\n')

                    if(file_type == 'forward'):
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
            print(i.seq, '\n', i.phread, sep = '')

        for i in self.pair2:
            print('Pair2 reads:')
            print(i.seq, '\n', i.phread, sep = '')

class Read:
    def __init__(self, idx, seq, phread, origin):
        self.__id = idx
        self.__seq = seq
        self.__phread = phread
        self.__origin = origin

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
        self.__chars = ['$','A','C','G','T']

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

    def index(self):
        s = self.__genome + '$'
        m = [s]
        for _ in range(len(s)-1):
            m.append(self.rot_left(s))

        idx = np.argsort(m)
        m = [m[i] for i in idx]
        L = [i[self.__size] for i in m]
        F = [i[0] for i in m]
        return(0)

    def bwt(self, string):
        return(0)

    def rot_left(self,s):
        return(s[1:] + s[0])

import numpy as np

G = Genome('file.txt')
print(G.get_genome())
print(G.rot_left('Hello'))
