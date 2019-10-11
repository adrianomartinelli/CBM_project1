class Reads:
    forward = []
    backward = []

    def read_file(self, file, file_type):
        with open(file, 'r') as f:
            line = f.readline()
            while line:
                if(line.startswith('@')):
                    #Process id line
                    line = line.strip('@').strip('\n')
                    id, origin = line.split('/')

                    #Process seq and skip next line which is a +
                    seq = f.readline().strip('\n')
                    f.readline()

                    #Process phread
                    phread = f.readline().strip('\n')

                    if(file_type == 'forward'):
                        self.forward.append(Read(id, seq, phread, origin))
                    else:
                        self.backward.append(Read(id, seq, phread, origin))
                    
                    #Read next line
                    line = f.readline()
                else:
                    line = f.readline()

    def print_reads(self):
        for i in self.forward:
            print(i.seq, '\n', i.phread, sep = '')

class Read:
    def __init__(self, id, seq, phread, origin):
        self.__id = id
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
