import sys
import random
import re


cdef class PopGenDefinitionReader:
    def __init__(self, str file, SequenceContainer sc):
        self.__filename = file
        self.__sc = sc
        self.__chasis = ""
        self.tuples = self.read_tuples()

    cdef __parse_chasis(self, str line):
        cdef:
            list a = line.split("=")
            str r = a[1].strip(" \"")
        self.__chasis = r

    cpdef str get_chasis(self):
        return self.__chasis

    cdef list read_tuples(self):
        cdef:
            list toret = [], sp
            str l
        
        for l in open(self.__filename):
            l = l.rstrip("\n")
            l = re.sub(r"#(.*)","",l) # remove comments
            l = l.strip()
            if l == "": # remove empty lines
                continue
            if(l.startswith("chassis")):
                self.__parse_chasis(l)
            elif "=" in l:
                self.__sc.addDefinition(l)
            else:
                sp = l.split(" ")
                toret.append(sp)
        
        self.insertions = len(toret)
        self.popsize = len(sp) - 1 # minus one since the first column is the position within the chasis
        return toret
    
    cpdef list read_transposed(self):
        if len(self.tuples) < 1:
            raise Exception("Invalid popgendef (population genome definition) file; at least one TE insertion site is required")
        
        cdef:
            list toret = [], tmp
            int i, k
            tuple toa
        for i in range(1, self.popsize+1):
            tmp = []
            for k in range(0, self.insertions):
                toa = (int(self.tuples[k][0]), self.tuples[k][i]) # position, tedefinition
                if toa[1] == "*":
                    continue
                tmp.append(toa)

            toret.append(tmp)
        
        return toret