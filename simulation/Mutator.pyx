import random
import math

cdef class PacBioMutator:
    def __init__(self, float errorrate, str errfrac):
        frac = errfrac.split(':')
        self.__er = errorrate
        self.__tr = {'A':['T','C','G'], 'T':['A','C','G'], 'C':['A','T','G'], 'G':['A','T','C'], 'N':['N']}
        self.__ins = ['A','T','C','G']
        self.__misf = float(frac[0])/100
        self.__delf = (float(frac[0]) + float(frac[1]))/100

    cdef str mutateseq(self, str seq):
        if(self.__er < 0.0000000001):
            return seq
        
        # err_f = open('ErrorAll.txt', 'a')
        cdef:
            list muts = []
            int i, j, length = len(seq)
            list lseq = list(seq)
            float prob
            list tr
            # int bin_len = len(seq)/50
            # list errs = [0,0,0]
        for i in range(0, length):
            prob = random.random()
            if prob < self.__er:
                muts.append(i)

        for i in reversed(muts):
            prob = random.random()
            if prob < self.__misf:
                # mismatch
                j = random.randint(0, 2)
                tr = self.__tr[lseq[i]]
                if len(tr) > 1:
                    lseq[i] = tr[j]
                # err_f.write('{0}M\n'.format(i/bin_len))
                # errs[0] += 1
            elif prob >= self.__misf and prob < self.__delf:
                # deletion
                del(lseq[i])
                # err_f.write('{0}D\n'.format(i/bin_len))
                # errs[1] += 1
            elif prob >= self.__delf:
                # insertion
                j = random.randint(0, 3)
                lseq.insert(i, self.__ins[j])
                # err_f.write('{0}I\n'.format(i/bin_len))
                # errs[2] += 1

        # err_f.close()
        # return "".join(lseq), errs
        return "".join(lseq)


cdef class PoisonSeqMutator:
    def __init__(self, float errorrate):
        self.__er = errorrate
        self.__tr = {'A':['T','C','G'], 'T':['A','C','G'], 'C':['A','T','G'], 'G':['A','T','C'], 'N':['N']}

    cdef int __getPoisson(self, float lam):
        cdef:
            float L = math.exp(-lam)
            float p = 1.0
            int k = 0
        while(True):
            k += 1
            p *= random.random()
            if(p < L):
                break
        return k-1
        
    cdef str mutateseq(self, str seq):
        if(self.__er < 0.0000000001):
            return seq
        
        cdef:
            list lseq = list(seq)
            int i, j, length = len(lseq)
            float prob
            list tr
        for i in range(0, length):
            # test every base if it should have an error
            prob = random.random()
            if prob < self.__er:
                j = random.randint(0, 2)
                tr = self.__tr[lseq[i]]
                if len(tr) > 1:
                    lseq[i] = tr[j]
        
        return "".join(lseq)


cdef class ExhaustiveSeqMutator:
    def __init__(self, float errorrate):
        self.__er = errorrate
        self.__tr = {'A':['T','C','G'], 'T':['A','C','G'], 'C':['A','T','G'], 'G':['A','T','C'], 'N':['N']}
    
    cdef str mutateseq(self, str seq):
        if(self.__er < 0.0000000001):
            return seq
        
        cdef:
            list lseq = list(seq)
            int i, j, length = len(lseq)
            float prob
            list tr
        for i in range(0, length):
            # test every base if it should have an error
            prob = random.random()
            if prob < self.__er:
                j = random.randint(0, 2)
                tr = self.__tr[lseq[i]]
                if len(tr) > 1:
                    lseq[i] = tr[j]
        
        return "".join(lseq)