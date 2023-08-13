import random
from scipy.stats import gamma, expon


cdef class RLDfactory:
    def __init__(self, str protocol):
        # ccs
        if protocol == "ccs":
            self.__alpha = 6.02924043340785
            self.__loc = 10682.156155519911
            self.__scale = 466.56937340727563
        # clr
        if protocol == "clr":
            self.__alpha = 2.1632592929739274
            self.__loc = -545.7711528901336
            self.__scale = 3865.8926027058233
        # ont
        if protocol == "ont":
            self.__loc = 213.989
            self.__scale = 6972.532
    
    cdef int gamma_nextl(self):
        cdef int rl = int(gamma.rvs(self.__alpha, self.__loc, self.__scale))
        return rl
    
    cdef int expon_nextl(self):
        cdef int rl = int(expon.rvs(self.__loc, self.__scale))
        return rl