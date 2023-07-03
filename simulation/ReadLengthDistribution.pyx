import random
from scipy.stats import gamma


# class RLDfactory_rldfile:
# 	def __init__(self,rldfile):
# 		self.__accuracy=1000000
# 		self.__lar=self.__load_rldfile(rldfile, self.__accuracy)
# 		self.__len=len(self.__lar)
	
# 	def __load_rldfile(self,file,accuracy):
# 		"""
# 		Readlength count
# 		5000	50
# 		7000	10
# 		8000	20
# 		"""
# 		totsum=0
# 		tsh={}
# 		for l in open(file):
# 			l=l.rstrip("\n")
# 			a=None
# 			if "\t" in l:
# 				a=l.split("\t")
# 			else:
# 				a=l.split(" ")
# 			rl=int(a[0])
# 			count=int(round(float(a[1])))
# 			if rl not in tsh:
# 				totsum+=count
# 				tsh[rl]=count
		
# 		lar=[]
# 		for rl,c in tsh.items():
# 			normcount=int(float(c*accuracy)/float(totsum))
# 			for i in range(0,normcount):
# 				lar.append(rl)
# 		return lar
	
# 	def next(self):
# 		index=random.randint(0,self.__len-1)
# 		return self.__lar[index];


cdef class RLDfactory_gamma:
    def __init__(self, float alpha, float loc, float beta):
        self.__alpha = alpha
        self.__loc = loc
        self.__beta = beta
    
    cdef int nextl(self):
        cdef int rl = int(gamma.rvs(self.__alpha, self.__loc, scale = self.__beta))
        return rl


cdef RLDfactory_gamma get_rld_factory(int mean, float alpha, float loc, float beta, str rldfile):
    if not rldfile:
        if alpha is None:
            raise Exception("Either provide a mean read length and standard deviation or a file with the read length distribution")
        print("Using gamma read length distribution with alpha:{0}, loc:{1} and beta:{2}".format(alpha, loc, beta))
        return RLDfactory_gamma(alpha, loc, beta)

    # else:
    #     print "Using read length distribution from file {0}".format(file)
    #     return RLDfactory_rldfile(file)