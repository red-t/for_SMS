#!/usr/bin/env python
import sys
import random
import collections
import numpy as np


def get_rld_factory(mean,p,file):
	if(file is None):
		if mean is None or p is None:
			raise Exception("Either provide a mean read length and standard deviation or a file with the read length distribution")
		print "Using gaussian read length distribution with mean {0} and standard deviation {1}".format(mean,p)
		return RLDfactory_nagative_binomial(mean,p)
	else:
		print "Using read length distribution from file {0}".format(file)
		return RLDfactory_rldfile(file)

	

class RLDfactory_rldfile:
	def __init__(self,rldfile):
		self.__accuracy=1000000
		self.__lar=self.__load_rldfile(rldfile, self.__accuracy)
		self.__len=len(self.__lar)
	
	def __load_rldfile(self,file,accuracy):
		"""
		Readlength count
		5000	50
		7000	10
		8000	20
		"""
		totsum=0
		tsh={}
		for l in open(file):
			l=l.rstrip("\n")
			a=None
			if "\t" in l:
				a=l.split("\t")
			else:
				a=l.split(" ")
			rl=int(a[0])
			count=int(round(float(a[1])))
			if rl not in tsh:
				totsum+=count
				tsh[rl]=count
		
		lar=[]
		for rl,c in tsh.items():
			normcount=int(float(c*accuracy)/float(totsum))
			for i in range(0,normcount):
				lar.append(rl)
		return lar
	
	def next(self):
		index=random.randint(0,self.__len-1)
		return self.__lar[index];
		
		

class RLDfactory_nagative_binomial:
	def __init__(self,mean,p):
		self.__m=mean*p/(1-p)
		self.__p=p
	
	def next(self):
		m=self.__m
		p=self.__p
		rl=int(np.random.negative_binomial(m, p))
		return rl