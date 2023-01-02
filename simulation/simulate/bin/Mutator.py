#!/usr/bin/env python
import random
import math
	
class PacBioMutator:

	def __init__(self,errorrate,errfrac):
		self.__er=errorrate
		self.__tr={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C'],'N':['N']}
		self.__ins=['A','T','C','G']
		
		errfrac = errfrac.split(':')
		self.__misf=float(errfrac[0])/100
		self.__delf=(float(errfrac[0])+float(errfrac[1]))/100
	
	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq
		
		err_f = open('ErrorAll.txt', 'a')
		muts=[]
		bin_len = len(seq)/50
		errs=[0,0,0]
		for i in range(0,len(seq)):
			if random.random()<self.__er:
				muts.append(i)
		
		ins=self.__ins
		lseq=list(seq)
		for i in reversed(muts):
			rand = random.random()
			if rand < self.__misf:
				#mismatch
				tr=self.__tr[lseq[i]]
				random.shuffle(tr)
				lseq[i]=tr[0]
				err_f.write('{0}M\n'.format(i/bin_len))
				errs[0] += 1
			elif rand >= self.__misf and rand < self.__delf:
				#deletion
				del(lseq[i])
				err_f.write('{0}D\n'.format(i/bin_len))
				errs[1] += 1
			elif rand >= self.__delf:
				#insertion
				random.shuffle(ins)
				lseq.insert(i,ins[0])
				err_f.write('{0}I\n'.format(i/bin_len))
				errs[2] += 1

		err_f.close()
		return "".join(lseq), errs






class PoisonSeqMutator:

	def __init__(self,errorrate):
		self.__er=errorrate
		self.__tr={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C'],'N':['N']}
	
	def __getPoisson(self,lam):
		L= math.exp(-lam)
		p=1.0
		k=0
		while(True):
			k+=1
			p*=random.random()
			if(p<L):
				break
		return k-1
	
	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq

		lseq=list(seq)
		for i in range(0,len(lseq)):
			# test every base if it should have an error
			if random.random()<self.__er:
				tr=self.__tr[lseq[i]]
				random.shuffle(tr)
				lseq[i]=tr[0]
		return "".join(lseq)
	
	
	
	
	

	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq
		aver=float(len(seq))*self.__er
		errors=self.__getPoisson(aver)
		lseq=list(seq) # sequence in a list
		pastpos=set([])  
		for i in range(0,errors):
			pos=random.randint(0,len(lseq)-1)
			if pos in pastpos:
				pos=random.randint(0,len(lseq)-1)
			pastpos.add(pos)
			tr=self.__tr[lseq[pos]]
			random.shuffle(tr)
			lseq[pos]=tr[0]
		return "".join(lseq)


class ExhaustiveSeqMutator:
	
	def __init__(self,errorrate):
		self.__er=errorrate
		self.__tr={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C'],'N':['N']}
	

	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq

		lseq=list(seq)
		for i in range(0,len(lseq)):
			# test every base if it should have an error
			if random.random()<self.__er:
				tr=self.__tr[lseq[i]]
				random.shuffle(tr)
				lseq[i]=tr[0]
		return "".join(lseq)
		
