#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import collections


 # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../../bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


class Measure:
	def __init__(self,chr,pos,measure):
		self.chr=chr
		self.pos=pos
		self.measure=measure


class MeasureReader:
	"""
	A light-weight pileup reader
	
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
		self.__buffer=[]

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()
		
	def buffer(self,tobuffer):
	    self.__buffer.append(tobuffer)


	def next(self):
		# first empty the buffer
		if len(self.__buffer)>0:
		    return self.__buffer.pop()
		
		# if empty read the next line    
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line.startswith("#"):
				continue
			if line != "":
				break
		
		a=line.split("\t")
		chr,pos,measure=a[0],a[1],a[-1]
		return Measure(chr,int(pos),float(measure))


parser = argparse.ArgumentParser(description="""           
Description
-----------
Identifies SNPs, true positives (simulated and identified), false positives (not simulated but identified) and false negatives (simulated but not identified)
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output

Authors
-------
    Robert Kofler
""")




parser.add_argument("--stat1", type=str, required=True, dest="stat1", default=None, help="A pileup file")
parser.add_argument("--stat2", type=str, required=True, dest="stat2", default=None, help="the minimum coverage")
parser.add_argument("--take", type=str, required=False, dest="take", default="take1", help="the minimum coverage [take1 | take2 | small | large]")
args = parser.parse_args()

take=args.take

mr1=MeasureReader(args.stat1)
mr2=MeasureReader(args.stat2)


activechr=None
while(True):
	try:
		m1=mr1.next()
		m2=mr2.next()
	except StopIteration:
		break
	if(activechr is None): # initialize the starting chromosomes
		if(m1.chr !=m2.chr):
			raise Exception("Invalid start of the two files")
		activechr=m1.chr
	if(m1.chr==m2.chr and m1.chr!=activechr):
		activechr=m1.chr
	
	
	if(m1.chr!=m2.chr):
		if(m1.chr==activechr):
			mr2.buffer(m2)
			continue
		elif(m2.chr==activechr):
			mr1.buffer(m1)
			continue
		else:
			raise Exception(m1.chr +" " +m2.chr)
	elif(m1.pos>m2.pos):
		mr1.buffer(m1)
		continue
	elif(m2.pos>m1.pos):
		mr2.buffer(m2)
		continue
	
	assert(m1.chr==m2.chr)
	assert(m1.pos==m2.pos)
	
	valtoprint=None
	# [take1 | take2 | small | large]
	if(take=="take1"):
		valtoprint=m1.measure
	elif(take=="take2"):
		valtoprint=m2.measure
	elif(take=="small"):
		if(m1.measure<m2.measure):
			valtoprint=m1.measure
		else:
			valtoprint=m2.measure
	elif(take=="large"):
		if(m1.measure> m2.measure):
			valtoprint=m1.measure
		else:
			valtoprint=m2.measure
	else:
		raise Exception("unknown")

	print "{0}\t{1}\t{2}".format(m1.chr,m1.pos,valtoprint)
				

    
        
    

