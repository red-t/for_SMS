#!/usr/bin/env python


        

class TESequence:
	def __init__(self,sequence,id,tsd):
		self.sequence=sequence
		self.id=id
		self.tsd=tsd

	


class SeqInserter:
	
	@classmethod
	def insertSequence(cls,ref,pos,toinsert,fout="",seqid=""):
		seq=toinsert.sequence
		tsd=toinsert.tsd
		tlp=pos-tsd
		if(tlp<0):
			raise Exception("Invalid position of TE; Insertion position minus TSD must be larger than zero")

		left=ref[:tlp] ## hmm guess we need tlp instead of pos
		trp=tlp+tsd
		tsdseq=ref[tlp:trp]
		right=ref[trp:] ## guess we need trp instead of pos

		# write out simulated insertion sequence
		if fout:
			p_left = 0
			p_right = len(ref)-1

			if tlp-2000 > 0:
				p_left = tlp-2000
			if trp+2000 < len(ref)-1:
				p_right = trp+2000
			
			t_left = ref[p_left:tlp]
			t_right = ref[trp:p_right]
			fout.write("{0}\t{1}\n".format(seqid, t_left+tsdseq+seq+tsdseq+t_right))

		return left+tsdseq+seq+tsdseq+right
	

	@classmethod
	def insertSequences(cls,ref,posinstuples,fout_name="",seqidtup=""):
		"""
		ref: reference genome
		posinstuples
		[(1,TESequence), (100, TESEquence)....]
		"""
		if fout_name:
			fout=open(fout_name, "a")
			tmp=sorted(posinstuples,key=lambda i:-i[0])
			idtmp=sorted(seqidtup,key=lambda i:-i[0])
			seq=ref
			i = 0
			for pos,toins in tmp:
				seq=SeqInserter.insertSequence(seq,pos,toins,fout,idtmp[i][1])
				i += 1

			fout.close()
		else:
			tmp=sorted(posinstuples,key=lambda i:-i[0])
			seq=ref
			for pos,toins in tmp:
				seq=SeqInserter.insertSequence(seq,pos,toins)

		return seq


"""
class TEDefinitionReader:

	A light-weight reader of TE definitions;
	returns a TEInsertionDefinition
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
		self.__prevheader=None

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()
	
	def next(self):
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
		a=line.split(" ")
		pos=int(a.pop(0))
		comment=a.pop(0)
		tedefs=[]
		for tmp in a:
			site=TEInsertionSite.parse_key(tmp)
			tedefs.append(site)	
		# position, comment
		e = TEInsertionDefinition(pos,comment,tedefs)
		return e
	
	@classmethod
	def readall(cls,file):
		entries=[]
		for e in TEDefinitionReader(file):
			entries.append(e)
		return entries

class TEDefinitionWriter:

	Write the content of a TE definition to a file
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"w")
		
	def write(self,tedef):
		position,comment,tedefinitions=tedef.position,tedef.comment,tedef.tesites
		topr=[str(position),comment]
		for ted in tedefinitions:
			topr.append(ted.getKey())
		form=" ".join(topr)
		self.__filehandle.write(form+"\n")
		return 1
	
	def close(self):
		self.__filehandle.close()
	
	@classmethod
	def writeall(cls,file,tedefentries):
		gw=TEDefinitionWriter(file)
		for e in tedefentries:
			gw.write(e)
		gw.close()
		return 1
"""