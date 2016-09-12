import sys

class NucSeq(object):
	""" A generic nucleotide sequence defined 5'-3'
	Attributes:
		seq: A string sequence
		doublestranded: True/False
		strand1: DNA/RNA 5'-3'
		strand2: DNA/RNA/None 3'-5'
	"""
	doublestranded = False 
	strand1 = None
	strand2 = None
	
	def __init__(self,seq):
		self.seq = seq.upper()
		
	def __str__(self):
		return self.seq
	
	@classmethod
	def fromfile(cls,file):
		with open(file, "r") as f:
			lines = f.readlines()
		lines = [x.strip() for x in lines]
		lines = [x for x in lines if x != '']
		lines = ["".join(x.split()) for x in lines]
		lines = [x.upper() for x in lines if x[0] != '>']
		return cls(''.join(lines))
		
	def toBNGSpecies(self):
		strand1_list = [build_base(x,self.strand1) for x in self.seq]
		strand2_list = []
		bondnum = 0;
		numbases = len(strand1_list)
		for i in range(numbases-1):
			bondnum = bondnum+1
			strand1_list[i] = strand1_list[i].replace(",p3,",",p3!"+str(bondnum)+",")
			strand1_list[i+1] = strand1_list[i+1].replace(",p5,",",p5!"+str(bondnum)+",")
		if self.strand2 is not None:
			strand2_list = [build_base(x,self.strand2) for x in sequence_complement(self.seq)]
			for i in range(numbases-1):
				bondnum = bondnum+1
				strand2_list[i] = strand2_list[i].replace(",p5,",",p5!"+str(bondnum)+",")
				strand2_list[i+1] = strand2_list[i+1].replace(",p3,",",p3!"+str(bondnum)+",")
			for i in range(numbases):
				bondnum = bondnum+1
				strand1_list[i] = strand1_list[i].replace(",c",",c!"+str(bondnum))
				strand2_list[i] = strand2_list[i].replace(",c",",c!"+str(bondnum))
			
		baselist = strand1_list+strand2_list
		bngspeciesstring = ".".join(baselist)
		return bngspeciesstring
		

def base_complement(base,basetype='DNA'):
	dna = {'A':'T','T':'A','U':'A','C':'G','G':'C'}
	rna = {'A':'U','T':'A','U':'A','C':'G','G':'C'}
	base = base.upper()
	if base not in ['A','T','U','C','G']:
		return None
	if basetype=='DNA':
		return dna[base]
	if basetype=='RNA':
		return rna[base]
		
def sequence_complement(seq,basetype='DNA'):
	baselist = [base_complement(x,basetype) for x in seq]
	return ''.join(baselist)
		
def build_base(base,basetype='DNA'):
	"""Nuc(b~A~T~C~G~U,na~d~r,p5,p3,c)"""
	str1 = "b~"+base
	letter = {'DNA':'d','RNA':'r'};
	str2 = "na~"+letter[basetype]
	basestring = "Nuc("+str1+","+str2+",p5,p3,c)"
	return basestring
	
class ssRNA(NucSeq):
	"""A single stranded RNA sequence. Subclass of NucSeq. """
	doublestranded = False
	strand1 = 'RNA'
	
class ssDNA(NucSeq):
	"""A single-stranded DNA sequence. Subclass of NucSeq. """
	doublestranded = False
	strand1 = 'DNA'
	
class dsDNA(NucSeq):
	"""A double-stranded RNA sequence. Subclass of NucSeq. """
	doublestranded = True
	strand1 = 'DNA'
	strand2 = 'DNA'
	
if __name__ == '__main__':
	file = sys.argv[1]
	a = dsDNA.fromfile(file)
	print a.toBNGSpecies()
	