##OverlapCheck.py
import optparse
from FNpool import  getSseq, seqonly, FastaTool, musclecall
from FNpool import OverlapCheck2 as OverlapCheck
import sys

class overlap(object):

	#def __init__(self,DBnameS,SeqS0_IDS,SeqS1_IDS):
	#	self.DBnameS = DBnameS
	#	self.SeqS0_IDS = SeqS0_IDS
	#	self.SeqS1_IDS = SeqS1_IDS
	#	self.SeqS0 = getSseq(self.DBnameS,self.SeqS0_IDS,'1')
	#	self.SeqS1 = getSseq(self.DBnameS,self.SeqS1_IDS,'1')

	def __init__(self,SeqS0,SeqS1):
		self.SeqS0 = SeqS0
		self.SeqS1 = SeqS1
		self.Seq0len = len(self.SeqS0)
		self.Seq1len = len(self.SeqS1)
		#OverlapL = self.overlapCall()
		

	def OverlapCall(self):
		SeqS0 = self.SeqS0
		SeqS1 = self.SeqS1		

		SeqS0F = seqonly(SeqS0)
		SeqS0R = seqonly(FastaTool(SeqS0).reversecomplement())
		SeqS1F = seqonly(SeqS1)
		SeqS1R = seqonly(FastaTool(SeqS1).reversecomplement())
	
		#Seq0 ---------->
		#Seq1      -------->
		T0_H1_OvI = OverlapCheck(SeqS0F,SeqS1F,1)
		#Seq1 <---------
		#Seq0      <----------
		H1_T0_OvI = OverlapCheck(SeqS1R, SeqS0R,1)


		#Seq0 ---------->
		#Seq1      <--------
		T0_T1_OvI = OverlapCheck(SeqS0F,SeqS1R,1)
		#Seq1 ---------->
		#Seq0     <-----------
		T1_T0_OvI = OverlapCheck(SeqS1F,SeqS0R,1)
	

		#Seq1 ---------->
		#Seq0     ------------>
		T1_H0_OvI = OverlapCheck(SeqS1F,SeqS0F,1)
		#Seq0 <----------
		#Seq1     <-----------	
		H0_T1_OvI = OverlapCheck(SeqS0R,SeqS1R,1)


		#Seq0 <---------
		#Seq1       ---------->
		H0_H1_OvI = OverlapCheck(SeqS0R, SeqS1F,1)
		#Seq1 <---------
		#Seq0       ---------->
		H1_H0_OvI = OverlapCheck(SeqS1R, SeqS0F,1)
	
		OverlapL = [T0_H1_OvI, H1_T0_OvI, T0_T1_OvI, T1_T0_OvI, T1_H0_OvI, H0_T1_OvI, H0_H1_OvI, H1_H0_OvI]
		return(OverlapL)

	def MSA(self, AlignNumI):
		SeqS0 = self.SeqS0
		SeqS1 = self.SeqS1		

		SeqS0F = seqonly(SeqS0)
		SeqS0R = seqonly(FastaTool(SeqS0).reversecomplement())
		SeqS1F = seqonly(SeqS1)
		SeqS1R = seqonly(FastaTool(SeqS1).reversecomplement())
		
		if AlignNumI == 0:
			MSAS = musclecall(">S0F\n"+SeqS0F+"\n>S1F\n"+SeqS1F,1)
		elif AlignNumI == 1:
			MSAS = musclecall(">S1R\n"+SeqS1R+"\n>S0R\n"+SeqS0R,1)

		elif AlignNumI == 2:
			MSAS = musclecall(">S0F\n"+SeqS0F+"\n>S1R\n"+SeqS1R,1)
		elif AlignNumI == 3:
			MSAS = musclecall(">S1F\n"+SeqS1F+"\n>S0R\n"+SeqS0R,1)
		
		elif AlignNumI == 4:
			MSAS = musclecall(">S1F\n"+SeqS1F+"\n>S0F\n"+SeqS0F,1)
		elif AlignNumI == 5:
			MSAS = musclecall(">S0R\n"+SeqS0R+"\n>S1R\n"+SeqS1R,1)

		elif AlignNumI == 6:
			MSAS = musclecall(">S0R\n"+SeqS0R+"\n>S1F\n"+SeqS1F,1)
		elif AlignNumI == 7:
			MSAS = musclecall(">S1R\n"+SeqS1R+"\n>S0F\n"+SeqS0F,1)

		return MSAS

def main():
	usage = "python OverlapCheck.py -d DBpath -1 1st_SeqID -2 2nd_SeqID"
	opt = optparse.OptionParser(usage)
	opt.add_option("-d",help="sequences database path")
	opt.add_option("-1",help="1st sequnce ID",dest="a")
	opt.add_option("-2",help="2nd sequnce ID",dest="b")
	opt.add_option("--align",help="2nd sequnce ID",dest="AlignNumI",default="-1")
	(options, args) = opt.parse_args()

	DBnameS, SeqS0_IDS, SeqS1_IDS, AlignNumI = options.d, options.a, options.b, int(options.AlignNumI)
	if None in [DBnameS, SeqS0_IDS, SeqS1_IDS]:
		print("missing variable")
		print(usage)
	else:
		
		SeqS0 = getSseq(DBnameS,SeqS0_IDS,'1')
		SeqS1 = getSseq(DBnameS,SeqS1_IDS,'1')
		print(DBnameS, SeqS0_IDS, SeqS1_IDS)
		ovelapL = overlap(SeqS0,SeqS1)
		print(ovelapL.OverlapCall())

	if AlignNumI != -1:
		print(ovelapL.MSA(AlignNumI))

if __name__ == "__main__":
	main()

