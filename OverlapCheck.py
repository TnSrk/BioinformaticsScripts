##OverlapCheck.py
import optparse
from FNpool import  getSseq, seqonly, FastaTool, musclecall, ConcurrentCall, consensusExtractKW, STDERR
from FNpool import OverlapCheck2 as OverlapCheck
import sys

class overlap(object):

	#def __init__(self,DBnameS,SeqS0_IDS,SeqS1_IDS):
	#	self.DBnameS = DBnameS
	#	self.SeqS0_IDS = SeqS0_IDS
	#	self.SeqS1_IDS = SeqS1_IDS
	#	self.SeqS0 = getSseq(self.DBnameS,self.SeqS0_IDS,'1')
	#	self.SeqS1 = getSseq(self.DBnameS,self.SeqS1_IDS,'1')

	def __init__(self,SeqS0,SeqS1,*arg):
		self.SeqS0 = FastaTool(SeqS0)
		self.SeqS1 = FastaTool(SeqS1)
		#OverlapL = self.overlapCall()
		if len(arg) == 1:
			self.OverlapI(int(arg[0]))
		else:
			self.OverlapI()
		#self.SeqS0NameS = FastaTool(SeqS0).name()
	
	def OverlapCall0(self):
		#SeqS0 = self.SeqS0
		#SeqS1 = self.SeqS1		

		SeqS0F = self.SeqS0.seqonly()
		SeqS0R = self.SeqS0.reversecomplement()
		SeqS1F = self.SeqS1.seqonly()
		SeqS1R = self.SeqS1.reversecomplement()
	
		#Seq0 ---------->
		#Seq1      -------->
		T0_H1_OvI = OverlapCheck(SeqS0F,SeqS1F,1)
		#T0_H1_OvI = [ SeqS0F, SeqS1F, 1 ] 
		#Seq1 <---------
		#Seq0      <----------
		H1_T0_OvI = OverlapCheck(SeqS1R, SeqS0R,1)
		#H1_T0_OvI = [ SeqS1R, SeqS0R, 1 ]

		#Seq0 ---------->
		#Seq1      <--------
		T0_T1_OvI = OverlapCheck(SeqS0F,SeqS1R,1)
		#T0_T1_OvI = [SeqS0F, SeqS1R, 1 ]
		#Seq1 ---------->
		#Seq0     <-----------
		T1_T0_OvI = OverlapCheck(SeqS1F,SeqS0R,1)
		#T1_T0_OvI = [SeqS1F, SeqS0R, 1 ]
	

		#Seq1 ---------->
		#Seq0     ------------>
		T1_H0_OvI = OverlapCheck(SeqS1F,SeqS0F,1)
		#T1_H0_OvI = [SeqS1F, SeqS0F, 1 ]

		#Seq0 <----------
		#Seq1     <-----------	
		H0_T1_OvI = OverlapCheck(SeqS0R,SeqS1R,1)
		#H0_T1_OvI = [SeqS0R, SeqS1R, 1 ]


		#Seq0 <---------
		#Seq1       ---------->
		H0_H1_OvI = OverlapCheck(SeqS0R, SeqS1F,1)
		#H0_H1_OvI = [SeqS0R, SeqS1F,1 ]

		#Seq1 <---------
		#Seq0       ---------->
		H1_H0_OvI = OverlapCheck(SeqS1R, SeqS0F,1)
		#H1_H0_OvI = [SeqS1R, SeqS0F, 1]

		
		OverlapL = [T0_H1_OvI, H1_T0_OvI, T0_T1_OvI, T1_T0_OvI, T1_H0_OvI, H0_T1_OvI, H0_H1_OvI, H1_H0_OvI]

		return(OverlapL)
		
		

	def OverlapCall1(self,cupNumI):
		#SeqS0 = self.SeqS0
		#SeqS1 = self.SeqS1
		EndLengthI = 20		

		SeqS0F = self.SeqS0.seqonly()
		SeqS0R = self.SeqS0.reversecomplement()
		SeqS1F = self.SeqS1.seqonly()
		SeqS1R = self.SeqS1.reversecomplement()
	
		#Seq0 ---------->
		#Seq1      -------->
		#T0_H1_OvI = OverlapCheck(SeqS0F,SeqS1F,1)
		T0_H1_OvI = [ SeqS0F, SeqS1F, 1, EndLengthI ] 
		#Seq1 <---------
		#Seq0      <----------
		#H1_T0_OvI = OverlapCheck(SeqS1R, SeqS0R,1)
		H1_T0_OvI = [ SeqS1R, SeqS0R, 1, EndLengthI ]

		#Seq0 ---------->
		#Seq1      <--------
		#T0_T1_OvI = OverlapCheck(SeqS0F,SeqS1R,1)
		T0_T1_OvI = [SeqS0F, SeqS1R, 1, EndLengthI ]
		#Seq1 ---------->
		#Seq0     <-----------
		#T1_T0_OvI = OverlapCheck(SeqS1F,SeqS0R,1)
		T1_T0_OvI = [SeqS1F, SeqS0R, 1, EndLengthI ]
	

		#Seq1 ---------->
		#Seq0     ------------>
		#T1_H0_OvI = OverlapCheck(SeqS1F,SeqS0F,1)
		T1_H0_OvI = [SeqS1F, SeqS0F, 1, EndLengthI ]

		#Seq0 <----------
		#Seq1     <-----------	
		#H0_T1_OvI = OverlapCheck(SeqS0R,SeqS1R,1)
		H0_T1_OvI = [SeqS0R, SeqS1R, 1, EndLengthI ]


		#Seq0 <---------
		#Seq1       ---------->
		#H0_H1_OvI = OverlapCheck(SeqS0R, SeqS1F,1)
		H0_H1_OvI = [SeqS0R, SeqS1F,1, EndLengthI ]

		#Seq1 <---------
		#Seq0       ---------->
		#H1_H0_OvI = OverlapCheck(SeqS1R, SeqS0F,1)
		H1_H0_OvI = [SeqS1R, SeqS0F, 1, EndLengthI]

		
		OverlapL = ConcurrentCall(OverlapCheck, [T0_H1_OvI, H1_T0_OvI, T0_T1_OvI, T1_T0_OvI, T1_H0_OvI, H0_T1_OvI, H0_H1_OvI, H1_H0_OvI], cupNumI)

		return(OverlapL)

	def OverlapI(self,*arg):
		if len(arg) == 0:
			OverlapL = self.OverlapCall0()
		elif len(arg) == 1:
			OverlapL = self.OverlapCall1(arg[0])

		matchI = -1
		matchF = 0.0
		AI = -1
		for i in (0,2,4,6):
			matchF0 = abs(OverlapL[i] - OverlapL[i+1])*1.0/(OverlapL[i] + OverlapL[i+1])*1.0/2.0
			if OverlapL[i] != -1 and (1.0 - matchF0) > matchF:
				matchI = OverlapL[i]
				matchF = 1.0 - matchF0
				AI = i

		self.OverlapI = matchI
		self.OverlapF = matchF
		self.OverlapL = OverlapL
		self.AI = AI

		return [matchI,matchF]

	def MSA(self, AlignNumI, *arg):
		SeqS0F = self.SeqS0.seqonly()
		SeqS0R = self.SeqS0.reversecomplement()
		SeqS1F = self.SeqS1.seqonly()
		SeqS1R = self.SeqS1.reversecomplement()
		AlignTagI = 0
		if len(arg) == 1:
			AlignTagI = arg[0]
		MSAS = ''
		if AlignNumI == 0:
			MSAS = musclecall(">S0_FWD\n"+SeqS0F+"\n>S1_FWD\n"+SeqS1F,AlignTagI)
		elif AlignNumI == 1:
			MSAS = musclecall(">S1_REV\n"+SeqS1R+"\n>S0_REV\n"+SeqS0R,AlignTagI)

		elif AlignNumI == 2:
			MSAS = musclecall(">S0_FWD\n"+SeqS0F+"\n>S1_REV\n"+SeqS1R,AlignTagI)
		elif AlignNumI == 3:
			MSAS = musclecall(">S1_FWD\n"+SeqS1F+"\n>S0_REV\n"+SeqS0R,AlignTagI)
		
		elif AlignNumI == 4:
			MSAS = musclecall(">S1_FWD\n"+SeqS1F+"\n>S0_FWD\n"+SeqS0F,AlignTagI)
		elif AlignNumI == 5:
			MSAS = musclecall(">S0_REV\n"+SeqS0R+"\n>S1_REV\n"+SeqS1R,AlignTagI)

		elif AlignNumI == 6:
			MSAS = musclecall(">S0_REV\n"+SeqS0R+"\n>S1_FWD\n"+SeqS1F,AlignTagI)
		elif AlignNumI == 7:
			MSAS = musclecall(">S1_REV\n"+SeqS1R+"\n>S0_FWD\n"+SeqS0F,AlignTagI)

		if len(arg) == 1 and arg[0] == 1:
			#STDERR("##########MSA=",MSAS) ##DEBUG
			MSASL = [">"+x for x in MSAS.split(">") if len(x) > 3]
			FirstSeq = [x for x in MSASL if x.find("S0_") != -1][0]
			SecondSeq = [x for x in MSASL if x.find("S1_") != -1][0]
			MSAS = FirstSeq + SecondSeq
			MSAS = MSAS.replace("S0",self.SeqS0.name()).replace("S1",self.SeqS1.name()).replace("_FWD","").replace("_REV","")
		return MSAS

def main():
	usage = "usage: python OverlapCheck.py -d DBpath -1 1st_SeqID -2 2nd_SeqID"
	opt = optparse.OptionParser(usage)
	opt.add_option("-d",help="sequences database path")
	opt.add_option("-1",help="1st sequnce ID",dest="a")
	opt.add_option("-2",help="2nd sequnce ID",dest="b")
	opt.add_option("-A","--align",help="output alignment in specific orientation",dest="AlignNumI",default="-1")
	opt.add_option("-t","--threads",help="cpu number to use",dest='t',default="1")
	opt.add_option("--consensus",help="cpu number to use",dest='consensus',default="F")	
	(options, args) = opt.parse_args()

	DBnameS, SeqS0_IDS, SeqS1_IDS, cpuNumI, AlignNumI, consensusS = options.d, options.a, options.b, int(options.t), int(options.AlignNumI), options.consensus
	if None in [DBnameS, SeqS0_IDS, SeqS1_IDS]:
		print("missing variable")
		print(usage) 
	else:
		
		SeqS0 = getSseq(DBnameS,SeqS0_IDS,'1')
		SeqS1 = getSseq(DBnameS,SeqS1_IDS,'1')
		print("#",DBnameS, SeqS0_IDS, SeqS1_IDS)
		ovelapL = overlap(SeqS0,SeqS1,cpuNumI)
		print("#","len(SeqS0)="+ str(len(SeqS0)))
		print("#","len(SeqS1)="+ str(len(SeqS1)))
		print("#",ovelapL.OverlapL)
		print("#",ovelapL.OverlapI)
		print("#",ovelapL.OverlapF)
		if ovelapL.AI >= 0 and AlignNumI == -1:
			if consensusS == 'T':
				print(consensusExtractKW(ovelapL.MSA(ovelapL.AI),0,0,0))
			else:
				print(ovelapL.MSA(ovelapL.AI))
			
 

	if AlignNumI > -1:
		print(ovelapL.MSA(AlignNumI))

	elif AlignNumI == -2:
		for i in (0,2,4,6):
			print("AlignNumI= "+str(i))
			print(ovelapL.MSA(i))
		

if __name__ == "__main__":
	main()

