#!/usr/bin/python
##TXblastX.py.py
##Takes input from blastx with -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qcovs' option
import OverlapCheck as OC
from SeqsMerger import SeqsMerger
from FNpool import ConcurrentCall, MSACheck as MCK, FastaTool, BlastHitGroup 
import sys,optparse, copy, subprocess, re
from concurrent.futures import ProcessPoolExecutor
from time import sleep
### Class And Function
def STDERR(*StrInS): #Function For Debugging
	outS = ' '.join([str(x) for x in StrInS])
	sys.stderr.write(str(outS)+'\n')

def ConcurrentCall0(Fnc, InputL, threadsI):
	QueD = {}
	pool = ProcessPoolExecutor(threadsI)
	numL = [x for x in range(len(InputL))]
	OutputL = []
	for i in numL:
		QueD["Q"+str(i)] = pool.submit(Fnc,*InputL[i])
		#STDERR("*InputL[i]=",*InputL[i])

	RoundI = 0
	while True:# and RoundI < waitI:
		xQueL =	[QueD["Q"+str(i)].done() for i in numL]
		STDERR(str(xQueL))
		if False not in xQueL:
			#OutputL.append(xQueL)
			break

		else:
			
			sleep(0.01)

		RoundI += 1

	OutputL = [QueD["Q"+str(i)].result() for i in numL if QueD["Q"+str(i)].done() == True]

	return OutputL
	

def chop(blastxS): #divide each hit into a list of attributes 
	if __TAG__ != '0':
		blastxL = [y for y in blastxS.splitlines() if y.find(__TAG__) != -1 and len(y) > 2]
	else:
		blastxL = [y for y in blastxS.splitlines() if len(y) > 2]		
	return blastxL

class hit(object):
	def __init__(self,blastxhitS):
		#self.AttributesL = [[y[0:2] + [float(y[2])] + [int(x) for x in y[3:10]] + [float(x) for x in y[10:12]] + [int(x) for x in y[12:]] for y in blastxhitS.split("\t")] ##split and cast type for each hit
		self.AttributesL = [y for y in blastxhitS.split("\t")] ##split and cast type for each hit
		##attributes qseqidS sseqidS pidentF lengthI mismatchI gapopenI qstartI qendI sstartI sendI evalueF bitscoreF slenI qlenI qcovsF
		self.qseqidS = self.AttributesL[0]	
		self.sseqidS = self.AttributesL[1]	
		self.pidentF = float(self.AttributesL[2])
		self.lengthI = int(self.AttributesL[3])
		self.mismatchI = int(self.AttributesL[4])
		self.gapopenI = int(self.AttributesL[5])
		self.qstartI = int(self.AttributesL[6])
		self.qendI = int(self.AttributesL[7])
		self.sstartI = int(self.AttributesL[8])
		self.sendI = int(self.AttributesL[9])
		self.evalueF = float(self.AttributesL[10])
		self.bitscoreF = float(self.AttributesL[11])
		self.slenI = int(self.AttributesL[12])
		self.qlenI = int(self.AttributesL[13])
		self.qcovsF = float(abs(self.qendI - self.qstartI))/float(self.qlenI)
		self.scovsF = float(abs(self.sendI - self.sstartI))/float(self.slenI)

def Pgroup(selectOJBL):
	nameL = [] ##create empty list to contain contigs names
	for OBJ in selectOJBL: 
		if OBJ.sseqidS not in nameL:
			nameL.append(OBJ.sseqidS)

	groupedL = [] ##create empty list to contain grouped hit
	while len(nameL) > 0:
		currentNameS = nameL.pop(0)
		#BestHitOBJ = sorted([x for x in selectOJBL if x.sseqidS == currentNameS], key=lambda LL:LL.bitscoreF)[::-1][0] ##select best hit for each 	
		OtherHitL = [x for x in selectOJBL if x.sseqidS == currentNameS]
		ExcNameL = [x.qseqidS for x in OtherHitL]
		#nameL = [x for x in nameL if x not in ExcNameL]
		STDERR("len(nameL)=",len(nameL))
		groupedL.append(OtherHitL)

	return groupedL


def merge5(Ori_inputlistOBJ,coverage):
	inputlistOBJ = copy.deepcopy(Ori_inputlistOBJ) #create new list of input blast results
	a = sorted(inputlistOBJ, key = lambda x:abs(x.sstartI - x.sendI) )[::-1] #sort by matched length 
	clustL = []

	while len(a) > 0:
	
		mark = a.pop(0)
		#clustL.append([])

		Mhead = min(mark.sstartI,mark.sendI)
		Mtail = max(mark.sstartI,mark.sendI)
		Mlength = abs(mark.sendI - mark.sstartI)
		
		## <------------------->
		##    <------------->
		temp = [x for x in a if min(x.sstartI,x.sendI) >= Mhead and max(x.sendI, x.sstartI) <= Mtail] #and abs(x.sendI - x.sstartI) > (coverage* Mlength)]

		##    <------------------->
		## <------------------>
		temp2 = [x for x in a if max(x.sendI, x.sstartI) < Mtail and min(x.sendI, x.sstartI) <= Mhead and ( max(x.sendI, x.sstartI) - Mhead ) >= coverage*(abs(x.sendI - x.sstartI))] 
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		
		## <------------------->
		##     <------------------>		
		temp2 = [x for x in a if min(x.sendI, x.sstartI) > Mhead and max(x.sendI, x.sstartI) > Mtail and (Mtail - min(x.sendI,x.sstartI)) >= coverage*(abs(x.sendI - x.sstartI))] 
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		temp = sorted(temp, key=lambda temp:temp.sstartI)
		
		a = [x for x in a if x not in temp]
		pick = [x for x in temp]
		pick.insert(0,mark) 
		clustL.append(pick) 
	del(inputlistOBJ)
	del(a)
	return clustL


def GroupScore(Eachgrouped1L): ## Scoring for each protein ID by match coverage calculated from protein regions matched by query sequences.
	groupedL =  merge5(Eachgrouped1L,0.75) ## Many query sequences may match to the same region. They will be merged before coverage calculation
	EachgroupedL = [] 
	for L in groupedL: ## for every merged groups
		CurrentL = sorted(L, key=lambda x:x.bitscoreF)[-1] ## Sort hitOBJs in each merged group by bitscore
		EachgroupedL.append(CurrentL) ## Select hits_OBJ with highest bitscore to represent the group match position
	
	OBJLsortedL = sorted(EachgroupedL, key=lambda x:min(x.sstartI, x.sendI) ) ## Sort by match position on SUBJECT sequence 

	HeadI = 0
	TailI = 0
	SLenI = OBJLsortedL[0].slenI

	CovNumI = 0
	CurrentTailI = 0
	for i in OBJLsortedL[:]:
		if i.sstartI < i.sendI:
			CurrentHeadI = i.sstartI
			CurrentTailI = i.sendI
		elif i.sstartI > i.sendI:
			CurrentHeadI = i.sendI
			CurrentTailI = i.sstartI

		if CurrentHeadI > TailI:
			CovNumI = CovNumI + abs(TailI - HeadI)
			HeadI = CurrentHeadI
			TailI = CurrentTailI
		else:
			
			TailI = CurrentTailI
	
	if TailI == CurrentTailI:
		CovNumI = CovNumI + abs(TailI - HeadI)

	return [[x.AttributesL for x in EachgroupedL],float(CovNumI)/float(SLenI)]
	

def group(selectOJBL):
	nameL = [] ##create empty list to contain contigs names
	for OBJ in selectOJBL: 
		if OBJ.qseqidS not in nameL:
			nameL.append(OBJ.qseqidS)

	groupedL = [] ##create empty list to contain grouped hit
	while len(nameL) > 0:
		currentNameS = nameL.pop(0)
		BestHitOBJ = sorted([x for x in selectOJBL if x.qseqidS == currentNameS], key=lambda LL:LL.bitscoreF)[::-1][0] ##select best hit for each 	
		OtherHitL = [x for x in selectOJBL if x.sseqidS == BestHitOBJ.sseqidS]
		ExcNameL = [x.qseqidS for x in OtherHitL]
		nameL = [x for x in nameL if x not in ExcNameL]
		STDERR("len(nameL)=",len(nameL))
		groupedL.append(OtherHitL)

	return groupedL
		
def TruncateCheck(AblastResultOBJ,EndNumI): ##filter_out blast hit, math by sub-sequence at middle of query.
	x = AblastResultOBJ	
	#if x.qstartI > min(x.qlenI*0.1, 10) and  x.qendI < max(x.qlenI*0.9, (x.qlenI - 10)): ##Truncate hit condidion
	if min(x.qstartI, x.qendI) > EndNumI and  max(x.qstartI, x.qendI) < (x.qlenI - EndNumI): ##Truncate hit condidion
		FlagI = 0
	else:
		FlagI = 1
	return FlagI

def QuerySeq(hitOBJ,__DBname__):
	QnameS = hitOBJ.qseqidS.replace("lcl|","").split()[0]
	DBnameS = __DBname__
	QueryStrandS = "plus"
	if hitOBJ.qstartI > hitOBJ.qendI:
			QueryStrandS = "minus"
	arg = "blastdbcmd -db "+DBnameS+" -strand " + QueryStrandS + " -entry \""+QnameS+"\";"
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()	
	TEMPseqs = x[0].strip() #.decode('utf-8')
	#STDERR("TEMPseqs=",TEMPseqs)
	del process
	return	TEMPseqs

def musclecallCLW(SeqS):
	#SeqS = self.SeqS
	#arg = "muscle -maxiters 32 -gapopen -1200 -quiet"
	ByteSeq = SeqS.encode('utf-8')
	#STDERR("musclecallCLW IS RUNNING")
	arg = "muscle -maxiters 32 -quiet -clw"
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(ByteSeq)[0]
	#STDERR("musclecallCLW out = ",x)
	if len(x) == 0:
		x = '0'.encode('utf-8')
	#STDERR("musclecallCLW IS RUNNING")
	return x.decode("utf-8")

def musclecall(SeqS):
	#SeqS = self.SeqS
	#arg = "muscle -maxiters 32 -gapopen -1200 -quiet"
	ByteSeq = SeqS.encode('utf-8')
	#STDERR("musclecallCLW IS RUNNING")
	arg = "muscle -maxiters 32 -quiet"
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(ByteSeq)[0]
	#STDERR("musclecallCLW out = ",x)
	if len(x) == 0:
		x = '0'.encode('utf-8')
	#STDERR("musclecallCLW IS RUNNING")
	return x.decode("utf-8")

class FastaTool0(object): ## class for manipulate sequences in fasta format
	def __init__(self,FastaInS):
		self.FastaInS = FastaInS
		self.seqLenI = self.seqlen()
		self.FastaSrev  = self.header() +"Rev\n" + self.reversecomplement()

	def seqonly(self):
		fastaL = self.FastaInS.splitlines()
		#fastaL = [x for x in fastaL if x.find(">") != 0]
		fastaL = [x for x in fastaL if x.find(">") == -1]
		fastaL = [x.strip() for x in fastaL]
		self.fastaS = ''.join(fastaL)
		return self.fastaS
	def FullSeqS(self):
		return self.FastaInS

	def SubSeqS(self,StartI,EndI):
		
		if StartI < EndI:
			SubSeqS = self.seqonly()[StartI:EndI]
			headerS = self.header() + ":"+str(StartI)+"-"+str(EndI)
			return headerS+"\n"+SubSeqS+"\n"
						
		elif StartI > EndI:
			SubSeqS = self.seqonly()[max(EndI-1,0):StartI]
			headerS = self.header() + ":"+str(EndI)+"-"+str(StartI)
			return headerS+"\n"+SubSeqS+"\n"

	def header(self):
		fastaL = self.FastaInS.splitlines()
		fastaL = [x for x in fastaL if x.find(">") != -1]
		fastaL = [x.strip() for x in fastaL]
		self.headerS = fastaL[0]
		return self.headerS

	def name(self):
		self.NameS = self.header().split(' ')[0].replace('>','').replace('lcl|','')
		return self.NameS

	def ManyLines(self,widthI):
		seqS = self.seqonly()
		seqLineI = len(seqS)/widthI
		remainsI = len(seqS)%widthI
		TEMPseqS = self.header()+"\n"
		for i in range(seqLineI):
			TEMPseqS = TEMPseqS + seqS[i*(widthI):(i+1)*(widthI)] + '\n'

		TEMPseqS = TEMPseqS + seqS[(-1)*(remainsI):] + '\n'
		return TEMPseqS


	def reversecomplement(self):
		string = self.seqonly()
		#header = self.header() + "Rev"
		reverse = string[::-1]
		reverse = reverse.replace("A","?")
		reverse = reverse.replace("T","A")
		reverse = reverse.replace("?","T")
		reverse = reverse.replace("G","?")
		reverse = reverse.replace("C","G")
		self.reverseS = reverse.replace("?","C")
		#print self.reverseS
		return self.reverseS ## complementary string ATCG --> CGAT

	def seqlen(self):
		
		self.seqLenI = len(self.seqonly())
		return self.seqLenI

def OverlapHit(HitOBJ1,HitOBJ2):#get two blast hit object return 0 if match position on subject not overlap else return 1
	#1 --------------        #1    ----------------
	#2    ------------------ #2  ----------------
	flagI = 0
	if min(HitOBJ1.sstartI,HitOBJ1.sendI) <= min(HitOBJ2.sstartI,HitOBJ2.sendI) and max(HitOBJ1.sstartI,HitOBJ1.sendI) > min(HitOBJ2.sstartI,HitOBJ2.sendI):
		flagI = 1
	elif min(HitOBJ1.sstartI,HitOBJ1.sendI) >= min(HitOBJ2.sstartI,HitOBJ2.sendI) and min(HitOBJ1.sstartI,HitOBJ1.sendI) < max(HitOBJ2.sstartI,HitOBJ2.sendI):
		flagI = 1
	return flagI

def AlignScore(SeqsL):

	scoredL = [] ## Empty list for score storing
	alignLenI = len(SeqsL[0][1])
	for seqS in SeqsL:
		AllGapI = seqS[1].count("-")
		HeadGapI = alignLenI - len(seqS[1].lstrip("-")) #count lead gap in alignment
		TailGapI = alignLenI - len(seqS[1].rstrip("-")) #count tail gap in alignment
		InnerGapI = AllGapI - (HeadGapI + TailGapI) #count inner gap in alignment
		InnerGapOpenI = len([x for x in seqS[1].strip("-").split("-") if len(x) > 0]) - 1 #count inner gap opening in alignment
		GapCountL = [seqS, AllGapI, HeadGapI, TailGapI, InnerGapI, InnerGapOpenI]
		scoredL.append(GapCountL)
	
	return 	scoredL

def MSAfastaSplit(MSAfasta):
	MSAfastaL = MSAfasta.split('>') ## Take input aligment as fasta format
	MSAfastaL = ['>' + x for x in MSAfastaL if len(x) > 0] ##Eliminate empty items
	MSAfastaL2 = []
	for x in MSAfastaL:## Separate each sequences
		xL = x.splitlines()
		MSAfastaL2.append( [ xL[0],''.join(xL[1:]) ] )

	return MSAfastaL2 ##Return 2D List [[NameS,seqS],[NameS,seqS]]

def MSACheck(MSAS): #Check consistency of MSA 
	MSAfastaL = [x for x in MSAfastaSplit(MSAS)]
	scoredL = AlignScore(MSAfastaL)
	#STDERR("MSACheck.scoredL",scoredL)
	return scoredL

def BlastHitJoiner3(HitL): ## Takes Input as Blast Hit, Return Joint Sequences and Names List of Each Group
	BigClustL = BlastHitGroup(HitL,0.8) ## Grouping blasthit
	BigOutL = [] 
	for ClustL in BigClustL:
		inL = [[x ,__DBname__] for x in ClustL] ##Prepare Input
		SeqSL = [x.decode("utf-8") for x in ConcurrentCall(QuerySeq, inL, 10)] ## Pull sequences from DB
		MergedL = SeqsMerger(SeqSL)  ## Merged Sequences
		BigOutL.append(MergedL)

	return(BigOutL)

def BlastHitJoiner2(HitL):
	SortedHitL = sorted(HitL, key = lambda x:x.bitscoreF)[::-1] ##sort by hit-bitscore
	GroupL = []
	UsedNameL = []
	NumI = 0
	LimitI = len(SortedHitL)
	MergedL = []
	while len(SortedHitL) > 0 and NumI < LimitI+1:
		CL = SortedHitL.pop(0) ##select hit with highest bit-score to use as anchor
		STDERR("ANCHOR=",CL.qseqidS.replace("lcl|",""))
		tmpL = [x for x in SortedHitL if OverlapHit(CL,x) > 0 and x.qseqidS.replace("lcl|","") not in UsedNameL]
		inL = [[x ,__DBname__] for x in tmpL]
		ToALignSL = [x.decode("utf-8") for x in ConcurrentCall(QuerySeq, inL, 10)]
		#STDERR("BlastHitJoiner.ToALignSL",ToALignSL)
		CLseqS = QuerySeq(CL,__DBname__).decode("utf-8")
		#STDERR("BlastHitJoiner.CLseqS",CLseqS)
		#ConInL = [[CLseqS , x, 8] for x in ToALignSL]
		ConInL = [[CLseqS , x] for x in ToALignSL]
		#STDERR("BlastHitJoiner.ConInL",ConInL)
		ChkL = [x for x in ConcurrentCall(OC.overlap, ConInL, 8)]
		PassedL = [x for x in ChkL if x.OverlapF > 0.2]
		#STDERR("BlasthitJoiner2.PassedL",[[x.SeqS0.name(),x.SeqS1.name(),x.OverlapI,x.OverlapF,x.AI] for x in PassedL]) ##DEBUG
		MSAInL = [x.MSA(x.AI,1) for x in PassedL]
		#STDERR("BlasthitJoiner2.MSAInL",MSAInL) ##DEBUG		
		ChkL2 = [x for x in ConcurrentCall(MCK, [[x] for x in MSAInL], 10)] ##Check alignment wether it good or not
		NumI += 1 ##endless while loop proof
		
	
		#STDERR("BlastHitJoiner.ChkL2",str(ChkL2).replace("]],","\n")) ##DEBUG
		#STDERR("BlastHitJoiner.ChkL2[0]",ChkL2[0][0].name,ChkL2[0].FragNumI,ChkL2[0].LargeRatioF)  ##DEBUG
		#break##DEBUG
		UsedNameL.append(CL.qseqidS.replace("lcl|",""))
		GroupL.append([CL.qseqidS.replace("lcl|","")])
		#STDERR("BlastHitJoiner2.len(ChkL2)=",len(ChkL2))##DEBUG
		TMPmergedSeqS = ''
		if len(ChkL2) > 0:
			TMPmergedSeqS = ">" + ChkL2[0][0].nameS +"\n"+ ChkL2[0][0].SeqS.replace("-","")
			#STDERR("BlastHitJoiner2.len(TMPmergedSeqS)=",len(TMPmergedSeqS))##DEBUG
			#STDERR("BlastHitJoiner2.TMPmergedSeqS ORIGINAL =",TMPmergedSeqS)##DEBUG
		else:
			TMPmergedSeqS = CLseqS
		for i in ChkL2: ##DEBUG
			#STDERR("############################################################ i ########################################################") ##DEBUG
			#STDERR(i)
			Seq0OBJ = i[0]
			Seq1OBJ = i[1]
			STDERR("Seq0OBJ.name=",Seq0OBJ.nameS,Seq0OBJ.FragNumI,Seq0OBJ.LargeRatioF,Seq0OBJ.SmallRatioF,Seq0OBJ.HeadGapI,Seq0OBJ.TailGapI,Seq0OBJ.GapRatioF,Seq0OBJ.AllGapI,Seq0OBJ.InnerGapPerSeqLenF) ##DEBUG
			STDERR("Seq1OBJ.name=",Seq1OBJ.nameS,Seq1OBJ.FragNumI,Seq1OBJ.LargeRatioF,Seq1OBJ.SmallRatioF,Seq1OBJ.HeadGapI,Seq1OBJ.TailGapI,Seq1OBJ.GapRatioF, Seq1OBJ.AllGapI, Seq1OBJ.InnerGapPerSeqLenF) ##DEBUG
			#if Seq0OBJ.FragNumI < 3 and Seq1OBJ.FragNumI < 3 and Seq0OBJ.GapRatioF < 0.1 and Seq1OBJ.GapRatioF < 0.1:
			#if Seq1OBJ.nameS not in UsedNameL and  Seq0OBJ.LargeRatioF > 0.9 and Seq1OBJ.LargeRatioF > 0.9 :
			if Seq1OBJ.nameS not in UsedNameL and ( Seq0OBJ.InnerGapPerSeqLenF < 0.1 and Seq1OBJ.InnerGapPerSeqLenF < 0.1  ) and ( ( Seq0OBJ.LargeRatioF > 0.25 and Seq0OBJ.N50NumI < 3 and Seq1OBJ.LargeRatioF > 0.25 and Seq1OBJ.N50NumI < 3 ) or ( Seq0OBJ.GapRatioF < __GAPLIMITF__  and Seq1OBJ.GapRatioF < __GAPLIMITF__ ) ):
				
									
				NextSeqS = ">" + Seq1OBJ.nameS +"\n"+ Seq1OBJ.SeqS.replace("-","")
				#STDERR("BlastHitJoiner2.FirstSeqS=",TMPmergedSeqS)##DEBUG
				#STDERR("BlastHitJoiner2.NextSeqS=",NextSeqS)##DEBUG
				TMPoverlap = OC.overlap(TMPmergedSeqS, NextSeqS)
				STDERR("BlastHitJoiner2.TMPoverlap=",TMPoverlap.OverlapL)##DEBUG
				ConsenseusSOBJ = FastaTool(TMPoverlap.conservedblock())
				if ConsenseusSOBJ.seqonly().count("N")/ConsenseusSOBJ.seqlen() < __NlimitF__:
					TMPmergedSeqS = TMPoverlap.merge()
					STDERR("BlastHitJoiner2.len(TMPmergedSeqS)=",len(TMPmergedSeqS))##DEBUG
					#STDERR("BlastHitJoiner2.TMPmergedSeqS=",TMPmergedSeqS)##DEBUG
					SeqNameS = Seq1OBJ.nameS
					if SeqNameS not in UsedNameL:
						GroupL[-1].append(SeqNameS)
						#TMPusednameL.append(SeqNameS)
						UsedNameL.append(SeqNameS)

			#STDERR("i[1][1:]=",i[1][1:]) ##DEBUG
		SortedHitL = [x for x in SortedHitL if x.qseqidS.replace("lcl|","") not in UsedNameL]
		SortedHitL = sorted( SortedHitL,  key = lambda x:x.bitscoreF)[::-1] ##sort by hit-bitscore
		MergedL.append(TMPmergedSeqS)
		#STDERR("BlastHitJoiner.TMPmergedSeqS",TMPmergedSeqS)

	STDERR("BlastHitJoiner.UsedNameL",UsedNameL)
	STDERR("BlastHitJoiner.GroupL",GroupL)
	return [GroupL, MergedL]

def FinalScore(OBJL,SeqS):
	SeqNameSL = [x for x in SeqS.splitlines()[0].replace(">","").replace("lcl|","").split()[0].split("_:") if x != "merged"  ]
	STDERR("FinalScore.SeqNameSL",SeqNameSL)
	sortedOBJL = sorted([x for x in OBJL if x.qseqidS.replace("lcl|","") in SeqNameSL], key=lambda x:x.bitscoreF )[::-1]
	TMPL = []
	UseNameSL = []
	for i in sortedOBJL:
		if i.qseqidS not in UseNameSL:
			TMPL.append(i)

	SumBitF = sum([x.bitscoreF for x in TMPL])
	SlenI = TMPL[0].slenI
	MemberNumI = len(TMPL)
	StartIOBJ = sorted(TMPL, key=lambda x:min(x.sstartI,x.sendI))[0]
	StartI = min(StartIOBJ.sstartI,StartIOBJ.sendI)
	EndIOBJ = sorted(TMPL, key=lambda x:max(x.sstartI,x.sendI))[-1]
	EndI =  max(EndIOBJ.sstartI, EndIOBJ.sendI)
	SumCovF = float(EndI - StartI)
	EvCovF = SumCovF/float(SlenI)
	EvBitF = sum([x.bitscoreF/float(x.slenI) for x in TMPL])/float(MemberNumI)
	EvBitLenF = EvBitF*EvCovF ##Everage bitscore multiply by Everage coverage
	return [StartI,EndI,SumCovF, EvCovF, SumBitF, EvBitF, EvBitLenF]


def main2(INS,TresholdF):
	INL = (chop(INS))
	hitOBJL = [hit(x) for x in INL]
	#hitOBJL = [x for x in hitOBJL if TruncateCheck(x,10) == 1] ##remove truncated hit before processing
	selectOJBL = [x for x in hitOBJL if x.qcovsF > __QCovF__ and x.evalueF < __EvalF__]

	#selectOJBL = selectOJBL[0:1000] ##DEBUG
	#STDERR([x.AttributesL for x in selectOJBL[0:2]]) ##DEBUG
	#STDERR("entering group function") ##DEBUG
	outS = ''
	#for OBJL in group(selectOJBL):
	for OBJL in Pgroup(selectOJBL):
		ScoreFL = GroupScore(OBJL)
		CrPnameS = OBJL[0].sseqidS
		STDERR("Subject ID =",CrPnameS)
		outS = outS + "\n#Grouping Sequences For "+ CrPnameS + "\n"
		
		sortedOBJL = sorted( OBJL, key=lambda x:( (x.sstartI + x.sendI)/2 ,min(x.sstartI, x.sendI) ) )
		TXnumS = str(len(ScoreFL[0]))
		outS = outS + "#" + '\n#'.join(['\t'.join((x.AttributesL)) for x in sortedOBJL]) + "\n#PoolCovScpre=\t" + str(ScoreFL[1]) +"\t"+ ScoreFL[0][0][1] + "\tTXnumS=\t"+TXnumS +  "\n############\n"		
		outS = outS + "#" + str(ScoreFL[0]).replace("""],""","""]\n""") + "\n############\n"
		#STDERR("TEMPSEQ") ##DEBUG
		#STDERR(''.join(MSAL)) ##DEBUG
		if ScoreFL[1] > TresholdF:
			#MergedGroupL = BlastHitJoiner2(OBJL)##DEBUG
			MergedGroupL = BlastHitJoiner3(OBJL)
			for g in MergedGroupL: ##DEBUG
				for gg in g:
					STDERR("main2.MergedGroupL.subGroup.eachGroup=",gg) ##DEBUG
			#GroupL = MergedGroupL[0]
			GroupL = [x[1] for x in MergedGroupL]
			#MergedSeqSL = sorted([x for x in MergedGroupL[1] if len(x) > 0], key=lambda x:len(x))[::-1]
			MergedSeqSL = sorted([x[0] for x in MergedGroupL if len(x[0]) > 0], key=lambda x:len(x))[::-1]
			#STDERR("main2.MergedSeqSL",MergedSeqSL) ##DEBUG
			outS = outS + "\n############# Align hit position \n"
			STDERR("main2.GroupL=",GroupL) ##DEBUG
			for i in GroupL:
				TMPobjL  = []
				STDERR("main2.GroupL.i=",i) ##DEBUG
				for ii in i:
					TMPobjL.append([x for x in sortedOBJL if x.qseqidS.replace("lcl|","") == ii][0])
				
				TMPstartI = min([min(x.sstartI, x.sendI) for x in TMPobjL])
				TMPendI = max([max(x.sstartI, x.sendI) for x in TMPobjL])
				outS = outS + "\n#"+ str(TMPstartI) + "-" + str(TMPendI) + "\t" + ' '.join(i)

			if __MergedFlagS__ != "0":

				LV2MergedSeqSL = SeqsMerger(MergedSeqSL)
				STDERR("main2.LV2MergedSeqSL[0]",LV2MergedSeqSL[0])##DEBUG

				MergedSeqsS = '\n'.join([x.splitlines()[0].split()[0].replace("lcl|","") +" "+ CrPnameS +" " + str(FinalScore(OBJL,x)) +'\n'+ '\n'.join(x.splitlines()[1:]) for x in LV2MergedSeqSL[1]]) + "\n"
				STDERR("main2.LV2MergedSeqSL",LV2MergedSeqSL[0])##DEBUG

			
				outS = outS + "\n###MERGED SEQ for "+ CrPnameS +"##\n"	+ MergedSeqsS

		else:
			
			outS = outS + "\n############# Overall Coverage not pass cut-off \n"
		
		if __OutTag__ == '0':
			print(outS)
			outS = ''
	#STDERR(group(selectOJBL)[0])
	return outS

def main3(INS,TresholdF):
	INL = (chop(INS))
	hitOBJL = [hit(x) for x in INL]
	#hitOBJL = [x for x in hitOBJL if TruncateCheck(x,10) == 1] ##remove truncated hit before processing
	selectOJBL = [x for x in hitOBJL if x.qcovsF > __QCovF__ and x.evalueF < __EvalF__]

	#selectOJBL = selectOJBL[0:1000] ##DEBUG
	#STDERR([x.AttributesL for x in selectOJBL[0:2]]) ##DEBUG
	#STDERR("entering group function") ##DEBUG
	outS = ''
	#for OBJL in group(selectOJBL):
	for OBJL in Pgroup(selectOJBL):
		ScoreFL = GroupScore(OBJL)
		CrPnameS = OBJL[0].sseqidS
		STDERR("Subject ID =",CrPnameS)
		outS = outS + "\n#Grouping Sequences For "+ CrPnameS + "\n"
		
		sortedOBJL = sorted( OBJL, key=lambda x:( (x.sstartI + x.sendI)/2 ,min(x.sstartI, x.sendI) ) )
		TXnumS = str(len(ScoreFL[0]))
		outS = outS + "#" + '\n#'.join(['\t'.join((x.AttributesL)) for x in sortedOBJL]) + "\n#PoolCovScpre=\t" + str(ScoreFL[1]) +"\t"+ ScoreFL[0][0][1] + "\tTXnumS=\t"+TXnumS +  "\n############\n"		
		outS = outS + "#" + str(ScoreFL[0]).replace("""],""","""]\n""") + "\n############\n"
		#STDERR("TEMPSEQ") ##DEBUG
		#STDERR(''.join(MSAL)) ##DEBUG
		if ScoreFL[1] > TresholdF:
			#MergedGroupL = BlastHitJoiner2(OBJL)##DEBUG
			MergedGroupL = BlastHitJoiner3(OBJL)
			MergedGroupL2D = []
			for g in MergedGroupL: 
				STDERR("main3.MergedGroupL.subGroup.g=",g) ##DEBUG
				for gg in g:
					STDERR("main3.MergedGroupL.subGroup.gg=",gg) ##DEBUG
					MergedGroupL2D.append(gg)
			#GroupL = MergedGroupL[0]
			GroupL = [x[1] for x in MergedGroupL2D]
			#MergedSeqSL = sorted([x for x in MergedGroupL[1] if len(x) > 0], key=lambda x:len(x))[::-1]
			MergedSeqSL = sorted([x[0] for x in MergedGroupL2D if len(x[0]) > 0], key=lambda x:len(x))[::-1]
			#STDERR("main2.MergedSeqSL",MergedSeqSL) ##DEBUG
			outS = outS + "\n############# Align hit position \n"
			#STDERR("main3.GroupL=",GroupL) ##DEBUG
			for i in GroupL:
				TMPobjL  = []
				#STDERR("main3.GroupL.i=",i) ##DEBUG
				for ii in i:
					TMPobjL.append([x for x in sortedOBJL if x.qseqidS.replace("lcl|","") == ii][0])
				
				TMPstartI = min([min(x.sstartI, x.sendI) for x in TMPobjL])
				TMPendI = max([max(x.sstartI, x.sendI) for x in TMPobjL])
				outS = outS + "\n#"+ str(TMPstartI) + "-" + str(TMPendI) + "\t" + ' '.join(i)

			if __MergedFlagS__ != "0":

				LV2MergedSeqSL = SeqsMerger(MergedSeqSL)
				LV3MergedSeqSL = [x[0] for x in LV2MergedSeqSL]
				#STDERR("main2.LV2MergedSeqSL",LV2MergedSeqSL)##DEBUG
				

				MergedSeqsS = '\n'.join([x.splitlines()[0].split()[0].replace("lcl|","") +" "+ CrPnameS +" " + str(FinalScore(OBJL,x)) +'\n'+ '\n'.join(x.splitlines()[1:]) for x in LV3MergedSeqSL]) + "\n"
				#STDERR("main2.LV2MergedSeqSL",LV2MergedSeqSL[0])##DEBUG

			
				outS = outS + "\n###MERGED SEQ for "+ CrPnameS +"##\n"	+ MergedSeqsS

		else:
			
			outS = outS + "\n############# Overall Coverage not pass cut-off \n"
		
		if __OutTag__ == '0':
			print(outS)
			outS = ''
	#STDERR(group(selectOJBL)[0])
	return outS




### Input and Option 
usage = "python TXblastX.py.py -i input.gff -o out_file"
opt = optparse.OptionParser(usage)
opt.add_option("-i",help="*input path, get input from stdin if ommit", default='0')
opt.add_option("-o",help="indicate output file name or print out as standard output",default="0")
opt.add_option("-e",help="E-value cutoff for selected hit",default="0.1")
opt.add_option("-c",help="query-length covery rate cutoff for selected hit",default="0.98")
opt.add_option("-t",help="cutoff for pool protein coverage treshold",default="0.9")
opt.add_option("--MergedFlagS",help="MergedFlagS to output merged sequences",default="1",dest="MergedFlagS")
opt.add_option("--NlimitF",help="Maximum N ratio allowed in Consensus 0.0 to 1.0",default="0.05",dest="NlimitF")
opt.add_option("--GAPLIMITF",help="Maximum GAP ratio allowed in alignment 0.0 to 1.0",default="0.1",dest="GAPLIMITF")
opt.add_option("--TAG",help="Spicies Tag in protein name",default="0",dest="TAG")
opt.add_option("-q",help="*query database path",dest='q',default="DB_PATH")

(options, args) = opt.parse_args()

##Documentation part
__QCovF__ = float(options.c) ;STDERR("__QCovF__ = ",__QCovF__)
__EvalF__ = float(options.e)	;STDERR("__EvalF__ = ",__EvalF__)
__TAG__ = options.TAG ;STDERR("__TAG__ = ",__TAG__)
__DBname__ = options.q ;STDERR("__DBname__ = ",__DBname__)
__TresholdF__ = float(options.t) ;STDERR("__TresholdF__ = ",__TresholdF__)
__NlimitF__ = float(options.NlimitF) ;STDERR("__NlimitF__ = ",__NlimitF__)
__OutTag__ = options.o ;STDERR("__OutTag__",__OutTag__)
__MergedFlagS__ = options.MergedFlagS ;STDERR("__MergedFlagS__ ",__MergedFlagS__ )
__GAPLIMITF__ = float(options.GAPLIMITF) ;STDERR("__GAPLIMITF__ ",__GAPLIMITF__ )

if options.i == '0': ##get input from pipe
	INS = sys.stdin.read()

else:#open file
	f=open(options.i,'r')
	INS = f.read()
	f.close()
if __OutTag__ == '0': ##print output to pipe
	print(main3(INS,__TresholdF__))
else:#write output to a flie
	f=open(options.o,'w')
	f.write(main3(INS,__TresholdF__))
	f.close()	



