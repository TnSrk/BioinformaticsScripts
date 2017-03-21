#!/usr/bin/python
#BlastOut2Consensus.py 2017 Feb 09
##Takes input from blastx with -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qcovs' option
import sys, subprocess, optparse, re, copy
from multiprocessing import Pool, cpu_count

### Class And Function
def STDERR(*StrInS): #Function For Debugging
	outS = ' '.join([str(x) for x in StrInS])
	sys.stderr.write(str(outS)+'\n')

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
		self.qcovsF = float(self.AttributesL[14])

def merge4(Ori_inputlistOBJ,coverage):
	inputlistOBJ = copy.deepcopy(Ori_inputlistOBJ)
	a = sorted(inputlistOBJ, key = lambda x:x.bitscoreF)[::-1] # sort blast hits by bitscore
	clustL = []

	while len(a) > 0:
	
		mark = a.pop(0)
		#clustL.append([])

		Mhead = min(mark.qstartI,mark.qendI)
		Mtail = max(mark.qstartI,mark.qendI)
		Mlength = abs(mark.qendI - mark.qstartI)
		
		## <------------------->
		##    <------------->
		temp = [x for x in a if min(x.qstartI,x.qendI) >= Mhead and max(x.qendI, x.qstartI) <= Mtail] #and abs(x.sendI - x.sstartI) > (coverage* Mlength)]

		##    <------------------->
		## <------------------>
		temp2 = [x for x in a if max(x.qendI, x.qstartI) <= Mtail and max(x.qendI, x.qstartI) >= Mhead] #and (Mhead - min(x.sstartI,x.sendI)) <= coverage*(abs(x.sendI - x.sstartI))] 
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		
		## <------------------->
		##     <------------------>		
		temp2 = [x for x in a if min(x.qendI, x.qstartI) >= Mhead and min(x.qendI, x.qstartI) >= Mtail ] #and (max(x.sendI, x.sstartI) - Mtail) <= coverage*(abs(x.sendI - x.sstartI))]
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		temp = sorted(temp, key=lambda temp:temp.qstartI)
		
		a = [x for x in a if x not in temp]
		pick = [x for x in temp]
		#mark.AttributesL[0] = "*"+mark.qseqidS
		pick.insert(0,mark) 
		clustL.append(pick) 
	del(inputlistOBJ)
	del(a)
	return clustL

def getSseq4(DBname,Scontigname,Shead,Stail):
#	print DBname,Scontigname,Shead,Stail ##DEBUG
	#arg1 = "fastacmd -p F -d "+DBname+" -s \""+Scontigname+"\" -L "+str(Shead)+","+str(Stail)+" ;"
	arg1 = "blastdbcmd -db "+DBname+" -entry \""+Scontigname+"\" -range "+str(Shead)+"-"+str(Stail)+" ;"
	#arg2 = "fastacmd -p F -S 2 -d "+DBname+" -s \""+Scontigname+"\" -L "+str(Stail)+","+str(Shead)+" ;"
	arg2 = "blastdbcmd -strand minus -db "+DBname+" -entry \""+Scontigname+"\" -range "+str(Shead)+"-"+str(Stail)+" ;"
		
	if int(Shead) < int(Stail):		
		arg = arg1
			
	elif int(Shead) > int(Stail):
		arg = arg2

	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	
	TEMPseqs = x[0]
	#sys.stderr.write(str(x[1])) ##DEBUG
	del process
	return	TEMPseqs

def getSeqFull(DBname,Scontigname):
#	print DBname,Scontigname,Shead,Stail ##DEBUG
	#arg1 = "fastacmd -p F -d "+DBname+" -s \""+Scontigname+"\" -L "+str(Shead)+","+str(Stail)+" ;"
	arg = "blastdbcmd -db "+DBname+" -entry \""+Scontigname+"\";"
		
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	
	TEMPseqs = x[0].strip()
	#STDERR("str(x[1])=",str(x[1])) ##DEBUG
	del process
	return	TEMPseqs

class FastaTool(object): ## class for manipulate sequences in fasta format
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

def musclecall(SeqS):
	#SeqS = self.SeqS
	#arg = "muscle -maxiters 32 -gapopen -1200 -quiet"
	arg = "muscle -maxiters 32 -quiet"
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(SeqS)[0]
	if len(x) == 0:
		x = '0'
	return x

def musclecallCLW(SeqS):
	#SeqS = self.SeqS
	#arg = "muscle -maxiters 32 -gapopen -1200 -quiet"
	arg = "muscle -maxiters 32 -quiet -clw"
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(SeqS)[0]
	if len(x) == 0:
		x = '0'
	return x

def AlignCheck(MSAfasta):
	SeqsL = MSAfasta.split('>') ## Take input aligment as fasta format
	SeqsL = ['>' + x for x in SeqsL if len(x) > 0] ##Eliminate empty items
	SeqsL = [FastaTool(x) for x in SeqsL] ## Create fasta object for each sequences
	conS  = ">CoNs\n" +  consensusExtractKW(MSAfasta,0,0.1) + "\n" ##Create consensus from input alignment 
	scoreL = [] ## Empty list for score storing
	for i in SeqSL:
		ToalignS = conS + i.FastaInS + "\n"
		NewAlign = musclecall(ToalignS)
		SeqsL = MSAfasta.split('>') ## Take input aligment as fasta format
		SeqsL = ['>' + x for x in SeqsL if len(x) > 0] ##Eliminate empty items
		NewAlignS = [''.join(x.splitlines()[1:]) for x in SeqsL if len(x) > 0 and x.find("CoNs") != 0][0] ##Eliminate empty items
		alignLenI = len(NewAlignS)
		NewHeadGap = alignLenI - len(NewAlignS.lstrip("-")) #count lead gap in alignment
		NewTailGap = alignLenI - len(NewAlignS.rstrip("-")) #count tail gap in alignment
		NewinnerGap = alignLenI - len(NewAlignS.strip("-")) #count inner gap in alignment
		GapCountL = [conS,NewHeadGap,NewTailGap,NewinnerGap]
		scoreL.append(GapCountL)
				
		del SeqsL

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
	MSAfastaL = [[x.splitlines()[0],''.join(x.splitlines()[1:])] for x in MSAfastaL] ## Separate each sequences
	return MSAfastaL ##Return 2D List [[NameS,seqS],[NameS,seqS]]

def UglyKicker(scoredL):
	EvrGapF = float(sum( [ x[5] for x in scoredL]))/float(len(scoredL))
	EvrInnerGapF = float(sum( [ x[4] for x in scoredL]))/float(len(scoredL))
	ExcludeSeqNameS = sorted(scoredL, key = lambda x:(abs(float(x[5]) - EvrGapF) + abs(float(x[4]) - EvrGapF)) )[-1][0][0].split()[0] ##get name of gap-prone sequence 
	ToAlignS = '\n'.join([x[0]+"\n"+x[1].replace('-','') for x in MSAfastaL if x[0].find(ExcludeSeqNameS) == -1])
	return snatch_scoredL

def SimilarPick(TargetSeqS,RelatedSeqSL,InnerGapLimitI): ## Aling Each match Sequences to target sequences then pick only one with no gap more than threshold
	TagetSeqNameS = TargetSeqS.splitlines()[0].split()[0]
	PickedL = [TargetSeqS]
	for i in RelatedSeqSL:
		ToALignS = TargetSeqS + "\n" + i
		MSAfastaS =  musclecall(ToALignS)
		#STDERR("SimilarPick.MSAfastaSplit(MSAfastaS)=",MSAfastaSplit(MSAfastaS))
		MSAfastaL = [x for x in MSAfastaSplit(MSAfastaS) if x[0].find(TagetSeqNameS) == -1]
		#STDERR("SimilarPick.MSAfastaL=",MSAfastaL)
		scoredL = AlignScore(MSAfastaL)
		for each in scoredL: 
			if each[5] <= InnerGapLimitI:
				PickedL.append(each[0][0]+"\n"+each[0][1].replace("-",""))

	return PickedL ## return [TargetSeqS, PickedSeqS1st, .... , PickedSeqSNth]

def AlignImprove(MSAfasta):
	MSAfastaL = MSAfastaSplit(MSAfasta)
	selectedNameL = []
	iterNUMI = 0
	SeqsL = MSAfastaL[:]
	STDERR("Before WHILE LOOP MSAfastaL",[x[0] for x in SeqsL])
	while len(selectedNameL) == 0 and iterNUMI < 20:
		scoredL = AlignScore(SeqsL)
		EvrGapF = float(sum( [ x[5] for x in scoredL]))/float(len(scoredL))
		EvrInnerGapF = float(sum( [ x[4] for x in scoredL]))/float(len(scoredL))
		#STDERR("scoredL=",scoredL) ##DEBUG
		#alignLenI = len(scoredL[0][0][1])
		kickL = [x[0] for x in scoredL if x[5] > 0 ]#or float(x[4])/float(alignLenI) > 0.1]
		selectedL = [x for x in scoredL if x[0] not in kickL ]
		selectedNameL = [x[0][0] for x in selectedL]
		kickNameL = [x[0] for x in kickL]
		#STDERR("selectedNameL=",selectedNameL) ##DEBUG		
		if len(selectedNameL) == 0:
			EvrGapF = float(sum( [ x[5] for x in scoredL]))/float(len(scoredL))
			EvrInnerGapF = float(sum( [ x[4] for x in scoredL]))/float(len(scoredL))
			ExcludeSeqNameS = sorted(scoredL, key = lambda x:(abs(float(x[5]) - EvrGapF) + abs(float(x[4]) - EvrGapF)) )[-1][0][0].split()[0] ##get name of gap-prone sequence 
			ToAlignS = '\n'.join([x[0]+"\n"+x[1].replace('-','') for x in MSAfastaL if x[0].find(ExcludeSeqNameS) == -1])
			STDERR("ToAlignS=",ToAlignS) ##DEBUG
			NewFastaAlignmentS = musclecall(ToAlignS)
			SeqsL = MSAfastaSplit(NewFastaAlignmentS) 
			#SeqsL = [[x[0],x[1].replace('-','')] for x in MSAfastaL if x[0].find(ExcludeSeqNameS) == -1] ##remove gap-prone sequence from alignment-to-be sequences set
			#STDERR("While_NewFastaAlignmentS=",NewFastaAlignmentS) ##DEBUG
			#STDERR("ExcludeSeqNameS=",ExcludeSeqNameS) ##DEBUG
		iterNUMI += 1
		

	if len(kickL) == 0: 
		NewFastaAlignmentS = MSAfasta

	if len(kickL) != 0 and len(selectedL) > 1:
		ToAlignS = '\n'.join(['\n'.join(x[0]) for x in selectedL])
		NewFastaAlignmentS = musclecall(ToAlignS)
	elif len(kickL) != 0 and len(selectedL) == 1:
		STDERR("selectedL=1 ::",selectedL) ##DEBUG

	#if len(selectedL) == 0:
	#	selectedL
	STDERR("NewFastaAlignmentS=",NewFastaAlignmentS) ##DEBUG
	STDERR("kickL=",kickL) ##DEBUG
	STDERR("selectedNameL=",selectedNameL) ##DEBUG
	STDERR("kickNameL=",kickNameL) ##DEBUG

	return [NewFastaAlignmentS,kickL,selectedNameL,kickNameL]


def consensusExtractKW(MSAfasta,GapWeightF,MinPF,GapIgnoreFlagI): ##MinPF = minimum percent for consensus

	SeqsL = MSAfasta.split('>') ## Take input aligment as fasta format
	SeqsL = ['>' + x for x in SeqsL if len(x) > 0] ##Eliminate empty items
	SeqsL = [FastaTool(x) for x in SeqsL] ## Create fasta object for each sequences
	SeqsL = [x.seqonly().upper() for x in SeqsL] ## extract sequences from each fasta
	SeqsL = [x for x in SeqsL if len(x) > 0]  ##Eliminate empty items

	Cons = '' ##creat variable for consensus string storing 
	AnumI = 0;TnumI = 0;GnumI=0;CnumI=0;GapNum=0 ##innitiate gap count variables
	GapS = '-' ## define gap character 
	columndepthF = float(len(SeqsL)) ## Define column depth for voting 

	#EndGapI = 10 ##to modify (lower) gap score at the end of MSA
	#SwitchFlag = 0
	for i in range(len(SeqsL[0])): ## loop over coluns of alignment
		anMSAcolumnL = [x[i] for x in SeqsL] ## extract each base from each sequences in column i_th
		NucD = {} ## create dictionary for keep bases occurence
		RAWGapF = float(anMSAcolumnL.count(GapS))/columndepthF ## calculate gap ratio
		NucD[GapS] = RAWGapF*GapWeightF ## calculate gap ratio
		NucD['A'] = float(anMSAcolumnL.count('A'))/columndepthF ## calculate base "A" ratio
		NucD['T'] = float(anMSAcolumnL.count('T'))/columndepthF ## calculate base "T" ratio
		NucD['G'] = float(anMSAcolumnL.count('G'))/columndepthF ## calculate base "G" ratio
		NucD['C'] = float(anMSAcolumnL.count('C'))/columndepthF ## calculate base "C" ratio		

		VotesNumF = max(NucD.values()) ## select max occurences
		VotesL = [x for x in NucD if NucD[x] == VotesNumF] ## find out which base is the most popular
		if GapIgnoreFlagI == 0:
			VotesL = [x for x in NucD if NucD[x] == VotesNumF and NucD[x] >= MinPF]  ##  if the base is not popular enought then vote to N
		elif GapIgnoreFlagI == 1:
			VotesL = [x for x in NucD if NucD[x] == VotesNumF and NucD[x] >= (MinPF - NucD[GapS])] ## To ignore gap, reduces gap weight by gap ratio
		if len(VotesL) == 0: ## If there is not a popular one then vtoe to N for the column
			VotesL = ['N']

		if len(VotesL) == 1: ## If there is only one popular base in a list then vtoe to that base for the column
			Cons = Cons + VotesL[0]

		elif len(VotesL) == 2 and GapS in VotesL:  ## If there are one popular base and gap in a list then select that base
			VotesL = [x for x in VotesL if x != GapS]
			Cons = Cons + VotesL[0]
		elif len(VotesL) > 1 and  GapS not in VotesL: ## If more than one popular base then define that base in the column as N 
			Cons = Cons + 'N'

	Consout = Cons.strip().strip('N') ## delete N characters at both ends

	return Consout

def HitJoin(seedOBJ,neighborL): ## join blast hit if they are overlap with each other 
	if seedOBJ.sstartI < seedOBJ.sendI or seedOBJ.qstartI < seedOBJ.qendI:
		direction = 0 ## 5' --> 3'
	else:
		direction = 1 ## 3' --> 5'
	if direction == 1:
		headJoinL = [x for x in neighborL if x.sstartI > seedOBJ.sstartI and x.sstartI < seedOBJ.sstartI  ]
	return None

def PreAlign(SubClusterL):
	PlusToAlignSeqSL = [getSeqFull(DBname,x.sseqidS) for x in SubClusterL if x.qstartI < x.qendI and x.sstartI < x.sendI ]
	#STDERR("PlusToAlignSeqSL=",PlusToAlignSeqSL) ##DEBUG
	PreMinusToAlignSeqSL = [getSeqFull(DBname,x.sseqidS) for x in SubClusterL if x.qstartI > x.qendI or x.sstartI > x.sendI ]
	#STDERR("PreMinusToAlignSeqSL=",PreMinusToAlignSeqSL) ##DEBUG
	OBJMinusToAlignSeqSL = [FastaTool(x) for x in PreMinusToAlignSeqSL]
	MinusToAlignSeqSL = [x.FastaSrev for x in OBJMinusToAlignSeqSL]
	del(OBJMinusToAlignSeqSL); del(PreMinusToAlignSeqSL); 
	ToAlignSeqSL = MinusToAlignSeqSL + PlusToAlignSeqSL
	ToAlignSeqSS = '\n'.join(ToAlignSeqSL)
	return ToAlignSeqSS
	

def main(SubcookedBlastL): ##Take input as PROCESSED blast(n,x) results LIST
	#cookedBlastL = [hit(x) for x in chop(INS)] #transform blast result into list of blast hit object
	ClusteredL = merge4(SubcookedBlastL,0.05)
	#STDERR("len(ClusteredL)=",len(ClusteredL))##DEBUG
	for i in ClusteredL: 
		STDERR('\n'.join([str(x.AttributesL) for x in i])) ##DEBUG
		STDERR("Each iter in ClusteredL ++++++++++++++++++++") ##DEBUG
		#STDERR("ToAlignSeqSL\n",ToAlignSeqSS)
		ToAlignSeqSS = PreAlign(i)
		FastaAlignmentS = musclecall(ToAlignSeqSS)
		FastaAlignmentSclw = musclecallCLW(ToAlignSeqSS)
		STDERR("FastaAlignmentS=",FastaAlignmentSclw)##DEBUG		
		
		AlignImproveL = AlignImprove(FastaAlignmentS)	
		NewAlignS = AlignImproveL[0]
		#STDERR("NewAlignS=",NewAlignS) ##DEBUG
		consensusS = consensusExtractKW(NewAlignS,0,0.01,1) ##Generate consensus with Minnimum base ration and ignore gap 
		STDERR("consensusS=",consensusS)
		#STDERR("AlignImproveL[1]=",AlignImproveL[1])
		STDERR("selectedNameL=",AlignImproveL[2])
		STDERR("kickNameL=",AlignImproveL[3])
		return AlignImproveL

def main1(CurrentNameS,SubcookedBlastL): ## Direct Checking Matched sequence by Performing Alignment
	ClusteredL = merge4(SubcookedBlastL,0.05)
	#STDERR("main1.ClusteredL=",ClusteredL)
	SelectedClustL = []
	for i in ClusteredL:
		if len([x.sseqidS.find(CurrentNameS) != -1 for x in i ]) > 0:
			SelectedClustL.append(i) 
		
	SeqSL = [">" + x for x in PreAlign(SelectedClustL[0]).split(">") if len(x) > 3]
	#STDERR("main1.SeqSL=",SeqSL)
	CurrentSeqS = [x for x in SeqSL if x.find(CurrentNameS+" ") != -1][0]
	CompareSeqSL = [x for x in SeqSL if x.find(CurrentNameS+" ") == -1]
	PickedL = SimilarPick(CurrentSeqS,CompareSeqSL,0)

	return PickedL ## [TargetSeqS, PickedSeqS1st, .... , PickedSeqSNth]
	


def main0(INS):
	#cookedBlastL = sorted([hit(x) for x in chop(INS)], key=lambda x: x.qlenI )[::-1]
	cookedBlastL = sorted([hit(x) for x in chop(INS)], key=lambda x: ( x.qlenI,x.qseqidS )  )[::-1]	
	STDERR("len(cookedBlastL)=",len(cookedBlastL))
	STDERR("cookedBlastL[:5]=",cookedBlastL[:5])
	ContigNamelistL = ['']
	for i in cookedBlastL:
		if i.qseqidS != ContigNamelistL[-1] or i.qseqidS.split()[0].replace('lcl|','') not in ContigNamelistL:
		#if i.qseqidS.split()[0].replace('lcl|','') not in ContigNamelistL:
			ContigNamelistL.append(i.qseqidS.split()[0].replace('lcl|',''))
		STDERR(i.qseqidS,i.qlenI) ##DEBUG

	ContigNamelistL.pop(0)
	#STDERR("ContigNamelistL[0:4]=",ContigNamelistL[0:4]) ##DEBUG	
	#STDERR("len(ContigNamelistL)=",len(ContigNamelistL)) ##DEBUG

	if options.o != '0':
		f=open(options.o,'w')
	
		
	while len(ContigNamelistL) > 0:
		CurrentNameS = ContigNamelistL.pop(0) ## Pop 1 of name in the list to use as seed 
		STDERR("CurrentNameS=",CurrentNameS) ##DEBUG
		SubcookedBlastL = [x for x in cookedBlastL if x.qseqidS.split()[0].replace('lcl|','') == CurrentNameS]
		#STDERR("SubcookedBlastL=",SubcookedBlastL) ##DEBUG
		selectedSeqSL = main1(CurrentNameS,SubcookedBlastL) ## Get Sequnces to create consensus and remove the member from ContigNamelistL
		PickedOBJL = [FastaTool(x) for x in selectedSeqSL] ## Transform each SeqS into Fasta object
		#STDERR("selectedSeqS=",selectedSeqS)
		if len(selectedSeqSL) == 1:
			consensusS = selectedSeqSL[0]
		elif selectedSeqSL > 1:
			FastaAlignmentS = musclecall('\n'.join(selectedSeqSL))
			consensusS = ">"+ CurrentNameS + " ConS\n" + consensusExtractKW(FastaAlignmentS,0,0.01,1) + "\n"
		STDERR("main0.consensusS=",consensusS) ##DEBUG
		ExcludeNameL = [x.name() for x in PickedOBJL] ## Get Name For Excluding
		STDERR("main0.ExcludeNameL=",ExcludeNameL) ##DEBUG

		ContigNamelistL = [x for x in ContigNamelistL if x not in ExcludeNameL]
		STDERR("######## len(ContigNamelistL)=",len(ContigNamelistL)) ##DEBUG
		
		MemberS = ";".join(ExcludeNameL)
		STDERR("main0.MemberS=",MemberS)
		
		if options.o != '0':
			f.write(consensusS)
		else:
			print consensusS.strip()

	if options.o != '0':		
		f.close()

### Input and Option 
usage = "python Blastout2cnsensus.py -i input.gff -o out_file"
opt = optparse.OptionParser(usage)
opt.add_option("-i",help="*input path, get input from stdin if ommit", default='0')
opt.add_option("-o",help="indicate output file name or print out as standard output",default="0")
opt.add_option("-q",help="*query database path",dest='q',default="DB_PATH")
opt.add_option("-d",help="subject database path",dest='d')
opt.add_option("-e",help="E-value cutoff for selected hit",dest="e",default="10.0")
opt.add_option("-c",help="query-length covery rate cutoff for selected hit",default="0.5")
opt.add_option("--TAG",help="Spicies Tag in protein name",default="0",dest="TAG")
(options, args) = opt.parse_args()

##Documentation part
__QCovF__ = float(options.c)*100.0 ;STDERR("__QCovF__ = ",__QCovF__)
__EvalF__ = float(options.e)	;STDERR("__EvalF__ = ",__EvalF__)
__TAG__ = options.TAG
DBname = options.q

if options.i == '0': ##get input from pipe
	INS = sys.stdin.read()

else:#open file
	f=open(options.i,'r')
	INS = f.read()
	f.close()
if options.o == '0': ##print output to pipe
	main0(INS)
	#print(main(INS))
else:#write output to a flie
	f=open(options.o,'w')
	f.write(main0(INS))
	f.close()	



