#!/usr/bin/python
##TXblastX.py.py
##Takes input from blastx with -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qcovs' option

import sys,optparse, copy
### Class And Function
def STDERR(*StrInS): #Function For Debugging
	outS = ' '.join([str(x) for x in StrInS])
	sys.stderr.write(str(outS)+'\n')

def chop(blastxS): #divide each hit into a list of attributes 
	if __TAG__ != '0':
		blastxL = [y for y in blastxS.splitlines() if y.find(__TAG__) != -1]
	else:
		blastxL = [y for y in blastxS.splitlines()]		
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

def merge4(Ori_inputlistOBJ,coverage):
	inputlistOBJ = copy.deepcopy(Ori_inputlistOBJ)
	a = sorted(inputlistOBJ, key = lambda x:abs(x.sstartI - x.sendI) )[::-1]
	clustL = []

	while len(a) > 0:
	
		mark = a.pop(0)
		#clustL.append([])

		Mhead = mark.sstartI
		Mtail = mark.sendI
		Mlength = abs(mark.sendI - mark.sstartI)
		
		## <------------------->
		##    <------------->
		temp = [x for x in a if x.sstartI >= Mhead and x.sendI <= Mtail] #and abs(x.sendI - x.sstartI) > (coverage* Mlength)]

		##    <------------------->
		## <------------------>
		temp2 = [x for x in a if x.sendI < Mtail and Mlength >= abs(x.sendI - x.sstartI) >= (coverage* Mlength) and (Mhead - x.sstartI) <= coverage*(abs(x.sendI - x.sstartI))] 
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		
		## <------------------->
		##     <------------------>		
		temp2 = [x for x in a if x.sstartI > Mhead and Mlength >= abs(x.sendI - x.sstartI) >= (coverage* Mlength) and (x.sendI - Mtail) <= coverage*(abs(x.sendI - x.sstartI))]
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		temp = sorted(temp, key=lambda temp:temp.sstartI)
		
		a = [x for x in a if x not in temp]
		
		pick = [x for x in temp]
		pick.insert(0,mark) 
		clustL.append(pick)
	return clustL

def GroupScore(Eachgrouped1L):
	groupedL =  merge4(Eachgrouped1L,0.75)
	EachgroupedL = []
	for L in groupedL:
		CurrentL = sorted(L, key=lambda x:x.bitscoreF)[-1]
		EachgroupedL.append(CurrentL)
	
	OBJLsortedL = sorted(EachgroupedL, key=lambda x:x.sstartI )
	if OBJLsortedL[0].sstartI < OBJLsortedL[0].sendI:
		HeadI = OBJLsortedL[0].sstartI
		TailI = OBJLsortedL[0].sendI
	elif OBJLsortedL[0].sstartI > OBJLsortedL[0].sendI:
		HeadI = OBJLsortedL[0].sendI
		TailI = OBJLsortedL[0].sstartI

	SLenI = OBJLsortedL[0].slenI

	CovNumI = 0
	CurrentTailI = 0
	for i in OBJLsortedL[1:]:
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
		
	

def main(INS):
	INL = (chop(INS))
	hitOBJL = [hit(x) for x in INL]
	if __TAG__ != '0':
		selectOJBL = [x for x in hitOBJL if x.qcovsF > __QCovF__ and x.evalueF < __EvalF__]
	selectOJBL = [x for x in hitOBJL if x.qcovsF > __QCovF__ and x.evalueF < __EvalF__]
	#selectOJBL = selectOJBL[0:1000] ##DEBUG
	
	STDERR([x.AttributesL for x in selectOJBL[0:2]]) ##DEBUG
	STDERR("entering group function") ##DEBUG
	outS = ''
	#for OBJL in group(selectOJBL):
	for OBJL in Pgroup(selectOJBL):
		ScoreFL = GroupScore(OBJL)
		sortedOBJL = sorted([x.AttributesL for x in OBJL], key=lambda x:(int(x[8]) + int(x[9]))/2 )
		#outS = outS + '\n'.join(['\t'.join(x.AttributesL) for x in sortedOBJL]) + "\n############\n"
		#outS = outS + '\n'.join(['\t'.join(x) for x in sortedOBJL]) + "\nPoolCovScpre=" + str(ScoreFL[1]) + "\nEachgroupedL=" + str(ScoreFL[0]) + "\n############\n"
		outS = outS + '\n'.join(['\t'.join(x) for x in sortedOBJL]) + "\nPoolCovScpre=\t" + str(ScoreFL[1]) +"\t"+ ScoreFL[0][0][1] + "\n############\n"		
		outS = outS + str(ScoreFL[0]).replace("""],""","""]\n""") + "\n############\n"
	#STDERR(group(selectOJBL)[0])
	return outS

### Input and Option 
usage = "python TXblastX.py.py -i input.gff -o out_file"
opt = optparse.OptionParser(usage)
opt.add_option("-i",help="*input path, get input from stdin if ommit", default='0')
opt.add_option("-o",help="indicate output file name or print out as standard output",default="0")
opt.add_option("-e",help="E-value cutoff for selected hit",default="0.1")
opt.add_option("-c",help="query-length covery rate cutoff for selected hit",default="0.98")
opt.add_option("--TAG",help="Spicies Tag in protein name",default="0",dest="TAG")
(options, args) = opt.parse_args()

##Documentation part
__QCovF__ = float(options.c)*100.0 ;STDERR("__QCovF__ = ",__QCovF__)
__EvalF__ = float(options.e)	;STDERR("__EvalF__ = ",__EvalF__)
__TAG__ = options.TAG

if options.i == '0': ##get input from pipe
	INS = sys.stdin.read()

else:#open file
	f=open(options.i,'r')
	INS = f.read()
	f.close()
if options.o == '0': ##print output to pipe
	print(main(INS))
else:#write output to a flie
	f=open(options.o,'w')
	f.write(Main(INS))
	f.close()	



