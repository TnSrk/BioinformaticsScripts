##FNpool
#!/usr/bin/python
#ContigDominoE10_batch5 to join contigs into putative longer contigs using Reads data 
#add Contig-over-reads priority and minimum reads number setting
## Read clustering for MSA
## resolve homopolymer with grep+regex
## Merge while loops for Read Grep  
## add Read select for ReadGrep
## regex module in longerSeq function
## adding realign module
## add ReadGrepNum option
## use simple match in realign
## 2014Oct27 exclude unitigname which used

import sys,  subprocess, optparse, re, time
from time import sleep
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor


def STDERR(*StrInS): #Function For Debugging
	outS = ' '.join([str(x) for x in StrInS])
	sys.stderr.write(str(outS)+'\n')


def ParallelLoop(FunctionS,ListL): ##Warning input Function should not print or write anything because output into a file at a time may cause problem
	npI = int(cpu_count())
	ChunkSizeI = npI*2
	chunkNumI = len(ListL)/ChunkSizeI
	BigChunkL = []
	for i in range(chunkNumI):
		
		BigChunkL.append(ListL[ChunkSizeI*i:ChunkSizeI*(i + 1)])

	remainsI = len(ListL)%ChunkSizeI
	if remainsI != 0:
		BigChunkL.append(ListL[(ChunkSizeI*chunkNumI):])

	BigClustL = []
	for chunkL in BigChunkL:
		chunkL = [[FunctionS] + [x] for x in chunkL]
		#STDERR("chunkL="+str(chunkL)) ##DEBUG
		chnkNPI = len(chunkL)
		pool = Pool(processes=chnkNPI)
		BigtempL = pool.map(multiCall, chunkL)	
		for j in BigtempL:
			BigClustL.append(j)

		pool.close()
		del BigtempL

	return BigClustL

def multiCall(FunctionSinputL):##for calling function in parallel mode
	Funct = FunctionSinputL[0]
	FunctionOut = Funct(*FunctionSinputL[1])
	#STDERR(str(Funct)+str(time.strftime("%H:%M:%S"))) ##DEBUG
	return FunctionOut

def ConcurrentCall(Fnc, InputL, threadsI, *arg):
	waitF = 0.0
	if len(arg) == 1:
		waitF = float(arg[0])
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
		#STDERR(str(xQueL))
		if False not in xQueL:
			#OutputL.append(xQueL)
			break
		else:
			if waitF > 0.0:
				sleep(waitF)

		RoundI += 1

	OutputL = [QueD["Q"+str(i)].result() for i in numL if QueD["Q"+str(i)].done() == True]

	return OutputL

class FastaTool(object): ## class for manipulate sequences in fasta format
	def __init__(self,FastaInS):
		self.FastaInS = FastaInS
		self.seqLenI = self.seqlen()
		self.FastaSrev  =  self.header() +"Rev\n" + self.reversecomplement()

	def seqonly(self):
		fastaL = self.FastaInS.splitlines()
		#fastaL = [x for x in fastaL if x.find(">") != 0]
		fastaL = [x for x in fastaL if x.find(">") == -1]
		fastaL = [x.strip() for x in fastaL]
		fastaS = ''.join(fastaL)
		return fastaS

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
		headerS = fastaL[0]
		return headerS

	def name(self):
		NameS = self.header().split(' ')[0].replace('>','').replace('lcl|','')
		return NameS

	def ManyLines(self,widthI):
		seqS = self.seqonly()
		seqLineI = len(seqS)/widthI
		remainsI = len(seqS)%widthI
		#TEMPseqS = self.header()+"\n"
		TEMPseqS = ''
		for i in range(seqLineI):
			TEMPseqS = TEMPseqS + seqS[i*(widthI):(i+1)*(widthI)] + '\n'

		TEMPseqS = TEMPseqS + seqS[(-1)*(remainsI):] + '\n'
		return TEMPseqS


	def reversecomplement(self):
		string = self.seqonly()
		#header = self.name() + "c"
		reverse = string[::-1].upper()

		reverse = reverse.replace("A","?")
		reverse = reverse.replace("T","A")
		reverse = reverse.replace("?","T")
		reverse = reverse.replace("G","?")
		reverse = reverse.replace("C","G")
		reverse = reverse.replace("?","C")
		if "U" in reverse:
			reverse = reverse.replace("A","?")
			reverse = reverse.replace("U","A")
			reverse = reverse.replace("?","U")
		reverseS = reverse	
		#print self.reverseS
		return reverseS ## complementary string ATCG --> CGAT

	def seqlen(self):
		
		seqLenI = len(self.seqonly())
		return seqLenI

def realign(fastaMsaS):
	SeqsL = fastaMsaS.split('>')
	SeqsL = ['>' + x.strip() for x in SeqsL if len(x) > 0]

	AlignLenI = len(seqonly(SeqsL[0]))
	SeqsL25HeadL = [x for x in SeqsL if len(seqonly(x)) - len(seqonly(x).lstrip('-')) <= AlignLenI*0.25]
	SeqsL25TailL = [x for x in SeqsL if len(seqonly(x)) - len(seqonly(x).rstrip('-')) <= AlignLenI*0.25]
	SeqsL25MidL = [x for x in SeqsL if x not in SeqsL25HeadL and x not in SeqsL25TailL] 

	allL = [SeqsL25HeadL,SeqsL25MidL,SeqsL25TailL]
	allL = ['\n'.join(x) for x in allL if len(x) > 0]
	
	#sys.stderr.write("THIS IS allL="+str(allL)+'\n')	
	STDERR("THIS IS allL="+str(allL)) ##DEBUG

	alignS = ''
	numI = 1
	for i in allL:
		#print "###############"
		if i.count('>') == 1:alignoutS = i.replace('-','')
		elif i.count('>') > 1:alignoutS = consensusExtract2(musclecall(i.replace('-','')),1.0,0.25)

		#print alignoutS
		alignS = alignS + '>' + str(numI) + '\n' + seqonly(alignoutS).strip()  + '\n'
		numI += 1

	sys.stderr.write("THIS IS alignS="+str(alignS)+'\n')	
	
	allL = alignS.split('>')
	allL = ['>'+x for x  in allL if x != '']

	#sys.stderr.write("THIS IS allL="+str(allL)+'\n')
	STDERR("THIS IS allL after consensus extract in  ="+str(allL)) ##DEBUG
	LongerSeq = allL.pop(0)
	while len(allL) > 0:
		popseqS = allL.pop(0) 
		if DirectionS == '3':
			PflagI = 0
		else:
			PflagI = 1
		#InS = LongerSeq + '\n' + popseqS
		sys.stderr.write("simplejoin2 is called in realign\n") ##DEBUG
		#LongerSeq = LongerSeqSmain(InS,DirectionS,PflagI)
		if SubAlignCheck(seqonly(LongerSeq),seqonly(popseqS)) == 1:
			sys.stderr.write("simplejoin2 is NOT called in realign\n") ##DEBUG
			sys.stderr.write("LongerSeq=\n"+LongerSeq+'\n') ##DEBUG
			sys.stderr.write("popseqS=\n"+popseqS+'\n') ##DEBUG

		elif len(seqonly(LongerSeq)) >= len(seqonly(popseqS)):
			LongerSeq = simplejoin2(seqonly(LongerSeq),seqonly(popseqS),0)
		else:
			LongerSeq = simplejoin2(seqonly(LongerSeq),seqonly(popseqS),1)
		#LongerSeqSmain(TEMPbridgeS1+'\n'+NextContigS,DirectionS,PflagI)		
		
	#print LongerSeq
	#alignoutS = consensusExtract2(musclecall(alignS),0.1,0.25)
	LongerSeq = seqonly(LongerSeq)
	return LongerSeq

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


def ConsMod2(seqS,WildCardS):
	ModSeqS = seqS.splitlines()
	SubS = ''.join(ModSeqS)
	SplitSubSL = SubS.split(WildCardS)
	#SubSL =  [x ofr x in SplitSubSL if len(x) > 0]

	ModSeqS = ''
	TempS = ''	
	GapcountI = 3

	#sys.stderr.write("this is SplitSubSL before loop\n"+str(SplitSubSL)+"\n") ##DEBUG
	if max([len(x) for x in SplitSubSL]) >= 5 and len(SplitSubSL) > 0:
		HeadI = 0
		TailI = 0
		i = -1
		while HeadI < 5:
			i += 1
			HeadI =  len(SplitSubSL[i])
			
		HeadI = i
		i = 0
		while TailI < 5:
			i += 1
			TailI =  len(SplitSubSL[-i])
			
		TailI = -i

		if TailI == -1:TailI = len(SplitSubSL) + 1
		else: TailI += 1
		SplitSubSL = SplitSubSL[HeadI:TailI]

		#sys.stderr.write("this is TailI HeadI\n"+str(HeadI)+" " + str(TailI) +"\n") ##DEBUG		
		#sys.stderr.write("this is SplitSubSL\n"+str(SplitSubSL)+"\n") ##DEBUG
		for i in SplitSubSL[::-1]:
			if len(i) > GapcountI:
				#TempS = WildCardS + i + TempS
				TempS =  i + TempS
				GapcountI += len(i)*-1
	
		ModSeqS = ModSeqS + TempS	
	else:ModSeqS = ''
	return ModSeqS


def consensusExtract2(MSAfasta,GapWeightF,MinPF): ##MinPF = minimum percent for consensus
	#sys.stderr.write("MSAfasta ="+MSAfasta+'\n') ##DEBUG

	SeqsL = MSAfasta.split('>')
	SeqsL = ['>' + x for x in SeqsL if len(x) > 0]
	SeqsL = [seqonly(x).upper() for x in SeqsL]
	SeqsL = [x for x in SeqsL if len(x) > 0]

	Cons = ''
	AnumI = 0;TnumI = 0;GnumI=0;CnumI=0;GapNum=0
	GapS = '-'
	columndepthF = float(len(SeqsL))
	#sys.stderr.write("columndepthF=\n"+str(columndepthF)+'\n') ##DEBUG

	#EndGapI = 10 ##to modify (lower) gap score at the end of MSA
	#SwitchFlag = 0
	for i in range(len(SeqsL[0])):
		anMSAcolumnL = [x[i] for x in SeqsL]
		NucD = {}
		NucD['A'] = float(anMSAcolumnL.count('A')/columndepthF)
		NucD['T'] = float(anMSAcolumnL.count('T')/columndepthF)
		NucD['G'] = float(anMSAcolumnL.count('G')/columndepthF)
		NucD['C'] = float(anMSAcolumnL.count('C')/columndepthF)
		NucD[GapS] = float((anMSAcolumnL.count(GapS)/columndepthF)*GapWeightF)

		#sys.stderr.write("max(NucD)=\n"+str(max(NucD.values()))+'\n') ##DEBUG
		VotesNumF = max(NucD.values())		
		#sys.stderr.write("VotesNumF=\n"+str(VotesNumF)+'\n') ##DEBUG
		#sys.stderr.write("NucD=\n"+str(NucD)+'\n') ##DEBUG
		VotesL = [x for x in NucD if NucD[x] == VotesNumF]
		VotesL = [x for x in NucD if NucD[x] == VotesNumF and NucD[x] >= MinPF]
		if len(VotesL) == 0:
			VotesL = ['N']
				
		#sys.stderr.write("VotesL="+str(VotesL)+ str(NucD) +'\n') ##DEBUG
		
		if len(VotesL) == 1:
			Cons = Cons + VotesL[0]

		elif len(VotesL) == 2 and GapS in VotesL:
			VotesL = [x for x in VotesL if x != GapS]
			Cons = Cons + VotesL[0]
		elif len(VotesL) > 1 and  GapS not in VotesL:
			Cons = Cons + 'N'

	#sys.stderr.write("this is Cons\n") ##DEBUG
	#sys.stderr.write(Cons+'\n') ##DEBUG

	#consensus2 = str(Cons).replace('-','')
	#Consout = consensus2.strip().strip('N')
	Consout = Cons.strip().strip('N')
	ConsoutS = ConsMod2(Consout,'N')
	ConsoutS = ConsoutS.strip().strip('-')
	ConsoutS = ConsMod2(ConsoutS,'-')

	#sys.stderr.write("this is ConsoutS\n") ##DEBUG
	#sys.stderr.write(ConsoutS+'\n') ##DEBUG

	return ConsoutS

def BlastCall2(seqS,CdbS):
	arg = "blastall -p blastn -m 8 -F F -d "+CdbS	
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(seqS)
	blastout = str(x[0])
	del process
	return blastout

def BlastCall3(seqS,CdbS,DirectionS):

	arg = "blastall -p blastn -W 7 -m 8 -F F -d "+CdbS	
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(seqS)
	blastout = str(x[0])
	del process
	return blastout

def blastTrans(blastResult,Eval):
	b = blastResult.splitlines()
	b = [x for x in b if len(x) > 0]
	b = [x.split('\t') for x in b]
	pL = [[x[0], x[1], float(x[2]), int(x[3]), int(x[4]), int(x[5]), int(x[6]), int(x[7]), int(x[8]), int(x[9]), float(x[10]), float(x[11])] for x in b]
	pL = [x for x in pL if float(x[-2]) < Eval]
	return pL

def unitigSelect(seqS,CdbS,blastoutnameS):
	arg2 = """Seqname=`head -n 1 __blastoutnameS__ | cut -f 1`;Seqlen=`grep -v "^>" __seqS__ | tr -d "\n" | wc -c` ; cat __blastoutnameS__ | sed "s/"$Seqname"/Query/" | sed 's/_/\t/' | sed 's/_/\t/' | awk -v Seqlen=$Seqlen '{if ($11 > $12 && $12 > 5 && $10 > (Seqlen - 5) && $11 > ($3 - 5)) print $0 "\t" $2 "_" $3 "_" $4 " -S 2";if ($11 < $12 && $10 > (Seqlen - 5) && $11 < 5 && $12 < ($3 - 25) ) print $0 "\t" $2 "_" $3 "_" $4 " -S 1" }' | sed 's/Query/'$Seqname'/' |  sort -k14,14gr | awk '{if ($6 > 15) print $0}' | grep -v "32850_1593_219178_5065" |head -n 1 """
	arg2 = arg2.replace("__seqS__",seqS)
	arg2 = arg2.replace("__blastoutnameS__",blastoutnameS)

	#print arg2 ##DEBUG
	process = subprocess.Popen(arg2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	blastout = x[0]
	#print blastout ##DEBUG
	del process
	return blastout.strip()

def unitigSelect2(seqS,blastResultL,DirectionS):
	OverlapI = 10

	SSeqS = seqonly(seqS)
	SSeqLenI = len(SSeqS)
	#print SSeqLenI ##DEBUG
	blastResultLmod	= [x[:] + [int(x[1].split('_')[1])] for x in blastResultL]
	plusResultL = [x for x in blastResultLmod if x[8] < x[9]]	
	minusResultL = [x for x in blastResultLmod if x[8] > x[9]]

	if DirectionS == '5':
		plusResultL = [x for x in plusResultL if x[6] < 3 and x[8] > OverlapI and x[9] > (x[12] - 3)]
		#plusResultL = [x for x in plusResultL if x[6] < 3 and x[9] > (x[12] - 3)]
		#minusResultL = [x for x in minusResultL if x[6] < 3 and x[8] < (x[12] - OverlapI) and x[9] < 3]
		minusResultL = [x for x in minusResultL if x[6] < 3 and x[9] < 3]

	if DirectionS == '3':
		plusResultL = [x for x in plusResultL if x[7] > (SSeqLenI - 3) and x[8] < 3 and (x[12] - x[9]) > OverlapI]
		minusResultL = [x for x in minusResultL if x[7] > (SSeqLenI - 3) and x[8] > (x[12] - 3) and x[9] > OverlapI]

	selectL = plusResultL + minusResultL
	if len(selectL) > 0:selectL =  sorted(selectL, key = lambda selectL:selectL[11])[-1]
	else: selectL = []
	return selectL

def FullOverlap(seqS,blastResultL):

	SSeqS = seqonly(seqS)
	SSeqLenI = len(SSeqS)
	#print SSeqLenI ##DEBUG
	blastResultLmod	= [x[:] + [int(x[1].split('_')[1])] for x in blastResultL]
	plusResultL = [x for x in blastResultLmod if x[8] < x[9] and x[8] < 3 and  x[9] >= (x[12] - 3)]	
	minusResultL = [x for x in blastResultLmod if x[8] > x[9] and x[9] < 3 and  x[8] >= (x[12] - 3)]

	selectL = plusResultL + minusResultL
	if len(selectL) > 0:selectL =  sorted(selectL, key = lambda selectL:selectL[11])
	else: selectL = []
	return selectL

def ReadsSelelct(seqS,dbpathS,blastResultL,DirectionS):
	OverlapI = 10
	mismatchI = 2

	SSeqS = seqonly(seqS)
	SSeqLenI = len(SSeqS)

	sys.stderr.write(str(SSeqLenI)+"\n") ##DEBUG
	blastResultLmod	= [x[:] + [SeqLenFastacmd(dbpathS,x[1])] for x in blastResultL]
	blastResultLmod	=  [x for x in blastResultLmod if x[5] < mismatchI] 

	for i in blastResultLmod: ##DEBUG
		sys.stderr.write(str(i)+"\n")
	plusResultL = [x for x in blastResultLmod if x[8] < x[9]]	
	minusResultL = [x for x in blastResultLmod if x[8] > x[9]]

	if DirectionS == '5':
		
		##READ		[8]--------->[9]
		##consensus	    [6]------------->[7]
		plusResultL = [x for x in plusResultL if x[6] < 3 and x[8] > OverlapI and x[9] > (x[12] - 3)]

		##READ		[8]<---------[9]
		##consensus	    [6]------------->[7]
		minusResultL = [x for x in minusResultL if x[6] < 3 and x[8] < (x[12] - OverlapI) and x[9] < 3]

	if DirectionS == '3':
		plusResultL = [x for x in plusResultL if x[7] > (SSeqLenI - 3) and x[8] < 3 and (x[12] - x[9]) > OverlapI]
		minusResultL = [x for x in minusResultL if x[7] > (SSeqLenI - 3) and x[8] > (x[12] - 3) and x[9] > OverlapI]

	selectL = plusResultL + minusResultL
	
	return selectL

def SeqLenFastacmd(dbpathS,SeqNameS):
	        
	arg = "fastacmd -d __dbpathS__ -s __SeqNameS__ | grep -v '>' | tr -d '\n' | wc -c """
	arg = arg.replace('__dbpathS__',dbpathS); arg = arg.replace('__SeqNameS__',SeqNameS)
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()[0]
	if len(x) < 1:
		x = '0'
	return int(x)


def F3primeReadCall(blastResult,OutputName,QueryFileName,DatabasesPath):
	arg = """./3primeReadV3.sh """+ blastResult +" "+ OutputName +" "+ QueryFileName +" "+ DatabasesPath
	
	#print arg ##DEBUG
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	TEMP = x[0]
	sys.stderr.write(str(x[1])) ##DEBUG
	del process      
	return TEMP

def getSseq(DBname,Scontigname,Direction):
	arg = "fastacmd -p F -d "+DBname+" -s "+Scontigname+" ;"  
	if Direction == '2':   
		arg = "fastacmd -p F -d "+DBname+" -s "+Scontigname+" -S 2;"    
                
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	TEMPseqs = x[0].decode()
	#sys.stderr.write(str(x[1])) ##DEBUG
	del process
	return  TEMPseqs

def musclecall(SeqS,*arg):
	AlignTagI = 0
	argS = "muscle -maxiters 32 -quiet"
	if len(arg) == 1:
		AlignTagI = arg[0]
	if AlignTagI == 0:
		argS = argS + " -clwstrict"
	elif AlignTagI == 1:
		argS = argS

	SeqSb = SeqS.encode()
	#arg = "muscle -maxiters 32 -gapopen -1200 -quiet"
	
	process = subprocess.Popen(argS, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(SeqSb)[0].decode()
	if len(x) == 0:
		x = '0'
	return x

def musclecallGap(SeqS):
	SeqSb = SeqS.encode()
	arg = "muscle -maxiters 32 -gapopen -1200 -quiet"
	#arg = "muscle -maxiters 32 -quiet"
	#print "this is alignInS",alignInS #DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(SeqSb)[0].decode()
	if len(x) == 0:
		x = '0'
	return x


def simplejoin2(seqIn1S,seqIn2S,PI):

	alignlenI = min([len(seqonly(seqIn1S)),len(seqonly(seqIn2S))])
	seq1S =seqonly(seqIn1S)[-alignlenI:];seq2S =seqonly(seqIn2S)[:alignlenI]

	alignInS = ">S1S\n"+seq1S+"\n>S2S\n"+seq2S+"\n"

	arg = "muscle -maxiters 32 -quiet"
	#print "this is alignInS",alignInS #DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(alignInS)[0]
	#print "this is x ",x #DEBUG
	xL = x.split('>')
	xL = [x for x in xL if len(x) > 0]
	#print "This is xL " ,xL #DEBUG
	del process

	HeadL = [x for x in xL if x.find("S1S") != -1]
	#print "This is HeadL " ,HeadL #DEBUG
	HeadL = HeadL[0].splitlines()[1:]
	HeadL = [x.strip() for x in HeadL]
	HeadS = ''.join(HeadL)
#	print "This is HeadS " ,HeadS #DEBUG
	TailL = [x for x in xL if x.find("S2S") != -1]
	TailL = TailL[0].splitlines()[1:]
	TailL = [x.strip() for x in TailL]
	TailS = ''.join(TailL)
#	print "This is TailS " ,TailS #DEBUG
	HGapNumI = HeadS.count('-')
	TGapNumI = TailS.count('-')
#	print "GapNumI=",GapNumI ##DEBUG
	#InnerGapI = HeadS[0:(-1)*(len(HeadS) - GapNumI)].count('-')
	HInnerGapI = HeadS.rstrip('-').count('-')
	TInnerGapI = TailS.lstrip('-').count('-')
	HGapNumI2 = HGapNumI - HInnerGapI 
	TGapNumI2 = TGapNumI - TInnerGapI 

	if PI == 0 and HGapNumI2 == 0:
		AlignoutS = seq1S  
	elif PI == 0 and HGapNumI2 > 0:
		AlignoutS = seq1S + TailS[(-1)*HGapNumI2:] 

	elif PI == 1 and TGapNumI2 == 0:
		AlignoutS = seq2S 
	elif PI == 1 and TGapNumI2 > 0:
		AlignoutS = HeadS[:HGapNumI2] + seq2S 
	else:
		print( "ERR at LongerSeq")

	AlignoutS = '>' + "simplejoin_OUT\n" + AlignoutS +"\n"
	return AlignoutS

def MSAfastaSplit(MSAfasta):
		MSAfastaL = MSAfasta.split('>') ## Take input aligment as fasta format
		MSAfastaL = ['>' + x for x in MSAfastaL if len(x) > 0] ##Eliminate empty items
		MSAfastaL2 = []
		for x in MSAfastaL:## Separate each sequences
			xL = x.splitlines()
			MSAfastaL2.append( [ xL[0],''.join(xL[1:]) ] )

		return MSAfastaL2 ##Return 2D List [[HeaderLine0S,seq0S],[HeaderLine1S,seq1S]]

def MSACheck(MSAS): #Check consistency of MSA 
	MSAfastaL = [x for x in MSAfastaSplit(MSAS)]
	scoredL = [MSAQual(x) for x in MSAfastaL] #create list of object of Alignment Score of each sequence in MSA
	#STDERR("MSACheck.scoredL",scoredL)
	return scoredL

class MSAQual(object): #Take a line of alignment in fasta format then return alignment quality parameters for a sequence
	def __init__(self,OneFastaMSAS):
		self.OneFastaMSAS = OneFastaMSAS
		self.nameS = OneFastaMSAS[0].replace(">","").split()[0].replace("lcl|","")
		self.SeqS = OneFastaMSAS[1]		
		self.alignLenI, self.AllGapI, self.HeadGapI, self.TailGapI, self.InnerGapI, self.GapOpenI, self.FragNumI, self.LargeRatioF, self.SmallRatioF = self.AlignScore(self.SeqS)
		self.GapRatioF = self.InnerGapI/len(self.SeqS)
	
	def AlignScore(self,Seqs):
		alignLenI = len(Seqs)
		AllGapI = Seqs.count("-")
		HeadGapI = alignLenI - len(Seqs.lstrip("-")) #count lead gap in alignment
		TailGapI = alignLenI - len(Seqs.rstrip("-")) #count tail gap in alignment
		InnerGapI = AllGapI - (HeadGapI + TailGapI) #count inner gap in alignment
		GapOpenI = len([x for x in Seqs.strip("-").split("-") if len(x) > 0]) - 1 #count inner gap opening in alignment
		FragmentL = sorted([x for x in re.sub('\-+','-',Seqs).split('-') if len(x) != 0], key=lambda x:len(x))
		FragNumI = len(FragmentL)
		LargeRatioF = len(FragmentL[-1])/alignLenI
		SmallRatioF = len(FragmentL[0])/alignLenI
		if SmallRatioF == LargeRatioF:
			SmallRatioF = 0.0
		scoredL = [alignLenI, AllGapI, HeadGapI, TailGapI, InnerGapI, GapOpenI, FragNumI, LargeRatioF, SmallRatioF]
		return 	scoredL


def OverlapCheck(seq1,seq2, *arg):## If seq1 tail overlap with seq2 head then return overlap length, else return 0
	if len(arg) > 0:
		relaxI = arg[0]
	else:
		relaxI = 5
	seq1S = str(seq1[:]);seq2S = str(seq2[:])
	Ls1I = len(seq1);Ls2I = len(seq2)

	SeqGrepNumI = 20
	EndS = str(seq1[-SeqGrepNumI:])
	#HeadS = str(seq1[:SeqGrepNumI])

	EndpatternS = homopolymerRegex(EndS,relaxI,1).replace('\\',"")
	
	
	slideLimI = Ls1I - SeqGrepNumI;slideI = 1
	MatchL = re.findall(EndpatternS,seq2S)
	if len(MatchL) > 0:MatchS = sorted(MatchL, key = lambda MatchL:len(MatchL))[::-1][0]
	else:MatchS = ''
	while (len(MatchL) == 0 or MatchL.count('') >= len(MatchL) or MatchS == '') and slideI < slideLimI:
		
		TEMPEndS = seq1[-(SeqGrepNumI+slideI):-slideI] 
		EndpatternS = homopolymerRegex(TEMPEndS,relaxI,1).replace('\\',"")
		#STDERR("OverlapCheck.EndpatternS=",EndpatternS)
		MatchL = re.findall(EndpatternS,seq2S)
		if len(MatchL) > 0:
			MatchS = sorted(MatchL, key = lambda MatchL:len(MatchL))[::-1][0]
		slideI += 1
		#STDERR("WhileLoopPlus1")

	DEBUG = "EndpatternS in OverlapCheck while loop "+EndpatternS + str(re.findall(EndpatternS,seq2S)) + '\n' ##DEBUG
	#STDERR(DEBUG) ##DEBUG
	#STDERR('MatchS='+MatchS+'\n') ##DEBUG	

	#if MatchS != '':OverlapI = seq2S.find(MatchS) + len(MatchS) 
	if MatchS != '':OverlapI = seq2S.find(MatchS) + len(MatchS) + (slideI - 1)
	else:OverlapI = -1
 
	#STDERR("OverlapI="+str(OverlapI)+'\n') ##DEBUG
	#STDERR("slideI="+str(slideI)+'\n') ##DEBUG
	if OverlapI < 0:OverlapI = 0
	return OverlapI

def OverlapCheck2(seq1,seq2, *arg):## If seq1 tail overlap with seq2 head then return overlap length, else return 0
	#OverlapCheck2(seq1,seq2, relaxI, SeqGrepNumI,StepSizeI) 
	relaxI = 5
	SeqGrepNumI = 20
	StepSizeI = 5
	if len(arg) > 0:
		relaxI = arg[0]
	if len(arg) > 1:
		SeqGrepNumI = arg[1]
	if len(arg) > 2:
		StepSizeI = arg[2]
		
	seq1S = str(seq1[:]);seq2S = str(seq2[:])
	Ls1I = len(seq1);Ls2I = len(seq2)

	
	EndS = str(seq1[-SeqGrepNumI:])
	#HeadS = str(seq1[:SeqGrepNumI])

	#EndpatternS = homopolymerRegex(EndS,relaxI,1).replace('\\',"")
	EndpatternS = EndS 
	
	
	slideLimI = Ls1I - SeqGrepNumI;slideI = 1
	MatchL = re.findall(EndpatternS,seq2S)
	if len(MatchL) > 0:MatchS = sorted(MatchL, key = lambda MatchL:len(MatchL))[::-1][0]
	else:MatchS = ''
	SwitchFlagI = 1
	EndFlagI = 0
	TMPmatchS = ''
	while (len(MatchL) == 0 or MatchL.count('') >= len(MatchL) or MatchS == '') and slideI < slideLimI and SwitchFlagI != 0:
		
		TEMPEndS = seq1[-(SeqGrepNumI+slideI):-slideI] 
		#EndpatternS = homopolymerRegex(TEMPEndS,relaxI,1).replace('\\',"")
		EndpatternS = TEMPEndS 
		#STDERR("OverlapCheck.EndpatternS=",EndpatternS)
		MatchL = re.findall(EndpatternS,seq2S)

		if len(MatchL) > 0 and SwitchFlagI > 0:
			SwitchFlagI = StepSizeI*-1
			TMPmatchS = sorted(MatchL, key = lambda MatchL:len(MatchL))[::-1][0]

		elif len(MatchL) > 0 and SwitchFlagI >= 0:
			MatchS = sorted(MatchL, key = lambda MatchL:len(MatchL))[::-1][0]
		
		if SwitchFlagI > 0:
			slideI += StepSizeI

		elif SwitchFlagI < 0:
			slideI += -1
			SwitchFlagI += 1
		elif SwitchFlagI == 0:
			slideI += 1
		
		#STDERR("WhileLoopPlus1")
		#STDERR("SwitchFlagI=",SwitchFlagI,"slideI=",slideI)
		#STDERR("MatchL=",MatchL)	
	#STDERR("WhileLoop EXIT")

	#DEBUG = "EndpatternS in OverlapCheck while loop "+EndpatternS + str(re.findall(EndpatternS,seq2S)) + '\n' ##DEBUG
	#STDERR(DEBUG) ##DEBUG
	#STDERR('MatchS='+MatchS+'\n') ##DEBUG	

	#if MatchS != '':OverlapI = seq2S.find(MatchS) + len(MatchS) 
	if MatchS != '':
		OverlapI = seq2S.find(MatchS) + len(MatchS) #- (slideI - 1)
		#STDERR("seq2S.find(MatchS)=",seq2S.find(MatchS))
		#STDERR("len(MatchS)=",len(MatchS))
		#STDERR("(slideI - 1)",(slideI - 1))
	elif TMPmatchS != '':
		OverlapI = seq2S.find(TMPmatchS) + len(TMPmatchS) #- (slideI - 1)
				
	else:OverlapI = -1
	 
 
	#STDERR("OverlapI="+str(OverlapI)+'\n') ##DEBUG
	#STDERR("slideI="+str(slideI)+'\n') ##DEBUG
	#if OverlapI < 0:OverlapI = 0
	return OverlapI

def LongerSeqS4(seq1,seq2,PI):
	#PI #piority 0=relies on first sequence,1 second sequence
	seq1S = seq1[:];seq2S = seq2[:]
	#DEBUG = "seq1S = "+seq1+"\nseq2S = "+ seq2 + "\n" ##DEBUG
	#sys.stderr.write(DEBUG) ##DEBUG
	
	Ls1I = len(seq1);Ls2I = len(seq2)

	OverlapI = OverlapCheck(seq1S,seq2S)

	sys.stderr.write("OverlapI="+str(OverlapI)+'\n') #DEBUG	

	GapPlus = 0
	seq2S_remains = ''
	seq1S_remains = ''

	if (Ls1I - OverlapI) < len(seq1S):
		seq1S = seq1[OverlapI*(-1):]
		#print "seq1 longer than OverlapI" ##DEBUG
		seq1S_remains = seq1[:OverlapI*(-1)]
				
	if  (Ls2I -  OverlapI) < len(seq2S):
		seq2S = seq2[:OverlapI]
		#GapPlus = len(seq2) - len(seq2S)
		seq2S_remains = seq2[OverlapI:]
		#print "seq2 longer than seq1" ##DEBUG
	
	alignInS = ">S1S\n"+seq1S+"\n>S2S\n"+seq2S+"\n"

	#print "alignInS\n",alignInS ##DEBUG
	arg = "muscle -maxiters 32 -gapopen -2400 -quiet"
	arg = "muscle -maxiters 32 -quiet"
	#print "this is alignInS",alignInS #DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(alignInS)[0]
	#print "this is x ",x #DEBUG
	xL = x.split('>')
	xL = [x for x in xL if len(x) > 0]
	#print "This is xL " ,xL #DEBUG
	del process
	HeadL = [x for x in xL if x.find("S1S") != -1]

	#print "This is HeadL " ,HeadL #DEBUG
	HeadL = HeadL[0].splitlines()[1:]
	HeadL = [x.strip() for x in HeadL]
	HeadS = ''.join(HeadL)
#	print "This is HeadS " ,HeadS #DEBUG
	TailL = [x for x in xL if x.find("S2S") != -1]
	TailL = TailL[0].splitlines()[1:]
	TailL = [x.strip() for x in TailL]
	TailS = ''.join(TailL)
#	print "This is TailS " ,TailS #DEBUG
	HGapNumI = HeadS.count('-')
	TGapNumI = TailS.count('-')
#	print "GapNumI=",GapNumI ##DEBUG
	#InnerGapI = HeadS[0:(-1)*(len(HeadS) - GapNumI)].count('-')
	HInnerGapI = HeadS.rstrip('-').count('-')
	TInnerGapI = TailS.lstrip('-').count('-')
	HGapNumI2 = HGapNumI - HInnerGapI 
	TGapNumI2 = TGapNumI - TInnerGapI 

	#print "GapNumI2=",GapNumI2 ##DEBUG
	#AlignoutS = HeadS[:len(HeadS) - GapNumI2] + TailS[len(TailS) - GapNumI2:]

	
	if PI == 0 and HGapNumI2 == 0:
		AlignoutS = seq1 + seq2S_remains +"\n"
	elif PI == 0 and HGapNumI2 > 0:
		AlignoutS = seq1 + TailS[(-1)*HGapNumI2:] + seq2S_remains +"\n"

	elif PI == 1 and TGapNumI2 == 0:
		AlignoutS = seq1S_remains + seq2 +"\n"
	elif PI == 1 and TGapNumI2 > 0:
		AlignoutS = seq1S_remains + HeadS[:HGapNumI2] + seq2 +"\n"
	#print "HeadS[:]\n",HeadS[:] ##DEBUG
	#print "TailS[(-1)*GapNumI2:]\n" ,TailS[(-1)*GapNumI2:] #DEBUG
	#print "TailS[:]\n" ,TailS[:] #DEBUG
	#print "AlignoutS\n" ,ManyLines(AlignoutS,80) #DEBUG
	else:
		print("ERR at LongerSeq")
	return AlignoutS
	

def LongerSeqSmain(seqs,directions,PflagI):
	#directions = 3

	seqsL = seqs.split('>')
	seqsL = ['>'+x for x in seqsL if len(x) > 0]#[:5]
	#print seqsL ##DEBUG

	seq1 = seqsL.pop(0)
	while len(seqsL) > 0:
		seq1S = seqonly(seq1)
		#print "seq1S###\n",seq1S  ##DEBUG
		seq2S =  seqonly(seqsL.pop(0))
		#print "seq2S###\n",seq2S  ##DEBUG
		if directions == '3':
			seq = LongerSeqS4(seq1S,seq2S,PflagI)
		else:
			seq = LongerSeqS4(seq2S,seq1S,PflagI)
		
		#print "###\n",seq  ##DEBUG
		seq1 = ">LongerSeq\n" + seq
	return seq1
	#return  ManyLines(seq1,80)

def LongerSeqS2(seq1,seq2,PI):
	#PI #piority 0=relies on first sequence,1 second sequence
	seq1S = seq1;seq2S = seq2
	DEBUG = "seq1S = "+seq1+"\nseq2S = "+ seq2 + "\n" ##DEBUG
	#sys.stderr.write(DEBUG) ##DEBUG
	
	Ls1I = len(seq1);Ls2I = len(seq2)

	SeqGrepNumI = 20
	EndS = seq1S[-(SeqGrepNumI-1):]
	EndS2 = seq2S[1:SeqGrepNumI]
	if seq2S.find(EndS) != -1:
		##                               | EndS  |
		##seq1  |------------------------|-------|
		##seq2                     |------------------------------
		OverlapI = seq2S.find(EndS) + SeqGrepNumI

	elif seq1S.rfind(EndS2) != -1: 
		##seq1  |--------------------------------|
		##seq2                     |--------|---------------------
		##                         | EndS2  |
		OverlapI = len(seq1) - seq1S.rfind(EndS2) 
	else:
		sys.stderr.write("OverlapI ERROR AT "+EndS+ " | " +EndS2+'\n') #DEBUG	

	#if OverlapI > len(seq1) or OverlapI > len(seq2):
	#	del OverlapI
	#	sys.stderr.write("OverlapI ERROR AT "+EndS+'\n') #DEBUG

	sys.stderr.write("OverlapI="+str(OverlapI)+'\n') #DEBUG	

	GapPlus = 0
	seq2S_remains = ''
	seq1S_remains = ''

	if (Ls1I - OverlapI) > 20:
		seq1S = seq1[OverlapI*(-1):]
		#print "seq1 longer than OverlapI" ##DEBUG
		seq1S_remains = seq1[:OverlapI*(-1)]
				
	if  (Ls2I -  OverlapI) > 20:
		seq2S = seq2[:OverlapI]
		#GapPlus = len(seq2) - len(seq2S)
		seq2S_remains = seq2[OverlapI:]
		#print "seq2 longer than seq1" ##DEBUG
	
	alignInS = ">S1S\n"+seq1S+"\n>S2S\n"+seq2S+"\n"

	#print "alignInS\n",alignInS ##DEBUG
	arg = "muscle -maxiters 32 -gapopen -2400 -quiet"
	arg = "muscle -maxiters 32 -quiet"
	#print "this is alignInS",alignInS #DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(alignInS)[0]
	#print "this is x ",x #DEBUG
	xL = x.split('>')
	xL = [x for x in xL if len(x) > 0]
	#print "This is xL " ,xL #DEBUG
	del process
	HeadL = [x for x in xL if x.find("S1S") != -1]

	#print "This is HeadL " ,HeadL #DEBUG
	HeadL = HeadL[0].splitlines()[1:]
	HeadL = [x.strip() for x in HeadL]
	HeadS = ''.join(HeadL)
#	print "This is HeadS " ,HeadS #DEBUG
	TailL = [x for x in xL if x.find("S2S") != -1]
	TailL = TailL[0].splitlines()[1:]
	TailL = [x.strip() for x in TailL]
	TailS = ''.join(TailL)
#	print "This is TailS " ,TailS #DEBUG
	HGapNumI = HeadS.count('-')
	TGapNumI = TailS.count('-')
#	print "GapNumI=",GapNumI ##DEBUG
	#InnerGapI = HeadS[0:(-1)*(len(HeadS) - GapNumI)].count('-')
	HInnerGapI = HeadS.rstrip('-').count('-')
	TInnerGapI = TailS.lstrip('-').count('-')
	HGapNumI2 = HGapNumI - HInnerGapI 
	TGapNumI2 = TGapNumI - TInnerGapI 

	#print "GapNumI2=",GapNumI2 ##DEBUG
	#AlignoutS = HeadS[:len(HeadS) - GapNumI2] + TailS[len(TailS) - GapNumI2:]

	
	if PI == 0 and HGapNumI2 == 0:
		AlignoutS = seq1 + seq2S_remains +"\n"
	elif PI == 0 and HGapNumI2 > 0:
		AlignoutS = seq1 + TailS[(-1)*HGapNumI2:] + seq2S_remains +"\n"

	elif PI == 1 and TGapNumI2 == 0:
		AlignoutS = seq1S_remains + seq2 +"\n"
	elif PI == 1 and TGapNumI2 > 0:
		AlignoutS = seq1S_remains + HeadS[:HGapNumI2] + seq2 +"\n"
	#print "HeadS[:]\n",HeadS[:] ##DEBUG
	#print "TailS[(-1)*GapNumI2:]\n" ,TailS[(-1)*GapNumI2:] #DEBUG
	#print "TailS[:]\n" ,TailS[:] #DEBUG
	#print "AlignoutS\n" ,ManyLines(AlignoutS,80) #DEBUG
	else:
		print("ERR at LongerSeq")
	return AlignoutS

def seqonly(fastaS):
	fastaL = fastaS.splitlines()
	fastaL = [x for x in fastaL if x.find(">") != 0]
	fastaL = [x.strip() for x in fastaL]
	fastaS2 = ''.join(fastaL)
	return fastaS2

def reversecomplement(string):
        reverse = string[::-1]
        
        reverse = reverse.replace("A","?")
        reverse = reverse.replace("T","A")
        reverse = reverse.replace("?","T")
        reverse = reverse.replace("G","?")
        reverse = reverse.replace("C","G")
        reverse = reverse.replace("?","C")
        return reverse ## complementary string ATCG --> CGAT


def ReadGrep3(seqS,ReadsDB,greplenI,DirectionS,RegexI):
	sys.stderr.write('Entering ReadGrep3 \n') ##DEBUG
	OriginalSeqS = seqS
	if len(seqS) > greplenI and DirectionS == '3':
		seqS = seqonly(seqS)[-greplenI:]

	elif len(seqS) > greplenI and DirectionS == '5':
		seqS = seqonly(seqS)[:greplenI]

	seqSmod = seqS.strip()[:-1]
	if RegexI > 0:seqSmod = homopolymerRegex(seqSmod,RegexI,2)
	grepout =  RegexSelect(OriginalSeqS,ReadsDB,DirectionS,seqSmod)
	
	seqS = reversecomplement(seqS)
	RseqSmod = seqS.strip()[1:]
	if RegexI > 0:RseqSmod = homopolymerRegex(RseqSmod,RegexI,2)
	minusgrepout = RegexSelect(OriginalSeqS,ReadsDB,DirectionS,RseqSmod)

	minusgrepoutreverse = ''
	for i in minusgrepout.splitlines():
		if i.find('>') != -1:
			minusgrepoutreverse = minusgrepoutreverse + i + '\n'
		else:
			minusgrepoutreverse = minusgrepoutreverse + reversecomplement(i.strip()) + '\n'

	grepout = grepout  + minusgrepoutreverse
 
	#sys.stderr.write("grepout "+grepout+'\n') ##DEBUG

	return grepout

def RegexSelect(OriginalSeqS,ReadsDB,DirectionS,seqSmod):
	seqS = seqonly(OriginalSeqS)
	sys.stderr.write('Entering RegexSelect \n') ##DEBUG
	#arg = """grep -o '""" + seqSmod + """' """ + ReadsDB + """ | grep -v '^-' | sort -u """
	arg = """grep -o '""" + seqSmod + """' | grep -v '^-' | sort -u """ ## USE READs data from memory to reduce disk IO 

	sys.stderr.write(arg + '\n')  ##DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(GB_MemReadS)
	patternS = str(x[0])
	sys.stderr.write('patternS=' + patternS + '\n') ##DEBUG

	del process
	patternL = patternS.splitlines()
	patternL = [x for x in patternL if len(x) > 0]
	
	if len(seqS) > 200:EndSeqSI = 200
	else:EndSeqSI = len(seqS)

	SelectgrepoutS = ''
	for SubSeqS in patternL:
		#grepoutS = RegexDB(ReadsDB,SubSeqS) ##Alter to parallel version
		grepoutS = grepCall(GB_MemReadS,SubSeqS)
		
		#STDERR("THISI PARALLEL INPUT"+str([[x,SubSeqS] for x in MemReadSL][:10]))
		#grepoutL = ParallelLoop(grepCall,[[x,SubSeqS] for x in MemReadSL]) 
		#grepoutS = '\n'.join(grepoutL)
		#sys.stderr.write('RegexSelect grepoutS= \n'+grepoutS+'\n') ##DEBUG
		
		grepoutL = grepoutS.split('>') 
		grepoutL = ['>'+x for x in grepoutL if x != '']
		
		#sys.stderr.write('RegexSelect grepoutL= \n'+str(grepoutL)+'\n') ##DEBUG
		if DirectionS == '3':
			grepoutL1 = [x for x in grepoutL if len(seqonly(x).split(SubSeqS)[1]) > 5 and len(seqonly(x).split(SubSeqS)[0]) < 20] 
			#grepoutL2 = [x for x in grepoutL if len(seqonly(x).split(SubSeqS)[1]) > 5 and len(seqonly(x).split(SubSeqS)[0]) > 20 and RegexCheck(seqS[-EndSeqSI:],seqonly(x).split(SubSeqS)[0][:20]) == 1]
			grepoutL2 = [x for x in grepoutL if len(seqonly(x).split(SubSeqS)[1]) > 5 and len(seqonly(x).split(SubSeqS)[0]) > 20]
			RegexCheckInputL = [ [seqS[-EndSeqSI:],seqonly(x).split(SubSeqS)[0][:20]] for x in grepoutL2 ]
			RegexCheckOutputL = ParallelLoop(RegexCheck,RegexCheckInputL)
			grepoutL2 = zip(grepoutL2,RegexCheckOutputL)
			grepoutL2 = [x[0] for x in grepoutL2 if x[1] == 1]

		if DirectionS == '5':
			grepoutL1 = [x for x in grepoutL if len(seqonly(x).split(SubSeqS)[0]) > 5 and len(seqonly(x).split(SubSeqS)[1]) < 20] 
			#grepoutL2 = [x for x in grepoutL if len(seqonly(x).split(SubSeqS)[0]) > 5 and len(seqonly(x).split(SubSeqS)[1]) > 20 and RegexCheck(seqS[:EndSeqSI],seqonly(x).split(SubSeqS)[1][-20:]) == 1]

			grepoutL2 = [x for x in grepoutL if len(seqonly(x).split(SubSeqS)[1]) > 5 and len(seqonly(x).split(SubSeqS)[0]) > 20]
			RegexCheckInputL = [ [seqS[:EndSeqSI],seqonly(x).split(SubSeqS)[1][-20:]] for x in grepoutL2 ]
			RegexCheckOutputL = ParallelLoop(RegexCheck,RegexCheckInputL)
			grepoutL2 = zip(grepoutL2,RegexCheckOutputL)
			grepoutL2 = [x[0] for x in grepoutL2 if x[1] == 1]


		#sys.stderr.write('RegexSelect grepoutL1 =' +str(grepoutL1)+ '\n') ##DEBUG
		#sys.stderr.write('RegexSelect grepoutL2 =' +str(grepoutL2)+ '\n') ##DEBUG
		SelectgrepoutL = grepoutL1 + grepoutL2
		SelectgrepoutS = SelectgrepoutS + '\n'.join(SelectgrepoutL) 

	#sys.stderr.write('RegexSelect SelectgrepoutS =' +SelectgrepoutS+ '\n') ##DEBUG
	#sys.stderr.write('Ending RegexSelect \n') ##DEBUG
	return SelectgrepoutS

def RegexCheck(seqS,grepoutS):
	#sys.stderr.write('Entering RegexCheck \n') ##DEBUG
	seqSregexS = homopolymerRegex(grepoutS,5,2)
	grepout = RegexCall(seqS,seqSregexS)

	if len(grepout) > 0:checkI = 1
	else:checkI = 0
	return checkI

def SubAlignCheck(seqIn1S,seqIn2S): #check if seqIn2S is a substring of seqIn1S
	#alignlenI = min([len(seqonly(seqIn1S)),len(seqonly(seqIn2S))])
	#seq1S =seqonly(seqIn1S)[-alignlenI:];seq2S =seqonly(seqIn2S)[:alignlenI]

	seq1S =seqonly(seqIn1S);seq2S =seqonly(seqIn2S)
	SeqS = ">S1S\n"+seq1S+"\n>S2S\n"+seq2S+"\n"
	sys.stderr.write('Entering SubAlignCheck \n') ##DEBUG

	AlignOutS = musclecallGap(SeqS)
	#AlignOutS = musclecall(SeqS)
	sys.stderr.write('AlignOutS in SubAlignCheck ='+AlignOutS+'\n') ##DEBUG

	if AlignOutS.count('>') < 2:
		flagI = 1

	elif AlignOutS.count('>') == 2:

		aling1S = seqonly('>'+AlignOutS.split('>')[1])
		aling2S = seqonly('>'+AlignOutS.split('>')[2])

		if aling1S[-3:] != '---' and len(aling1S[:-3].rstrip('-')) < len(aling1S.rstrip('-')):
			aling1S = seqonly('>'+AlignOutS.split('>')[1])[:-3]
		if aling2S[:3] != '---' and len(aling2S[3:].lstrip('-')) < len(aling2S.lstrip('-')):
			aling2S = seqonly('>'+AlignOutS.split('>')[2])[3:]

		sys.stderr.write('aling1S in SubAlignCheck ='+aling1S+'\n') ##DEBUG	
		sys.stderr.write('aling2S in SubAlignCheck ='+aling2S+'\n') ##DEBUG	

		if len(aling1S) - len(aling1S.lstrip('-')) > 0 and len(aling2S) - len(aling2S.rstrip('-')) > 0:flagI = -1
		elif len(aling1S) - len(aling1S.rstrip('-')) > 0 and len(aling2S) - len(aling2S.lstrip('-')) > 0:flagI = -1
		else:flagI = 1

	return flagI
	

def AlignCheck(seqIn1S,seqIn2S,OverlapLI):

	alignlenI = min([len(seqonly(seqIn1S)),len(seqonly(seqIn2S))])
	seq1S =seqonly(seqIn1S)[-alignlenI:];seq2S =seqonly(seqIn2S)[:alignlenI]

	SeqS = ">S1S\n"+seq1S+"\n>S2S\n"+seq2S+"\n"
	#sys.stderr.write('Entering AlignCheck \n') ##DEBUG

	AlignOutS = musclecall(SeqS)
	#sys.stderr.write('AlignOutS in AlignCheck ='+AlignOutS+'\n') ##DEBUG
	
	if AlignOutS.count('>') < 2:
		flagI = -1
	elif AlignOutS.count('>') == 2 and seqonly('>'+AlignOutS.split('>')[0]).strip('-').find('-----') == -1 and seqonly('>'+AlignOutS.split('>')[1]).strip('-').find('-----') == -1 and seqonly('>'+AlignOutS.split('>')[1]).strip('-').count('-') < len(seq2S)*0.10 and seqonly('>'+AlignOutS.split('>')[0]).strip('-').count('-') < len(seq1S)*0.10:

		if len(consensusExtract2(AlignOutS,1.0,1.0)) >= OverlapLI:flagI = 1
		else:flagI = -1
	else:
		flagI = -1
	#sys.stderr.write('flagI in AlignCheck ='+str(flagI)+'\n') ##DEBUG
	return flagI
	

def RegexCall(seqS,SubSeqS):
	#sys.stderr.write('Entering RegexCall \n') ##DEBUG
	arg = """grep '""" + SubSeqS + """' | grep -v '^-'	"""
	#arg = "echo " + seqS
	#sys.stderr.write(arg + '\n')  ##DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(seqS)
	grepoutS = str(x[0])
	del process
	#sys.stderr.write('grepoutS RegexCall =' +grepoutS+ '\n') ##DEBUG
	return grepoutS

def RegexDB(DB,SubSeqS):
	sys.stderr.write('Entering RegexDB \n') ##DEBUG
	arg = """grep -B 1 '""" + SubSeqS + """' """ +DB+ """ | grep -v '^-'	"""

	sys.stderr.write(arg + '\n')  ##DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	grepoutS = str(x[0])
	del process
	#sys.stderr.write('grepoutS RegexDB =' +grepoutS+ '\n') ##DEBUG
	return grepoutS

def grepCall(SeqS,SubSeqS):
	#sys.stderr.write('Entering RegexDB \n') ##DEBUG
	arg = """ grep -B 1 '""" + SubSeqS + """' | grep -v '^-' """

	#sys.stderr.write(arg + '\n')  ##DEBUG
	process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate(SeqS)
	grepoutS = str(x[0])
	del process
	#sys.stderr.write('grepoutS RegexDB =' +grepoutS+ '\n') ##DEBUG
	return grepoutS


def homopolymerRegex(OriseqS,HomoNumI,minimumrepeatI):
	seqS = OriseqS
	BaseL = [] 	
	for j in ('A','T','C','G'):
		i=30
		foundIndexI = '0'
		while i > 0:
			i -=1
			if seqS.find(i*j) != -1 and seqS.find(i*j) != foundIndexI:
			
				BaseL.append([j, i, seqS.find(i*j), seqS.count(i*j)])
				
				foundIndexI = seqS.find(i*j)
				seqS = seqS.replace(i*j,i*"?")
	
	PatternL = [x for x in BaseL if x[1] > minimumrepeatI]	
	#DEBUG = "this is PatternL in homopolymerRegex ="+str(PatternL)+'\n' ##DEBUG
	#sys.stderr.write(DEBUG) ##DEBUG
	if HomoNumI > len(PatternL):HomoNumI = len(PatternL)
	#print "this is HomoNumI",str(HomoNumI) ##DEBUG

	selectL = sorted(PatternL, key = lambda PatternL:PatternL[1] )[::-1]
	selectL = selectL[:HomoNumI]
	#DEBUG = "this is selectL in homopolymerRegex ="+str(selectL)+'\n' ##DEBUG
	#sys.stderr.write(DEBUG) ##DEBUG

	PatternS = OriseqS
	for i in selectL:
		regaex = i[0]+'\{1,\}'
		PatternS = PatternS.replace((i[0]*i[1]),regaex)

	return PatternS

def ManyLines(seqS,widthI):
	
	seqLineI = len(seqS)/widthI
	remainsI = len(seqS)%widthI
	TEMPseqS = '>LongSeq\n'
	for i in range(seqLineI):
		TEMPseqS = TEMPseqS + seqS[i*(widthI):(i+1)*(widthI)] + '\n'

	TEMPseqS = TEMPseqS + seqS[(-1)*(remainsI):] + '\n'
	return TEMPseqS

def mainBatch(SseqS,DirectionS,IterationI,RdbS,CdbS,OutS,NameS,StartNum,ExcludedL):
   
	f = open(OutS,'a')
	f.close()

	IterationIcount = StartNum
	#SseqF = open(SseqS)
	TEMPbridgeS1 = SseqS
	#SseqF.close()

	PflagI = 0 #0=contig 1=Rreads consensus
	ReadPflagI = 0
	consensus = ''
	PreviousConsensusS = ''
	SubTempS = ''
	FLAG = 0
	NextContigS = ''
	#banReadL = []
	
	SubExcludedL = []	
	while IterationIcount < IterationI+1:
		FLAG = 0
		sys.stderr.write("Entering iter "+str(IterationIcount)+'\n') ##DEBUG	
		blastResultS = BlastCall2(TEMPbridgeS1,CdbS)
		#blastResultS = BlastCall3(TEMPbridgeS1,CdbS,DirectionS)
		blastResultL = blastTrans(blastResultS,1e-08)
		for i in blastResultL:sys.stderr.write(str(i)+'\n') ##DEBUG

		blastResultL = [x for x in blastResultL if x[1] not in ExcludedL]
		SubExcludedL = [x[1] for x in FullOverlap(TEMPbridgeS1,blastResultL)]
		ExcludedL = ExcludedL + SubExcludedL
		
		unitigBlast = unitigSelect2(TEMPbridgeS1,blastResultL,DirectionS)
		sys.stderr.write("len(seqonly(TEMPbridgeS1))="+str(len(seqonly(TEMPbridgeS1)))+'\n') ##DEBUG
		sys.stderr.write("UNITIG BLAST\n"+str(unitigBlast)+'\n') ##DEBUG

		
		SelectedReadsL = []
		if len(unitigBlast) > 0:
			f = open(OutS,'a')
			f.write("#"+'\t'.join([str(x) for x in unitigBlast])+'\n')
			f.close()

			FastacmdSeqName = unitigBlast[1]
			if unitigBlast[8] < unitigBlast[9]:
				NextContigS = getSseq(CdbS,FastacmdSeqName,'1')
			else:
				NextContigS = getSseq(CdbS,FastacmdSeqName,'2')
			#print NextContigS ##DEBUG
			#TempSeqDomino = 'TemSeqDomino'+IterationIcount
			f = open(OutS,'a')
			f.write(NextContigS)
			f.close()
			if ReadPflagI == 0: 
				PflagI = 0 #piorityflag

			else:
				PflagI = 1
			ReadPflagI = 0
			#LongerNextContigS = LongerSeqS2(seqonly(TEMPbridgeS1),seqonly(NextContigS),PflagI)
			sys.stderr.write("LongerSeqSmain is called in FIRST IF\n") ##DEBUG
			LongerNextContigS = LongerSeqSmain(TEMPbridgeS1+'\n'+NextContigS,DirectionS,PflagI)

			#sys.stderr.write("TEMPbridgeS1\n"+TEMPbridgeS1+'\n') ##DEBUG
			#sys.stderr.write("NextContigS\n"+NextContigS+'\n') ##DEBUG
			#print "LongerNextContigS\n",LongerNextContigS ##DEBUG
			TEMPbridgeS1 = LongerNextContigS
			ExcludedL.append(str(FastacmdSeqName))
						
		elif len(unitigBlast) < 1:
			blastResultS = BlastCall2(TEMPbridgeS1,RdbS)
			#blastResultS = BlastCall3(TEMPbridgeS1,RdbS,DirectionS)

			blastResultL = blastTrans(blastResultS,1e-08)
			sys.stderr.write("THIS IS blastResultL\n"+str(blastResultL)+'\n') ##DEBUG			
			#sys.stderr.write("TEMPbridgeS1\n"+TEMPbridgeS1+'\n') ##DEBUG
			SelectedReadsL = ReadsSelelct(TEMPbridgeS1,RdbS,blastResultL,DirectionS) 
			#sys.stderr.write("THIS IS ReadsSelelct(TEMPbridgeS1,RdbS,blastResultL,DirectionS) \n"+ str(SelectedReadsL) + '\n') ##DEBUG
			
			if len(SelectedReadsL) >= ReadGrepNum:
				f = open(OutS,'a')
				f.write("#"+'\n#'.join([str(x) for x in SelectedReadsL])+'\n')
				f.close()
					
				Reads = ''	
				for i in SelectedReadsL:
					if i[8] > i[9]:Reads = Reads + getSseq(RdbS,i[1],'2')
					else:Reads = Reads + getSseq(RdbS,i[1],'1')
			
				#print "READ for musclecall",Reads ##DEBUG	
				if Reads.count('>') > 1:	
					AlingINS = musclecall(Reads)
					#sys.stderr.write("AlingINS ="+ AlingINS + '\n')##DEBUG
					consensus =  consensusExtract2(AlingINS,1.0,0.25)
				elif Reads.count('>') == 1:
					consensus = seqonly(Reads)

				#consensus = F3primeReadCall(blastoutnameS,Tag,seqS,RdbS)
				if len(consensus) > 1:

					#Br = Tag+".fa"
					ReadsnumI = str(Reads.count('>'))
					f = open(OutS,'a')
					f.write("#ReadsNumber="+ReadsnumI+"\n")
					
					#f.write(">"+Tag)
					NextContigS = ">"+ NameS + str(IterationIcount) + "\n" + consensus +  "\n" 
					f.write(NextContigS)
					f.close()
					
					if ReadPflagI == 0: 
						PflagI = 0 #piorityflag

					elif ReadPflagI != 0 and ReadPflagI < ReadsnumI:
						PflagI = 0
				
					elif ReadPflagI != 0 and ReadPflagI > ReadsnumI:
						PflagI = 1
					ReadPflagI = ReadsnumI

					#LongerNextContigS = LongerSeqS2(seqonly(TEMPbridgeS1),seqonly(NextContigS),PflagI)
					sys.stderr.write("LongerSeqSmain is called in elif unitigblast\n") ##DEBUG
					LongerNextContigS = LongerSeqSmain(TEMPbridgeS1+'\n'+NextContigS,DirectionS,PflagI)
					TEMPbridgeS1 = LongerNextContigS


		if len(unitigBlast) < 1 and len(SelectedReadsL) < ReadGrepNum:
			#sys.stderr.write("if len(unitigBlast) < 1 and len(SelectedReadsL) < " + str(ReadGrepNum) + '\n') ##DEBUG
			#sys.stderr.write(str(len(unitigBlast))+ " " + str(len(SelectedReadsL)) + '\n') ##DEBUG

			#ReadS = ReadGrep(TEMPbridgeS1,RdbS,greplenI,DirectionS)
			#ReadS = ReadGrep3(TEMPbridgeS1,RdbS,50,DirectionS,0)
			#sys.stderr.write("ReadS.count('>')=" + str(ReadS.count('>')) + '\n')	##DEBUG	
			#if ReadS.count('>') < minimumReadI and greplenI < 25:
			#	TrimlenI = 1
			#	#SubTempS = seqonly(TEMPbridgeS1)
								
			#else:TrimlenI = 12 
			AlingINS = ''
			#minimumReadI = 5
			minimumReadI = ReadGrepNum

			ReadS = ''
			#while TrimlenI < 11 :
			SubTempS = TEMPbridgeS1
			TrimlenI = 1
			while ReadS.count('>') < minimumReadI and TrimlenI < 31 :  ##Trim last base when grep cannot find any read 
				
				greplenI = 50
				while ReadS.count('>') < minimumReadI and greplenI > 31:

					RegexI = 0
					while ReadS.count('>') < minimumReadI and RegexI < 3:
						
						ReadS = ReadGrep3(SubTempS,RdbS,greplenI,DirectionS,RegexI) 
						EachReadL = [x for x in ReadS.split('>') if len(x) > 0]
						InputL = [[SubTempS,'>'+x,greplenI] for x in EachReadL] 
						AlignCheckL = ParallelLoop(AlignCheck,InputL)
						CheckedReadL = zip(EachReadL,AlignCheckL)
						
						#ReadL = ['>'+x for x in ReadS.split('>') if len(x) > 0 and AlignCheck(SubTempS,'>'+x,greplenI) == 1] 
						ReadL = ['>'+x[0] for x in CheckedReadL if x[1] == 1] 

			
						
						sys.stderr.write("minimumReadI="+str(minimumReadI) + "TrimlenI="+str(TrimlenI) + " greplenI="+  str(greplenI)+ " RegexI=" + str(RegexI) + "\n")	##DEBUG	
						#sys.stderr.write("SubTempS= \n"+SubTempS + '\n')	##DEBUG
						#sys.stderr.write("grep seq = SubTempS[-greplenI:]= \n"+SubTempS + '\n')	##DEBUG

						RegexI = RegexI + 1

					greplenI = greplenI - 2 ##Reduce greplenI by 2 if grep cannot find any read

				if DirectionS == '3':
					SubTempS = seqonly(TEMPbridgeS1)[:-TrimlenI]
				if DirectionS == '5':
					SubTempS = seqonly(TEMPbridgeS1)[TrimlenI:]
	
				TrimlenI += 5

			if TrimlenI < 31 and ReadS.count('>') >= minimumReadI: ##check if SubTempS can be used as TEMPbridgeS1
			#if ReadS.count('>') >= minimumReadI:
				TEMPbridgeS1 = ">LonSeq\n"+SubTempS
	
			#AlingIN2S = ''
			if ReadS.count('>') >= minimumReadI:
				
				AlingINS = musclecall(ReadS)
				#AlingIN2S = AlingMod(AlingINS,DirectionS,5)

				if len(AlingINS) < 10:sys.stderr.write("len(AlingINS) < 10\n") ##DEBUG
				#sys.stderr.write("Entering Alignmod \n") ##DEBUG
				sys.stderr.write("this is AlingINS \n"+AlingINS+'\n') ##DEBUG
				

			if len(AlingINS) > 10:
				#consensus =  consensusExtract2(AlingINS,1.0,0.25)
				consensus = realign(AlingINS)
				sys.stderr.write("this is consensus \n"+consensus+'\n') ##DEBUG
				
			else:
				consensus = ''			
				
			if len(consensus) > 10:
		
				LongerNextContigS = 0
				
				NextContigS = ">"+ NameS + str(IterationIcount) + "GrepRead\n" + consensus +  "\n" 

				#sys.stderr.write("this is NextContigS in ReadGrep loop "+NextContigS+'\n') ##DEBUG
				#sys.stderr.write("LongerSeqSmain is called After ReadGrep3\n") ##DEBUG
				#sys.stderr.write("LongerSeqSmain seq1 ReadGrep3 "+TEMPbridgeS1+"\n") ##DEBUG
				#sys.stderr.write("LongerSeqSmain seq2 ReadGrep3 "+NextContigS+"\n") ##DEBUG
				
				LongerNextContigS = LongerSeqSmain(TEMPbridgeS1+'\n'+NextContigS,DirectionS,0)

				#sys.stderr.write("this is LongerNextContigS in ReadGrep loop "+LongerNextContigS+'\n') ##DEBUG
				if IterationIcount == 0:
					LongerCheckI = 10
				else:
					FlankingS = LongerSeqSmain(PreviousConsensusS+NextContigS,DirectionS,0)
					LongerCheckI = len(seqonly(FlankingS)) - len(seqonly(NextContigS)) 
					STDERR('FlankingS=' + str(FlankingS)) ##DEBUG
					STDERR('NextContigS=' + str(NextContigS)) ##DEBUG
					STDERR('PreviousConsensusS=' + str(PreviousConsensusS)) ##DEBUG
					STDERR('LongerCheckI=' + str(LongerCheckI)) ##DEBUG

				if (len(seqonly(TEMPbridgeS1)) + 5  < len(seqonly(LongerNextContigS))) and (LongerCheckI > 5):
					f = open(OutS,'a')
					f.write(NextContigS)
					f.close()

					TEMPbridgeS1 = LongerNextContigS
				#TEMPbridgeS1 = LongerNextContigS
				else:
					FLAG = 1

		if len(unitigBlast) < 1 and len(SelectedReadsL) < ReadGrepNum and consensus == '' or  FLAG == 1 :

			sys.stderr.write("end at iter "+str(IterationIcount)) ##DEBUG
			IterationIcount = IterationI+1
				
		sys.stderr.write("TEMPbridgeS1 for IterationIcount="+str(IterationIcount)+'\n') ##DEBUG
		if len(TEMPbridgeS1) < 200:sys.stderr.write(TEMPbridgeS1+'\n') ##DEBUG
		else:sys.stderr.write(TEMPbridgeS1[:100]+"..."+ TEMPbridgeS1[-50:] +'\n') ##DEBUG
		IterationIcount += 1
		#FLAG += 1
		PreviousConsensusS = str(NextContigS)
	LastCONTIG = ManyLines(seqonly(TEMPbridgeS1),80)

	#f = open(OutS,'a')
	#f.write(LastCONTIG)
	#f.close()
	return [LastCONTIG,ExcludedL]

def batchLoop(ContignamesF):
	f = open(ContignamesF,'r')
	nameL = [x.split()[0].replace('>','') for x in f.readlines() if x.find('>') == 0]
	f.close()

	nameL = [x.split('_') for x in nameL]
	#sys.stderr.write("nameL[:10] before sorted="+str(nameL[:10])+'\n') ##DEBUG
	nameL2 = sorted(nameL, key = lambda nameL:int(nameL[1]))[::-1]
	nameL = ['_'.join(x) for x in nameL2]
		
	#sys.stderr.write("nameL[:10] After sorted="+str(nameL[:10])+'\n') ##DEBUG
	if ExcludedS == '0':ExcludedL = []
	else:
		ExcludedF = open(ExcludedS,'r')
		ExcludedL = ExcludedF.readlines()
		ExcludedF.close()
		ExcludedL = [x.strip() for x  in ExcludedL]		
		ExcludedL = [x for x  in ExcludedL if len(x) > 1] 

	f = open(OutS,'w')
	f.close()

	for i in nameL[:]:

		if i not in ExcludedL:
			#ExcludedL = []
			sys.stderr.write("EXTENDING CONTIG "+ i +'\n') ##DEBUG
			sys.stderr.write("ExcludedL ="+ str(ExcludedL) +'\n') ##DEBUG

			#SseqS =  getSseq(ContignamesF,i,'3')
			f2 = open(ContignamesF,'r')
			SseqS = f2.read()
			f2.close()
			EachOUTL = mainBatch(SseqS,'3',IterationI,RdbS,CdbS,OutS,NameS,StartNum,ExcludedL) 
			ExcludedL = EachOUTL[1] #+ [i]

			f = open(OutS,'a')
			f.write("### Extended to 3 prime done \n" )
			f.close()

			SseqS =   EachOUTL[0].replace('LongSeq',i+"EXTENDED")
			EachOUTL = mainBatch(SseqS,'5',IterationI,RdbS,CdbS,OutS,NameS,StartNum,ExcludedL) 
			ExcludedL = EachOUTL[1] #+ [i]

			LoopTag = "### Extende From "+i+" BOTH SIDE\n" 

			f = open(OutS,'a')
			f.write(LoopTag)
			f.write(EachOUTL[0]+'\n')
			f.write("### Extended to 5 prime done NEXT contig\n" )
			f.close()
		
			print( EachOUTL[0])
		




