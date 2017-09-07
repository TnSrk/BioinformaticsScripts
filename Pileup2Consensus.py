#Pileup2Consensus.py
import sys, re


def STDERR(*StrInS): #Function For Debugging
	if __ERRFLAG__ == 1:
		outS = ' '.join([str(x) for x in StrInS])
		sys.stderr.write(str(outS)+'\n')

def patternTrans(pS):
	#print(" pS = "+pS)
	presceedS = pS[0]
	signS = pS[1]
	numI = int(pS[2:])
	modsignS = ''
	if signS == "+":
		modsignS = "\+"
	elif signS == "-":
		modsignS = "\-"
	a = "\\"
	patternOutS =  "(" + presceedS +  modsignS +  str(numI) + "[a-zA-Z]"  + "{" + str(numI) + "})"
	return patternOutS

def lineCons(lineS):
	#print(lineS) ##DEBUG
	lineL = [x for x in lineS.split()]
	ContigNameS = lineL[0]; STDERR("ContigNameS",ContigNameS)
	PositionS = lineL[1]; STDERR("PositionS",PositionS)
	RefSeqS = lineL[2]; STDERR("RefSeqS",RefSeqS)
	BamNumI = len(lineL)
	DepthI = 0
	CharS = ''
	QS = ''
	STDERR(ContigNameS + " " + PositionS + RefSeqS)##DEBUG
	for chunkI in range(3, BamNumI, 3):
		DepthI = DepthI + int(lineL[chunkI])
		CharS = CharS + lineL[chunkI+1]
		QS = QS + lineL[chunkI+2]

	patternS0 = "\^.|\$"
	#patternSA = "(\+|\-)"
	patternS1 = "(.\+\d{1,}|.\-\d{1,})"
	#patternS2 = "(.+[0-9].{1,3})"
	
	#NegativePattern1 = "(?:\.+\d.{1,}|\.-\d{1,})"
	CharS_1 = re.sub(patternS0,"",CharS)
	#CharS_2 = re.sub(patternSA," \\1",CharS)
	STDERR(CharS)##DEBUG
	STDERR(CharS_1)##DEBUG
	#STDERR("CharS_2="+CharS_2)##DEBUG
	p = re.compile(patternS1)
	#newString = re.split(patternS2,CharS)
	#STDERR("newString="+str(newString))##DEBUG
	matchL = p.findall(CharS_1)
	STDERR(matchL) ##DEBUG
	patternL = [patternTrans(x) for x in matchL]
	STDERR("len(patternL)="+str(len(patternL)))##DEBUG
	STDERR("patternL="+str(patternL))##DEBUG
	setpatternL = list(set(patternL))
	STDERR("len(setpatternL)="+str(len(setpatternL)))##DEBUG
	STDERR("patternL="+str(setpatternL))##DEBUG
	
	TMPseq2 = CharS_1
	countL = []

	for i in setpatternL:
		p = re.compile(i)
		m = p.search(TMPseq2)
		matchS = ''
		countI = 0
		if m:
			matchS = m.group()
			TMPseq2 = p.sub("_",TMPseq2)
			countI = TMPseq2.count("_")
			TMPseq2 =  TMPseq2.replace("_","")

			countL.append([matchS,countI])

	CharL = [x for x in TMPseq2]
	SetCharL = list(set(CharL))
	STDERR("CharL=",CharL) ##DEBUG
	STDERR("SetCharL=",SetCharL) ##DEBUG
	for i in SetCharL:
		STDERR("each i=",i)
		countI = TMPseq2.count(i)
		countL.append([i,countI])
	NumRe0 = re.compile("\+\d{1,}")
	ModCountL = [[ NumRe0.sub("",x[0].replace(",",".").upper().replace(".",RefSeqS) ), x[1]] for x in  countL]	
	STDERR("countL="+str(countL))##DEBUG
	STDERR("ModCountL="+str(ModCountL))##DEBUG
	DelCharL = [x for x in ModCountL if x[0].find("-") != -1 ]
	NumRe1 = re.compile("\-\d.{1,}")
	ModCountL = [[ NumRe1.sub("",x[0] ), x[1]] for x in  ModCountL]
	STDERR("ModCountL="+str(ModCountL))##DEBUG
	STDERR("DelCharL="+str(DelCharL))##DEBUG
	FinalD = {}
	for i in ModCountL:
		if i[0] not in FinalD:
			FinalD[i[0]] = i[1]
		else:
			NewValueI = FinalD[i[0]] + i[1]
			FinalD[i[0]] = NewValueI

	STDERR("FinalD =",FinalD)
	FinalL = sorted([[ x,FinalD[x]] for x in FinalD], key=lambda x:x[1])[::-1]
	STDERR("FinalL =",FinalL)
		
		
	#print(CharS_2)
	STDERR("QS=",QS)
	STDERR("#")

	return [ContigNameS,PositionS,FinalL,DelCharL]

	
def main(INL,LineWidhtI):
	ConsFasta = ''
	previousDelNum = 0
	ConS = ''
	CurrentContigNameS = ''
	CurrentPositionI = 0
	outS = ''
	NumI = 0
	for lineS in INL:
		lineCons_OUT = lineCons(lineS)
		STDERR("lineCons_OUT=",lineCons_OUT)
		if lineCons_OUT[2][0][0] == "*" and  previousDelNum > 0:
			ConS = ''

		elif lineCons_OUT[2][0][0] == "*" and  previousDelNum == 0:
			if len(lineCons_OUT[2]) > 1:			
				ConS = lineCons_OUT[2][1][0]
			else:
				ConS = '' 
		else:
			ConS = lineCons_OUT[2][0][0]

		if lineCons_OUT[3] != []: 
			previousDelNum = lineCons_OUT[3][0][1]

		#STDERR("ConL=",ConL)

		STDERR("ConS=",ConS)
		if CurrentContigNameS != lineCons_OUT[0]:
			outS = outS + ">"+lineCons_OUT[0] + "\n"
			CurrentContigNameS = lineCons_OUT[0]

		outS = outS + ConS
		NumI += len(ConS)
		if NumI == LineWidhtI:
			outS = outS + "\n"
			NumI = 0

	print(outS)

def main2(fname,LineWidhtI,NlimitI,LinePrintNumI):
	## NlimitI = allowed number of consecutive "N" or "Gap" in output if exceed split to new contigs
	ConsFasta = ''
	previousDelNum = 0
	ConS = ''
	CurrentContigNameS = ''
	CurrentPositionI = 0
	CollectiveGapI = 0
	NewCollectiveGapI = 0
	GapFlag = -1
	outS = ''
	NumI = 0
	L_numI = 0
	with open(fname, 'r') as f:
		for line in f:
			lineCons_OUT = lineCons(line)
			#sys.stderr.write("lineCons_OUT="+str(lineCons_OUT)) ##DEBUG
			CurrentPositionI = int(lineCons_OUT[1])
			if lineCons_OUT[2][0][0] == "*" and  previousDelNum > 0: ##If consensus base is * and DELETION presents in previous line 
				ConS = '' ## then this base should be empty.

			elif lineCons_OUT[2][0][0] == "*" and  previousDelNum == 0:##If consensus base is * and DELETION presents is not in previous line 
				if len(lineCons_OUT[2]) > 1: ## ,if There are more than one Character in list
					ConS = lineCons_OUT[2][1][0] ## then use base from Reference sequence
				else: ## ,if There are only one Character in list
					ConS = 'n' ## then this base should be n because of ambiguities.
			else: ## If no Deletion Count in previous line 
				ConS = lineCons_OUT[2][0][0] ## then this base is the highest occurrence. 

			if lineCons_OUT[3] != []: 
				previousDelNum = lineCons_OUT[3][0][1]
			elif lineCons_OUT[3] == [] and previousDelNum > 0:
				previousDelNum = previousDelNum - 1
			else:
				previousDelNum = 0

			#STDERR("ConL=",ConL)

			#STDERR("ConS=",ConS)
			if CurrentContigNameS != lineCons_OUT[0] and GapFlag < 1:
				HeaderS = ">"+lineCons_OUT[0] + ":" + str(CurrentPositionI) + "\n"
				L_numI += 1
				CurrentContigNameS = lineCons_OUT[0]
				outS = outS + "\n" + HeaderS
							
			
			if ConS != '' and ConS != "N" and ConS != "n":
				NumI += len(ConS)				
				GapFlag = 0
				NewCollectiveGapI = CollectiveGapI
				CollectiveGapI = 0	
				sys.stderr.write("ConS ="+str(ConS)+"\n")##DEBUG
				
				
			else:
				#ConS = ''
				CollectiveGapI += 1
				GapFlag = 1
				sys.stderr.write("CollectiveGapI ="+str(CollectiveGapI)+"\n")##DEBUG


			if NewCollectiveGapI >= NlimitI and GapFlag == 0:
				outS = outS + "\n>"+lineCons_OUT[0] + ":" + str(CurrentPositionI) + "\n"				
				L_numI += 1
				sys.stderr.write("Another Gap Break \n") ##DEBUG

			outS = outS + ConS
			#sys.stderr.write("ConS="+ConS+"\n") ##DEBUG
			if NumI >= LineWidhtI:
				outS = outS + "\n"
				NumI = 0
				#sys.stderr.write("N_numI + G_numI= "+ str(N_numI + G_numI) + " \n") ##DEBUG
				L_numI += 1

			

			if L_numI >= LinePrintNumI:
				sys.stdout.write(outS)
				outS = ''
				L_numI = 0
				sys.stderr.write("Another Print Passed \n") ##DEBUG

			#

		if L_numI > 0:
			sys.stdout.write(outS)

			
		
		
__ERRFLAG__ = 0
#INL = [x for x in sys.stdin.readlines() if x[0] != "#"]
#main(INL,60)
fname = sys.argv[1]

main2(fname,100,1,1000)

