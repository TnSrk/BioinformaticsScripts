#!/usr/bin/python
##augustus2TBL.py

import sys,optparse,copy
def STDERR(*StrInS):
	outS = ' '.join([str(x) for x in StrInS])
	sys.stderr.write(str(outS)+'\n')

def chop(Gffs): #divide each prediction into a list of Features 
	GffL = [[y.split("\t") for y in x.splitlines()] for x in Gffs.split("###\n")]	
	return GffL

def cut(GffL_elementL): ##separate each feature of each genes 
	MiscL = [x for x in GffL_elementL if x[0][0] == "#"]
	EdibleL = [x for x in GffL_elementL if x[0][0] != "#"]
	PretasteL = [x[2] for x in EdibleL]
	TasteL = []
	for i in PretasteL:
		if i not in TasteL:
			TasteL.append(i)
	EatL = []
	for i in TasteL:
		DishL = [x for x in EdibleL if x[2] == i]
		EatL.append(DishL)	

	return [MiscL,TasteL,EatL]

def pick(CutL,SelectedTasteL):
	ToChewL = []
	for i in SelectedTasteL:
		DeliciousPieceL = [x for x in CutL if x[0][2] == i] 
		ToChewL.append(DeliciousPieceL)

	return ToChewL


class Fruit(object):	
	def __init__(self,GffL_elementL):
		self.GffL_elementL = GffL_elementL
		self.CutL = cut(self.GffL_elementL)
		self.PeelL = self.CutL[0] ## augustus commnet part
		self.TypeL = self.CutL[1] ## feature types 
		self.PulpL = []
		for gffrowL in self.CutL[2]:
			for gffrow in gffrowL:
				eachPieceOBJ = Piece(gffrow)
				self.PulpL.append(eachPieceOBJ.pieceL) ## GFF row part

		self.DishL = []
		for TypeS in self.TypeL:
			groupL = [x for x in self.PulpL if x[0][2] == TypeS]
			groupL.append(TypeS)
			#groupL = [[x[0],x[1][TypeS]] for x in groupL]				
			self.DishL.append(groupL)

class Piece(object):
	def __init__(self,gffRow):
		self.gffRowL = gffRow
		self.pieceL = [self.gffRowL[0:8],self.chew(self.gffRowL)]
		
	def chew(self,Pulp):		
		LastColumnL = [x.strip()+'"' for x in Pulp[8].split('";')]
		if len(LastColumnL) > 1:
			TEMPL = [x.split() for x in  LastColumnL if len(x.split()) > 1]
			KL,VL = [x[0] for x in TEMPL ],[' '.join(x[1:]) for x in TEMPL] #generate Key and Value
			tasteD =  { k:v for k,v in zip(KL,VL)} ##gerate dictionary between ID and value
		else:
			tasteD = {Pulp[2]+"_id":LastColumnL[0]}
		return tasteD


def chew(DishL,SelectedTasteL):
	OutS = ''
	geneL = [x for x in DishL if x[-1] == 'gene'][0][0] ##select quallifier from gene record
	STDERR("THIS IS geneL") ##DEBUG
	STDERR(geneL) ##DEBUG
	STDERR("geneL END") ##DEBUG
	tagD = geneL[1] ##get dictionary from list
	tagD.update([x for x in DishL if x[-1] == 'transcript'][0][0][1]) ##select quallifier from gene record

	if geneL[0][6] == '-': ##Create gene position
		OutS = geneL[0][4] + "\t" + geneL[0][3] + "\t" + "gene" + "\n" + "\t"*3 + "locus_tag\t" + tagD['gene_id'] + "\n" ##minus strand feature
	elif geneL[0][6] == '+':
		OutS = geneL[0][3] + "\t" + geneL[0][4] + "\t" + "gene" + "\n" + "\t"*3 + "locus_tag\t" + tagD['gene_id'] + "\n" ##plus  strand feature
	
	for i in SelectedTasteL[1:]:
		PreCurrentL = [x[:-1] for x in DishL if x[-1] == i]
		if len(PreCurrentL) > 0:
			STDERR("THIS IS CurrentL") ##DEBUG
			CurrentL = PreCurrentL[0]
			STDERR(str(CurrentL)) ##DEBUG
			STDERR("str(CurrentL[0][0])") ##DEBUG			
			STDERR(str(CurrentL[0][0])) ##DEBUG			
			if geneL[0][6] == '-': ##Create other attributes positions MINUS
				STDERR("THIS IS CurrentL->"+str(CurrentL)+"<<<<<<<<<<<<<<<") ##DEBUG
				PositionL = [x[0][4]+"\t"+x[0][3] for x in sorted(CurrentL, key = lambda x:(int(x[0][3])+int(x[0][4]))/2)[::-1]]
				#PositionL = [x[0][4]+"\t"+x[0][3]+"MINUS" for x in CurrentL]
			elif geneL[0][6] == '+': ##Create other attributes positions PLUS
				STDERR("THIS IS CurrentL->"+str(CurrentL)+"<<<<<<<<<<<<<<<") ##DEBUG
				PositionL = [x[0][3]+"\t"+x[0][4] for x in CurrentL]
				
			NextLinesS = PositionL[0] + "\t" + i + "\n"
			for PositionS in PositionL[1:]:
				NextLinesS = NextLinesS + PositionS + "\t" + "\n"
			for att in ['product', 'gene_id', 'transcript_id']:
				if i in ['gene', 'transcript', 'CDS']:
					NextLinesS = NextLinesS + "\t"*3 + att + "\t" + tagD[att] + "\n"

			OutS = OutS + NextLinesS 
	return OutS

def SmoothChew(DishL,SelectedTasteL):
	SelectedTasteL = ['gene', 'transcript', 'CDS'] 
	OutS = ''
	geneL = [x for x in DishL if x[-1] == 'gene'][0][0] ##select quallifier from gene record
	STDERR("THIS IS geneL") ##DEBUG
	STDERR(geneL) ##DEBUG
	STDERR("geneL END") ##DEBUG
	tagD = geneL[1] ##get dictionary from list
	tagD.update([x for x in DishL if x[-1] == 'transcript'][0][0][1]) ##select quallifier from gene record
	##Partial feature Check
	StartCodonS = "0"
	StopCodonS = "0"
	StartFlagI = 0 ##Flag for partial feature editing
	StopFlagI = 0 ##Flag for partial feature editing
	STDERR("THIS IS DishL"+str(DishL)) ##DEBUG
	if geneL[0][6] == '+': #change last position of CDS to be stop codon
		StartCodonL = [x[:-1] for x in DishL if x[-1] == 'start_codon']
		if len(StartCodonL) > 0:
			StartCodonS = StartCodonL[0][0][0][3]
		else:
			StartCodonL = [x[:-1] for x in DishL if x[-1] == 'gene']
			StartCodonS = "<" + StartCodonL[0][0][0][3]
			StartFlagI = 1

		StopCodonL = [x[:-1] for x in DishL if x[-1] == 'stop_codon']
		if len(StopCodonL) > 0:
			StopCodonS = StopCodonL[0][0][0][4]
		else:
			StopCodonL = [x[:-1] for x in DishL if x[-1] == 'gene']
			StopCodonS = ">" + StopCodonL[0][0][0][4]
			StopFlagI = 1

	elif geneL[0][6] == '-': #change last position of CDS to be stop codon
		StartCodonL = [x[:-1] for x in DishL if x[-1] == 'start_codon']
		if len(StartCodonL) > 0:
			StartCodonS = StartCodonL[0][0][0][4]
		else:
			StartCodonL = [x[:-1] for x in DishL if x[-1] == 'gene']
			StartCodonS = "<" + StartCodonL[0][0][0][4]
			StartFlagI = 1

		StopCodonL = [x[:-1] for x in DishL if x[-1] == 'stop_codon']
		if len(StopCodonL) > 0:
			StopCodonS = StopCodonL[0][0][0][3]
		else:
			StopCodonL = [x[:-1] for x in DishL if x[-1] == 'gene']
			StopCodonS = ">" + StopCodonL[0][0][0][3]
			StopFlagI = 1

	STDERR("THIS IS START CODON"+str(StartCodonS)) ##DEBUG
	STDERR("THIS IS STOP CODON"+str(StopCodonS)) ##DEBUG

	if geneL[0][6] == '-': ##Create gene position
		#OutS = geneL[0][4] + "\t" + geneL[0][3] + "\t" + "gene" + "\n" + "\t"*3 + "locus_tag\t" + tagD['gene_id'].replace('"','') + "\n" ##minus strand feature
		OutS = StartCodonS + "\t" + StopCodonS + "\t" + "gene" + "\n" + "\t"*3 + "locus_tag\t" + tagD['gene_id'].replace('"','') + "\n" ##minus strand feature
	elif geneL[0][6] == '+':
		#OutS = geneL[0][3] + "\t" + geneL[0][4] + "\t" + "gene" + "\n" + "\t"*3 + "locus_tag\t" + tagD['gene_id'].replace('"','') + "\n" ##plus  strand feature
		OutS = StartCodonS + "\t" + StopCodonS + "\t" + "gene" + "\n" + "\t"*3 + "locus_tag\t" + tagD['gene_id'].replace('"','') + "\n" ##plus  strand feature
		
	##define mRNA position

	PreCurrentL = [x[:-1] for x in DishL if x[-1] == 'CDS']
	if len(PreCurrentL) > 0:
		CurrentL = PreCurrentL[0]
		if geneL[0][6] == '-': ##Create other attributes positions MINUS
			CurrentL2 = copy.deepcopy(sorted(CurrentL, key = lambda x:(int(x[0][3])+int(x[0][4]))/2)[::-1]) ##Create Minus Strand Feature
			CurrentL2[0][0][4] =  StartCodonS
			CurrentL2[-1][0][3] =  StopCodonS
			PositionL = [x[0][4]+"\t"+x[0][3] for x in CurrentL2]

		elif geneL[0][6] == '+': ##Create other attributes positions PLUS
			CurrentL2 = copy.deepcopy(CurrentL)
			CurrentL2[0][0][3] =  StartCodonS
			CurrentL2[-1][0][4] =  StopCodonS

			PositionL = [x[0][3]+"\t"+x[0][4] for x in CurrentL2]
			
		NextLinesS = PositionL[0] + "\t" + 'mRNA' + "\n"
		for PositionS in PositionL[1:]:
			NextLinesS = NextLinesS + PositionS + "\t" + "\n"
		for att in ['product', 'gene_id', 'transcript_id']:
			NextLinesS = NextLinesS + "\t"*3 + att + "\t" + tagD[att].replace('"','') + "\n"
		
		OutS = OutS + NextLinesS.replace('gene_id','protein_id') 	

	##define CDS position
	
	PreCurrentL = [x[:-1] for x in DishL if x[-1] == 'CDS']
	if len(PreCurrentL) > 0:
		CurrentL = PreCurrentL[0]
		if geneL[0][6] == '-': ##Create other attributes positions MINUS
			PositionL = [x[0][4]+"\t"+x[0][3] for x in sorted(CurrentL, key = lambda x:(int(x[0][3])+int(x[0][4]))/2)[::-1]]
			if StartFlagI == 1 and PositionL[0].split('\t')[0] == StartCodonS.replace("<",""):
				PositionL[0] = "<" + PositionL[0] ##change to Partial feature
			#OutS = OutS + str(StartFlagI) +" StartFlagI " + str(PositionL[0].split('\t')[0]) +" PositionL[0].split('\t')[0] " + StartCodonS.replace("<","") + ' StartCodonS.replace("<","")\n' ##DEBUG 
			if StopFlagI == 1 and PositionL[-1].split('\t')[1] == StopCodonS.replace(">",""):
				PositionL[0] = PositionL[-1].replace("\t","\t>") ##change to Partial feature


		elif geneL[0][6] == '+': ##Create other attributes positions PLUS
			PositionL = [x[0][3]+"\t"+x[0][4] for x in CurrentL]
			if StartFlagI == 1 and PositionL[0].split('\t')[0] == StartCodonS.replace("<",""):
				PositionL[0] = "<" + PositionL[0] ##change to Partial feature
			if StopFlagI == 1 and PositionL[-1].split('\t')[1] == StopCodonS.replace(">",""):
				PositionL[-1] = PositionL[-1].replace("\t","\t>") ##change to Partial feature
			#OutS = OutS + str(StartFlagI) +" StartFlagI " + str(PositionL[-1].split('\t')[1]) +" PositionL[-1].split('\t')[1] " + StartCodonS.replace(">","") + ' StartCodonS.replace(">","")\n' ##DEBUG 
		NextLinesS = PositionL[0] + "\t" + 'CDS' + "\n"
		for PositionS in PositionL[1:]:
			NextLinesS = NextLinesS + PositionS + "\t" + "\n"
		for att in ['product', 'gene_id', 'transcript_id']:
			NextLinesS = NextLinesS + "\t"*3 + att + "\t" + tagD[att].replace('"','') + "\n"
		
		OutS = OutS + NextLinesS.replace('gene_id','protein_id') 	

	
	##define intron position
	PreCurrentL = [x[:-1] for x in DishL if x[-1] == 'intron']
	if len(PreCurrentL) > 0:
		CurrentL = PreCurrentL[0]
		if geneL[0][6] == '-': ##Create other attributes positions MINUS
			PositionL = [x[0][4]+"\t"+x[0][3] for x in sorted(CurrentL, key = lambda x:(int(x[0][3])+int(x[0][4]))/2)[::-1]]
			if StartFlagI == 1 and PositionL[0].split('\t')[0] == StartCodonS.replace("<",""):
				PositionL[0] = "<" + PositionL[0] ##change to Partial feature
			#OutS = OutS + str(StartFlagI) +" StartFlagI " + str(PositionL[0].split('\t')[0]) +" PositionL[0].split('\t')[0] " + StartCodonS.replace("<","") + ' StartCodonS.replace("<","")\n' ##DEBUG 
			if StopFlagI == 1 and PositionL[-1].split('\t')[1] == StopCodonS.replace(">",""):
				PositionL[0] = PositionL[-1].replace("\t","\t>") ##change to Partial feature
			
		elif geneL[0][6] == '+': ##Create other attributes positions PLUS
			PositionL = [x[0][3]+"\t"+x[0][4] for x in CurrentL]
			if StartFlagI == 1 and PositionL[0].split('\t')[0] == StartCodonS.replace("<",""):
				PositionL[0] = "<" + PositionL[0] ##change to Partial feature
			if StopFlagI == 1 and PositionL[-1].split('\t')[1] == StopCodonS.replace(">",""):
				PositionL[-1] = PositionL[-1].replace("\t","\t>") ##change to Partial feature
			#OutS = OutS + str(StartFlagI) +" StartFlagI " + str(PositionL[-1].split('\t')[1]) +" PositionL[-1].split('\t')[1] " + StartCodonS.replace(">","") + ' StartCodonS.replace(">","")\n' ##DEBUG 

		IntronNumberI = 1 ##intron number increment
		NextLinesS = '' ## innitiate Next Line Text
		for PositionS in PositionL:
			NextLinesS = NextLinesS + PositionS + "\t" + 'intron' + "\n" + "\t"*3 + 'number' + "\t" + str(IntronNumberI) + "\n"
			IntronNumberI += 1

		OutS = OutS + NextLinesS 
	return OutS

def Main(Gffs):
	chopedL = chop(Gffs)[:] 
	#for i in chopedL:##DEBUG
	#	print len(i) ##DEBUG
	chopedL = [x for x  in chopedL if len(x) > 3]
	SelectedTasteL = ['gene', 'transcript', 'CDS', 'intron'] 
	
	ContigTagS = ''
	OutStringS = ''
	STDERR("THIS IS chopedL")
	STDERR(len(chopedL[0]),len(chopedL[10]),len(chopedL[-1])) ##DEBUG
	STDERR(str(chopedL[0]),"#####\n",str(chopedL[10]),"#########\n",str(chopedL[-1])) ##DEBUG
	for pieceE in chopedL: 
		FruitOBJ = Fruit(pieceE)

		STDERR("FruitOBJ.DishL") ##DEBUG
		STDERR(str(FruitOBJ.DishL).replace(", [[[","\n, [[[")) ##DEBUG
		
		#print(FruitOBJ.DishL[0][0][0][0]) ##DEBUG
		if FruitOBJ.DishL[0][0][0][0] != ContigTagS:
			ContigTagS = FruitOBJ.DishL[0][0][0][0]
			OutStringS = OutStringS + ">Features " + ContigTagS  + "\n"
		OutStringS = OutStringS + SmoothChew(FruitOBJ.DishL,SelectedTasteL)

	return(OutStringS) ##DEBUG


def TestMain():

	chopedL = chop(Gffs)[:] 

	SelectedTasteL = ['gene', 'transcript', 'intron', 'CDS','start_codon','stop_codon'] 
	
	ContigTagS = ''
	for pieceE in chopedL: 
		FruitOBJ = Fruit(pieceE)
		#for i in FruitOBJ.CutL: ##DEBUG
		#	print(str(i)) ##DEBUG
		ToChew = pick(FruitOBJ.CutL[2],SelectedTasteL)
		#print("###THIS#IS#TO#SHEW#")
		#print(str(ToChew)) ##DEBUG
		geneL = FirstRowL[0]
		
		if ContigTagS != geneL[0]:
			print(">Features "+geneL[0])
		ContigTagS = geneL[0]
		print(geneL[3]+"\t"+geneL[4]+"\t"+geneL[2])
		print("\t"*3 + "locus_tag" + "\t" + FirstRowL[1]['gene_id'])

		TranscriptL = [x for x in FruitOBJ.DishL if x[0][0][2] == "transcript"][0]
		print (TranscriptL[0][0][3]+"\t"+TranscriptL[0][0][4] + "\tmRNA")
		for i in TranscriptL[1:]:
			print(i[0][3]+"\t"+i[0][4])

		productS = [x for x in FruitOBJ.DishL if x[0][0][2] == "gene"][0][0][1]['product']
		qualifierS = ''
		qualifierS = qualifierS + "\t"*3+"product"+"\t"+ productS + "\n"
		protein_idS = FirstRowL[1]['gene_id'].replace('"','') + "_pid"
		qualifierS = qualifierS + "\t"*3 + "protein_id" + "\t" + protein_idS  + "\n"
		transcript_idS = TranscriptL[0][1]['transcript_id']
		qualifierS = qualifierS + "\t"*3 + "transcript_id" + "\t" + transcript_idS  + "\n"
		print(qualifierS)

		CDS_L = [x for x in FruitOBJ.DishL if x[0][0][2] == "CDS"][0]
		print (CDS_L[0][0][3]+"\t"+CDS_L[0][0][4] + "\tCDS")
		for i in CDS_L[1:]:
			print(i[0][3]+"\t"+i[0][4])
		print(qualifierS)
		

### Input and Option 
usage = "augustus2TBL.py -i input.gff -o out_file"
opt = optparse.OptionParser(usage)
opt.add_option("-i",help="*input path, get input from stdin if ommit", default='0')
opt.add_option("-o",help="indicate output file name or print out as standard output",default="0")
opt.add_option("-x",help="*exclude gene_id file", default='0')
opt.add_option("--id-fill",help="*exclude gene_id file", default='0')
(options, args) = opt.parse_args()


if options.i == '0': ##get input from pipe
	Gffs = sys.stdin.read()
else:#open file
	f=open(options.i,'r')
	Gffs = f.read()
	f.close()

ExcludeL = []
if options.x != '0': ##get exclude list from pipe
	f=open(options.x,'r')
	ExcludeL = [x.strip() for x in f.readlines() if len(x) > 0]
	f.close()
	GffL = [x.strip() for x in Gffs.split("###\n") if x.find(' gene_id ') != -1]
	STDERR(ExcludeL)
	NewGffS = ''
	for line in GffL:
		#STDERR(line.split(';'),"***********") ##DEBUG
		gene_idS = [xx for xx in list(set([x for x in line.split(';') if x.find(' gene_id ') == 0])) if len(xx) > 0][0].split()[-1].replace('"','')
		STDERR(gene_idS) ##DEBUG
		if gene_idS not in ExcludeL:
			NewGffS = NewGffS + line + "###\n"
	Gffs = NewGffS
	

if options.o == '0': ##print output to pipe
	print(Main(Gffs))
else:#write output to a flie
	f=open(options.o,'w')
	f.write(Main(Gffs))
	f.close()







