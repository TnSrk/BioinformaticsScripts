#!/usr/bin/python
##TXblastX.py.py
##Takes input from blastx with -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qcovs' option

import sys,optparse
### Class And Function
def STDERR(*StrInS): #Function For Debugging
	outS = ' '.join([str(x) for x in StrInS])
	sys.stderr.write(str(outS)+'\n')

def chop(blastxS): #divide each hit into a list of attributes 
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
		if OBJ.qseqidS not in nameL:
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
	selectOJBL = [x for x in hitOBJL if x.qcovsF > __QCovF__ and x.evalueF < __EvalF__]
	#selectOJBL = selectOJBL[0:1000] ##DEBUG
	
	STDERR([x.AttributesL for x in selectOJBL[0:2]]) ##DEBUG
	STDERR("entering group function") ##DEBUG
	outS = ''
	#for OBJL in group(selectOJBL):
	for OBJL in Pgroup(selectOJBL):
		sortedOBJL = sorted([x.AttributesL for x in OBJL], key=lambda x:(int(x[8]) + int(x[9]))/2 )
		#outS = outS + '\n'.join(['\t'.join(x.AttributesL) for x in sortedOBJL]) + "\n############\n"
		outS = outS + '\n'.join(['\t'.join(x) for x in sortedOBJL]) + "\n############\n"
	#STDERR(group(selectOJBL)[0])
	return outS

### Input and Option 
usage = "python TXblastX.py.py -i input.gff -o out_file"
opt = optparse.OptionParser(usage)
opt.add_option("-i",help="*input path, get input from stdin if ommit", default='0')
opt.add_option("-o",help="indicate output file name or print out as standard output",default="0")
opt.add_option("-e",help="E-value cutoff for selected hit",default="0.1")
opt.add_option("-c",help="query-length covery rate cutoff for selected hit",default="0.98")
(options, args) = opt.parse_args()

##Documentation part
__QCovF__ = float(options.c)*100.0 ;STDERR("__QCovF__ = ",__QCovF__)
__EvalF__ = float(options.e)	;STDERR("__EvalF__ = ",__EvalF__)

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



