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

		Mhead = mark.sstartI
		Mtail = mark.sendI
		Mlength = abs(mark.sendI - mark.sstartI)
		
		## <------------------->
		##    <------------->
		temp = [x for x in a if min(x.sstartI,x.sendI) >= Mhead and max(x.sendI, x.sstartI) <= Mtail] #and abs(x.sendI - x.sstartI) > (coverage* Mlength)]

		##    <------------------->
		## <------------------>
		temp2 = [x for x in a if max(x.sendI, x.sstartI) < Mtail and Mlength >= abs(x.sendI - x.sstartI) >= (coverage* Mlength) and (Mhead - min(x.sstartI,x.sendI)) <= coverage*(abs(x.sendI - x.sstartI))] 
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		
		## <------------------->
		##     <------------------>		
		temp2 = [x for x in a if min(x.sendI, x.sstartI) > Mhead and Mlength >= abs(x.sendI - x.sstartI) >= (coverage* Mlength) and (max(x.sendI, x.sstartI) - Mtail) <= coverage*(abs(x.sendI - x.sstartI))]
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

class FastaTool(object): ## class for manipulate sequences in fasta format
	def __init__(self,FastaInS):
		self.FastaInS = FastaInS
		self.seqLenI = self.seqlen()

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


def main(INS): ##Take input as blast(n,x) result string
	cookedBlastL = [hit(x) for x in chop(INS)] #transform blast result into list of blast hit object
	ClusteredL = merge4(cookedBlastL,0.1)
	for i in ClusteredL: 
		STDERR('\n'.join([str(x.AttributesL) for x in i])) ##DEBUG
		STDERR("Each iter in ClusteredL ++++++++++++++++++++") ##DEBUG


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
	print(main(INS))
else:#write output to a flie
	f=open(options.o,'w')
	f.write(Main(INS))
	f.close()	

main(INS)

