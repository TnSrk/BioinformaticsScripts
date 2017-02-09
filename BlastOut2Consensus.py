#!/usr/bin/python
#BlastOut2Consensus.py 2017 Feb 09
##Takes input from blastx with -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qcovs' option
import sys, subprocess, optparse, re
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

def main(INS): ##Take input as blast(n,x) result string
	cookedBlastL = [hit(x) for x in chop(INS)] #transform blast result into list of blast hit object
	
	


### Input and Option 
usage = "python Blastout2cnsensus.py -i input.gff -o out_file"
opt = optparse.OptionParser(usage)
opt.add_option("-i",help="*input path, get input from stdin if ommit", default='0')
opt.add_option("-o",help="indicate output file name or print out as standard output",default="0")
opt.add_option("-e",help="E-value cutoff for selected hit",default="0.1")
opt.add_option("-c",help="query-length covery rate cutoff for selected hit",default="0.5")
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

