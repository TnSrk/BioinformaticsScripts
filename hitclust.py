#!/usr/bin/python
## hitclust10_j.py

##add clustalw alignment program option [-a 2]
##20121128:change fastacmd arg to get lenght from 100000000 to 1000000000
##20130811:Add AlignReOrder() to fix "-stable" option missing  

import sys, subprocess, os, optparse

## +++++++ contigName(a) ++++++++++ 
def contigsName(a): ## Extract contig names from blasthits
	
	a = a.split('\n')
	#a = [x for x in a if x != '']
	a = [x.split('\t') for x in a][0:]
	a = [x for x in a if len(x) == 12]
	#a = [x for x in a if float(x[10]) <= Eval]
	contignum = 0
	
	contignameL = []
	contignameL.append(a[0][0])
	p = -1
	for i in a[1:]:
		p += 1
		if i[0] not in contignameL:
			contignum += 1 
			contignameL.append(i[0])
	contignameL = [x for x in contignameL if x != '']		
	return contignameL

## ++++++++++++ contigjuice(blastfilename,contigname) +++++++++
def contigJuice(a,contigname):

	a = a.split('\n')
	contig = contigname+'\t'

	b = []
	for i in a:
		if i.find(contig) != -1:
			b.append(i[0:])
	
	return b

## ++++++++ blastTrans(contig) +++++++++
def blastTrans(b):

	b = [x.split('\t') for x in b]
	pL = [[x[0], x[1], float(x[2]), int(x[3]), int(x[4]), int(x[5]), int(x[6]), int(x[7]), int(x[8]), int(x[9]), x[10], float(x[11])] for x in b]
	pL = [x for x in pL if float(x[-2]) < Eval]
	return pL

## ++++ merge() ++++++
def merge4(inputlist,coverage):
	a = sorted(inputlist, key = lambda inputlist:abs(inputlist[7] - inputlist[6]) )[::-1]
	clustL = []

	while len(a) > 0:
	
		mark = a.pop(0)
		#clustL.append([])

		Mhead = mark[6]
		Mtail = mark[7]
		Mlength = abs(mark[7] - mark[6])
		
		## <------------------->
		##    <------------->
		temp = [x for x in a if x[6] >= Mhead and x[7] <= Mtail and abs(x[7] - x[6]) > (coverage* Mlength)]

		##    <------------------->
		## <------------------>
		temp2 = [x for x in a if x[7] < Mtail and Mlength >= abs(x[7] - x[6]) >= (coverage* Mlength) and (Mhead - x[6]) <= coverage*(abs(x[7] - x[6]))] 
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		
		## <------------------->
		##     <------------------>		
		temp2 = [x for x in a if x[6] > Mhead and Mlength >= abs(x[7] - x[6]) >= (coverage* Mlength) and (x[7] - Mtail) <= coverage*(abs(x[7] - x[6]))]
		temp2 = [x for x in temp2 if x not in temp]

		temp = temp + temp2
		temp = sorted(temp, key=lambda temp:temp[6])
		
		a = [x for x in a if x not in temp]
		
		pick = [x for x in temp]
		pick.insert(0,mark) 
		clustL.append(pick)
	return clustL
	
def Rmerge(inputlist,coverage):
	a = sorted(inputlist, key = lambda inputlist:inputlist[3])
	clustL = []

	while len(a) > 0:
	
		mark = a.pop(0)
		clustL.append([])

		Mhead = mark[6]
		Mtail = mark[7]
		Mlength = mark[3]
		
		##Ref    <------------->
		##Sub <------------------->
		temp = [x for x in a if x[6] <= Mhead and x[7] >= Mtail]
		
		temp = trim(temp,Mhead,Mtail,Mlength,coverage)		
			
		tL = [x for x in temp]
		##Ref    <------------------->
		##Sub <------------------>
		temp2 = [x for x in a if x[6] <= Mhead and x[7] <= Mtail and x[6] <= Mhead <= x[7] and Mlength*coverage <= (x[7] - Mhead)] 
		temp2 = trim(temp2,Mhead,Mtail,Mlength,coverage)
		tL += [x for x in temp2 if x not in tL] 

		##Ref <------------------->
		##Sub     <------------------>		
		temp3 = [x for x in a if x[6] > Mhead and x[7] >= Mtail and x[6] <= Mtail <= x[7] and Mlength*coverage <= (Mtail - x[6])]
		temp3 = trim(temp3,Mhead,Mtail,Mlength,coverage)
		tL += [x for x in temp3 if x not in tL]

		#temp = temp + temp2
		temp = sorted(tL, key=lambda tL:tL[6])
		
		a = [x for x in a if x not in temp]
		
		pick = [x for x in temp]
		pick.insert(0,mark) 
		clustL[-1] = pick
	return clustL

def Rmerge2(inputlist,coverage): ##include trim and remains
	a = [x for x in inputlist if abs(x[9] - x[8] + 1) >= 50]
	a = sorted(a, key = lambda a:a[3])
	clustL = []
	
	while len(a) > 0:
	
		mark = a[0]
		#clustL.append([])

		Mhead = mark[6]
		Mtail = mark[7]
		Mlength = mark[3]
		
		##Ref    <------------->
		##Sub <------------------->
		temp = [x for x in a if x[6] <= Mhead and x[7] >= Mtail]
		
		#temp = trim(temp,Mhead,Mtail,Mlength,coverage)		
			
		tL = [x for x in temp]
		##Ref    <------------------->
		##Sub <------------------>
		temp2 = [x for x in a if x[6] <= Mhead and x[7] <= Mtail and x[6] <= Mhead <= x[7] and Mlength*coverage <= (x[7] - Mhead)] 
		#temp2 = trim(temp2,Mhead,Mtail,Mlength,coverage)
		tL += [x for x in temp2 if x not in tL] 

		#temp = temp + temp2
		
		##Ref <------------------->
		##Sub     <------------------>		
		temp3 = [x for x in a if x[6] > Mhead and x[7] >= Mtail and x[6] <= Mtail <= x[7] and Mlength*coverage <= (Mtail - x[6])]
		#temp3 = trim(temp3,Mhead,Mtail,Mlength,coverage)
		tL += [x for x in temp3 if x not in tL]

		#temp = temp + temp2
		#temp = sorted(temp, key=lambda temp:temp[6])
		temp = sorted(tL, key=lambda tL:tL[6])
		
		a = [x for x in a if x not in temp]
		tempL = chop2(temp)
		
		if len(tempL[0]) > 0:			
			pick = [x for x in tempL[0]]
			clustL.append(pick)
		if len(tempL[1]) > 0:
			a += [x for x in tempL[1]]
		#pick.insert(0,mark) 
		
	return clustL

#### what chop() do 
##REF  <----------------------------------------->     <----------------->  <------------> <------->
##RE1  <----------------->                             <----------------->
##RE2                       <------------>         =>                       <------------>
##RE3  <----------------->  <-------------------->     <----------------->  <------------> <------->
def chop2(tempL):
	chopL = tempL[:]	
	chopedL = []
	midL,median,medL,medp = mid(chopL)
	#chopL = sorted(midL,key = lambda midL:abs((midL[7] + midL[6])/2 - medp))
	chopL = sorted(midL,key = lambda midL:abs(abs(midL[7] - midL[6]) - medL))
#	print 'chopL=',chopL ##DEBUG
	mark = chopL[0]
	Mhead = mark[6]; Mtail = mark[7]
#	print 'Mhead=',Mhead,' Mrtail=',Mtail  ##DEBUG
	remainL = remains(chopL,Mhead,Mtail)
#	print 'remainL=',remainL ##DEBUG
	remainL = [x for x in remainL if abs(x[7] - x[6] + 1) >= 50]
	trimL = trim(chopL,Mhead,Mtail)
#	print 'trimL=',trimL ##DEBUG
	trimL = [x for x in trimL if abs(x[7] - x[6] + 1) >= 50]
	chop2L = [trimL,remainL]
#	print 'chop2L=',chop2L ##DEBUG
	return chop2L

def trim(tempL,Mhead,Mtail):
	temp = list(tempL)
	trimL = []
	for i in temp:
		hremainI = Mhead - i[6] + 1
		tremainI = i[7] - Mtail + 1
#		print i,"###BEFORE TRIMING###DEBUG" ##DEBUG
#		print 'hremainI=',hremainI,'tremainI=',tremainI#DEBUG
		trim = i[:]
		if hremainI >= 50:
			trim[6] = Mhead
			if i[8] < i[9]:	trim[8] = i[8] + hremainI				
			else:           trim[8] = i[8] - hremainI
		if tremainI >= 50:
			trim[7] = Mtail
			if i[8] < i[9]:	trim[9] = i[9] - tremainI
			else:           trim[9] = i[9] + tremainI
		
		trimL.append(trim)
		
#		print i,">>>>>> AFTER TRIMING <<<< DEBUG" ##DEBUG

	trimL = [x for x in trimL if (abs(x[9] - x[8])+1) >= 50]
	return trimL

def remains(tempL,Mhead,Mtail):
	keepL = list(tempL)
	remainL = []
	
	for i in keepL:
		hremainI = Mhead - i[6] + 1
		tremainI = i[7] - Mtail + 1
#		print 'i=',i ##DEBUG
#		print 'hremain=',hremainI,'tremain=',tremainI ##DEBUG
		if hremainI >= 50:
			#temp = list(i)
			temp = i[0:7] + [[]] + [i[8]] + [[]] + i[10:]
			temp[7] = Mhead
			if i[8] < i[9]:
				#      <--------->
				# <-------------->
				temp[9] = i[8] + hremainI				
			else:
				temp[9] = i[8] - hremainI
			remainL.append(temp)
#			print "hremainI >= 50:" ##DEBUG
		if tremainI >= 50:
			#temp = list(i)
			temp = i[0:6] + [[]] + [i[7]] + [[]] + i[9:]
			temp[6] = Mtail
			if i[8] < i[9]:
				temp[8] = i[9] - tremainI
			else:
				temp[8] = i[9] + tremainI
			remainL.append(temp)
#			print "tremainI >= 50:" ##DEBUG
	sys.stderr.write(str(remainL)+'\n')
	remainL = [x for x in remainL if (abs(x[9] - x[8]) + 1) >= 50]
	return remainL


def bitcount(fname):
	Dbit = {}
	f = open(fname,'r')
	for line in f:
		x = line.split()
		if len(x) == 12:
			if x[1] in Dbit:
				Dbit[x[1]] += float(x[11])
			else:
				Dbit[x[1]] = float(x[11])
	f.close()
	return Dbit

def bitcount2(fname):
	Dbit = {}
	f = open(fname,'r')
	eachf = f.readlines()
	f.close()
	for line in eachf:
		x = line.strip().split()
		if len(x) == 12:
			if x[1] in Dbit:
				Dbit[x[1]] += float(x[11])
			else:
				Dbit[x[1]] = float(x[11])

	return Dbit

def bitselect(clustL,mode):

	midL,median,medL,medp = mid(clustL)
	if len(midL) > 5:
		if   mode == "1":midL = sorted(midL,key = lambda midL:midL[-2])[0:5] #Top five least bitscore
		elif mode == "2":midL = bitpick(midL) #average bit score 
		elif mode == "3":midL = sorted(midL,key = lambda midL:abs(abs(midL[7] - midL[6]) - medL))[0:5] #selecting by median length
		elif mode == "4":midL = sorted(midL,key = lambda midL:float(midL[11]))[::-1][0:5] #select five most similar sequences (tpo 5 maximum bit score)
		elif mode == "5":midL = sorted(midL,key = lambda midL:abs(midL[13]-medp))[0:5] #select five middle sequences, sorted by middle position on target sequence 
	return midL

def mid(clustL):
	
        a = sorted(clustL,key = lambda clustL:abs(clustL[7] - clustL[6]))#[::-1] ## ORDER BY match length
	midL = [x[:12]+[Dbit[x[1]]]+[(x[6]+x[7])/2] for x in a] ## to add bitcount of each genome and middle of match position
        midL = sorted(midL,key = lambda midL:midL[-1])[::-1] ##ORDER BY middle of match position
	al = len(midL)

 	if al%2 == 1:
 		median = (midL[(al/2) - 1][-1] + midL[(al/2)][-1])/2 ##ORDER BY bit-score 
 		medL = (abs(midL[(al/2) - 1][7] - midL[(al/2) - 1][6]) + abs(midL[(al/2)][7] - midL[(al/2)][6]))/2 ##ORDER BY match length
		medp = (midL[(al/2) - 1][13] + midL[(al/2)][13])/2 ##ORDER BY middle of match position
 	else:
 		median = midL[(al/2)][-1] ##ORDER BY bit-score
 		medL = abs(midL[(al/2)][7] - midL[(al/2)][6]) ##ORDER BY match length
		medp = midL[(al/2)][13] ##ORDER BY middle of match position
	
	midset = [midL,median,medL,medp]	
	return midset
	
def bitpick(bitL):
	bitsum = [x[-2] for x in bitL]
	bitsum = sum(bitsum)/len(bitsum)

	bitL = sorted(bitL, key=lambda bitL:bitL[-2])
	aL = [bitL.pop(0)]
	lastL = bitL.pop(-1)
	#sum
	for i in (0.25,0.50,0.75):		
		bitL = sorted(bitL, key=lambda bitL:abs(bitL[-2] - bitsum*i))
		aL.append(bitL.pop(0))
	aL.append(lastL)  
	return aL

## ++++ theone() +++++++++++++++++++++++++++++++++++++++
def theone(clustL):

	a = sorted(clustL,key = lambda clustL:clustL[11])[::-1]
	b = [x[1] for x in a]
	nameL = relatename(a)
	num = []

	for i in nameL:
		num.append(b.index(i))

	uniqueL = []

	for i in num:
		uniqueL.append(a[i])
		
	if logname != '0':
		x = open(logname,'a')
		x.write(str(uniqueL)+'\n')
		x.close()
#	sys.stderr.write("\n"+str(uniqueL)+"\n<<<<<<<THE ONE<<<<<<<<<<<<<DEBUG") ##DEBUG
	return uniqueL	

def theone3(clustL):

	a = sorted(clustL,key = lambda clustL:clustL[11])[::-1]
	b = [x[1] for x in a]
	nameL = relatename(a)
	num = []

	for i in nameL:
		num.append(b.index(i))
		hitjoin(a,i)
		
	uniqueL = []

	for i in num:
		uniqueL.append(a[i])
		
	if logname != '0':
		x = open(logname,'a')
		x.write(str(uniqueL)+'\n')
		x.close()

	return uniqueL	
	
def hitjoin(hitL,name):
	newhitL = [x for x in hitL if x[1].find(name) != -1][:]
	newhitL.sort(key=lambda newhitL:(min(newhitL[8],newhitL[9])))
	newhitL = newhitL[::-1]
	tempL = []
	mark = list(newhitL[0])
	for i in newhitL[1:]:
		xl = [min(i[6],i[7]) - max(mark[6],mark[7]),i[8] - mark[9]]
		mark = list(i)
		xl = i + xl
		tempL.append(xl)
	return tempL	 
	
def theone2(clustL):

	a = sorted(clustL,key = lambda clustL:clustL[11])[::-1]
	b = [x[1] for x in a]
	nameL = relatename(a)
	num = []
	uniqueL = []
	
	for i in nameL:
		eachL = [x for x in a if x[1] == i]
		eachL = sorted(eachL,key = lambda eachL:eachL[11])[::-1]
		if len(eachL) > 1:GPL = GPjoin(eachL);print GPL
		else:GPL = i
		uniqueL.append(GPL)	
		
	if logname != '0':
		x = open(logname,'a')
		x.write(str(uniqueL)+'\n')
		x.close()

	return uniqueL		
	
## ++++++ relatename(clustL) +++++++	
def relatename(clustL):
#	print "RELATED NAME" ##DEBUG
#	print clustL,"\n############sep" ##DEBUG
	a = sorted(clustL,key = lambda clustL:clustL[1])
	a = [x[1] for x in a]
#	print a ##DEBUG
#	print "END OF RELATED NAME" ##DEBUG
	nameL = [a.pop(0)]
	
	for i in a:
		if i not in nameL:
			nameL.append(i)
	return nameL

def GPjoin(eachL):
	tL = eachL[:]
	#GPL = []
	GPL = [tL.pop(0)]
	while len(tL) > 0:
		TEMP = GPL[-1]
		exL = [x for x in tL if TEMP[6] <= (x[6]+x[7])/2.0 <= TEMP[7] or x[6] <= (TEMP[6]+TEMP[7])/2.0 <= x[7]]
		tL = [x for x in tL if x not in exL]
		if len(tL) > 0:GPL.append(tL.pop(0))
#	print GPL ## DEBUG
#	print '##GPL < #################' DEBUG
	return GPL
	
		

## +++++++ getAlign(DB,qeury,clustlist) ++++++++++++++++++++++++++
## use to get reference sequence and related sequences for using as MSA input 


## +++++++ getAlign() ++++++++++++++++++++++++
def getAlign3(DB,qeury,clust):
	qStartL = [int(x[6]) for x in clust]
	#qStartL.sort()
	#Qhead = qStartL[0]; del qStartL #old
	Qhead = min(qStartL)#; #del qStartL

	qStopL = [int(x[7]) for x in clust]
	#qStopL.sort()
	#Qtail = qStopL[-1]; del qStopL #old
	Qtail = max(qStopL)#; del qStopL

	qname = clust[0][0] ### clean query name before get string
	align = Qseqjuice5(qeury,qname,Qhead,Qtail); #del qname, Qhead, Qtail

	SL = clust
	SL = [[x[i] for i in (1,8,9)] for x in SL]

	for i in  SL:
		
		Sseq = getSseq4(DB,i[0],i[1],i[2])
		align += Sseq
	return align
	
def  Qseqjuice5(query,name,Qhead,Qtail):
	if name.find('|') != -1:
		pre_head = name.split('|')[1]
		pre_head = int(pre_head.split(':')[0])
	else:
		pre_head = 0

	head = pre_head + int(Qhead) - 1
	tail = pre_head + int(Qtail) - 1

	if head > tail:
		head,tail = tail,head

	if name.find('|') != -1:
		name = name.split('|')[0:-1]
		name = ''.join(name)
		name = name.replace('>','')
	
	seqlen = Dseqlim[str(name)]

	if tail > seqlen:
		tail = seqlen

	arg = """fastacmd -p F -d __DBNAME__ -s "__name__" -L __head__,__tail__ ;"""
	arg = arg.replace('__DBNAME__',query)
	arg = arg.replace('__name__',name)
	arg = arg.replace('__head__',str(head))
	arg = arg.replace('__tail__',str(tail))

	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()
	
	TEMPseqs = x[0]
	#sys.stderr.write(str(x[1])) ##DEBUG
	del process
	return	TEMPseqs


## ++++++++++ clean(text) +++++++++
def clean(text):
	text = text.upper()
	text = [x for x in text if x == 'A' or x == 'T' or x == 'C' or x == 'G']
	text = ''.join(text)
	return text

## ++++++++++ reversecomplement(string) ++++++++++
def reversecomplement(string):
	reverse = string[::-1]
	
	reverse = reverse.replace("A","?")
	reverse = reverse.replace("T","A")
	reverse = reverse.replace("?","T")
	reverse = reverse.replace("G","?")
	reverse = reverse.replace("C","G")
	reverse = reverse.replace("?","C")

	return reverse ## complementary string ATCG --> CGAT

## getSseq4(DBname,Scontigname,Shead,Stail)
def getSseq4(DBname,Scontigname,Shead,Stail):
#	print DBname,Scontigname,Shead,Stail ##DEBUG
	arg1 = "fastacmd -p F -d "+DBname+" -s \""+Scontigname+"\" -L "+str(Shead)+","+str(Stail)+" ;"
	arg2 = "fastacmd -p F -S 2 -d "+DBname+" -s \""+Scontigname+"\" -L "+str(Stail)+","+str(Shead)+" ;"
		
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

## +++++ MSAgen() +++++
def MSAgen4(align,Ip):

	if Ip == '0':
	        arg = "echo \'##########\n\' ; muscle_old -clwstrict -gapopen __GAP__ -stable -maxiters 32 -quiet"
		arg = arg.replace("__GAP__",GAP)

		process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process.communicate(align)
		#x = x[0]

	elif Ip == '1':
		Dn = shortname(align)
		for n in Dn: ##to short name in alignment
			label = len(Dn[n]) - len(n)
			label = n+' '*label
			align = align.replace(Dn[n],label)
		#print align
		f = open('TEMPseqs','w')
		f.write(align)
		f.close()
		arg = """echo \'##########\n\nCLUSTAL W (1.81) multiple sequence alignment' ; mafft --clustalout --quiet --maxiterate 1000 --localpair TEMPseqs | grep -v "CLUSTAL format alignment by MAFFT" | tr "[:lower:]" "[:upper:]" """

		process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process.communicate()
		#x = x[0]

	elif Ip == '2':
		Dn = shortname(align)
		for n in Dn: ##to short name in alignment
			label = len(Dn[n]) - len(n)
			label = n+' '*label
			align = align.replace(Dn[n],label)
		#print align
		f = open('TEMPseqs','w')
		f.write(align)
		f.close()

		arg = """clustalw -infile=TEMPseqs -outorder=INPUT -quiet -outfile=TEMPseqs.aln"""

		process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process.communicate()
		del x

		arg = """echo \'##########\n\nCLUSTAL W (1.81) multiple sequence alignment' ;cat TEMPseqs.aln | grep -v "CLUSTAL 2.1 multiple sequence alignment"  """

		process2 = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process2.communicate()
		
		

	alignout = x[0]
	if Ip == '1' or Ip == '2':
		
		for keys in Dn.keys():
			Dn[keys] = Dn[keys].split(' ')[0][:]

		Ilngth = max([len(x) for x in Dn.values()])
#		print Ilngth  ##DEBUG
		Klngth = max([len(x) for x in Dn.keys()])
#		print Klngth  ##DEBUG
		for n in Dn: ##to short name in alignment
			label = Ilngth - len(Dn[n])
			label = Dn[n]+' '*label
			alignout = alignout.replace(n,label)
		Saster = '\n '+(' '*(Ilngth - Klngth))
		alignout = alignout.replace('\n ',Saster)
		#alignout = Salignout.replace('#',' ')
		
#	sys.stderr.write(str(x[1])) ##DEBUG
	del process
	return	alignout

## +++++ MSAgen5() +++++
## rename for all alignment program
def MSAgen5(align,Ip):

	Dn = shortname(align)
	for n in Dn: ##to short name in alignment
		label = len(Dn[n]) - len(n)
		label = n+' '*label
		align = align.replace(Dn[n],label)
		#print align
	f = open('TEMPseqs','w')
	f.write(align)
	f.close()

	if Ip == '0':
	        arg = "echo \'##########\n\' ; muscle -clwstrict -gapopen __GAP__ -maxiters 32 -quiet"
		arg = arg.replace("__GAP__",GAP)

		process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process.communicate(align)
		#x = x[0]

	elif Ip == '1':

		arg = """echo \'##########\n\nCLUSTAL W (1.81) multiple sequence alignment' ; mafft --clustalout --quiet --maxiterate 1000 --localpair TEMPseqs | grep -v "CLUSTAL format alignment by MAFFT" | tr "[:lower:]" "[:upper:]" """

		process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process.communicate()
		#x = x[0]

	elif Ip == '2':
		arg = """clustalw -infile=TEMPseqs -outorder=INPUT -quiet -outfile=TEMPseqs.aln"""

		process = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process.communicate()
		del x

		arg = """echo \'##########\n\nCLUSTAL W (1.81) multiple sequence alignment' ;cat TEMPseqs.aln | grep -v "CLUSTAL 2.1 multiple sequence alignment"  """

		process2 = subprocess.Popen(arg, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		x = process2.communicate()
		
		

	alignout = x[0]
	if Ip == '0':
		name = "S1S"
		alignout = AlignReOrder(alignout,name)
		alignout = alignout + '\n\n'

	if Ip == '0' or Ip == '1' or Ip == '2':
		
		for keys in Dn.keys():
			Dn[keys] = Dn[keys].split(' ')[0][:]

		Ilngth = max([len(x) for x in Dn.values()])
#		print Ilngth  ##DEBUG
		Klngth = max([len(x) for x in Dn.keys()])
#		print Klngth  ##DEBUG
		for n in Dn: ##to short name in alignment
			label = Ilngth - len(Dn[n])
			label = Dn[n]+' '*label
			alignout = alignout.replace(n,label)
		Saster = '\n '+(' '*(Ilngth - Klngth))
		alignout = alignout.replace('\n ',Saster)
		#alignout = Salignout.replace('#',' ')
		
#	sys.stderr.write(str(x[1])) ##DEBUG
	del process
	return	alignout
        
## ++++++++++++++++++++++

def seqLIM2(dbpath):

	arg = """grep '>' __dbpath__ | cut -f 1 -d ' ' | sed 's/>//' """
	arg = arg.replace('__dbpath__',dbpath)
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()[0]

	Lname = x.split('\n')
	Lname = [x.strip() for x in Lname if x.strip() != '']

	Tseqlen = {}
	for i in Lname:
		
		Tseqlen[i.replace('lcl|','')] = fastacmdcall2(dbpath,i)
	del process
	return Tseqlen

def seqLIM3(dbpath,col):
	
	arg = """sort -u -k """ + col + "," + col + " " + blast +  """ | awk '{print $"""+ col +"""}' """
	#arg = arg.replace('__dbpath__',dbpath)
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()[0]

	if col == '1':
		Lname = x.split('\n')
		#Lname = [x.replace('lcl|','').split('|')[0] for x in Lname if x.strip() != '']
		Lname = [x.strip().replace('lcl|','') for x in Lname if x.strip() != '']
		Lname = [x.split('|')[0] for x in Lname]
		

	else:
		Lname = x.split('\n')
		Lname = [x.strip() for x in Lname if x.strip() != '']

	
	Lname2 = []
	for i in Lname:
		if i not in Lname2:
			Lname2.append(i)

		#sys.stderr.write(str(i)+"\n") ##DEBUG
	Tseqlen = {}
	for i in Lname2:
		#ti = i.replace('lcl|','').split('|')[0]
		ti = i
		#sys.stderr.write(ti+'\n') # DEBUG
		Tseqlen[ti] = fastacmdcall2(dbpath,i)
	del process
	return Tseqlen

def fastacmdcall2(dbpath,name):

	sys.stderr.write(dbpath+" "+name+"\n")  ##DEBUG
	arg = """fastacmd -d __dbpath__ -s "__name__" -L 0,1000000000 > /dev/null"""
	arg = arg.replace('__dbpath__',dbpath); arg = arg.replace('__name__',name)
	process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	x = process.communicate()[1]

	seq_lim = x[:]
        seq_lim = seq_lim.split("ERROR: From location cannot be greater than ")[1]
        seq_lim = int(seq_lim.split(". Ignoring sequence")[0])

	del process
	return seq_lim

def fastacmdcall(dbpath,name):

	arg = """fastacmd -d __dbpath__ -s "__name__" -L 0,1000000000 2> TEMPseq | wc > /dev/null"""
	arg = arg.replace('__dbpath__',dbpath); arg = arg.replace('__name__',name)
	os.system(arg)

	sf = open('TEMPseq','r')
	seq_lim = sf.read()
	sf.close()

	sys.stderr.write(seq_lim)
	seq_lim = seq_lim.split("ERROR: From location cannot be greater than ")[1]
	seq_lim = int(seq_lim.split(". Ignoring sequence")[0])
	
	return seq_lim


def shortname(align):
	N = (align.split('\n'))
	N = [x.split()[0][1:] for x in N if len(x) > 0 and x[0] == '>']
	#N = [x.split()[0][1:] for x in N if len(x) > 0 and x[0] == '>']
	#N = [str(x.strip()) for x in N if len(x) > 0 and x[0] == '>']


	num = 1
	Dn = {}
	for x in N:
		label = 'S'+str(num)+'S'
		Dn[label] = x
		num += 1

	return Dn

def AlignReOrder(align,name): ##to fix muscle -stable order option
	StAL = align.split('\n\n')
	NewL = [StAL[0]] 
	for i in StAL[1:]:
		xi = i.splitlines()
		Target = [x for x in xi if x.find(name) == 0]
		Related = [x for x in xi if x not in Target] 
		OrderedL = Target + Related
		OrderedL = '\n'.join(OrderedL)
		NewL.append(OrderedL)
		#print NewL
	NewL = '\n\n'.join(NewL)
	return NewL

###############################


def main2(DB,query):

	f = open(blast,'r')
	blasthit = f.read()
	f.close()

	nameL = contigsName(blasthit)
	if logname != '0':
		k = open(logname,'w')
		k.write('\n'.join(nameL))
		k.close()
		
#	nameL = nameL[0:2] ## for testing only ##DEBUG
	if outname != "0":
		fout = open(outname,'w')
#	print "nameL=",[x for x in nameL]##DEBUG
	sys.stderr.write("enter MAIN loop\n")
	sys.stderr.write("nameL length ="+str(len(nameL))+" \n")
	for i in nameL[:]:

		hits = contigJuice(blasthit,i)
		hits = blastTrans(hits)
		hits = [x for x  in hits if x[1] != ban] ### exclude ban sequences
		if grouping == "1":clustL = merge4(hits,coverage)
		elif grouping == "2":clustL = Rmerge2(hits,coverage)
		
		#clustL = clustL[0:20] ##<<<<<<<<<< set for testing <<<<<<< ##DEBUG
		x = 0
		sys.stderr.write("clustL length ="+str(len(clustL))+" \n")
		for i in clustL[:]:
			#sys.stderr.write("#FOR I IN CLUSTL #####\n") ##DEBUG
			#sys.stderr.write("each clust="+str(i)+"\n<<<<<<<<<<<<<<<<<<<<<<<<<DEBUG") ##DEBUG
			xi = theone(i)
			#sys.stderr.write("\n"+str(xi)+"\n<<<<<<<<<<<<<<<<<<<<<<<<<DEBUG") ##DEBUG
			#CHECK DIRECTION IN TRIM FUNCTION#
			if mode != "0":					
			#	xi = theone(i)
				xi = bitselect(xi,mode)
				#sys.stderr.write("\n"+str(xi)+"\n<<<<<<<<<<<<<<<<<<<<<<<<<DEBUG") ##DEBUG
				#xi = chop2(xi)[0]
				#print len(i) ##DEBUG
			xi = [x for x in xi if abs(x[7] - x[6] ) > 48]	
			Isnum = len(xi)	
			sys.stderr.write(str('\n'.join([str(x) for x in xi]))+'\n')
			sys.stderr.write(str('\n'.join([str(x[6:10]+[x[7] - x[6]]) for x in xi]))+'\n')
			if Isnum >= (Iseqnem - 1):
				#sys.stderr.write("$$$$$$$$$$$$$$$$$$$$SEP"+'\n') ##DEBUG
				#for j in xi:sys.stderr.write(str(j)+'\n') ##DEBUG
				align = getAlign3(DB,query,xi)
				Saln = MSAgen5(align,Ip)
				if outname == "0":print Saln
				else:fout.write(Saln)

	if outname != "0":
		fout.close()


#####################################################################
####  Variable define part      #####################################

usage = """-i blast_result -r related_database -t target_database"""
opt = optparse.OptionParser(usage)
opt.add_option("-r",help="*indicate related database path")
opt.add_option("-t",help="*indicate target database path")
opt.add_option("-i",help="*indicate input file name")
opt.add_option("-l",help="indicate logfile path", default='0')
opt.add_option("-g",help="indicate muscle Gap penalties as float [default = -400.0]")
opt.add_option("-c",help="indicate coverage percentage (0.0 to 1.0) as float")
opt.add_option("-o",help="indicate output file name or print out as standard output",dest="o",default="0")
opt.add_option("-e",help="indicate E-value cutoff")
opt.add_option("--grouping",help="grouping method 1=grouping by longest sequences; 2=grouping by shortest sequence then trim default=2",default="2",dest="grouping")
opt.add_option("-m","--mode",help="bit selecting mode 0=only rnazWindow 1=least bit, 2=average bit, 3=median length, 4=most similar, 5=middle of sequences, 6=similar length", dest="m", default='2')
opt.add_option("-n","--seqnum",help="minimum number of sequenceses in alignments ", dest="n", default='2')
opt.add_option("-a","--align",help="select MSA program; default = 0 [MUSCLE]  1 = MAFFT local 2 = clustalw", dest="a", default='0')
opt.add_option("-b","--ban",help="excluded sequence name in database to prevent selecting itself", dest="b", default='0')
(options, args) = opt.parse_args()

DB = options.r
query = options.t
blast = options.i
mode = options.m
Ip = options.a
outname = options.o
logname = options.l
Iseqnem = int(options.n)
grouping = options.grouping
ban = options.b
ban = (ban.strip('"')).strip("'")

if None in (DB,query,blast):
	opt.print_help()
	sys.stderr.write('missing input(s)\n')
	sys.exit()

if options.c != None: coverage = float(options.c) ## coverage parameter
else:coverage=0.5 
if coverage > 1.0 or coverage < 0.5:
	sys.stderr.write('"coverage" between 0 to 1\n')
	sys.exit()

if options.e != None and float(options.e) <= 10.0 : Eval = float(options.e) ## E-value parameter
elif options.e == None:Eval = 10.0
else:
	sys.stderr.write('"E" must be float\n')
	sys.exit()

if options.g != None: GAP = str(float(GAP)) ## output name
elif options.g == None:GAP= "-400.0"
if float(GAP) > 0:
	sys.stderr.write('"gap cost" must be negative\n')
	sys.exit()

if mode not in ("0","1","2","3","4","5","6"):
	sys.stderr.write("selecting mode incorrect\n")
	sys.exit()

if Ip not in ('0','1','2'):
	sys.stderr.write("incorrect selecting MSA program\n")
	sys.exit()
if grouping not in ('1','2'):
	sys.stderr.write("incorrect grouping mode\n")
	sys.exit()
#sys.stderr.write("enter Dseqlim loop\n") # DEBUG
Dseqlim = seqLIM3(query,'1')
Rseqlim = seqLIM3(DB,'2')
#sys.stderr.write("exit Dseqlim loop\n") # DEBUG

#sys.stderr.write(str(Dseqlim)+'\n') # DEBUG
if logname != '0':
	arg = "date >> " + logname
	os.system(arg)

#sys.stderr.write("enter Dbit loop\n") # DEBUG
Dbit = bitcount2(blast)
#sys.stderr.write("exit Dbit loop\n") # DEBUG

main2(DB,query)

