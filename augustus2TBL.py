#!/usr/bin/python
##augustus2TBL.py

import sys,optparse

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
	ToShewL = []
	for i in SelectedTasteL:
		DeliciousPieceL = [x for x in CutL if x[0][2] == i] 
		ToShewL.append(DeliciousPieceL)

	return ToShewL


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
			self.DishL.append(groupL)

class Piece(object):
	def __init__(self,gffRow):
		self.gffRowL = gffRow
		self.pieceL = [self.gffRowL[0:8],self.shew(self.gffRowL)]
		
	def shew(self,Pulp):		
		LastColumnL = [x.strip()+'"' for x in Pulp[8].split('";')]
		if len(LastColumnL) > 1:
			TEMPL = [x.split() for x in  LastColumnL if len(x.split()) > 1]
			KL,VL = [x[0] for x in TEMPL ],[' '.join(x[1:]) for x in TEMPL] #generate Key and Value
			tasteD =  { k:v for k,v in zip(KL,VL)} ##gerate dictionary between ID and value
		else:
			tasteD = {Pulp[2]+"_id":LastColumnL[0]}
		return tasteD


def TestMain():

	chopedL = chop(Gffs)[:] 

	SelectedTasteL = ['gene', 'transcript', 'intron', 'CDS'] 
	ContigTagS = ''
	for pieceE in chopedL: 
		FruitOBJ = Fruit(pieceE)
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
		

Gffs = sys.stdin.read()
TestMain()

