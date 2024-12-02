#python3

from Bio import SeqIO
from Bio.Seq import Seq
#from io import StringIO
from Bio import pairwise2
import re
import glob
import os
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
#from Bio.Align import PairwiseAligner 
    
def alignedSeqQual(AlignedSeq, QL):
    exL = ["",""] #seq, qual
    i = 0
    for letter in AlignedSeq:
       if (letter != "-"):
           exL[0].append(letter)
           exL[1].append(QL.pop(0))
    return exL

def primerSearch(seqA_str,seqB_str, shortlimit_I=3):
    query_S = seqA_str
    print("query_S",query_S)
    found_I = -1
    matchlen_I = 0
    while found_I == -1:
        found_I = seqB_str.find(query_S)        
        if (found_I != -1 and len(query_S) > shortlimit_I):
            print(found_I ,query_S)            
            if ( len(query_S) > matchlen_I):
                matchlen_I = len(query_S)

        query_S = query_S[1:]

    found_I = -1
    query_S = seqA_str    
    while found_I == -1:
        found_I = seqB_str.find(query_S)        
        if (found_I != -1 and len(query_S) > shortlimit_I):
            print(found_I ,query_S)
            if ( len(query_S) > matchlen_I):
                matchlen_I = len(query_S)

        query_S = query_S[:-1]
    
    found_I = -1
    query_S = str(Seq(seqA_str).reverse_complement())    
    print("query_S Reverse Complemented ",query_S)
    while found_I == -1:
        found_I = seqB_str.find(query_S)        
        if (found_I != -1 and len(query_S) > shortlimit_I):
            print(found_I ,query_S)
            if ( len(query_S) > matchlen_I):
                matchlen_I = len(query_S)

        query_S = query_S[1:]
    
    found_I = -1
    query_S = str(Seq(seqA_str).reverse_complement())    
    while found_I == -1:
        found_I = seqB_str.find(query_S)        
        if (found_I != -1 and len(query_S) > shortlimit_I):
            print(found_I ,query_S)
            if ( len(query_S) > matchlen_I):
                matchlen_I = len(query_S)

        query_S = query_S[:-1]
    return found_I

def cutoff(intI,ct=0):
    strX = "-"
    if (intI < ct):
        strX = "!"
    return strX

output_str = ""

def Oprint(*InStr):
    global output_str 
    output_str += "".join(InStr) + "\n"
## Variable definition
#input_file = "seqF.ab1"
ab1FL = glob.glob("*.ab1")
for f in ab1FL:
    print(f)
    Oprint(f)

#seqF_file = "seqF.ab1"
output_fname = "SangerTrimOUT.txt"

if len(ab1FL) < 2:
    output_str = "INCORRECT INPUT OR MISSING FILE"
    print(output_str)
    f = open(output_fname,'w')
    f.write(output_str)
    f.close()
    exit()

seqF_file = ab1FL[0]
seqR_file = ab1FL[1]


seqF_obj = SeqIO.parse(seqF_file, "abi")  
seqF_seq = str([x for x in seqF_obj][0].seq)
seqF_score_L = [x for x in SeqIO.parse(seqF_file, "abi")  ][0].letter_annotations["phred_quality"]


seqR_obj = SeqIO.parse(seqR_file, "abi")  
seqR_seq = str([x for x in seqR_obj][0].seq.reverse_complement())
seqR_score_L = [x for x in SeqIO.parse(seqR_file, "abi")  ][0].letter_annotations["phred_quality"][::-1]

primerFL = glob.glob("*_PRIMER.fasta")
for f in primerFL:
    print(f)
    Oprint(f)

#primer_file = "adapterF.fasta"
primer_file = primerFL[0]
#primer_fileB = "adapterR.fasta"
primer_fileB = primerFL[1]

primer_obj = SeqIO.parse(primer_file, "fasta")
primer_seq = str([x for x in primer_obj][0].seq)
#for o in adapter_obj:
#    adapter_seq = str( o.seq )

primer_objB = SeqIO.parse(primer_fileB, "fasta")
primer_seqB = str([x for x in primer_objB][0].seq)
print("+++++++ INPUT ++++++ ")
print("Forward Sequence")
print(str([x for x in SeqIO.parse(seqF_file, "abi")][0].seq))
print("Reverse Sequence")
print(str([x for x in SeqIO.parse(seqR_file, "abi")][0].seq))
print("Forward Primer")
print(primer_seq)
print("Reverse Primer")
print(primer_seqB)

Oprint("+++++++ INPUT ++++++ ")
Oprint("Forward Sequence")
Oprint(str([x for x in SeqIO.parse(seqF_file, "abi")][0].seq))
Oprint("Reverse Sequence")
Oprint(str([x for x in SeqIO.parse(seqR_file, "abi")][0].seq))
Oprint("Forward Primer")
Oprint(primer_seq)
Oprint("Reverse Primer")
Oprint(primer_seqB)

"""
print("++++++ seF_seq")
#for i in range(len(seqF_seq))[50]:
for i in range(0,50):
    print(seqF_seq[i] + "  ",end="")

print('')
print('+')

#for i in range(len(seqF_seq))[50]:
for i in range(0,50):
    if (seqF_score_L[i] < 10):
        print( str(seqF_score_L[i]) + "  ", end ="")
    else:
        print( str(seqF_score_L[i]),"", end ="")
print('')

print(seqF_seq)
#print("\033[1;32m "+ seqF_seq +"  \n")
print("++++++ seqR_seq")
print(seqR_seq)
"""
#seqR_len_I = len(seqR_seq)
#print(seqR_len_I/2)

#search_part_str = seqR_seq[:int(seqR_len_I/2)]
#print(search_part_str)

#alignments = pairwise2.align.globalms(seqF_seq,seqR_seq,5,-4,-5,-.1)
alignments = pairwise2.align.localms(seqF_seq,seqR_seq,5,-4,-5,-.1)

#for align in alignments:
    #print(pairwise2.format_alignment(*align))


#print(pairwise2.format_alignment(*alignments[0]))
#print(alignments[0])
#print(alignments[0].score)

"""
print(alignments[0].seqA[0:alignments[0].start + 0])
gap_F_numI = alignments[0].seqA[0:alignments[0].start + 0].count("-")
print(gap_F_numI)

seqB_str = alignments[0].seqB[0:alignments[0].start + 0]
print(alignments[0].seqB[0:alignments[0].start + 0])
gap_R_numI = alignments[0].seqB[0:alignments[0].start + 0].count("-")
print(gap_R_numI)
"""

print("+++++++ INPUT END ++++++ ")
Oprint("+++++++ INPUT END ++++++ ")

alignLenI = int(alignments[0].end)
half_I = int(alignLenI/2)
#print(alignments[0].seqA[(half_I - 5 ):(half_I + 5 )])
#print(alignments[0].seqB[0:half_I],end="")
#print(alignments[0].seqA[half_I:])

seqA_score_L = []
Lnum = 0
for c in str(alignments[0].seqA):
    if (c != "-"):
        seqA_score_L.append(seqF_score_L[Lnum])
        Lnum += 1

    else:
       seqA_score_L.append(0)

seqB_score_L = []
Lnum = 0
for c in str(alignments[0].seqB):
    if (c != "-"):
        seqB_score_L.append(seqR_score_L[Lnum])
        Lnum += 1
    else:
        seqB_score_L.append(0)

print("+++++++ ALIGNED ++++++++++++++  ")
Oprint("+++++++ ALIGNED ++++++++++++++  ")

#print("score A:")
print( "".join([ str(int(x/10)) for x in seqA_score_L]))
Oprint( "".join([ str(int(x/10)) for x in seqA_score_L]))

#print()
#print(pairwise2.format_alignment(*alignments[0]))
print(alignments[0].seqA)
Oprint(alignments[0].seqA)
#print()
print(alignments[0].seqB)
Oprint(alignments[0].seqB)
#print()
print( "".join([ str(int(x/10)) for x in seqB_score_L]))
Oprint( "".join([ str(int(x/10)) for x in seqB_score_L]))
#print("score B")

concat_seq_str_OBJ = Seq(alignments[0].seqB[0:half_I] + alignments[0].seqA[half_I:])
print("+++++++ JOINED ++++++++++++++  ")
Oprint("+++++++ JOINED ++++++++++++++  ")

print("+++++ Q-Score Select ++++")
print("len(alignments[0].seqA[a])")
print(len(alignments[0].seqA))
print("len(alignments[0].seqB[a])")
print(len(alignments[0].seqB))
print("alignLenI")
print(alignLenI)

gapped_align_ken_I = len(alignments[0].seqB)
BestScoreFullSeq_str = ""
BestScoreFullScore_L = []
for a in range(gapped_align_ken_I):
    curr_alphabet = str(alignments[0].seqB[a])
    curr_score = seqB_score_L[a]
    if (seqA_score_L[a] > seqB_score_L[a]):
        curr_alphabet = str(alignments[0].seqA[a])
        curr_score = seqA_score_L[a]
    
    BestScoreFullSeq_str += curr_alphabet
    BestScoreFullScore_L.append(curr_score)

print( "".join([ str(int(x/10)) for x in BestScoreFullScore_L]))
print(BestScoreFullSeq_str)
cutoff_F = 40.0
print( "".join([ cutoff(x,cutoff_F) for x in  BestScoreFullScore_L]))
print( "Low Quality Bases (" , cutoff_F ,")=", len([ x for x in  BestScoreFullScore_L if ( x < cutoff_F ) ]) )

target_score = 50.0
window_len = 3
curr_pos_I = 0
trimStartPos_I = 0
trimEndPos_I = 0
trim_score_F = 0.0
while( (trim_score_F < target_score) and (curr_pos_I <= (gapped_align_ken_I - window_len)) ):
    scoresL = [x for x in BestScoreFullScore_L[curr_pos_I:( curr_pos_I + window_len ) ]]
    trim_score_F = float( sum(scoresL) )/(window_len)
    curr_pos_I += 1

trimStartPos_I = curr_pos_I - 1
print("5' Trim site",trimStartPos_I, " Mean Score=", trim_score_F)

curr_pos_I = (gapped_align_ken_I - window_len)
trim_score_F = 0
while( (trim_score_F < target_score) and (curr_pos_I > window_len) ):
    scoresL = [x for x in BestScoreFullScore_L[curr_pos_I:( curr_pos_I + window_len ) ]]    
    trim_score_F = float( sum(scoresL) )/(window_len)    
    curr_pos_I -= 1

trimEndPos_I = curr_pos_I + window_len + 1
print("3' Trim site",trimEndPos_I, " Mean Score=", trim_score_F)

print("trimmed_seq")
print("".join(BestScoreFullSeq_str[trimStartPos_I:trimEndPos_I]))
print( "".join([ str(int(x/10)) for x in  BestScoreFullScore_L[trimStartPos_I:trimEndPos_I]]))

cutoff_F = 40.0
print( "".join([ cutoff(x,cutoff_F) for x in  BestScoreFullScore_L[trimStartPos_I:trimEndPos_I]]))
print( "Low Quality Bases (" , cutoff_F ,")=", len([ x for x in  BestScoreFullScore_L[trimStartPos_I:trimEndPos_I] if ( x < cutoff_F ) ]) )


msa_obj = MultipleSeqAlignment(
        [  SeqRecord( Seq(BestScoreFullSeq_str), id="consensus"),
            SeqRecord( Seq(alignments[0].seqA), id="SeqA"),
            SeqRecord( Seq(alignments[0].seqB), id="SeqB"),
        
    ])

#print(msa_obj.alignment)


print("concatenated_seq")
Oprint("concatenated_seq")
print(str(concat_seq_str_OBJ))
Oprint(str(concat_seq_str_OBJ))

target_score = 5.0
ali_len = len(str(concat_seq_str_OBJ))
window_len = 3
start_pos = 0
end_pos = 0

#trimim 5prime end by calculate mean of 3 bases phredscore/10 if greater than target then mark the position
fiveprimePos = 0
score_F = 0.0
while( (score_F < target_score) and (fiveprimePos <= (ali_len - window_len)) ):
    score_F = (( seqB_score_L[fiveprimePos] + seqB_score_L[fiveprimePos + 1] + seqB_score_L[fiveprimePos + 2]  )/window_len)/10
    fiveprimePos += 1

start_pos = fiveprimePos - 1
print("5' Trim site",fiveprimePos, " Mean Score=", score_F)
Oprint("5' Trim site", str(fiveprimePos), " Mean Score=", str(score_F))

fiveprimePos = (ali_len - window_len)
score_F = 0.0
while( (score_F < target_score) and (fiveprimePos >= 0) ):
    score_F = (( seqA_score_L[fiveprimePos] + seqA_score_L[fiveprimePos + 1] + seqA_score_L[fiveprimePos + 2]  )/window_len)/10
    fiveprimePos -= 1
end_pos = fiveprimePos + window_len + 1
print("3' Trim site",fiveprimePos, " Mean Score=", score_F)
Oprint("3' Trim site",str(fiveprimePos), " Mean Score=", str(score_F))
print("++++++++++++++")
Oprint("++++++++++++++")
print("Trimmed sequence (Q score >= 50)")
Oprint("Trimmed sequence (Q score >= 50)")
print(str(concat_seq_str_OBJ)[start_pos:end_pos])
Oprint(str(concat_seq_str_OBJ)[start_pos:end_pos])

#Write output to a file
f = open(output_fname,'w')
f.write(output_str)
f.close()

"""
for i in range(len(concat_seq_str_OBJ)):
    print(str(concat_seq_str_OBJ)[i] + "  ",end="")
    

print("")
print("+++++++")
for i in range(0,half_I):
    curr_score_I = seqR_score_L[i]
    if (curr_score_I < 10):
        print(str(curr_score_I) + "  ",end="")
    else:
        print(str(curr_score_I) + " ",end="")

SeqF_GapCount_I = alignments[0].seqA.count('-')
partB_seq = alignments[0].seqA[half_I:]
sub_Tail_seqF_score_L = seqF_score_L[ (-1 * half_I):]
#print("###",end="")
for i in range(0,len(sub_Tail_seqF_score_L)  ):
    curr_I = sub_Tail_seqF_score_L[i]
    if (curr_I < 10):
        print(str(curr_I) + "  ",end="")
    else:
        print(str(curr_I) + " ",end="")
#print(concat_seq_str_OBJ)

print("+++++++")
print(len(alignments[0].seqA))
print(len(seqF_score_L))
print("B +++++++")
print(len(alignments[0].seqB))
print(len(seqF_score_L))
"""

#primerSearch(adapter_seq,str(concat_seq_str_OBJ))
#primerSearch(adapter_seqB,str(concat_seq_str_OBJ))

"""
pwaligner = PairwiseAligner("blastn")
pwaligner.mode = "global"


pwAlingmentG = pwaligner.align(seqF_seq,seqR_seq)
print(pwAlingmentG.score)
print("len(pwAlingmentG)")
print(len(pwAlingmentG))
print(pwAlingmentG[0])

print(pwAlingmentG[0].aligned)
print(len(pwAlingmentG[0].frequencies))
print(pwAlingmentG[0].frequencies.keys())



find_I = -1
query_str = adapter_seq
print(query_str)
while ( find_I == -1 and len(query_str) > 3):
    if (concat_seq_str_OBJ.find(query_str) != -1 ) :        
        print(query_str)

    query_str = query_str[1:]

query_str = adapter_seq
print(query_str)
while ( find_I == -1 and len(query_str) > 3):
    if (concat_seq_str_OBJ.find(query_str) != -1 ) :        
        print(query_str)
    
    query_str = query_str[:-1]

query_str = str(Seq(adapter_seq).reverse_complement())
print(query_str)
while ( find_I == -1 and len(query_str) > 3):
    if (concat_seq_str_OBJ.find(query_str) != -1 ) :        
        print(query_str)
    
    query_str = query_str[1:]

query_str = str(Seq(adapter_seq).reverse_complement())
print(query_str)
while ( find_I == -1 and len(query_str) > 3):    
    if (concat_seq_str_OBJ.find(query_str) != -1 ) :        
        print(query_str)

    query_str = query_str[:-1]

print("+++++++")

adapter_seq_obj = Seq.Seq(adapter_seq)
adapter_seq_RC_seq = adapter_seq_obj.reverse_complement()
print(adapter_seq_RC_seq)
adapter_obj_aln = pairwise2.align.localms( adapter_seq, seqB_str,5,-4,-5,-.1)
print(adapter_obj_aln[0])
print(pairwise2.format_alignment(*adapter_obj_aln[0]))
"""








