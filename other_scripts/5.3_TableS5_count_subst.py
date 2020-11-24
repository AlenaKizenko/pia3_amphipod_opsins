from Bio import AlignIO
#import ete3
#from Bio import SeqIO
#from Bio.SeqIO.FastaIO import SimpleFastaParser
#from Bio.SeqRecord import


alignment = "Baikal_and_MWS_Gammaridae_all_opsins_Bta_Tpa.aln" 
ref_seq_name = "RHO_Bos_taurus_AAA30674.1"


aln_list = AlignIO.read(open(alignment), "fasta")


## Sanity check: divide all sequences into groups (manually)
idlist = []
for record in aln_list:
    idlist.append(record.id)
print("0:2 = Outgroup:", idlist[0:2])
print("2:103 = Baikal; ", "\n from:", idlist[2:3] ,"\n to:", idlist[103:104])
print("104:117 = LWS from MWS+; ", "\n from:", idlist[104:105] ,"\n to:", idlist[117:118])
print("118:123 = MWS; ", "\n from:", idlist[118:119] ,"\n to:", idlist[123:124])


## Count positions in B. taurus for easier numbering
letter_counter = 0; gap_counter = 0; isGap = False

refPos = []
refSeq = aln_list[0].seq
for letter in refSeq:
    if letter == "-":
        gap_counter += 1
        isGap = "Gap"
    else:
        letter_counter += 1
        isGap = ""
    refPos.append(str(letter_counter) + isGap)
#print(refPos)

## And in squid also count
letter_counter = 0; gap_counter = 0; isGap = False
refPosSquid = []
refSeqSquid = aln_list[1].seq
for letter in refSeqSquid:
    if letter == "-":
        gap_counter += 1
        isGap = "Gap"
    else:
        letter_counter += 1
        isGap = ""
    refPosSquid.append(str(letter_counter) + isGap)



### Iterate over columns and check amino acids
for pos in range(len(aln_list[0])):
    #print(pos)
    column = aln_list[:, pos]
    BaikalSet = set(column[2:104])
    MWSplusSet = set(column[104:118])
    MWSSet = set(column[118:124])
    for aa in BaikalSet:
        aacnt = column[2:104].count(aa)
        if not aa == "-" and aa in MWSSet and not aa in MWSplusSet:
            print("MWS-like in Baikal only:", aa, ";", pos+1,
                  "in alignment; or", refPos[pos], "in B. taurus; or",
                  refPosSquid[pos], "in Todarodes pacificus; ", 
                  "found ", aacnt, " times in Baikal seqs",
                  BaikalSet, MWSplusSet, MWSSet)

cnt = 0

## Small switches to check co-occurrence of different substitutions

#for record in aln_list:
#    if record.seq[255] == "A" and record.seq[295] == "A":
#        print(record.id)
#        cnt += 1


#for record in aln_list:
#    if record.seq[306] == "T" and record.seq[295] == "A":
#        print(record.id)
#        cnt += 1

#for record in aln_list:
#    if record.seq[306] == "T" and record.seq[255] == "A":
#        print(record.id)
#        cnt += 1

#for record in aln_list:
#    if record.seq[255] == "A" and record.seq[295] == "A" and record.seq[306] == "T":
#        print(record.id)
#        cnt += 1


#for record in aln_list:
#    if record.seq[255] == "A" and record.seq[339] == "W" and record.seq[306] == "T":
#        print(record.id)
#        cnt += 1

### 3 most common
for record in aln_list:
    if record.seq[255] == "A" and record.seq[295] == "A" and record.seq[306] == "T":
        print(record.id)
        cnt += 1


print(cnt)
