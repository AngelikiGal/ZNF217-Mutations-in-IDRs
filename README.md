# ZNF217-Mutations-in-IDRs

Authors: Panagiota-Aggeliki Galliou (ag.gal.work@gmail.com)

This tool identifies point mutations and short linear motifs regions (SLiMs) in the disordered regions of a protein (IDRs). Firstly, it predicts potential SLiMs and residues in a protein using given protein binding motifs (like phosphorylation motifs for kinases, ligation motifs for ligases, etc). Then, it reads the protein's IDRs from an input file, keeps only the the predicted SLiMs that are located inside of IDRs and categorize them according to the protein and its binding motif that used for the prediction of each SLiM. Then, it reads the mutated residues found in the protein from mass spectrometry, identifies those that are located inside of IDRs and categorize them according to the predicted SLiM they could potentially obstruct. 


# Requirements

This tool is written in Python 2 and therefore, it can run in most operating systems.

# How to run

## Needed files

PhosphoKin uses as input two .txt files and one .xlsx file; one containing the protein sequence (.txt), one containing the protein’s IDRs (.txt) and one containing the mutated residues found in the protein (.xlsx).


### 1.) The Protein Sequence file:
The file must have a faste format (https://en.wikipedia.org/wiki/FASTA_format)

e.g. 

>sp|O75362|ZN217_HUMAN Zinc finger protein 217 OS=Homo sapiens OX=9606 GN=ZNF217 PE=1 SV=1
MQSKVTGNMPTQSLLMYMDGPEVIGSSLGSPMEMEDALSMKGTAVVPFRATQEKNVIQIE
GYMPLDCMFCSQTFTHSEDLNKHVLMQHRPTLCEPAVLRVEAEYLSPLDKSQVRTEPPKE
KNCKENEFSCEVCGQTFRVAFDVEIHMRTHKDSFTYGCNMCGRRFKEPWFLKNHMRTHNG

### 2.) The IDRs file:
The file must have the following format: [Start_of_IDR_1]-[End_of_IDR_1], [Start_of_IDR_2]-[End_of_IDR_2], etc.

e.g.

>[1]-[51], [107]-[122]

### 2.) The Mutations file:
The file must have two columns. The first will contain the mutated residues and the second the Functional impact. In the first row the names of each column: cell [0,1] should be "Mutations" and cell [0,2] should be "Functional Impact"

e.g.

> |Mutations| Functional Impact|
> |---------| -----------------|
> |N8H| deleterious|
> |Y17F| tolerated|

where the 8th residue was N and was mutated to H and this mutation had a deleterious impact for the protein.


### Note:
It should be noted that the proteins and their motifs are already included in the tool, and in order to use this tool for different kinases the code should me altered. 
The kinases and their motifs included in this tool are displayed below:

>MAPK: [KR]{0,2}[KR].{0,2}[KR].{2,4}[ILVM].[ILVF], ...[ST]P..  
>MAPKAPK1: Kinase, "[RK].R..S",mod  
>MAPKAPK2: Kinase, "S...[ST]",mod  
>ERK1: Kinase, "[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]",doc,"[RK].{2,4}[LIVMP].[LIV].>[LIVMF]",doc,"SP",mod,"..SP",mod,".[ST]P",mod  
>ERK2: Kinase, "[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]",doc,"[RK].{2,4}[LIVMP].[LIV].>[LIVMF]",doc,"SP",mod,"..SP",mod,".[ST]P",mod  
>GSK-3: Kinase, "..SP",mod,".[ST]P",mod,"...[ST]...[ST]",mod  
>P38: Kinase, "[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]",doc,"[RK].{2,4}[LIVMP].[LIV].[LIVMF]",doc))
>JNK: Kinase, "[RK]P[^P][^P]L.[LIVMF]",doc  
>CDK: Kinase, "...[ST]P[RK]",mod,"...[ST]P..[RK]",mod,"SP.[RK].",mod,"[ST]P.[RK]",mod  
>CKI: Kinase, "S..[ST]...",mod,"[ST]..[ST]",mod,"[ED]..[ST]",mod,"[ST]...[ST][M/L/V/I/F]",mod,"SP..[ST]",mod  
>CKII: Kinase,"S..[EST]",mod,"[ST]..[EDSY]",mod,"S..[ED]",mod,"[ST]..[ED]",mod,"S.[EST]",mod  
>NEK-2: Kinase, "[FLM][^P][^P][ST][^DEP][^DE]",mod  
>PIKK: Kinase,"...[ST]Q..",mod  
>ATM: Kinase,"SQ",mod,"LSQE",mod  
>DNA-PK: Kinase, "P[ST].",mod,".SQ",mod  
>PKA: Kinase,"[RK][RK].[ST][^P]..",mod,".R.[ST][^P]..",mod,"[RK].[ST]",mod,"[ST].[RK]",mod,"K...[ST]",mod,"K..>[ST]",mod,"R..S",mod,"[RK][RK].[ST]",mod,"R.S",mod,"KR..S",mod  
>PKC Epsilon: Kinase,"R[KER].S",mod  
>PKC: Kinase,"[RK].[ST]",mod,"[ST].[RK]",mod,"[RK]..[ST]",mod,"[RK]..[ST].[RK]",mod,"[RK].[ST].[RK]",mod  
>PLK-1: Kinase,".[DNE][^PG][ST][FYILMVW]...",mod, ".[DNE][^PG][ST][^PEDGKN][FWYLIVM].",mod,"S[ST].",lig  
>PLK-2: Kinase,"[DE]..[ST][EDILMVFWY][DE].",mod, "[DE]..[ST][EDILMVFWY].[DE]",mod  
>PLK-3: Kinase,"[DE]..[ST][EDILMVFWY][DE].",mod, "[DE]..[ST][EDILMVFWY].[DE]",mod  
>PLK-4: Kinase,"..[^IRFW][ST][ILMVFWY][ILMVFWY]",mod  
>CaMKII: Kinase,"R..S",mod,"R..[ST]",mod,"[MVLIF].[RK]..[ST]..",mod  
>CaMKIV: Kinase,"[MILVFY].R..[ST]",mod  
>Chk1: Kinase,"[MILV].[RK]..[ST]",mod  
>CLK-1: Kinase,"R..[ST]..R",mod  
>CLKB-1: Kinase,"LRT",mod  
>mTOR: Kinase,"FTY",mod  
>NIMA: Kinase,"FR.[ST]",mod  
>EGFRK: Kinase,".[ED]Y.",mod,".[ED]Y[ILV]",mod  
>JAK2: Kinase,"Y..[LIV]",mod  
>SRC: Kinase,"[IVLS].Y..[LI]",mod,"Y[AGSTED]",mod  
>Itk: Kinase,"Y[AEV][YFESNV][PFIH]",mod  
>BRCA1: Ligase, ".S..F",doc,".S..F.K	",doc  
>APC/C: Ligase,".KEN.",deg  
>SPOP: Ligase,"[AVP].[ST][ST][ST]",deg  
>PP2B: Phosphatase, "L.[LIVAPM]P",doc  
>PP2C-delta: Phosphatase, ".T.Y.",mod  
>SHP-1: Phosphatase,"[ED].Y",mod, "[IV].Y..[LV]",lig  
>TC-PTP: Phosphatase,"[EDY]Y",mod  
>CtBP1: Dehydrogenase, "P[LVIPME][DENS][LM][VASTRG]", lig, "G[LVIPME][DENS][LM][VASTRG]K", lig, "G[LVIPME][DENS][LM][VASTRG].[KR]", lig  
>Caspase-3: Protease,"[DSTE][^P][^DEWHFYC]D[GSAN]",clv  
>Caspase-7: Protease,"[DSTE][^P][^DEWHFYC]D[GSAN]",clv  
>USP7: Deubiquitinase, "[PA][^P][^FYWIL]S[^P]",doc,"	K...K",doc  
>Pin1: Isomerase,"...[ST]P.",doc  
>SUMO: ETC,"[DEST]{0,5}.[VILPTM][VIL][DESTVILMA][VIL].{0,1}[DEST]{1,10}",lig  
>SUMO-1: ETC,"[SDE].{0,5}[DE].K.{0,1}[AIFLMPSTV]",mod,"[DEST]{0,5}.[VILPTM][VIL][DESTVILMA][VIL].{0,1}[DEST]{1,10}",lig  
>FHA Domain: ETC,"..T..[ILV].",lig,"..T..[DE].",lig  
>MYND: ETC,"P.L.P",lig  
>TRF1-2 Domain: ETC,"[FY].L.P",lig  
>NAE1-UBA3: ETC,"[ILM][ILMF].{1,2}[ILM].{0,4}K",lig  
>Importin-A: ETC,"[PKR].{0,1}[^DE]K[RK][^DE][KR][^DE]", mod, "[PKR].{0,1}[^DE]K[RK][KR][^DE][^DE]", mod, "[PKR].{0,1}[^DE]RK[^DE][KR][^DE]",mod, "[PKR].{0,1}[^DE]RK[KR][^DE][^DE]", mod, "[PKR]K[RK][^DE][KR][^DE]", mod, "[PKR]K[RK][KR][^DE][^DE]", mod,  "[PKR]RK[^DE][KR][^DE]",mod, "[PKR]RK[KR][^DE][^DE]", mod  
>MDC-1: ETC,"S[ST].",lig  
>Cyclin A: ETC,".[^EDWNSG][^D][RK][^D]L.{0,1}[FLMP].{0,3}[EDST]",doc, "[KRH].{0,3}[^EDWNSG][^D][RK][^D]L.{0,1}[FLMP].{0,3}[EDST]", doc  
>14-3-3: ETC,"R[^DE]{0,2}[^DEPG][ST][FWYLMV].",doc, "R[^DE]{0,2}[^DEPG][^PRIKGN]P", doc, "R[^DE]{0,2}[^DEPG][^PRIKGN].{2,4}[VILMFWYP]", doc, "R..S", doc  
>CDC20: ETC,"[ILVMF].[ILMVP][FHY].[DE]",doc  
>GRB-1: ETC,"Y.N",mod,"Y[MILV].[MILV]",mod  

First is the name of the protein, after the ":" character is the type of the protein (ETC means that the protein didn't fit in the other categories or that it is a domain), and inside the characters "" are the protein motifs divided by commas.  After each motif there is a description of the motifs function (doc: docking, clv: cleavage, mod: modification, lig: ligation, deg: degradation). 


# In Terminal

`python2 ZNF-217_ALL.py`

This tool will not ask from the user to separately give each file needed as an input. The needed input files should be present in the same folder as the tool before the running. The names of the input files should be:  "ZNF271_Sequence.txt",  "ZNF217_disordered_regions.txt" and "ZNF217_mutations.xlsx"

# The output

As output the tool produces four .txt files; one containing the predicted SLiMs in the protein, one containing those SLiMs that are located inside the IDRs of the protein, one containing the mutated residues of the protein that are located inside of the IDRs of the protein, and one containing the mutated residues inside of the protein's IDRs that actually could obstruct the predicted SLiMs of the protein. 

# Got a Question?

Please do not hesitate to ask any question via email: Panagiota-Angeliki Galliou (ag.gal.work@gmail.com).
