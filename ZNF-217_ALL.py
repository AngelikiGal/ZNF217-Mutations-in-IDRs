#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 14:35:46 2018

@author: angel
"""
import re
import pandas as pd

# ~ ~ classes ~ ~ 
class slim:
    def __init__(self):
        self.start = ""
        self.end = ""
        self.sequence = ""
class result_motif:
    def __init__(self):
        self.regex = ""
        self.Aslims = []
        self.interaction_type = ""  #what type of interaction this motif is. It takes doc (docking), mod(modification), lig (ligation), deg (degradation), clv(cleavage)

class result_protein:
    def __init__(self):
        self.name = ""
        self.type="" #kinase,ligase, phosphotase, dehydrogenase, protease, deubiquitinase, isomerase, etc
        self.Amotifs = []
        self.total_slims= ""
        

class protein:
    def __init__(self):
        self.name = ""  #name of the protein
        self.type = ""  #kinase,ligase, phosphotase, dehydrogenase, protease, deubiquitinase, isomerase, etc
        self.motifs = [] # list of motifs the protein has. it takes motif elements

class motif:
    def __init__(self):
        self.reg_expr = ""  #regular expression of the motif
        self.interaction_type = ""  #what type of interaction this motif is. It takes doc (docking), mod(modification), lig (ligation), deg (degradation), clv(cleavage)


class disordered_region:
    def __init__(self):
        self.start = 0  #start
        self.end = 0  #end
        

class mutation:
    def __init__(self):
        self.before = ""  #original amino acid
        self.after = ""  #mutated amino acid
        self.residue = "" #position of the mutated amino acid
        self.impact = "" #the impact of this mutation according to SIFT
        self.count = 1  #counts the times of the same mutation in the array
# ~ ~ classes ~ ~   
        

# ~ ~ functions ~ ~
def create_disordered_region_entry(array,region):
    ## ~ ~ ~ creates a disordered_region class ~ ~ ~ ##
    start,end=region.split("-",1)
    start=start.replace("[","")
    start=start.replace("]","")
    end=end.replace("[","")
    end=end.replace("]","")
            
    dire=disordered_region()
    dire.start=int(start)
    dire.end=int(end)
    
    array.append(dire)
    return
          
def create_slim(st,en,seq):
    s=slim()
    s.start=st     #shows the location of residue in the start of the slim
    s.end=en       #shows the location of residue in the end of the slim
    s.sequence=seq  #shows the sequence of the slim
    return s

def create_result_motif(regex,Aslim,interaction_type):
    rm=result_motif()
    rm.regex=regex    #shows the sequence of the motif
    rm.Aslims=Aslim   # an array of reults of slims identified by this motif. it holds slim objects
    rm.interaction_type=interaction_type
    return rm

def create_result_protein(name,Amotifs,total_slims,thetype):
    rp=result_protein()       
    rp.name=name         #name of a kinase
    rp.Amotifs=Amotifs   #an array of results of motifs for identified for this kinase. it holds result_motif objects.
    rp.total_slims=total_slims
    rp.type=thetype
    return rp
    
def create_motif(reg_expr,it):
    m=motif()
    m.reg_expr=reg_expr
    m.interaction_type=it
    return m
 
def create_protein (name,the_type,m1="",t1="",m2="",t2="",m3="",t3="",m4="",t4="",m5="",t5="",m6="",t6="",m7="",t7="",m8="",t8="",m9="",t9="",m10="",t10=""):
    #m1->motif 1. Maximum number of motifs = 10 for a protein.
    protein_motifs=[]
    
    #check how many motifs where given and for every motif given create an motif element and append it on an array
    if((m1 != "") and (t1 != "")):        
        protein_motifs.append(create_motif(m1,t1))       
    if((m2 != "") and (t2 != "")):
        protein_motifs.append(create_motif(m2,t2))   
    if((m3 != "") and (t3 != "")):
        protein_motifs.append(create_motif(m3,t3)) 
    if((m4 != "") and (t4 != "")):
        protein_motifs.append(create_motif(m4,t4)) 
    if((m5 != "") and (t5 != "")):
        protein_motifs.append(create_motif(m5,t5)) 
    if((m6 != "") and (t6 != "")):
        protein_motifs.append(create_motif(m6,t6))    
    if((m7 != "") and (t7 != "")):
        protein_motifs.append(create_motif(m7,t7)) 
    if((m8 != "") and (t8 != "")):
        protein_motifs.append(create_motif(m8,t8)) 
    if((m9 != "") and (t9 != "")):
        protein_motifs.append(create_motif(m9,t9)) 
    if((m10 != "") and (t10 != "")):
        protein_motifs.append(create_motif(m10,t10)) 
    
    """    
    for i in protein_motifs:
        print(i.reg_expr, i.interaction_type)
    """
    
    #create a protein element with a name, a type and the array with its motifs
    p=protein()
    p.name=name
    p.type=the_type
    p.motifs=protein_motifs
    #return the protein element
    return p

def create_mutation_entry(bef,aft,pos,imp,co):
    ## ~ ~ ~ creates a mutation class ~ ~ ~ ##
    mut=mutation()
    mut.before=bef
    mut.after=aft
    mut.residue=pos
    mut.impact=imp
    mut.count=co
    return mut
 

def get_sequence (filename):
    ## ~ ~ ~ Opens the Protein.txt file and gets protein sequence ~ ~ ~ ##
    seq=""
    protein_file= open(filename) # here you can put any protein sequence in fasta format
    if protein_file:
        print("Input file ",filename," opened.\n")
        for line in protein_file:
            line=line.strip()
            if ">" in line:
                continue
            else:
                seq+=line        
    else:
        raise IOError("Sequence file ",filename," could not open")
    
    protein_file.close()    
    #print(seq)        
    return seq

def get_disordered_regions(filename):
    ## ~ ~ ~ Opens the Protein_disordered_regions.txt file and gets the disordered regions ~ ~ ~ ##
    disorderded_regions=[]
    region=""
    seq=""
    protein_file= open(filename) 
    if protein_file:
        print("\nInput file ",filename," opened.\n")
        for line in protein_file:
            line=line.strip()
            seq+=line    
        
        if("," in seq):  #then there are more than one disoredered regions
            region,rest=seq.split(",",1)
            create_disordered_region_entry(disorderded_regions,region)            
                        
            while(len(rest)>2):            
                if("," in rest): #there is another disordered region                 
                    region,rest=rest.split(",",1)
                    create_disordered_region_entry(disorderded_regions,region)
                                
                else: #this is the last disordered region
                    create_disordered_region_entry(disorderded_regions,rest)
                    break
                                              
        else: #there is only one disordered region
            create_disordered_region_entry(disorderded_regions,seq)
        
    else:
        raise IOError("Sequence file ",filename," could not open")
    protein_file.close()                        
    return disorderded_regions

def get_mutations(output,outputname, filename):
    ## ~ ~ ~ reads an excel file and its columns(Mutations and Functional Impact), and compines them to create an mutation class for each row on the document ~ ~ ~ ##
    df = pd.read_excel(filename, sheetname='Sheet1')
    
    print("\nInput file ",filename," opened.\n")    
    
    mutations_E=df["Mutation"]
    impact_E=df["Functional Impact"]
    #print(mutations_E)
    #print(impact_E)

    
    all_mutations=[]
    
    num=0 
    for i in mutations_E:  # as arrays start from zero the first entry will be mutations_E[0] and will correspond to impact_E[0]
        before=i[0]
        after=i[-1]
        position=i[1:-1]
        impact=impact_E[num]
           
        all_mutations.append(create_mutation_entry(before,after,position,impact,1))
        num=num+1

    count_duplicates(all_mutations)
    
    print("Output file ", outputname, " opened for writting.")       
    protein_name=outputname.split("-",1)
    protein_name=str(protein_name[0])
    output.write("\n{:s} {:s} (Total={:d})\n".format(protein_name,"All Mutations:",len(all_mutations)))
        
    for x in all_mutations:
        output.write("{:s} {:s} {:s} {:s} {:d}\n".format(x.before,x.residue,x.after,x.impact,x.count))    

    all_deleterious, all_tolerated=count_impact(all_mutations)
    all_deleterious_percentage=int((float(all_deleterious)/len(all_mutations))*100)
    all_tolerated_percentage=int((float(all_tolerated)/len(all_mutations))*100)
    output.write("From All mutations: {:d} are deleterious ({:d}/{:d}={:d}%) and {:d} are tolerated ({:d}/{:d}={:d}%)\n".format(all_deleterious,all_deleterious,len(all_mutations),all_deleterious_percentage,all_tolerated,all_tolerated,len(all_mutations), all_tolerated_percentage))
     
    del df #release memory
    return all_mutations

# returns yes if slim inside in a disordered region, else it returns no.
def is_in_dire(slim,dire):
    status="No"
    if ((dire.start<=slim.start) and (slim.end<=dire.end)):
        status="Yes"
    """
    if((slim.start<dire.start) and (slim.end<=dire.start)):
        status="Yes"
    if((slim.start<=dire.end)and(dire.end<slim.end)):
        status="Yes"
    """
    return status

def count_duplicates(array):
    ## ~ ~ ~ counts the same mutations in different samples and changes the count attribute of the mutation (object from mutation class) to show the number of times each mutation appears in the all_mutations array ~ ~ ~ ##
    for i in range(0,len(array)): 
        current=array[i]
        for x in range(i+1, len(array)):
            check=array[x]
            if ((current.before==check.before) and (current.residue==check.residue)):
                current.count=current.count+1
                check.count=check.count+1
    return
            
def count_impact(array):
    deleterious=0
    tolerated=0
    for i in array:
        if (i.impact=="deleterious"):
            deleterious=deleterious+1
        if(i.impact=="tolerated"):
            tolerated=tolerated+1
    return deleterious, tolerated       


def find_motif_for_protein(kinase, protein_seq):
    ## ~ ~ ~ Finds slims on the given protein based on the motifs of a single protein from the known ones~ ~ ~ ##
    slim_array=[]
    motifs_array=[]
    
    #print(kinase.name,"\n")    
    
    motifs_buffer=kinase.motifs
    total_slims=0
    for i in motifs_buffer:
        
        #print (i.reg_expr)                        
        slim_buffer=re.findall(i.reg_expr, protein_seq)    
        #print (slim_buffer)
        
        slim_buffer=set(slim_buffer) #convert it to a set to keep only the unique entries
        slim_buffer=list(slim_buffer) # convert it to a list for easier manipulation
        total_slims=total_slims+len(slim_buffer) #to count all slims found for a protein
        for x in slim_buffer:                
            for m in (re.finditer(x,protein_seq)):
                st=m.start()+1  #start of the slim. +1 because list start from 0 but residues from 1
                en=m.end()+0    #end if the slim.
                
            new_slim=create_slim(st,en,x) #creates an new slim object
            slim_array.append(new_slim) #appends the new slim object in a buffer slim array that has all the slims for a given motif.
            
        new_motif=create_result_motif(i.reg_expr,slim_array,i.interaction_type) #creates new result motif object          
        motifs_array.append(new_motif)   #appends the new result motif object into a buffer motifs array that has all the results for the motifs of a given kinase 
        slim_array=[] # clearing the buffer array so it can have only the slims of the next motif and not also the slims of the previous motifs
        
    new_protein=create_result_protein(kinase.name,motifs_array,total_slims,kinase.type) #creates an new result_kinase object     
    motifs_array=[] #clear the array
    total_slims=0
    return new_protein

def find_motifs_for_all_proteins (kinases_array, protein_seq):
    ## ~ ~ ~ Finds slims in a given protein based on all the motifs of all of the known proteins ~ ~ ~ ##
    results_kinases_array=[]
    
    for i in kinases_array:
        result_protein=find_motif_for_protein(i, protein_seq)
        results_kinases_array.append(result_protein)  #appends the new result_protein object to an array that holds the results for all the kinases examined. 
    return results_kinases_array

# return an array that has only the slims of the motifs of kinases that are inside disordered regions 
def search_in_dire(results_array,direfile):
    
    disordered_regions=get_disordered_regions(direfile)
    #print (disordered_regions)
    
    results_dire_array=[]
    slims_dire_array=[]
    motifs_dire_array=[]
    total_dire_slims=0
    
    for i in results_array:
        for a in i.Amotifs:
            for b in a.Aslims:
                for c in disordered_regions:
                    if(is_in_dire(b,c)=="Yes"):
                        slims_dire_array.append(b)
            
            slims_dire_array=set(slims_dire_array) #convert it to a set to keep only the unique entries
            slims_dire_array=list(slims_dire_array) # convert it to a list for easier manipulation    
            
            total_dire_slims=total_dire_slims+len(slims_dire_array) #to count all slims found for a protein
            #check if slims_dire_array is empty. if empty there is no point of creating new object
            if(len(slims_dire_array)>0):
                new_motif=create_result_motif(a.regex,slims_dire_array,a.interaction_type) #creates new result motif object          
                motifs_dire_array.append(new_motif)   #appends the new result motif object into a buffer motifs array that has all the results for the motifs of a given kinase 
                slims_dire_array=[] #clears array
        #check if motifs_dire_array is empty. if empty there is no point of creating new object
        if(len(motifs_dire_array)>0):
            new_protein=create_result_protein(i.name,motifs_dire_array,total_dire_slims,i.type) #creates an new result_kinase object     
            results_dire_array.append(new_protein)
            motifs_dire_array=[] #clears array
            total_dire_slims=0
    return results_dire_array

def find_mutations_in_disorded_regions(output,outputname,all_mutations,disordered_regions):
    #finds the mutations that are in the disordered regions and appends them in an array
    mutations_in_dire=[]
    for i in all_mutations:
        for x in disordered_regions: 
            if (int(x.start)<=int(i.residue)<=int(x.end)):
                mutations_in_dire.append(i)
            
                
    
    print("Output file ", outputname, " opened for writting.")       
    protein_name=outputname.split("-",1)
    protein_name=str(protein_name[0])
    output.write("\n{:s} {:s} (Total={:d})\n".format(protein_name,"Mutations In Disordered Regions:",len(mutations_in_dire)))
        
    for i in mutations_in_dire:
        output.write("{:s} {:s} {:s} {:s} {:d}\n".format(i.before,i.residue,i.after,i.impact,i.count))    
    
    in_dire_deleterious, in_dire_tolerated=count_impact(mutations_in_dire)
    in_dire_deleterious_percentage=int((float(in_dire_deleterious)/len(mutations_in_dire))*100)
    in_dire_tolerated_percentage=int((float(in_dire_tolerated)/len(mutations_in_dire))*100)
    
    output.write("From mutations IN Disordered Regions: {:d} are deleterious ({:d}/{:d}={:d}%) and {:d} are tolerated ({:d}/{:d}={:d}%)\n".format(in_dire_deleterious,in_dire_deleterious,len(mutations_in_dire),in_dire_deleterious_percentage,in_dire_tolerated,in_dire_tolerated,len(mutations_in_dire),in_dire_tolerated_percentage))
    
    #finds the mutations that are in ordered regions and appends them in another array
    mutations_in_ordered_regions=[]
    for i in all_mutations:
        if (i in mutations_in_dire):
            pass
        else:
            mutations_in_ordered_regions.append(i)
    
    in_orderedre_deleterious, in_orderedre_tolerated=count_impact(mutations_in_ordered_regions)
    in_orderedre_deleterious_percentage=int((float(in_orderedre_deleterious)/len(mutations_in_ordered_regions))*100)
    in_orderedre_tolerated_percentage=int((float(in_orderedre_tolerated)/len(mutations_in_ordered_regions))*100)
    
    output.write("From mutations Outside of Disordered Regions (In Ordered regions): {:d} are deleterious ({:d}/{:d}={:d}%) and {:d} are tolerated ({:d}/{:d}={:d}%)\n".format(in_orderedre_deleterious,in_orderedre_deleterious,len(mutations_in_ordered_regions),in_orderedre_deleterious_percentage,in_orderedre_tolerated,in_orderedre_tolerated,len(mutations_in_ordered_regions),in_orderedre_tolerated_percentage))
    
    return mutations_in_dire

def find_mutations_in_a_slim(slim, mutations_in_dire):
    #for a given slim search all the mutations in the mutations_in_dire array and returns an array with the mutations inside the slim
    mutations_in_slim_array=[]
    for i in mutations_in_dire:
        if int(slim.start)<=int(i.residue)<=int(slim.end): #mutation in a slim
            mutations_in_slim_array.append(i)  
    return mutations_in_slim_array

def print_motifs_for_kinase(outputname,output,results_array):
    ## ~ ~ ~ prints the results on a file ~ ~ ~ ##
    underscores="_____________________________________________"
    
    count=0
    for a in results_array:
        output.write("{:s}\n{:s} [{:s}] :  (Total slims found= {:d})\n\n".format(underscores,a.name,a.type,a.total_slims))
        for b in a.Amotifs:
            count=count+1
            output.write(" {:d}) {:s}: ({:s})\n".format(count,b.regex,b.interaction_type))
            for c in b.Aslims:
                output.write("       {:d} - {:d} {:s}\n".format(c.start,c.end,c.sequence))
            output.write("\n")
        count=0 # set count to zero again to start counting the motifs of the next kinase in the list
    output.write(underscores)
    
    protein_name=outputname.split("-",1)
    protein_name=str(protein_name[0])
    print(protein_name, "Results written in ",outputname, "file.\n")  
    return

def print_In_DIRE_motifs_for_kinase(outputname,output,results_dire_array,mutations_dire_array):
    ## ~ ~ ~ prints the results on a file ~ ~ ~ ##
    underscores="_____________________________________________"
    
    count=0
    for a in results_dire_array:
        output.write("{:s}\n{:s} [{:s}] :  (Total slims found= {:d})\n\n".format(underscores,a.name,a.type,a.total_slims))
        for b in a.Amotifs:
            count=count+1
            output.write(" {:d}) {:s}: ({:s})\n".format(count,b.regex,b.interaction_type))
            for c in b.Aslims:
                output.write("       {:d} - {:d} {:s}\n".format(c.start,c.end,c.sequence))
                slim_mutations=find_mutations_in_a_slim(c,mutations_dire_array)
                output.write("           [Mutations: (Total = {:d})".format(len(slim_mutations)))
                string_buffer=""
                for i in slim_mutations:
                    string_buffer=string_buffer+"{:s} {:s} {:s} {:s} {:d}, ".format(i.before,i.residue,i.after,i.impact,i.count)   
                output.write(" {:s}]\n\n".format(string_buffer))    
            output.write("\n")
        count=0 # set count to zero again to start counting the motifs of the next kinase in the list
    output.write(underscores)
    
    protein_name=outputname.split("-",1)
    protein_name=str(protein_name[0])
    print(protein_name, "Results written in ",outputname, "file.\n")  
    return
    return

def search_protein(sequencefile,direfile):
    ## ~ ~ ~ The last function that does it all. Finds slims on a given protein for all known protein motifs and write them on a file ~ ~ ~ ##
    
    protein_seq=get_sequence(sequencefile)
    disordered_regions=get_disordered_regions(direfile)
    
    results_array=find_motifs_for_all_proteins(All_proteins,protein_seq)
    results_dire_array=search_in_dire(results_array,direfile)
    #print(results_dire_array)
    
    protein_name=sequencefile.split("_",1)
    protein_name=str(protein_name[0])
    
    interactors_and_slims_str="Possible Interactors and Slims"
    
    outputname=protein_name+"-SLiMs.txt"
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")
            
            dashes="-------------------------------------------------------------"
            output.write("{:s}\n{:s}\n".format(protein_name,dashes))
            
            
            output.write("{:s} Sequence:\n".format(protein_name))
            output.write("{:s} \n".format(protein_seq))
            
            print(protein_name, "Sequence written in ", outputname, "file.\n")                            
            output.write("{:s}\n\n".format(dashes))
                       
            output.write("{:s} {:s}:\n".format(protein_name,interactors_and_slims_str))
            
            print_motifs_for_kinase(outputname,output,results_array)
               
        else:
            raise IOError("Output file ",outputname," not created")
    output.close()
    

    
    outputname=protein_name+"-Mutations-in-disordered-regions.txt"
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")

            output.write("{:s}\n{:s}\n".format(protein_name,dashes))
        
            output.write("{:s} Sequence:\n".format(protein_name))
            output.write("{:s} \n".format(protein_seq))
        
            print(protein_name, "Sequence written in ", outputname, "file.\n")                            
            output.write("{:s}\n\n".format(dashes))
        
            output.write("\n{:s} {:s} (Total={:d})\n".format(protein_name,"Disordered Regions:",len(disordered_regions)))
            for i in disordered_regions:
                output.write("{:d} - {:d}\n".format(i.start,i.end))            
            print(protein_name, "Disordered regions written in ", outputname, "file.\n")
            output.write("{:s}\n\n".format(dashes)) 
        
            protein_mutations=get_mutations(output, outputname, "ZNF217_mutations.xlsx")
            #print(mutations)
        
            protein_mutations_in_dire=find_mutations_in_disorded_regions(output, outputname,protein_mutations,disordered_regions)

            output.write("{:s}\n\n".format(dashes))
        else:
            raise IOError("Output file ",outputname," not created")
        output.close()  
     
       
    #to find interactions in the possible new domain or subdomain of znf-217
    for i in results_array:
        for a in i.Amotifs:
            for b in a.Aslims:
                if ((285<=int(b.start)) and (int(b.end)<=310)):
                    print("{:d} - {:d} {:s} protein: {:s}".format(b.start, b.end, b.sequence, i.name))
    for i in protein_mutations_in_dire:
        if(285<=int(i.residue)<=310):
           print("{:s} {:s} {:s} {:s} {:d}\n".format(i.before,i.residue,i.after,i.impact,i.count))
      
         
    
    outputname=protein_name+"-DIRE_SLiMs_Mutations_2.txt"
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")
            
            output.write("{:s}\n{:s}\n".format(protein_name,dashes))
            
            
            output.write("{:s} Sequence:\n".format(protein_name))
            output.write("{:s} \n".format(protein_seq))
            
            print(protein_name, "Sequence written in ", outputname, "file.\n")                            
            output.write("{:s}\n\n".format(dashes))
            
            output.write("\n{:s} {:s} (Total={:d})\n".format(protein_name,"Disordered Regions:",len(disordered_regions)))
            for i in disordered_regions:
                output.write("{:d} - {:d}\n".format(i.start,i.end))            
            print(protein_name, "Disordered regions written in ", outputname, "file.\n")
            output.write("{:s}\n\n".format(dashes)) 
            
            output.write("{:s} in Disordered Regions of {:s}:\n".format(interactors_and_slims_str,protein_name))
            
            print_In_DIRE_motifs_for_kinase(outputname,output,results_dire_array,protein_mutations_in_dire)
               
        else:
            raise IOError("Output file ",outputname," not created")
    output.close()
             

     
    return 


# ~ ~ functions ~ ~  



# ~ ~ type of interaction for the motifs of each protein ~ ~ 
mod="modification"
doc="docking"
lig="ligation"
deg="degradation"
clv="cleavage"
# ~ ~ 

All_proteins=[]
# create a protein element and directly appent it to an array

#  ~ ~ all proteins and their motifs ~ ~ Data
All_proteins.append(create_protein("MAPK", "Kinase", "[KR]{0,2}[KR].{0,2}[KR].{2,4}[ILVM].[ILVF]",doc, "...[ST]P..",mod))
All_proteins.append(create_protein("MAPKAPK1", "Kinase", "[RK].R..S",mod))
All_proteins.append(create_protein("MAPKAPK2", "Kinase", "S...[ST]",mod))
All_proteins.append(create_protein("ERK1", "Kinase", "[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]",doc,"[RK].{2,4}[LIVMP].[LIV].[LIVMF]",doc,"SP",mod,"..SP",mod,".[ST]P",mod))
All_proteins.append(create_protein("ERK2", "Kinase", "[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]",doc,"[RK].{2,4}[LIVMP].[LIV].[LIVMF]",doc,"SP",mod,"..SP",mod,".[ST]P",mod))
All_proteins.append(create_protein("GSK-3", "Kinase", "..SP",mod,".[ST]P",mod,"...[ST]...[ST]",mod))
All_proteins.append(create_protein("P38", "Kinase", "[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]",doc,"[RK].{2,4}[LIVMP].[LIV].[LIVMF]",doc))
All_proteins.append(create_protein("JNK", "Kinase", "[RK]P[^P][^P]L.[LIVMF]",doc))
All_proteins.append(create_protein("CDK", "Kinase", "...[ST]P[RK]",mod,"...[ST]P..[RK]",mod,"SP.[RK].",mod,"[ST]P.[RK]",mod))
All_proteins.append(create_protein("CKI", "Kinase", "S..[ST]...",mod,"[ST]..[ST]",mod,"[ED]..[ST]",mod,"[ST]...[ST][M/L/V/I/F]",mod,"SP..[ST]",mod))
All_proteins.append(create_protein("CKII", "Kinase","S..[EST]",mod,"[ST]..[EDSY]",mod,"S..[ED]",mod,"[ST]..[ED]",mod,"S.[EST]",mod))
All_proteins.append(create_protein("NEK-2", "Kinase", "[FLM][^P][^P][ST][^DEP][^DE]",mod))
All_proteins.append(create_protein("PIKK", "Kinase","...[ST]Q..",mod))
All_proteins.append(create_protein("ATM", "Kinase","SQ",mod,"LSQE",mod))
All_proteins.append(create_protein("DNA-PK", "Kinase", "P[ST].",mod,".SQ",mod))
All_proteins.append(create_protein("PKA", "Kinase","[RK][RK].[ST][^P]..",mod,".R.[ST][^P]..",mod,"[RK].[ST]",mod,"[ST].[RK]",mod,"K...[ST]",mod,"K..[ST]",mod,"R..S",mod,"[RK][RK].[ST]",mod,"R.S",mod,"KR..S",mod))
All_proteins.append(create_protein("PKC Epsilon", "Kinase","R[KER].S",mod))
All_proteins.append(create_protein("PKC ", "Kinase","[RK].[ST]",mod,"[ST].[RK]",mod,"[RK]..[ST]",mod,"[RK]..[ST].[RK]",mod,"[RK].[ST].[RK]",mod))
All_proteins.append(create_protein("PLK-1", "Kinase",".[DNE][^PG][ST][FYILMVW]...",mod, ".[DNE][^PG][ST][^PEDGKN][FWYLIVM].",mod,"S[ST].",lig))
All_proteins.append(create_protein("PLK-2", "Kinase","[DE]..[ST][EDILMVFWY][DE].",mod, "[DE]..[ST][EDILMVFWY].[DE]",mod))
All_proteins.append(create_protein("PLK-3", "Kinase","[DE]..[ST][EDILMVFWY][DE].",mod, "[DE]..[ST][EDILMVFWY].[DE]",mod))
All_proteins.append(create_protein("PLK-4", "Kinase","..[^IRFW][ST][ILMVFWY][ILMVFWY]",mod))
All_proteins.append(create_protein("CaMKII", "Kinase","R..S",mod,"R..[ST]",mod,"[MVLIF].[RK]..[ST]..",mod))
All_proteins.append(create_protein("CaMKIV", "Kinase","[MILVFY].R..[ST]",mod))
All_proteins.append(create_protein("Chk1", "Kinase","[MILV].[RK]..[ST]",mod))
All_proteins.append(create_protein("CLK-1", "Kinase","R..[ST]..R",mod))
All_proteins.append(create_protein("CLKB-1", "Kinase","LRT",mod))
All_proteins.append(create_protein("mTOR", "Kinase","FTY",mod))
All_proteins.append(create_protein("NIMA", "Kinase","FR.[ST]",mod))
All_proteins.append(create_protein("EGFRK", "Kinase",".[ED]Y.",mod,".[ED]Y[ILV]",mod))
All_proteins.append(create_protein("JAK2", "Kinase","Y..[LIV]",mod))
All_proteins.append(create_protein("SRC", "Kinase","[IVLS].Y..[LI]",mod,"Y[AGSTED]",mod))
All_proteins.append(create_protein("Itk", "Kinase","Y[AEV][YFESNV][PFIH]",mod))
All_proteins.append(create_protein("BRCA1", "Ligase", ".S..F",doc,".S..F.K	",doc))
All_proteins.append(create_protein("APC/C", "Ligase",".KEN.",deg))
All_proteins.append(create_protein("SPOP", "Ligase","[AVP].[ST][ST][ST]",deg))
All_proteins.append(create_protein("PP2B", "Phosphatase", "L.[LIVAPM]P",doc))
All_proteins.append(create_protein("PP2C-delta", "Phosphatase", ".T.Y.",mod))
All_proteins.append(create_protein("SHP-1", "Phosphatase","[ED].Y",mod, "[IV].Y..[LV]",lig))
All_proteins.append(create_protein("TC-PTP", "Phosphatase","[EDY]Y",mod))
All_proteins.append(create_protein("CtBP1", "Dehydrogenase", "P[LVIPME][DENS][LM][VASTRG]", lig, "G[LVIPME][DENS][LM][VASTRG]K", lig, "G[LVIPME][DENS][LM][VASTRG].[KR]", lig))
All_proteins.append(create_protein("Caspase-3", "Protease","[DSTE][^P][^DEWHFYC]D[GSAN]",clv))
All_proteins.append(create_protein("Caspase-7", "Protease","[DSTE][^P][^DEWHFYC]D[GSAN]",clv))
All_proteins.append(create_protein("USP7", "Deubiquitinase", "[PA][^P][^FYWIL]S[^P]",doc,"	K...K",doc))
All_proteins.append(create_protein("Pin1", "Isomerase","...[ST]P.",doc))
All_proteins.append(create_protein("SUMO", "ETC","[DEST]{0,5}.[VILPTM][VIL][DESTVILMA][VIL].{0,1}[DEST]{1,10}",lig))
All_proteins.append(create_protein("SUMO-1", "ETC","[SDE].{0,5}[DE].K.{0,1}[AIFLMPSTV]",mod,"[DEST]{0,5}.[VILPTM][VIL][DESTVILMA][VIL].{0,1}[DEST]{1,10}",lig))
All_proteins.append(create_protein("FHA Domain", "ETC","..T..[ILV].",lig,"..T..[DE].",lig))
All_proteins.append(create_protein("MYND", "ETC","P.L.P",lig))
All_proteins.append(create_protein("TRF1-2 Domain", "ETC","[FY].L.P",lig))
All_proteins.append(create_protein("NAE1-UBA3", "ETC","[ILM][ILMF].{1,2}[ILM].{0,4}K",lig))
All_proteins.append(create_protein("Importin-A", "ETC","[PKR].{0,1}[^DE]K[RK][^DE][KR][^DE]", mod, "[PKR].{0,1}[^DE]K[RK][KR][^DE][^DE]", mod, "[PKR].{0,1}[^DE]RK[^DE][KR][^DE]",mod, "[PKR].{0,1}[^DE]RK[KR][^DE][^DE]", mod, "[PKR]K[RK][^DE][KR][^DE]", mod, "[PKR]K[RK][KR][^DE][^DE]", mod,  "[PKR]RK[^DE][KR][^DE]",mod, "[PKR]RK[KR][^DE][^DE]", mod))
All_proteins.append(create_protein("MDC-1", "ETC","S[ST].",lig))
All_proteins.append(create_protein("Cyclin A", "ETC",".[^EDWNSG][^D][RK][^D]L.{0,1}[FLMP].{0,3}[EDST]",doc, "[KRH].{0,3}[^EDWNSG][^D][RK][^D]L.{0,1}[FLMP].{0,3}[EDST]", doc))
All_proteins.append(create_protein("14-3-3", "ETC","R[^DE]{0,2}[^DEPG][ST][FWYLMV].",doc, "R[^DE]{0,2}[^DEPG][^PRIKGN]P", doc, "R[^DE]{0,2}[^DEPG][^PRIKGN].{2,4}[VILMFWYP]", doc, "R..S", doc))
All_proteins.append(create_protein("CDC20", "ETC","[ILVMF].[ILMVP][FHY].[DE]",doc))
All_proteins.append(create_protein("GRB-1", "ETC","Y.N",mod,"Y[MILV].[MILV]",mod))
# ~ ~ all proteins and their motifs ~ ~ 

# ~ ~ main ~ ~ 
search_protein("ZNF271_Sequence.txt","ZNF217_disordered_regions.txt")
# ~ ~ main ~ ~ 

