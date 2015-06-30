from __future__ import division
from os import getcwd  , path 
from biom.table import Table
import numpy as np

def make_partial_table(table_fp):

    if not path.exists(table_fp):
        raise IOError("File doesnt exist ! "+getcwd())
    
    current_ko = None
    otu_list = []
    
    #Fixing the full and midial path management
    table_p = "/".join(table_fp.split("/")[:-1])
    table_name = table_fp.split("/")[-1]
    output_name = table_p+"/ko_summary."+table_name
    
    int_result = open(output_name , "w")
    
    for line in open(table_fp , "r"):
        #parse the column
        ko,sample,otu = line.strip("\n").split("\t")[:3]
        
        if 'Gene' in ko :
            continue
                   
        #case #1 first entry on the file 
        elif current_ko == None :
            current_ko = ko
            otu_list.append(otu)
            continue
        #case #2 change on the ko 
        elif not current_ko == ko:
                int_result.write("%s|%s\n" % ( ko , "\t".join(otu_list) ))
                current_ko = ko
                otu_list = [otu]
                continue
        #case #3 new otu realtion for the ko 
        else :
            otu_list.append(otu)
    
    int_result.close()
    return output_name

def construct_biom_KO_contrib_tables(table_fp):
    KO = []
    Otus = set()
    
    #Fixing for file path 
    file_fp = make_partial_table(table_fp)
    output_p = "/".join(file_fp.split("/")[:-1])
    table_name = file_fp.split("/")[-1]
    
    #something odd
    # determine the Otus and have the KO in list 
    for line in open(file_fp , "r"):
        ko , otu_list = line.strip("\n").split("|")
        otu_list = otu_list.split("\t")
        KO.append(ko)
        Otus = Otus.union(set(otu_list))
    
    Otus = list(Otus)
    data = np.zeros((len(Otus) , len(KO)))
    
    #Populate the adj matrix with the data parsed from the file 
    for line in open(file_fp , "r"):
       ko , otu_list = line.strip("\n").split("|")
       otu_list = otu_list.split("\t")
       k_index = KO.index(ko)
       
       for o in otu_list:
           o_index = Otus.index(o)
           data[o_index , k_index ] = 1
    
    
    #write down the tables  
    ##sample = ko
    ko_as_sample = Table(data,Otus,KO)
    doc = open(output_p+"/ko_as_sample."+table_name , "w")
    doc.write(ko_as_sample.to_json("KO_CONTRIB_SCRIPT"))
    doc.close()
    
    
    ##sample = otu
    Otu_as_sample = Table(data.T , KO , Otus)
    doc = open(output_p+"/otu_as_sample."+table_name, "w")
    doc.write(Otu_as_sample.to_json("KO_CONTRIB_SCRIPT"))
    doc.close()
         
            

