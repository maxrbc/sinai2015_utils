#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:26:31 2015

@author: xaradrim
"""
GET_QUANT_BY_GROUP = lambda x : (x[0] , len(x[1]))


def extract_otus_with_taxa(biom_table ,taxa_level = 1):
    result = {}
    
    for otu in biom_table.ids(axis='observation'):
        phylum = biom_table.metadata(otu , axis='observation')['taxonomy'][taxa_level]
        
        if phylum in result.keys():
            result[phylum].append(otu)
            
        else:
            result[phylum] = [otu,]
    
    return result

def extract_kos_with_taxa(biom_table , path_level = 0):
    result = {}
    
    for ko in biom_table.ids(axis='sample'):
        pathway_list = biom_table.metadata(ko)['KEGG_Pathways']
        
        for pathway in pathway_list: 
            
            if 'None' in pathway:
                continue
            
            elif pathway[path_level] in result.keys():
                result[pathway[path_level]].append(ko)
                
            else:
                result[pathway[path_level]] = [ko,]
                
    return result
                

## chr - ID LABEL START END COLOR
chrm_tmpl = "chr - %(ID)s %(LABEL)s %(START)s %(END)s %(COLOR)s\n"
band_tmpl = "band %(PARENT)s %(ID)s %(LABEL)s %(START)s %(END)s %(COLOR)s\n"

def making_karyotype(biom_table , taxa_level = 1 , path_level = 0):
    
    #extracting the groups and names 
    otus = extract_otus_with_taxa(biom_table , taxa_level)
    kos = extract_kos_with_taxa(biom_table , path_level)
    
    # Acq chrm header band information for the groups
    otus_chrm = map(GET_QUANT_BY_GROUP , otus.iteritems())
    kos_chrm = map(GET_QUANT_BY_GROUP , kos.iteritems())        
    
    #Makig otu karyotype
    tmpl = {}
    otu_karyo = open('otu_karyotype.txt' , "w")
    for chrm , size  in otus_chrm:
        tmpl['ID'] = chrm 
        tmpl['LABEL'] = chrm
        tmpl['START'] = 0
        tmpl['END'] = (size*2)-1 
        tmpl['COLOR'] = chrm
        line  = chrm_tmpl % tmpl
        otu_karyo.write(line)
    
    #add the bands 
    for chrm , size in otus_chrm:
        otu_list = otus[chrm]
        index = 0
        count = 0
        while index < size:
            tmpl['PARENT'] = chrm
            tmpl['ID'] =  otu_list[index]
            tmpl['LABEL'] = otu_list[index]
            tmpl['START'] = count
            tmpl['END'] = count+1
            tmpl['COLOR'] = chrm
            
            count = count + 2
            index = index + 1
            
            line  = band_tmpl % tmpl
            otu_karyo.write(line)
    otu_karyo.close()
    
    #Making ko karyotype 
    ko_karyo = open('ko_karyotype.txt', 'w')
    for chrm , size  in kos_chrm:
        tmpl['ID'] = chrm 
        tmpl['LABEL'] = chrm
        tmpl['START'] = 0
        tmpl['END'] = (size*2)-1 
        tmpl['COLOR'] = chrm
        line  = chrm_tmpl % tmpl
        ko_karyo.write(line)
    
    #add bands now
    for chrm , size in kos_chrm:
        ko_list = kos[chrm]
        index = 0
        count = 0
        while index < size:
            tmpl['PARENT'] = chrm
            tmpl['ID'] =  ko_list[index]
            tmpl['LABEL'] = ".".join([ko_list[index] , chrm])
            tmpl['START'] = count
            tmpl['END'] = count+1
            tmpl['COLOR'] = chrm
            
            count = count + 2
            index = index + 1
            
            line  = band_tmpl % tmpl
            ko_karyo.write(line)
    ko_karyo.close()
    return

    
    
    
    