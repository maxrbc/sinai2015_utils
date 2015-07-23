#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:26:31 2015

@author: xaradrim
"""



def extract_otus_with_taxa(biom_table):
    result = {}
    
    for otu in biom_table.ids(axis='observation'):
        phylum = biom_table.metadata(otu , axis='observation')['taxonomy'][1]
        
        if phylum in result.keys():
            result[phylum].append(otu)
        else:
            result[phylum] = [otu,]
    return result
    