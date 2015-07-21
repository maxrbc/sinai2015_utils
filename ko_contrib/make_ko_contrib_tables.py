#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 09:02:43 2015

@author: xaradrim
"""
from ko_contrib import make_ko_contrib_table as construct
import argparse as args
import os

def take_args():
    parse = args.ArgumentParser()
    parse.add_argument("-i" , help="Biom table with the otus from the analysis" , required= True , dest='otu_table' , type = str )
    parse.add_argument("-p" , help='prediction table resulting from the picrust analysis' , required=True , dest='predicted_table' , type=str )
    parse.add_argument("-c" , help='KO contributions map from the analysis done with picrust' , required=True , dest='ko_contrib' , type=str)
    parse.add_argument("-f" , help="Force write the output dir" , dest="force_write" , action="store_true" ,  type = bool)
    parse.add_argument("-o" , help="Output directory" , required=True , dest='output_dir' , type=str )
    
    return parse.parse_args()




def main():
    opt = take_args()
    
    if not os.path.exists(opt.otu_table):
        raise IOError("Otu table doesnt exist!")
    
    if not os.path.exists(opt.predicted_table):
        raise IOError("Prediction table doesnt exits!")
    
    if not os.path.exists(opt.ko_contrib):
        raise IOError("KO contributions file doesnt exist")
    
    if os.path.exists(opt.output_dir) and opt.force_write or not os.path.exists(opt.output_dir):
        os.mkdir(opt.output_dir)
    
    
    
                        
    construct(opt.ko_contrib , opt.otu_table , opt.predicted_table , opt.output_dir)
    pass






if __name__ == "__main__":
    main()