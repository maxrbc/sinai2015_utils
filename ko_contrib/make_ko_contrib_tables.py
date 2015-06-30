#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 09:02:43 2015

@author: xaradrim
"""
from ko_contrib import construct_biom_KO_contrib_tables as construct
import argparse as args

def take_args():
    parse = args.ArgumentParser()
    parse.add_argument("-i" ,help="KO contribution map " , required= True , dest='ko_table' , type = str )
    
    
    return parse.parse_args()




def main():
    opt = take_args()
    construct(opt.ko_table)
    pass






if __name__ == "__main__":
    main()