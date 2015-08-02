#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:26:31 2015

@author: xaradrim
"""


import seaborn as sns


GET_QUANT_BY_GROUP = lambda x : (x[0] , len(x[1]))
TO_RGB = lambda x : ( int(x[0]*255) , int(x[1]*255) , int(x[2]*255) )



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


def get_values_associated_with_id(table , id , axis = 'observation'):
        axis_class = 'sample'
        if axis_class == axis :
                axis_class = 'observation'
        print axis_class , axis , "\n"

        temp = {id : [] }
        for  i in table.ids(axis_class):
                if axis_class == 'sample' :
                        sample = i
                        obs = id
                else :
                        sample = id
                        obs = i
                if table.get_value_by_ids(obs , sample) > 0:
                        temp[id].append(i)
        return temp




def generate_band_location(biom_table  , taxa_level = 1 , path_level = 0):
    #extracting the groups and names
    otus = extract_otus_with_taxa(biom_table , taxa_level)
    kos = extract_kos_with_taxa(biom_table , path_level)

    # Acq chrm header band information for the groups
    otu_taxa_group = map(GET_QUANT_BY_GROUP , otus.iteritems())
    ko_path_group = map(GET_QUANT_BY_GROUP , kos.iteritems())

    #Structure for storing the band positioning

    '''
    result {
                otus : {
                           OTU_ID : ( id.band_parent , start , end )

                        }

                Keggs: {
                            KEGG_ID : [ ( id.band_parent , start , end ) , ( id.band_parent , start , end ) ]
                        }

    '''

    result = {'otus' : {} , 'keggs' : {} }
    for taxa , size in otu_taxa_group:
        otu_list = otus[taxa]
        index = 0
        count = 0
        while index < size :
            result['otus'][otu_list[index]] = (".".join([otu_list[index] , taxa]), count , count+1)
            count = count + 2
            index = index + 1

    for path , size in ko_path_group :
        ko_list = kos[path]
        index = 0
        count = 0
        while index < size :
            #make the tuple
            temp = (".".join([ko_list[index] , path.replace(" " , "_")]) , count , count + 1)
            kegg = ko_list[index]
            if not kegg in result['keggs'].keys() :
                result['keggs'][kegg] = [temp,]
            else:
                result['keggs'][kegg].append(temp)
            count = count +2
            index = index +1

    return result





def making_karyotype(biom_table ,output_dir=".", taxa_level = 1 , path_level = 0):

    ## chr - ID LABEL START END COLOR
    chrm_tmpl = "chr - %(ID)s %(LABEL)s %(START)s %(END)s %(COLOR)s\n"
    band_tmpl = "band %(PARENT)s %(ID)s %(LABEL)s %(START)s %(END)s %(COLOR)s\n"
    
    #extracting the groups and names 
    otus = extract_otus_with_taxa(biom_table , taxa_level)
    kos = extract_kos_with_taxa(biom_table , path_level)
    
    # Acq chrm header band information for the groups
    otus_chrm = map(GET_QUANT_BY_GROUP , otus.iteritems())
    kos_chrm = map(GET_QUANT_BY_GROUP , kos.iteritems())        
    
    #Makig otu karyotype
    tmpl = {}
    otu_karyo = open(output_dir+'/otu_karyotype.txt' , "w")
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
    ko_karyo = open(output_dir+'/ko_karyotype.txt', 'w')
    for chrm , size  in kos_chrm:
        chrm = chrm.replace(" " , "_")
        tmpl['ID'] = chrm 
        tmpl['LABEL'] = chrm
        tmpl['START'] = 0
        tmpl['END'] = (size*2)-1 
        tmpl['COLOR'] = chrm
        line  = chrm_tmpl % tmpl
        ko_karyo.write(line)
    
    #add bands now
    for chrm , size in kos_chrm:
        # carefull with order in here
        #chrm needs to be replace but reference in the dict
        #must be found first then change the reference name
        ko_list = kos[chrm]
        chrm = chrm.replace(" " , "_")

        index = 0
        count = 0
        while index < size:
            tmpl['PARENT'] = chrm
            tmpl['ID'] = ".".join([ko_list[index] , chrm])
            tmpl['LABEL'] = ko_list[index]

            tmpl['START'] = count
            tmpl['END'] = count+1
            tmpl['COLOR'] = chrm
            
            count = count + 2
            index = index + 1
            
            line  = band_tmpl % tmpl
            ko_karyo.write(line)
    ko_karyo.close()
    return




def make_color_reference(biom_table , taxa_level = 1 , path_level = 0):
    ko_quant = map(GET_QUANT_BY_GROUP , extract_kos_with_taxa(biom_table , path_level).iteritems())
    otu_quant = map(GET_QUANT_BY_GROUP , extract_otus_with_taxa(biom_table , taxa_level).iteritems())

    ko_name , ko_size = zip(*ko_quant)
    otu_name , otu_size = zip(*otu_quant)

    group = []
    group.extend(ko_name)
    group.extend(otu_name)


    color_reference = {}
    colors = sns.color_palette('Set2' , len(group))

    for member , color  in zip(group , colors):
        #Names cant have whitespaces between
        #they must be separated using something else
        #this case an underscore

        member = member.replace(" " , "_")
        color_reference[member] = TO_RGB(color)

    return color_reference


    

    pass






def make_circos_conf(biom_table , ko_path , out_path , output_dir):
    templ = '''

################################################################
# The remaining content is standard and required. It is imported
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files,
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>

<ideogram>


<spacing>
default = 0.005r
</spacing>

radius    = 0.9r
thickness = 20p
fill      = yes


# ##########################################################################
# label
#
# show_label       = yes
# label_font       = light
#
# # 50 pixels outside the outer ideogram radius
# label_radius = dims(ideogram,radius_outer) + 200p
#
# # 5% of inner radius outside outer ideogram radius
# # label_radius = dims(ideogram,radius_outer) + 0.05r
#
# # inside ideogram
# # label_radius = (dims(ideogram,radius_inner)+dims(ideogram,radius_outer))/2-24
#
# # 100 pixels inside the ideogram radius
# # label_radius = dims(ideogram,radius_inner) - 100p
#
# label_with_tag   = yes
# label_size       = 5
# label_parallel   = no
# label_case       = upper
#
# ##########################################################################


</ideogram>


##########################################################################
#LINKS

<links>
<link>
ribbon        = yes
radius        = 0.95r
bezier_radius = 0.1r
crest         = 0.85
file          = links.txt
color         = black
</link>
</links>



##########################################################################


'''
    conf_fp = "/".join([output_dir , "circos.conf"])

    tmp = ["\n<colors>"]

    for name , color in make_color_reference(biom_table).iteritems():
        n_color = "{0[0]},{0[1]},{0[2]}".format(color)
        t = "{} = {} ".format(name , n_color)
        tmp.append(t)
    tmp.append("</colors>\n")
    coloring = "\n".join(tmp)
    karyo = "karyotype = {},{}\n\n".format(ko_path , out_path)




    doc = open(conf_fp , "w")
    doc.write(karyo)
    doc.write(templ)
    doc.write(coloring)
    doc.close()

    pass

