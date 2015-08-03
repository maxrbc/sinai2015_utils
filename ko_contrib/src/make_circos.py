#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:26:31 2015

@author: xaradrim
"""

import os
from utils import extract_kos_with_taxa , extract_otus_with_taxa , GET_QUANT_BY_GROUP , \
    make_color_reference , get_values_associated_with_id




def generate_band_location(biom_table, taxa_level=1, path_level=0):

    # extracting the groups and names
    otus = extract_otus_with_taxa(biom_table, taxa_level)
    kos = extract_kos_with_taxa(biom_table, path_level)

    # Acq chrm header band information for the groups
    otu_taxa_group = map(GET_QUANT_BY_GROUP, otus.iteritems())
    ko_path_group = map(GET_QUANT_BY_GROUP, kos.iteritems())

    # Structure for storing the band positioning

    '''
    result {
                otu_chrm : [ ( id.band , start , end ) , ( id.band , start , end ) ]

                kegg_chrm : [ ( id.band , start , end ) , ( id.band , start , end ) ]

                otus : {
                           OTU_ID : ( id.band_parent , start , end ) ,

                        }

                keggs: {
                            KEGG_ID : [ ( id.band_parent , start , end ) , ( id.band_parent , start , end ) ] ,

                        }

    '''

    result = { 'otu_chrm' : [] , 'otus': {}, 'kegg_chrm' : [] , 'keggs': {} }


    #otu bands

    for taxa, size in otu_taxa_group:
        result['otu_chrm'].append((taxa , 0 , (size*2) - 1))
        otu_list = otus[taxa]
        index = 0
        count = 0
        while index < size:
            result['otus'][otu_list[index]] = (".".join([otu_list[index], taxa]), count, count + 1)
            count = count + 2
            index = index + 1

    for path, size in ko_path_group:
        result['kegg_chrm'].append((path , 0 , (size*2) -1 ))
        ko_list = kos[path]
        index = 0
        count = 0
        while index < size:
            # make the tuple
            temp = (".".join([ko_list[index], path.replace(" ", "_")]), count, count + 1)
            kegg = ko_list[index]
            if not kegg in result['keggs'].keys():
                result['keggs'][kegg] = [temp, ]
            else:
                result['keggs'][kegg].append(temp)
            count = count + 2
            index = index + 1

    return result


def make_links(generated_bands ,table, output_dir='.'):

    link_doc = open(os.path.join(output_dir , 'links.txt') , 'w')

    for otu in generated_bands['otus'].keys():

        # getting related keggs to the specify otu
        kegg_related = get_values_associated_with_id(table , otu)
        otu_band    = generated_bands['otus'][otu][0].split('.')[1]
        otu_start   = str(generated_bands['otus'][otu][1])
        otu_end     = str(generated_bands['otus'][otu][2])
        otu_link    = " ".join([otu_band , str(otu_start) , str(otu_end) ])

        for kegg in kegg_related[otu]:
            for kegg_info in generated_bands['keggs'][kegg]:
                kegg_band   = kegg_info[0].split('.')[1]
                kegg_start  = str(kegg_info[1])
                kegg_end    = str(kegg_info[2])
                kegg_link   = " ".join([kegg_band , kegg_start , kegg_end])

                link  = " ".join([otu_link , kegg_link , "\n"])
                link_doc.write(link)

    link_doc.close()

    pass


def making_karyotype(generated_bands, output_dir="."):

    ## chr - ID LABEL START END COLOR
    chrm_tmpl = "chr - %(ID)s %(LABEL)s %(START)s %(END)s %(COLOR)s\n"
    band_tmpl = "band %(PARENT)s %(ID)s %(LABEL)s %(START)s %(END)s %(COLOR)s\n"
    tmpl = {}

   #############################################################################
   ##### Making the Otu karyotype

    otu_karyo = open(output_dir + '/otu_karyotype.txt', "w")
    for taxa, start , end  in generated_bands['otu_chrm']:
        tmpl['ID'] = taxa
        tmpl['LABEL'] = taxa
        tmpl['START'] = start
        tmpl['END'] = end
        tmpl['COLOR'] = taxa
        line = chrm_tmpl % tmpl
        otu_karyo.write(line)

    # add the bands
    for id_parent , start, end  in generated_bands['otus'].values():

        tmpl['PARENT'] = id_parent.split(".")[1]
        tmpl['ID'] = id_parent.split(".")[0]
        tmpl['LABEL'] = id_parent.split(".")[0]
        tmpl['START'] = start
        tmpl['END'] = end
        tmpl['COLOR'] = id_parent.split(".")[1]

        line = band_tmpl % tmpl
        otu_karyo.write(line)

    otu_karyo.close()

    #############################################################################

    #############################################################################
    #### Making KEGG karyotype

    ko_karyo = open(output_dir + '/ko_karyotype.txt', 'w')
    for path, start , end  in generated_bands['kegg_chrm']:
        tmpl['ID'] = path
        tmpl['LABEL'] = path
        tmpl['START'] = start
        tmpl['END'] = end
        tmpl['COLOR'] = path
        line = chrm_tmpl % tmpl
        ko_karyo.write(line)

    # add bands now
    for t in generated_bands['keggs'].values():

        for id_parent ,start ,end in t:

            tmpl['PARENT'] = id_parent.split(".")[1]
            tmpl['ID'] = id_parent
            tmpl['LABEL'] = id_parent.split(".")[0]
            tmpl['START'] = start
            tmpl['END'] = end
            tmpl['COLOR'] = id_parent.split(".")[1]
            line = band_tmpl % tmpl
            ko_karyo.write(line)

    ko_karyo.close()

    #############################################################################
    pass



def make_circos_conf(biom_table, ko_path, out_path, output_dir):
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
    conf_fp = "/".join([output_dir, "circos.conf"])

    tmp = ["\n<colors>"]

    for name, color in make_color_reference(biom_table).iteritems():
        n_color = "{0[0]},{0[1]},{0[2]}".format(color)
        t = "{} = {} ".format(name, n_color)
        tmp.append(t)
    tmp.append("</colors>\n")
    coloring = "\n".join(tmp)
    karyo = "karyotype = {},{}\n\n".format(ko_path, out_path)

    doc = open(conf_fp, "w")
    doc.write(karyo)
    doc.write(templ)
    doc.write(coloring)
    doc.close()

    pass
