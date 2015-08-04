__author__ = 'xaradrim'

from seaborn import color_palette

GET_QUANT_BY_GROUP = lambda x: (x[0], len(x[1]))
TO_RGB = lambda x: (int(x[0] * 255), int(x[1] * 255), int(x[2] * 255))


def extract_otus_with_taxa(biom_table, taxa_level=1):
    result = {}

    for otu in biom_table.ids(axis='observation'):
        phylum = biom_table.metadata(otu, axis='observation')['taxonomy'][taxa_level]

        if phylum in result.keys():
            result[phylum].append(otu)

        else:
            result[phylum] = [otu, ]

    return result


def extract_kos_with_taxa(biom_table, path_level=0):
    result = {}

    for ko in biom_table.ids(axis='sample'):
        pathway_list = biom_table.metadata(ko)['KEGG_Pathways']

        for pathway in pathway_list:

            # Need to check this with jose , but for now
            # So it works
            if 'None' in pathway:
                if pathway[0] in result.keys():
                    result[pathway[0]].append(ko)
                else:
                    result[pathway[0]] = [ko, ]

            elif pathway[path_level] in result.keys():
                result[pathway[path_level]].append(ko)

            else:
                result[pathway[path_level]] = [ko, ]

    return result


def make_color_reference(group):

    color_reference = {}
    colors = color_palette('Set2', len(group))

    for member, color in zip(group, colors):
        # Names cant have whitespaces between
        # they must be separated using something else
        # this case an underscore

        member = member.replace(" ", "_")
        color_reference[member] = TO_RGB(color)

    return color_reference


def get_values_associated_with_id(table, value_id):
    '''
    Get all values that are associated with the given id

    :rtype : dict
    :param table: biom table object
    :param value_id: id to be searched { ex. '15478' if OTU or 'KO5487' if KEGG}
    :return: { id : [link1 , link2 , link 3]
    '''
    result = {}
    temp = []

    if value_id in table.ids():
        # we got a sample ( KEGGS )
        for obs in table.ids(axis='observation'):
            if table.get_value_by_ids(obs, value_id) > 0:
                temp.append(obs)

    elif value_id in table.ids(axis='observation'):
        # we got a observation ( OTU )
        for sample in table.ids():
            if table.get_value_by_ids(value_id, sample) > 0:
                temp.append(sample)
    else:
        raise AttributeError('The entered ID : ' + value_id + ', Doesnt exist in biom table ')

    result[str(value_id)] = temp
    return result
