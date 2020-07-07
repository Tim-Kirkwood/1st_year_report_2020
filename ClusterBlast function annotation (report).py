# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:17:23 2020

@author: u03132tk
"""

import json
import sys

def cleaner(string):
    '''this takes in a string of numbers and chars, and returns a list of numbers split by :'''
    lis_out_final=[]
    lis_in=[]
    for i in string:
        if i==':' or i in str(list(range(10))).replace('[','').replace(']',''):
            lis_in.append(i) 
    lis_out_object=''.join(lis_in).split(':')
    for i in lis_out_object:
        lis_out_final.append(int(i))
    return lis_out_final


def region_bool(nested_region_index_list,query_start,query_end):
    '''return region if indexes are in region, or false if they are not.  region list format: '[[start region1, end region 1],[start region 2, end region 2]...]'''
    for i in nested_region_index_list:
        for ii in i:
            if ii<0:
                sys.exit('start index too low')
    inregion=False
    for i in nested_region_index_list:
        if i[0]<=start<=i[1] and i[0]<=end<=i[1]:
            return nested_region_index_list.index(i)+1
    return inregion
        

'''this will find all genes and their smcog classification from antismash file.  use to give coordinates of dna seqs to send to multigene blast'''
#load the dem genome json file into dictionary with looooots of nests
test='NC_003888.json'
sample='DEMgenome.json'
data = json.load(open(sample,'r'))


            
#go through dict to relevant nested index
record_link=data['records']
#for i in range(len(record_link)):
vicky_3=record_link[2]
features=vicky_3['features']
genome_dic=vicky_3['seq']
genome = genome_dic['data']

 
#get your regions to filter non bgc other genes
lis=[]
regions=[]
cand_cluster=[]
lis_whole_genome=[]
for i in features:
    if i['type']=='region':
        locus_cleaned_list=cleaner(i['location'])
        regions.append(locus_cleaned_list) 
    if i['type']=='cand_cluster':
        locus_cleaned_list=cleaner(i['location'])
        cand_cluster.append(locus_cleaned_list) 
               
for i in features:
    qualifiers = i['qualifiers']    
    if i['type']=='CDS':#**doesnt seem to have join when i tested it, or the other characters**some also have '>' and '<'?  Cut these out too**join is in some of the locations.  maybe important, but cut out and see if it affects final output of gene grouos** i added this in so only cds are added to loci.  otherwise the locus cleaned list may include stuff for protoclusters etc.
        if len(qualifiers['locus_tag'])>1:
            sys.exit('Error 2')
        locus_cleaned_list=cleaner(i['location'])
        start=locus_cleaned_list[0]
        end=locus_cleaned_list[1]
        #if the gene has a kind assigned by smcog (additional, core or transport) then it is appended to lis
        region=region_bool(regions,start,end)
        if region: 
            if 'gene_kind' in qualifiers.keys():
                lis.append((region,qualifiers['locus_tag'],qualifiers['gene_kind'],qualifiers['translation'],start,end))#nb smcog annotataes some other genes but not all of them.    MAYBE THE ONES IT ANNOTATES ARE MORE IMPORTANT?  CHECK WITH GENMTAMICIN CLUSTER
                lis_whole_genome.append((region,qualifiers['locus_tag'],qualifiers['gene_kind'],qualifiers['translation'],start,end))
            else:
                test=[]
                for item in lis:
                    test.append(item[1])
                if qualifiers['locus_tag'] not in test:
                    lis.append((region,qualifiers['locus_tag'],'OTHER',qualifiers['translation'],start,end))
                    lis_whole_genome.append((region,qualifiers['locus_tag'],'OTHER',qualifiers['translation'],start,end))
        else:
            lis_whole_genome.append((qualifiers['locus_tag'],'Non-cluster',qualifiers['translation'],start,end))

#####################################################
##find cluster blast annotations for protein function
clusterblast_proteins_allregions=data['records'][2]['modules']['antismash.modules.clusterblast']['general']['results']

region_number=list(range(1,31))
region_results=[]
for i in region_number:
    #make variable of the ranked clusterblast clusters for the region you are on 
    ranked_clusterblast_hits= clusterblast_proteins_allregions[i-1]['ranking']
    #build dict of orfs in region under analysis, update values as you loop through known clusters
    clusterorfs={}
    for ii in lis:
        if ii[0]==i:
            name=str(ii[1]).strip('[]\'')
            clusterorfs[name]=None
    #go to top ranked cluster homolog from clusterblast (ie number[0]).  
    for number in range(len(ranked_clusterblast_hits)):
        topclusterdata=ranked_clusterblast_hits[number]#abstract this to loop until you have a hoit for all proteins
        protein_blast_summary=topclusterdata[1]['pairings'] #[0] is just info on homolog cluster, in [1] the other non-pairing keys have scoring values eg synteny
        protein_homolog=[]
        #update orf dict value
        for iii in protein_blast_summary:
            #build tuple variables for value.  iii[0] is hit title, [1] is int - maybe number of matches?, [2] is protein-specific blast data
            name_start=iii[0].index('ctg')
            name_end=iii[0].rindex('|')#last instance
            dem_cluster=iii[0][name_start:name_end]
            test_annot=iii[2]['annotation']
            evalue=iii[2]['evalue']
            source=iii[2]['name']
            #if key in orf dict is none (ie hasnt been picked up by a higher ranked cluster) update value
            if dem_cluster in clusterorfs.keys() and clusterorfs[dem_cluster] is None: #if value=none?
                clusterorfs[dem_cluster]=(source,test_annot, 'Cluster rank:'+str(number),'Total clusters compared: '+str(len(ranked_clusterblast_hits)),evalue)
            if None not in clusterorfs.values(): #if there are orf keys without a proper value cycle through the rest of the cluster hits until they are all filled or there are no more homolog clusters
                break
    region_results.append((i,clusterorfs))