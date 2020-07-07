# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:17:23 2020

@author: u03132tk
"""


import json
#load the dem genome json file into dictionary with looooots of nests
sample='sequence.json'
data = json.load(open(sample,'r'))
       
#for i in range(len(record_link)):
base=data['records'][0]
features=base['features']
genome_dic=base['seq']
genome = genome_dic['data']
######################################################
##find number of predicted and known pbtA/X/AX clusters 

def list_checker(sourcelist,checklist):
    shared_items=[]
    for item in sourcelist:
        if item in checklist:
            shared_items.append(item)
    return shared_items
    

region_number=[4]


clusterblast_proteins_allregions=base['modules']['antismash.modules.clusterblast']['general']['results']
pbtA_clusters=[]
pbtA_clusters_more=[]

id_pbta='ctg1_782'
id_pbtx='ctg1_779'
id_pbtr='ctg1_775'
id_pbtg1='ctg1_776'
id_pbtb1='ctg1_777'
id_pbto='ctg1_778'
id_pbtm1='ctg1_780'
id_pbtm2='ctg1_781'
id_pbtb='ctg1_783'
id_pbtc='ctg1_784'
id_pbtd='ctg1_785'
id_pbte='ctg1_786'
id_pbtf='ctg1_787'
id_pbtg='ctg1_788'
id_pbth='ctg1_789'
id_pbtm3='ctg1_790'
id_pbtm4='ctg1_791'

pbta=[]
pbtx=[]
pbtr=[]
pbtg1=[]
pbtb1=[]
pbto=[]
pbtm1=[]
pbtm2=[]
pbtb=[]
pbtc=[]
pbtd=[]
pbte=[]
pbtf=[]
pbtg=[]
pbth=[]
pbtm3=[]
pbtm4=[]


for region in region_number:
    #make variable of the ranked clusterblast clusters for the region you are on 
    ranked_clusterblast_hits= clusterblast_proteins_allregions[region-1]['ranking']
    #go to top ranked cluster homolog from clusterblast (ie number[0]).  
    for number in range(len(ranked_clusterblast_hits)):
        topclusterdata=ranked_clusterblast_hits[number]#abstract this to loop until you have a hoit for all proteins
        protein_blast_summary=topclusterdata[1]['pairings'] #[0] is just info on homolog cluster, in [1] the other non-pairing keys have scoring values eg synteny
        #update orf dict value
        for iii in protein_blast_summary:
            #build tuple variables for value.  iii[0] is hit title, [1] is int - maybe number of matches?, [2] is protein-specific blast data
            name_start=iii[0].index('ctg')
            name_end=iii[0].rindex('|')#last instance
            dem_cluster=iii[0][name_start:name_end]
            test_annot=iii[2]['annotation']
            evalue=iii[2]['evalue']
            source_protein=iii[2]['name']
            source_cluster=iii[2]['genecluster']
            if dem_cluster == id_pbta:
                pbta.append(source_cluster)#, source_protein,test_annot))
                pbtA_clusters_more.append((source_cluster, evalue))
            if dem_cluster==id_pbtx:
                pbtx.append(source_cluster)
            if dem_cluster==id_pbtr:
                pbtr.append(source_cluster)
            if dem_cluster==id_pbtg1:
                pbtg1.append(source_cluster)
            if dem_cluster==id_pbtb1:
                pbtb1.append(source_cluster)
            if dem_cluster==id_pbto:
                pbto.append(source_cluster)
            if dem_cluster==id_pbtm1:
                pbtm1.append(source_cluster)
            if dem_cluster==id_pbtm2:
                pbtm2.append(source_cluster)
            if dem_cluster==id_pbtb:
                pbtb.append(source_cluster)
            if dem_cluster==id_pbtc:
                pbtc.append(source_cluster)
            if dem_cluster==id_pbtd:
                pbtd.append(source_cluster)
            if dem_cluster==id_pbte:
                pbte.append(source_cluster)
            if dem_cluster==id_pbtf:
                pbtf.append(source_cluster)
            if dem_cluster==id_pbtg:
                pbtg.append(source_cluster)
            if dem_cluster==id_pbth:
                pbth.append(source_cluster)
            if dem_cluster==id_pbtm3:
                pbtm3.append(source_cluster)
            if dem_cluster==id_pbtm4:
                pbtm4.append(source_cluster)

pbtax=list_checker(pbta,pbtx)
pbtar=list_checker(pbta,pbtr)
pbtag1=list_checker(pbta,pbtg1)
pbtab1=list_checker(pbta,pbtb1)
pbtao=list_checker(pbta,pbto)
pbtam1=list_checker(pbta,pbtm1)
pbtam2=list_checker(pbta,pbtm2)
pbtab=list_checker(pbta,pbtb)
pbtac=list_checker(pbta,pbtc)
pbtad=list_checker(pbta,pbtd)
pbtae=list_checker(pbta,pbte)
pbtaf=list_checker(pbta,pbtf)
pbtag=list_checker(pbta,pbtg)
pbtah=list_checker(pbta,pbth)
pbtam3=list_checker(pbta,pbtm3)
pbtam4=list_checker(pbta,pbtm4)
                

                         
                
                

            
