# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:59:32 2020

@author: u03132tk
"""

import pandas as pd
import primer3 as p3
import primers_functions as pf
def seqinmain(primer5to3,substring_length,genome,complement):
    seq=primer5to3[len(primer5to3)-substring_length:len(primer5to3)]
    test1=genome.count(seq)>1
    test4=complement.count(seq[::-1])>1
    out=test1 or test4
    return out

def hetero_homo_hairpin_test(primer_seq,primer_list, threshold_tm):
    for primerA in primer_list:
        if primerA==primer_seq:
            continue
        else:
            heterodimer_tm=p3.calcHeterodimer(primerA,primer_seq).tm
            if heterodimer_tm>threshold_tm:
                return ('heterodimer_tm>threshold_tm',
                        ('test primer', primer_seq), 
                        ('hit primer',primerA))
            
            homodimer_tm=p3.calcHomodimer(primer_seq).tm
            if homodimer_tm>threshold_tm:
                return ('heterodimer_tm>threshold_tm',
                        ('test primer', primer_seq), 
                        ('hit primer',primerA))
            hairpin_tm=p3.calcHairpin(primer_seq).tm
            if hairpin_tm>threshold_tm:
                return ('hairpin_tm>threshold_tm',
                        ('test primer', primer_seq), 
                        ('hit primer',primerA))
    return 'This is a non-boolean test'#not boolean as i want the tm if it fails

with open ("DEM32671genome.txt", "r") as myfile:
    genome = myfile.read()
print ('genome loaded')
complement_genome=pf.complement(genome)
print ('complement genome loaded')

url_cluster5_20bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588100119&job_key=6-E04D8TMrsVhTeAOuATskD7AoBt6BmdbA'
url_cluster5_25bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key=trxpvfcX-r_dhWqAZ-BOsh37X4Aw6ESdMQ'
url_cluster5_30bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key=vrRhtf8X8r_VhWKAb-BGshX7V4A46EydOQ+'
url_cluster8_20bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588100988&job_key=KiD1If8S8rrVhPeB-uHTs4D6woGt6dmcrA'
url_cluster8_25bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588100999&job_key=wcseyoEWjL6rhByBEeE4s2v6KYFG6TKcRw'
url_cluster8_30bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588102627&job_key=z8UQNoqdhzWgC50OkG65POp1qA7HZrMTxg'
url_cluster17_20bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588101128&job_key=b2WwliqdJzUACz0OMG4ZPEp1CA5nZhMTZg'
url_cluster17_25bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588101143&job_key=y8EUwIsWhr6hhBaBG-Eys2H6I4FM6TicTQ'
url_cluster17_30bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588101149&job_key=GhDFEc8SwrrlhMeByuHjs7D68oGd6emcnA'
url_cluster22_20bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588100843&job_key=09kM2JMWnr65hA6BA-Eqs3n6O4FU6SCcVQ'
url_cluster22_25bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588102891&job_key=i4FUgF4SU7p0hFaBW-FysyH6Y4EM6XicDQ'
url_cluster22_30bp='https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1588100868&job_key=LCbzJ_kS9LrThPGB_OHVs4b6xIGr6d-cqg'


cluster_urls=[(5,url_cluster5_20bp),
              (8,url_cluster8_20bp),
              (17,url_cluster17_20bp),
              (22,url_cluster22_20bp),
              (5,url_cluster5_25bp),
              (8,url_cluster8_25bp),
              (17,url_cluster17_25bp),
              (22,url_cluster22_25bp),
              (5,url_cluster5_30bp),
              (8,url_cluster8_30bp),
              (17,url_cluster17_30bp),
              (22,url_cluster22_30bp)]

primer_blast_params_list=[]
dataframe_out= pd.DataFrame()
primer_pair=0
print ('Empty dataframe built')

for cluster in cluster_urls:
    cluster_number=cluster[0]
    url=cluster[1]#output of primer blast
    tables = pd.read_html(url)
    #get running params
    primer_blast_params_list.append(('cluster number:  '+str(cluster_number), tables[0])) 
    for dataframe in tables[1::]:
        productlength=dataframe['Length'][3]
        newdata=dataframe.drop([2,3,4,5,6,7])
        newdata['primer_pair']=primer_pair#group forward and reverse primer pairs
        newdata['product length']=productlength#allow sorting by product length
        newdata['cluster']=cluster_number
        dataframe_out=dataframe_out.append(newdata, ignore_index=True)
        primer_pair+=1
dataframe_out.to_excel("rawdataframe_out.xlsx", sheet_name='Sheet_name_1')
print ('Full dataframe made')
sublenlist=[12,11,10]
for i in sublenlist:
    print ("filtering for unique 3' substrings of length:  "+str(i))
    dataframe_out_copy=dataframe_out.copy()#make a copy because you are running different filters on the same dataframe so dont alter in place.  could probs just put inplace=true instead tbf
    #do this on list of primers for all clusters:
    primer_seq_list=list(dataframe_out_copy["Sequence (5'->3')"])
    no_hetero_homo_hairpin_primers=[]
    three_prime_unique=[]
    substringlen=i
    thresholdtm=50
    for seq in primer_seq_list:
        test=hetero_homo_hairpin_test(seq,primer_seq_list,thresholdtm)#increase primer generation tm so this tm is ok 
        if test=='This is a non-boolean test':
            no_hetero_homo_hairpin_primers.append(seq)
        threeprime_test=seqinmain(seq,substringlen,genome,complement_genome)
        if threeprime_test==False:
            three_prime_unique.append(seq)
            
    print ('Primer sequences tested - remaining non-heterodimerising primers:  ', len(primer_seq_list) )   
     
    for primer_seq in no_hetero_homo_hairpin_primers:
        dataframe_out_copy.loc[dataframe_out_copy["Sequence (5'->3')"] == primer_seq, 'sec_test'] = 'OK'
    for primer_seq_prime in three_prime_unique:
        dataframe_out_copy.loc[dataframe_out_copy["Sequence (5'->3')"] == primer_seq_prime, 'three_prime'] = 'unique'
    test_dataframe=dataframe_out_copy
    dataframe_out_copy=dataframe_out_copy.dropna(subset=['sec_test'])#remove primers that form heterodimers
    dataframe_out_copy=dataframe_out_copy.dropna(subset=['three_prime'])#remove primers that have common 3 prime ends 
    test2=dataframe_out_copy
    primerPairs=list(dataframe_out_copy["primer_pair"])
    not_primerPairs=[]
    
    for pairID in primerPairs:
        if primerPairs.count(pairID)<2:
            if pairID not in not_primerPairs:
                not_primerPairs.append(pairID)
    
    dataframe_out_whole_pairs=dataframe_out_copy[~dataframe_out_copy['primer_pair'].isin(not_primerPairs)]#~ is reverse in general - 'is not in'
    
    with pd.ExcelWriter('primerBLASTanalyser_output.xlsx', engine="openpyxl", mode='a') as writer:  
        dataframe_out_whole_pairs.to_excel(writer, sheet_name='3 prime unique substring_'+substringlen+'_threshold tm_'+thresholdtm )
print ('Dataframes cleaned and exported')
    

