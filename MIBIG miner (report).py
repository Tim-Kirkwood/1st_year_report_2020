# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 19:55:39 2020

@author: u03132tk
"""

import json
import os

path='C:/Users/u03132tk/.spyder-py3/mibig_json_2.0'

json_list=[]
types=[]
types_dict={}

weird=[]

polycluster=[]
nrpcluster=[]
othercluster=[]
saccharidecluster=[]
rippcluster=[]

polycluster_pop=0
nrpcluster_pop=0
othercluster_pop=0
saccharidecluster_pop=0
rippcluster_pop=0

directory = os.fsencode(path)
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".json"): 
         filename_total=path+'/'+filename
         data = json.load(open(filename_total,'r'))
         json_list.append(data)
         continue
     else:
         continue
hybrid_count=0
hybrid_2=0   
hybrid_3=0   
hybrid_4=0   
max_hybrid=0
types_in_mibig=0

test_sub=0
test_main=0
types_in_mibig_sub=0
types_in_mibig_main=0
types_in_mibig_sub1=0
types_in_mibig_sub2=0
types_in_mibig_sub11=0
wtf=[]
non_min_sac=[]
non_min_nrp=[]
non_min_other=[]
non_min_ripp=[]
non_min_poly1=[]
non_min_poly2=[]
non_min_poly3=[]
non_min_terp=[]
non_min_alk=[]

polytypeI=[]
polytypeII=[]
polytypeIII=[]
polytypeAT=[]
for json_object in json_list:
    dict_object=json_object['cluster']
    types_in_mibig+=len(set(dict_object['biosyn_class']))
    if len(dict_object['biosyn_class'])>1:
        hybrid_count+=1
    if len(dict_object['biosyn_class'])==2:
        hybrid_2+=1  
    if len(dict_object['biosyn_class'])==3:
        hybrid_3+=1   
    if len(dict_object['biosyn_class'])==4:
        hybrid_4+=1   
    if len(dict_object['biosyn_class'])>max_hybrid:
        max_hybrid=len(dict_object['biosyn_class'])
    clustertype=[]
    if 'minimal' not in dict_object.keys():
        print ('oh no')
        break
    #issues are coming here.  you have same number of types, but some are getting lost in the logic.  its def from here because number of clusters main class > num subclass
    if dict_object['minimal']==False:
        types_in_mibig_sub+=len(dict_object['biosyn_class'])
        testripp='ripp' in dict_object.keys()
        testpoly='polyketide' in dict_object.keys()
        testsac='saccharide' in dict_object.keys()
        testnrp='nrp' in dict_object.keys()
        testother='other' in dict_object.keys()
        testterpene='terpene' in dict_object.keys()
        testalk='alkaloid' in dict_object.keys()
        if testpoly or testsac or testnrp or testother or testripp or testterpene or testalk:
            types_in_mibig_sub1+=len(dict_object['biosyn_class'])
            test_sub+=1
            if testalk:
                if 'subclass' in dict_object['alkaloid'].keys():
                    types_in_mibig_sub11+=1
                    name_type=['Alkaloid']
                    add_list_alk=[dict_object['alkaloid']['subclass']]
                    types+=[':'.join(name_type+add_list_alk)]
                else:
                    types_in_mibig_sub11+=1
                    non_min_alk+=[dict_object]
                    types+=['Alkaloid']
            if testterpene:
                if 'carbon_count_subclass' in dict_object['terpene'].keys() or 'structural_subclass' in dict_object['terpene'].keys():
                    if 'carbon_count_subclass' in dict_object['terpene'].keys():
                        carbon_count_subclass = dict_object['terpene']['carbon_count_subclass']
                    if 'structural_subclass' in dict_object['terpene'].keys():
                        structural_subclass = dict_object['terpene']['structural_subclass']
                    types_in_mibig_sub11+=1
                    name_type=['Terpene']
                    types+=[':'.join(name_type+[carbon_count_subclass] + [structural_subclass])]
                    terp=[dict_object]
                else:
                    types_in_mibig_sub11+=1
                    non_min_terp+=[dict_object]
                    types+=['Terpene']
            if testpoly: #'''keep going through this layer by layer and see wher data is leaking.  leak = where sub1.1 < sub 1.  its in one of these tests'''
                if 'synthases' in dict_object['polyketide'].keys():
                    if 'subclass' in dict_object['polyketide']['synthases'][0].keys():
                        count_typeI=0
                        count_typeII=0
                        count_typeIII=0
                        name_type=['Polyketide']
                        add_list_polyketide=dict_object['polyketide']['synthases'][0]['subclass']
                        #types_in_mibig_sub11+=1
                        for a_type in add_list_polyketide:
                            testtypeI='type I' in a_type and 'type II' not in a_type or 'Type I' in a_type and 'Type II' not in a_type
                            testtypeII='type II' in a_type and 'type III' not in a_type or 'Type II' in a_type and 'Type III' not in a_type
                            testtypeIII='type III' in a_type or 'Type III' in a_type
                            if testtypeI:
                                count_typeI+=1
                            if testtypeII:
                                count_typeII+=1
                            if testtypeIII:
                                count_typeIII+=1
                        if count_typeI == len(add_list_polyketide):#dont artificially inflate typeX clusters by separateing from their subclass
                            types_in_mibig_sub11+=1
                            types+=[':'.join(name_type+add_list_polyketide)]
                            polytypeI+=[dict_object]
                        elif count_typeII == len(add_list_polyketide):
                            types_in_mibig_sub11+=1
                            types+=[':'.join(name_type+add_list_polyketide)]
                            polytypeII+=[dict_object]
                        elif count_typeIII == len(add_list_polyketide):
                            types_in_mibig_sub11+=1
                            types+=[':'.join(name_type+add_list_polyketide)]
                            polytypeIII+=[dict_object]
                        else:       
                            for a_type in add_list_polyketide:
                                print (a_type)
                                types_in_mibig_sub11+=1
                                types+=[':'.join(name_type+[a_type])]
                            print('   ')
                            non_min_poly2+=[dict_object]
                        polycluster+=[':'.join(name_type+add_list_polyketide)]
                    else:
                        types_in_mibig_sub11+=1
                        types+=['Polyketide']
                        non_min_poly1+=[dict_object]
                elif 'subclasses' in dict_object['polyketide'].keys():#not here
                    name_type=['Polyketide']
                    add_list_polyketide=dict_object['polyketide']['subclasses']
                    if 'Trans' in ':'.join(name_type+add_list_polyketide) and 'ype I' not in ':'.join(name_type+add_list_polyketide):
                        polytypeAT=[dict_object]
                    if dict_object['polyketide']['subclasses']==['Tetracycline']:
                        print (dict_object['mibig_accession'])
                    for a_type in add_list_polyketide:
                        types_in_mibig_sub11+=1
                        types+=[':'.join(name_type+[a_type])]
                    non_min_poly3+=[dict_object]
                    print (add_list_polyketide)
                    print ('done')
                else:
                    types_in_mibig_sub11+=1
                    types+=['Polyketide']
                    
            if testsac:
                if 'subclass' in dict_object['saccharide'].keys():
                    types_in_mibig_sub11+=1
                    name_type=['Saccharide']
                    add_list_saccharide=[dict_object['saccharide']['subclass']]
                    types+=[':'.join(name_type+add_list_saccharide)]
                    saccharidecluster+=add_list_saccharide
                    saccharidecluster_pop+=1
                else:
                    types_in_mibig_sub11+=1
                    non_min_sac+=[dict_object]
                    types+=['Saccharide']
            if testnrp:
                if 'subclass' in dict_object['nrp'].keys():
                    types_in_mibig_sub11+=1
                    name_type=['NRP']
                    add_list_nrp=[dict_object['nrp']['subclass']]
                    types+=[':'.join(name_type+add_list_nrp)]
                    nrpcluster+=add_list_nrp
                    nrpcluster_pop+=1
                else:
                    types_in_mibig_sub11+=1
                    non_min_nrp+=[dict_object]
                    types+=['NRP']
            if testother:
                if 'subclass' in dict_object['other'].keys():
                    types_in_mibig_sub11+=1
                    name_type=['Other']
                    add_list_other=[dict_object['other']['subclass']]
                    types+=[':'.join(name_type+add_list_other)]
                    othercluster+=add_list_other
                    othercluster_pop+=1
                    #print(clustertype)
                else:
                    types_in_mibig_sub11+=1
                    non_min_other+=[dict_object]
                    types+=['Other']
            if testripp:
                rippcluster_pop+=1
                rippcluster+=dict_object
                if 'subclass' in dict_object['ripp'].keys():
                    types_in_mibig_sub11+=1
                    name_type=['RiPP']
                    add_list_other=[dict_object['ripp']['subclass']]
                    types+=[':'.join(name_type+add_list_other)]
                    rippcluster+=add_list_other
                else:
                    types_in_mibig_sub11+=1
                    non_min_ripp+=[dict_object]
                    types+=['RiPP']
        else:
            types_in_mibig_sub2+=len(dict_object['biosyn_class'])
            weird.append(dict_object)
            types+=dict_object['biosyn_class']
    elif dict_object['minimal']==True:
        types_in_mibig_main+=len(dict_object['biosyn_class'])
        test_main+=1
        if type(dict_object['biosyn_class']) is list:
            types+=dict_object['biosyn_class']
        else:
            print ('oh no')
            break
    else:
        wtf.append(dict_object)

for i in types:
    types_dict[i]=0

for type_object in types:
    if type_object in types_dict.keys():
        types_dict[type_object]+=1
    else:
        print (type_object, types_dict.keys())
        print ('oh no')
        
hybrid_percent=100*(hybrid_count/len(json_list))
hybrid_percent2=100*(hybrid_2/len(json_list))
hybrid_percent3=100*(hybrid_3/len(json_list))
hybrid_percent4=100*(hybrid_4/len(json_list))
