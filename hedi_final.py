
import os
import sys
import pandas as pd
from pandas import Series, DataFrame, rolling_median
import numpy as np
import re
#import matplotlib.pyplot as plt
import math
from collections import OrderedDict
from datetime import date, datetime
#from ipykernel import kernelapp as app
#get_ipython().magic('matplotlib inline')
import csv

infile = sys.argv[1]

outfolder = "/Volumes/Pathology/Molecular\ Pathology/Test\ Development/Bioinformatics_Projects/Databases/Parsing_HEDI/testing_hedi/"
#print(outfolder.replace(" ", "\ "))
print("infile is : %s" %infile)
xl = pd.read_excel(sys.argv[1], sheet_name='Detail')

print(xl.columns)


# In[4]:

#xl[['CASE','CASE TYPE', 'TEST', 'SOURCE','TEST.1','COMPONENT NAME', 'RESULT','#.1', 'COMPONENT NAME.1', 'RESULT.1', '#.2', 'COMPONENT NAME.2', 'RESULT.2', '#.3', 'COMPONENT NAME.3',
#      'RESULT.3', '#.4', 'COMPONENT NAME.4', 'RESULT.4', '#.5','COMPONENT NAME.5', 'RESULT.5','COMPONENT NAME.6','RESULT.6']]


#import nltk
#nltk.download()


#s = xl[['TEST','RESULT.1','RESULT.2','RESULT.3']]
new = xl.loc[xl['TEST.1'] == "CANCER MUTATION PROFILING"]
new1 = new.reset_index(drop = True)


#new2
new11 = new1[['SOURCE','CASE','PAT MRN','TEST','RESULT','RESULT.1','RESULT.2','RESULT.3','RESULT.4','RESULT.5','RESULT.6']]
new2 = new11.loc[(new11['RESULT']== "Detected")|(new11['RESULT'] == "See Text")]


#s['BACKGROUND'] = s['RESULT.3'].str.split('BACKGROUND:', expand=True)[1]
#s
#get text after INTERPRETATION word
new2['Inter'] = new2['RESULT.5'].str.split('INTERPRETATION', expand=True)[1]
new2
# get text upto performed word
new2['INTERPRETATION_1'] = new2['Inter'].str.rsplit('REGIONS WITH INADEQUATE COVERAGE').str[0]
new2
new2['INTERPRETATION_2'] = new2['Inter'].str.rsplit('REGIONS OF INADEQUATE COVERAGE').str[0]
new2
new2['INTERPRETATION'] = new2['INTERPRETATION_1'] + new2['INTERPRETATION_2']
#get text after References word
new2['References'] = new2['RESULT.5'].str.split('REFERENCES', expand=True)[1]
new2
new2['References_test'] = new2['RESULT.5'].str.split('Variants and clinical interpretation', expand=True)[1]
new2
new2['chr'] = new2['RESULT.5'].str.split('REFSEQ ID', expand=True)[1]
new2
new2['chr1'] = new2['chr'].str.rsplit('variants with potential clinical significance').str[0]
new2
new2['chr3'] = new2['chr'].str.rsplit('chr$').str[0]
new2
#s['c_position'] = s[2].str.split('=', expand=True)[1]
#s['p_position'] = s[3].str.split('=', expand=True)[1]
#s.columns
#s


# In[9]:

new2 ['Variant1'] = new2['INTERPRETATION'].str.split('\s{4,}', expand=True)[0]
new2 ['Variant2'] = new2['INTERPRETATION'].str.split('\s{4,}', expand=True)[1]
new2 ['Variant3'] = new2['INTERPRETATION'].str.split('\s{4,}', expand=True)[2]
new2


# In[10]:

new3 = new2[['SOURCE','CASE','TEST','RESULT','RESULT.1','RESULT.2','RESULT.3','RESULT.4','RESULT.5','INTERPRETATION','Variant1','Variant2','Variant3','References_test']]
new3


# In[11]:

#new3 = new2.drop(new2.index[4])
#new3.reset_index(drop = True)

#import nltk
#nltk.download()
#from nltk.tokenize import sent_tokenize, word_tokenize
#from nltk.corpus import stopwords
#from nltk.tokenize import word_tokenize


#new3['tokenized_sents'] = new3.apply(lambda row: nltk.sent_tokenize(row['INTERPRETATION']), axis=1)
#new3

#new3.to_csv('interpret.csv')

#new3.tokenized_sents


# In[11]:

new3['gene.1'] = new3['RESULT.1'].str.rsplit(' ').str[0]
new3
new3['gene.2'] = new3['RESULT.2'].str.rsplit(' ').str[0]
new3
new3['gene.3'] = new3['RESULT.3'].str.rsplit(' ').str[0]
new3
new3['1_c'] = new3['RESULT.1'].str.split('c.', expand=True)[1]
new3['1_p'] = new3['RESULT.1'].str.split('p.', expand=True)[1]

new3['gene.1_c'] = new3['1_c'].str.rsplit('p.').str[0]
new3
new3['1.1_p'] = new3['1_p'].str.rsplit('p.').str[0]
new3
new3['gene.1_p'] = new3['1.1_p'].str.rsplit(' ').str[0]
new3
new3['2_c'] = new3['RESULT.2'].str.split('c.', expand=True)[1]
new3['2_p'] = new3['RESULT.2'].str.split('p.', expand=True)[1]
new3['gene.2_c'] = new3['2_c'].str.rsplit('p.').str[0]
new3
new3['gene.2_p'] = new3['2_p'].str.rsplit(' ').str[0]
new3
new3['3_c'] = new3['RESULT.3'].str.split('c.', expand=True)[1]
new3['3_p'] = new3['RESULT.3'].str.split('p.', expand=True)[1]
new3['gene.3_c'] = new3['3_c'].str.rsplit('p.').str[0]
new3
new3['gene.3_p'] = new3['3_p'].str.rsplit(' ').str[0]
new3
new3['Run'] = new3['RESULT.5'].str.split('======================================', expand=True)[1]
new3
new3['Run_BC'] = new3['Run'].str.split('________________________', expand=True)[0]
new3
#s
#get text after INTERPRETATION word
#new2['Inter'] = new2['RESULT.5'].str.split('INTERPRETATION', expand=True)[1]
#new2
# get text upto performed word
#new2['INTERPRETATION'] = new2['Inter'].str.rsplit('REGIONS WITH INADEQUATE COVERAGE').str[0]
#new2
#s['c_position'] = s[2].str.split('=', expand=True)[1]
#s['p_position'] = s[3].str.split('=', expand=True)[1]
#s.columns
#s


new4 = new3[['SOURCE','CASE','TEST','RESULT','gene.1','gene.1_c','gene.1_p','gene.2','gene.2_c','gene.2_p','gene.3','gene.3_c','gene.3_p','Variant1','Variant2','Variant3','References_test','Run_BC']]
#new4

result1 = new4[['SOURCE','CASE','TEST','RESULT','gene.1','gene.1_c','gene.1_p','Variant1','References_test','Run_BC']]
result1

result2 = new4[['SOURCE','CASE','RESULT','TEST','gene.2','gene.2_c','gene.2_p','Variant2','References_test','Run_BC']]
result2

result3 = new4[['SOURCE','CASE','RESULT','TEST','gene.3','gene.3_c','gene.3_p','Variant3','References_test','Run_BC']]
result3


result_1 = result1.rename(columns = {'gene.1':'gene','gene.1_c':'gene_c','gene.1_p':'gene_p', 'Variant1':'Variant_information'})
result_1

result_2 = result2.rename(columns = {'gene.2':'gene','gene.2_c':'gene_c','gene.2_p':'gene_p','Variant2':'Variant_information'})
result_2

result_3 = result3.rename(columns = {'gene.3':'gene','gene.3_c':'gene_c','gene.3_p':'gene_p','Variant3':'Variant_information'})
result_3


result = pd.concat([result_1,result_2,result_3], axis =0)
result


result_final = result.sort_values(by='CASE', ascending=1)


result_finals = result_final[['SOURCE','CASE','RESULT','TEST','gene','gene_c','gene_p','Variant_information','References_test','Run_BC']]
result_finals


#getting c. p.codes from variant information column
#result_final['c_code'] = result_final['Variant_information'].str.extract('(c.[0-9]\w{0,})', expand=True)
#result_final['p_code'] = result_final['Variant_information'].str.extract('(p.[A-Z][0-9]\w{0,})', expand=True)
##getting proper Run and barcodefor sample from Run_BC columns
#result_final['Run_details'] = result_final['Run_BC'].str.extract('(R[0-9]|MP#: _\)', expand=True)

#result_final['Run_details'] = result_final['Run_BC'].str.extract('(R[0-9]+)', expand=False).str.rstrip()
result_finals.to_csv(infile+'_VariantsList_1.csv',index=False)

