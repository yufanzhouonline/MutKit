########################################################################
##### A python library for reading and analysis of mutation data  ######
#####                by Yufan (Harry) Zhou at 2024/09/13          ######
########################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Mutlib(object):
    '''
    Library to 
    1) read VCF files
    2) extract fastq file from genome
    3) get mutational context
    '''
    def __init__(self):
        self.chrlist = ['chr' + str(i) for i in range(1,20)] + ['chrX', 'chrY']
        
    def read_vcf(self, filename='', rename_columns = True):
        '''
        read vcf file
        startline: line number of vcf data header line
        lines: list of vcf files in lines
        return:
        self.header: list of VCF header, each line is the element of the list
        self.vcf: data frame of mutations
        '''
        with open(filename, 'rt') as f:
            lines = [line[:-1] for line in f]
        for i in range(len(lines)):
            #print(i)
            if lines[i][:2] == '#C':
                startline = i
                break
        colnames = lines[startline][1:].split('\t')
        vcf = [i.split('\t') for i in lines[startline+1:]]
        vcfdf = pd.DataFrame(vcf)
        if rename_columns:
            vcfdf.columns = colnames
        self.header = lines[:startline+1]
        self.vcf = vcfdf
        
    def pretreat(self, vcfdf, posstr = True):
        '''
        Extract chrosome in chrlist
        Sort in the order of chromosome and position
        '''
        ### sort vcf
        vcfdf = vcfdf[vcfdf.CHROM.isin(self.chrlist)].copy()
        vcfdf['POS'] = vcfdf['POS'].astype(int)
        vcfdf.sort_values(by=['CHROM', 'POS'], ascending=True, inplace=True)
        vcfdf.reset_index(drop = True, inplace = True)
        if posstr:
            vcfdf['POS'] = vcfdf['POS'].astype(str)
        return vcfdf
    
    def sbs2dbs_with_loop(self, vcfdf):
        '''
        transfer continuous SBS to DBS with loop
        '''
        vcfdf['dbstag'] = 0
        for i in vcfdf.index[1:]:
            print(i + 1, '/', vcfdf.shape[0])
            if vcfdf.loc[i, 'CHROM'] == vcfdf.loc[i-1, 'CHROM'] and vcfdf.loc[i, 'POS'] == vcfdf.loc[i-1, 'POS'] + 1:
                vcfdf.loc[i - 1, 'dbstag'] = 1
                vcfdf.loc[i, 'POS'] = vcfdf.loc[i - 1, 'POS']
                vcfdf.loc[i, 'REF'] = vcfdf.loc[i - 1, 'REF'] + vcfdf.loc[i, 'REF']
                vcfdf.loc[i, 'ALT'] = vcfdf.loc[i - 1, 'ALT'] + vcfdf.loc[i, 'ALT']
        ### update vcf
        vcfdf = vcfdf[vcfdf.dbstag==0].copy()
        vcfdf.reset_index(drop = True, inplace = True)
        vcfdf['POS'] = vcfdf['POS'].astype(str)
        vcfdf.drop('dbstag', inplace=True, axis=1)
        return vcfdf
        
    def sbs2dbs(self,vcfdf):
        '''
        transfer continuous SBS to DBS with apply
        '''
        ### label removed first sbs
        ### label removed first sbs
        vcfdf['postchrom'] =  list(vcfdf.CHROM[1:]) + ['chr0']
        vcfdf['postpos'] =  list(vcfdf.POS[1:]) + [0]
        vcfdf['postref'] =  list(vcfdf.REF[1:]) + ['X']
        vcfdf['postalt'] = list(vcfdf.ALT[1:]) + ['X']
        vcfdf['dbstag'] = vcfdf.apply(lambda row : 1 if row['CHROM'] == row['postchrom'] and row['POS'] == row['postpos'] - 1 else 0, axis=1)
        ### updated value of second sbs to dbs
        vcfdf['prechrom'] = ['chr0'] + list(vcfdf.CHROM[:-1])
        vcfdf['prepos'] = [0] + list(vcfdf.POS[:-1])
        vcfdf['preref'] = ['X'] + list(vcfdf.REF[:-1])
        vcfdf['prealt'] = ['X'] + list(vcfdf.ALT[:-1])
        vcfdf['newpos'] =  vcfdf.apply(lambda row : row['prepos'] if row['CHROM'] == row['prechrom'] and row['POS'] == row['prepos'] + 1 else row['POS'], axis=1)
        vcfdf['newref'] =  vcfdf.apply(lambda row : row['preref'] + row['REF'] if row['CHROM'] == row['prechrom'] and row['POS'] == row['prepos'] + 1 else row['REF'], axis=1)
        vcfdf['newalt'] =  vcfdf.apply(lambda row : row['prealt'] + row['ALT'] if row['CHROM'] == row['prechrom'] and row['POS'] == row['prepos'] + 1 else row['ALT'], axis=1)
        vcfdf['POS'] = vcfdf['newpos']
        vcfdf['REF'] = vcfdf['newref']
        vcfdf['ALT'] = vcfdf['newalt']
        ### update vcf
        vcfdf = vcfdf[vcfdf.dbstag==0].copy()
        vcfdf.reset_index(drop = True, inplace = True)
        vcfdf['POS'] = vcfdf['POS'].astype(str)
        vcfdf.drop(['dbstag', 'postchrom', 'postpos', 'postref', 'postalt', 'prechrom', 'prepos', 'preref', 'prealt', 'newpos', 'newref', 'newalt'], inplace=True, axis=1)
        return vcfdf
    
    def muse(self, vcfdf):
        '''
        extract AD, DP and AF in MuSE
        '''
        vcfdf['ad'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[2].split(',')[1]), axis=1)
        vcfdf['dp'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[1]), axis=1)
        vcfdf['af'] = [round(i, 4) for i in vcfdf['ad'] / vcfdf['dp']]
        
        return vcfdf.iloc[:,[0, 1, 2, 3, 4, 5, 6, 11, 12, 13]].copy()
    
    def mutect2(self, vcfdf):
        '''
        extract AD, DP and AF in Mutect2
        '''
        vcfdf['ad'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[1].split(',')[1]), axis=1)
        vcfdf['dp'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[3]), axis=1)
        vcfdf['af'] = vcfdf.apply(lambda row : float(row[vcfdf.columns[10]].split(':')[2]), axis=1)
        
        return vcfdf.iloc[:,[0, 1, 2, 3, 4, 5, 6, 11, 12, 13]]
        
    def strelka(self, vcfdf):
        '''
        extract AD, DP and AF in Strelka
        '''
        ### Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
        
        vcfdf['ad'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[0]), axis=1)
        vcfdf['dp'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[0]), axis=1)
        vcfdf['af'] = 1.0
        
        ### for SNPs
        vcfdf['au1'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[4].split(',')[0]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['cu1'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[5].split(',')[0]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['gu1'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[6].split(',')[0]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['tu1'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[7].split(',')[0]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['au2'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[4].split(',')[1]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['cu2'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[5].split(',')[1]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['gu2'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[6].split(',')[1]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['tu2'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[7].split(',')[1]) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        
        vcfdf['altcounts'] = vcfdf.apply(lambda row : int(row['au2'] if row['ALT'] == 'A' else row['cu2'] if row['ALT'] == 'C' else row['gu2'] if row['ALT'] == 'G' else row['tu2'] if row['ALT'] == 'T' else 0) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        vcfdf['refcounts'] = vcfdf.apply(lambda row : int(row['au1'] if row['REF'] == 'A' else row['cu1'] if row['REF'] == 'C' else row['gu1'] if row['REF'] == 'G' else row['tu1'] if row['REF'] == 'T' else 0) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else 0, axis=1)
        
        ### for INDELs
        vcfdf['indelaltcounts'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[10]].split(':')[3].split(',')[0]) if row[vcfdf.columns[8]].split(':')[1] == 'DP2' else 0, axis=1)
        vcfdf['indelrefcounts'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[9]].split(':')[2].split(',')[0]) if row[vcfdf.columns[8]].split(':')[1] == 'DP2' else 0, axis=1)
        
        ### for SNPs and INDELs
        #vcfdf['af'] = vcfdf['altcounts'] / (vcfdf['altcounts'] + vcfdf['refcounts'])
        #vcfdf['af'] = vcfdf.apply(lambda row : row['altcounts'] / (row['altcounts'] + row['refcounts']) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else row['indelaltcounts'] / (row['indelaltcounts'] + row['indelrefcounts']) if (row['indelaltcounts'] + row['indelrefcounts']) > 0 else 0, axis=1)
        vcfdf['af'] = vcfdf.apply(lambda row : row['altcounts'] / (row['altcounts'] + row['refcounts']) if row[vcfdf.columns[8]].split(':')[1] == 'FDP' else row['indelaltcounts'] / (row['indelaltcounts'] + row['indelrefcounts']), axis=1)
        vcfdf['ad'] = vcfdf.apply(lambda row : int(np.ceil(row['dp'] * row['af'])), axis=1)
        
        return vcfdf.iloc[:,[0, 1, 2, 3, 4, 5, 6, 11, 12, 13]]
    
    def varscan(self, vcfdf):
        '''
        extract AD, DP and AF in Varscan
        '''
        vcfdf['ad'] = 0
        vcfdf['dp'] = vcfdf.apply(lambda row : int(row[vcfdf.columns[7]].split(';')[1].split('=')[1]), axis=1)
        vcfdf['af'] = vcfdf.apply(lambda row : float(row[vcfdf.columns[7]].split(';')[0].split('=')[1]), axis=1)
        vcfdf['ad'] = [int(i) for i in vcfdf['dp'] * vcfdf['af']]
        
        return vcfdf.iloc[:,[0, 1, 2, 3, 4, 5, 6, 8, 9, 10]]
        
    def twoplus(self, allvcf):
        '''
        Integrate 4 callers to 2+ mutations
        '''
        ### sort and label duplicates
        allvcf.sort_values(by=['CHROM', 'POS'], ascending=True, inplace=True)
        allvcf.reset_index(drop = True, inplace = True)
        
        ### prepare for count duplicates
        allvcf['mut'] = allvcf.apply(lambda row : row['CHROM'] + '_' + row['POS'] + '_' + row['REF'] + '_' + row['ALT'], axis=1)
        allvcf['premut'] = ['firstmut'] + list(allvcf.mut[:-1])
        
        ### count duplicates
        allvcf['equal'] = 1
        allvcf['preequal'] = 1
        allvcf['equal'] = allvcf.apply(lambda row : 2 if row['mut'] == row['premut'] and row['preequal'] == 1 else row['equal'], axis=1)
        allvcf['preequal'] = [1] + list(allvcf.equal[:-1])
        allvcf['equal'] = allvcf.apply(lambda row : 3 if row['mut'] == row['premut'] and row['preequal'] == 2 else row['equal'], axis=1)
        allvcf['preequal'] = [1] + list(allvcf.equal[:-1])
        allvcf['equal'] = allvcf.apply(lambda row : 4 if row['mut'] == row['premut'] and row['preequal'] == 3 else row['equal'], axis=1)
        
        ### extract unique and save
        savevcf = allvcf[allvcf.equal==2].copy()
        savevcf.reset_index(drop = True, inplace = True)
        return savevcf
    
    def remove_dup(self, savevcf):
        '''
        Remove duplicated mutations
        '''
        savevcf['prechr'] = ['chr0'] + list(savevcf.CHROM[:-1])
        savevcf['prepos'] = ['0'] + list(savevcf.POS[:-1])
        savevcf['duptag'] = savevcf.apply(lambda row : 1 if row['CHROM'] == row['prechr'] and row['POS'] == row['prepos'] else 0, axis=1)
        savevcf = savevcf[savevcf.duptag == 0].copy()
        savevcf.reset_index(drop = True, inplace = True)
        return savevcf
        
    def savevcf(self, savevcf, header, savepath, filename, postfix = '.call2plus.vcf'):
        '''
        Save vcf to savepath filename
        '''
        savevcflist = []
        for i in savevcf.index:
            savevcflist.append(savevcf.loc[i, 'CHROM'] + '\t' + str(savevcf.loc[i, 'POS']) + '\t' + savevcf.loc[i, 'ID'] + '\t' + savevcf.loc[i, 'REF'] + '\t' + savevcf.loc[i, 'ALT'] + '\t' + savevcf.loc[i, 'QUAL'] + '\t' + savevcf.loc[i, 'FILTER'] + '\t' + str(savevcf.loc[i, 'ad']) + '\t' + str(savevcf.loc[i, 'dp']) + '\t' + str(round(savevcf.loc[i, 'af'], 4)) + '\n')
        
        header = [i + '\n' for i in header]
        header = header + ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAD\tDP\tAF\n']
        savevcflist = header + savevcflist
        
        savefilename = savepath + filename.split('.')[0] + postfix
        with open(savefilename, 'wt') as f:
            f.writelines(savevcflist)
    
    def fourcallers(self, sampleno, filenamelist, filepath, savepath):
        '''
        Run four callers
        '''
        ### muse
        filename = filenamelist[sampleno * 4]
        print('Treat sample: ' + filename.split('.')[0], 'MuSE')
        self.read_vcf(filepath + filename)
        self.vcf = self.pretreat(self.vcf, posstr=False)
        self.vcf = self.sbs2dbs(self.vcf)
        #self.vcf = self.sbs2dbs_with_loop(self.vcf)
        allvcf = self.muse(self.vcf)
        
        ### mutect2
        print('Treat sample: ' + filename.split('.')[0], 'Mutect2')
        filename = filenamelist[sampleno * 4 + 1]
        self.read_vcf(filepath + filename)
        header = self.header[:-1]
        self.vcf = self.pretreat(self.vcf, posstr=True)
        allvcf = pd.concat([allvcf, self.mutect2(self.vcf)], axis = 0)
        
        ### strelka
        print('Treat sample: ' + filename.split('.')[0], 'Strelka')
        filename = filenamelist[sampleno * 4 + 2]
        self.read_vcf(filepath + filename)
        self.vcf = self.pretreat(self.vcf, posstr=True)
        allvcf = pd.concat([allvcf, self.strelka(self.vcf)], axis = 0)
        
        ### varscan
        print('Treat sample: ' + filename.split('.')[0], 'Varscan')
        filename = filenamelist[sampleno * 4 + 3]
        self.read_vcf(filepath + filename)
        self.vcf = self.pretreat(self.vcf, posstr=True)
        allvcf = pd.concat([allvcf, self.varscan(self.vcf)], axis = 0)
        
        ### integrate by caller 2+
        print('Treat sample: ' + filename.split('.')[0], 'Integration')
        savevcf = self.twoplus(allvcf)
        ### transfer SBS to DBS
        savevcf = self.pretreat(savevcf, posstr=False)
        savevcf = self.sbs2dbs(savevcf)
        #savevcf = self.sbs2dbs_with_loop(savevcf)
        ### remove duplicates
        savevcf = self.remove_dup(savevcf)
        ### save savevcf
        self.savevcf(savevcf = savevcf, header = header, savepath = savepath, filename = filename, postfix = '.call2plus.vcf')
