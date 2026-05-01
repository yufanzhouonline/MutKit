########################################################################
##### A python library for reading and analysis of mutation data  ######
#####                by Yufan (Harry) Zhou at 2024/12/20          ######
########################################################################

from mutpattern import *

class Mutsv(Mutpattern):
    '''
    Integration of four SV callers: manta, gridss, svaba and delly
    '''
    def manta(self, filepath = '', filenamelist = '', somatic_score_cutoff = 39):
        '''
        Read manta SV
        '''
        super().read_sv(filepath = filepath, filenamelist = filenamelist, somatic_score_cutoff = somatic_score_cutoff)
    
    def svaba(self, filepath = '', filenamelist = ''):
        '''
        Read svaba SV
        '''
        self.svvcflist = []
        for i in filenamelist:
            print(i)
            self.read_vcf(filepath + i, rename_columns = False)
            columns = self.header[-1].split('\t')
            #columns = ['CHROM'] + columns[1:9] + ['c10', 'c11', 'c12'] + columns[9:]
            columns = ['CHROM'] + columns[1:]
            if self.vcf.shape[0]:
                self.vcf.columns = columns
                self.vcf['svtype'] = self.vcf.apply(lambda row : re.search(r'SVTYPE=([A-Z]+)', row[self.vcf.columns[7]]).group(1), axis=1)
                self.vcf = self.vcf[self.vcf.CHROM.isin(self.chrlist)].copy()
                self.vcf.reset_index(drop = True, inplace = True)
            self.svvcflist.append(self.vcf)
    
    def delly(self, filepath = '', filenamelist = ''):
        '''
        Read delly SV
        '''
        self.svvcflist = []
        for i in filenamelist:
            print(i)
            self.read_vcf(filepath + i, rename_columns = False)
            columns = self.header[-1].split('\t')
            columns = ['CHROM'] + columns[1:]
            if self.vcf.shape[0]:
                self.vcf.columns = columns
                self.vcf['svtype'] = self.vcf.apply(lambda row : re.search(r'SVTYPE=([A-Z]+)', row[self.vcf.columns[7]]).group(1), axis=1)
                self.vcf = self.vcf[self.vcf.CHROM.isin(self.chrlist)].copy()
                self.vcf.reset_index(drop = True, inplace = True)
            self.svvcflist.append(self.vcf)
    
    def gridss(self, filepath = '', filenamelist = '', qual_cutoff = 1000):
        '''
        Read gridss SV
        '''
        self.svvcflist = []
        for i in filenamelist:
            print(i)
            self.read_vcf(filepath + i, rename_columns = False)
            columns = self.header[-1].split('\t')
            columns = ['CHROM'] + columns[1:]
            if self.vcf.shape[0]:
                self.vcf.columns = columns
                self.vcf['QUAL'] = self.vcf['QUAL'].astype(float)
                #self.vcf['svtype'] = self.vcf.apply(lambda row : re.search(r'SVTYPE=([A-Z]+)', row[self.vcf.columns[7]]).group(1), axis=1)
                ### get simple type of sv by R
                self.vcf['svtype'] = self.vcf.apply(lambda row : re.search(r'SIMPLE_TYPE=([A-Z]+)', row[self.vcf.columns[7]]).group(1) if re.search(r'SIMPLE_TYPE=([A-Z]+)', row[self.vcf.columns[7]]) else 'BND', axis=1)
                self.vcf = self.vcf[(self.vcf.CHROM.isin(self.chrlist)) & (self.vcf.QUAL >=qual_cutoff)].copy()
                self.vcf.reset_index(drop = True, inplace = True)
            self.svvcflist.append(self.vcf)
    
    def pre_two_plus(self):
        '''
        Treatment before two caller plus
        '''
        callerlist = []
        for i in range(len(self.svvcflist)):
            print(i)
            pdcaller = pd.DataFrame({'CHROM':[], 'POS':[], 'newchr':[], 'newpos':[], 'svtype':[]})
            pdcaller = pdcaller.astype({'CHROM':str, 'POS':int, 'newchr':str, 'newpos':int, 'svtype':str})
            if self.svvcflist[i].shape[0]:
                svi = self.svvcflist[i][self.svvcflist[i].svtype=='BND'].copy()
                if svi.shape[0]:
                    svi.reset_index(drop = True, inplace = True)
                    svi['newchr'] = svi.apply(lambda row : re.search(r'chr(\d{1,2}|X|Y)', row[svi.columns[4]]).group(0) if re.search(r'chr(\d{1,2}|X|Y)', row[svi.columns[4]]) else 'chr-', axis=1)
                    svi['newpos'] = svi.apply(lambda row : int(re.search(r':(\d+)', row[svi.columns[4]]).group(1)) if re.search(r':(\d+)', row[svi.columns[4]]) else 0, axis=1)
                    pdcaller = pd.concat([pdcaller.loc[:, ['CHROM', 'POS', 'newchr', 'newpos', 'svtype']], svi.loc[:, ['CHROM', 'POS', 'newchr', 'newpos', 'svtype']]], axis=0)
                svj = self.svvcflist[i][self.svvcflist[i].svtype!='BND'].copy()
                if svj.shape[0]:
                    svj['newchr'] = svj['CHROM']
                    svj['newpos'] = svj['POS']
                    pdcaller = pd.concat([pdcaller.loc[:, ['CHROM', 'POS', 'newchr', 'newpos', 'svtype']], svj.loc[:, ['CHROM', 'POS', 'newchr', 'newpos', 'svtype']]], axis=0)
            pdcaller = pdcaller.astype({'CHROM':str, 'POS':int, 'newchr':str, 'newpos':int, 'svtype':str})
            pdcaller = pdcaller[(pdcaller.newchr != 'chr-') | (pdcaller.newpos != 0)].copy()
            pdcaller.reset_index(drop = True, inplace = True)
            callerlist.append(pdcaller)
        return callerlist
    
    def sv_four_callers(self, sampleno = 10, filepath = '', filenamelist1 = '', filenamelist2 = '', filenamelist3 = '', filenamelist4 = '', savepath= '', somatic_score_cutoff = 39,  qual_cutoff = 500):
        '''
        Get SVs in at least present in two callers
        '''
        ### gridss
        print('Treat sample: gridss')
        self.gridss(filepath = filepath, filenamelist = filenamelist1, qual_cutoff = qual_cutoff)
        gridsslist = self.pre_two_plus()
        ### manta
        print('Treat sample: manta')
        self.manta(filepath = filepath, filenamelist = filenamelist2, somatic_score_cutoff = somatic_score_cutoff)
        mantalist = self.pre_two_plus()
        ### svaba
        print('Treat sample: svaba')
        self.svaba(filepath = filepath, filenamelist = filenamelist3)
        svabalist = self.pre_two_plus()
        ### delly
        print('Treat sample: delly')
        self.delly(filepath = filepath, filenamelist = filenamelist4)
        dellylist = self.pre_two_plus()
        ### integrate
        twocallerlist = []
        ### integrate
        for iii in range(sampleno):
            gridsslist[iii]['caller'] = 1
            mantalist[iii]['caller'] = 2
            svabalist[iii]['caller'] = 3
            dellylist[iii]['caller'] = 4
            alllist = pd.concat([gridsslist[iii], mantalist[iii], svabalist[iii], dellylist[iii]], axis=0)
            alllist.sort_values(by=['svtype', 'CHROM', 'POS', 'newchr', 'newpos', 'caller'], ascending=True, inplace=True)
            alllist.reset_index(drop = True, inplace = True)
            alllist['overlap'] = 1
            for i in alllist.index[1:]:
                ### conditions which is overlapped SV
                if alllist.loc[i, 'svtype'] == alllist.loc[i-1, 'svtype'] and alllist.loc[i, 'CHROM'] == alllist.loc[i-1, 'CHROM'] and alllist.loc[i, 'newchr'] == alllist.loc[i-1, 'newchr'] and abs(alllist.loc[i, 'POS'] - alllist.loc[i-1, 'POS']) <= 100 and abs(alllist.loc[i, 'newpos'] == alllist.loc[i-1, 'newpos']) <= 100 and alllist.loc[i, 'caller'] != alllist.loc[i-1, 'caller']:
                    alllist.loc[i, 'overlap'] = alllist.loc[i-1, 'overlap'] + 1
            twocaller = alllist[alllist.overlap==2].copy()
            twocaller.reset_index(drop = True, inplace = True)
            twocallerlist.append(twocaller)
            print(twocaller.shape[0])
            fileprefix = filenamelist1[iii].split('.')[0]
            twocaller.iloc[:,:5].to_csv(savepath + fileprefix + '.combined.txt', index=False, header=True, sep='\t')
        return twocallerlist

