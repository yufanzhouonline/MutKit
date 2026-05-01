########################################################################
#####     A python library for Analysis of mutational patterns     #####
#####                by Yufan (Harry) Zhou at 2024/10/08          ######
########################################################################

from mutlib import *
import re
from scipy.stats import pearsonr
from scipy.stats import linregress

class Mutpattern(Mutlib):
    '''
    Analysis of mutational patterns
    '''
    def get_seq_dict(self, filenamepath = '/data/harry/mm10chr/mm10'):
        '''
        Get whole genome sequence
        '''
        self.seqdict = {}
        for chrno in self.chrlist:
            filename = filenamepath + chrno + '.fa'
            print('Read', filename)
            with open(filename, 'rt') as f:
                lines = [line[:-1] for line in f]
            geneseq = ''.join(lines[1:])
            self.seqdict.update({chrno:geneseq})
    
    def add_value(self):
        '''
        Add SNV types and mutational context of VCF files
        '''
        self.vcf['snvtype'] =  self.vcf.apply(lambda row : 'SBS' if len(row['REF']) == 1 and len(row['ALT']) == 1 else 'DBS' if len(row['REF']) == len(row['ALT']) else 'INDEL', axis=1)
        self.vcf['context'] =  self.vcf.apply(lambda row : self.seqdict[row['CHROM']][(int(row['POS'])-2):(int(row['POS'])+1)].lower(), axis=1)
    
    def read_context(self, filepath = '', filenamelist = '', AF_cutoff = 0.05, DP_cutoff = 19):
        '''
        Read mutational motifs and context of VCF files.
        '''
        self.mutvcflist = []
        self.unfilterlist = []
        self.filterlist = []
        self.sbslist = []
        self.dbslist = []
        self.indellist = []
        self.apobeclist = []
        self.apobecsbslist = []
        self.apobecdbslist = []
        self.apobecindellist = []
        for filename in filenamelist:
            print('Read', filename)
            self.read_vcf(filepath + filename)
            self.add_value()
            self.vcf['DP'] = self.vcf['DP'].astype(int)
            self.vcf['AF'] = self.vcf['AF'].astype(float)
            self.unfilterlist.append(self.vcf.shape[0])
            self.vcf = self.vcf[(self.vcf.AF > AF_cutoff) & (self.vcf.DP > DP_cutoff)].copy()
            self.vcf.reset_index(drop = True, inplace = True)
            apobec = self.vcf[(self.vcf.context == 'tca') | (self.vcf.context == 'tct')].copy()
            self.filterlist.append(self.vcf.shape[0])
            self.sbslist.append(self.vcf[self.vcf.snvtype=='SBS'].shape[0])
            self.dbslist.append(self.vcf[self.vcf.snvtype=='DBS'].shape[0])
            self.indellist.append(self.vcf[self.vcf.snvtype=='INDEL'].shape[0])
            self.apobeclist.append(apobec.shape[0])
            self.apobecsbslist.append(apobec[apobec.snvtype=='SBS'].shape[0])
            self.apobecdbslist.append(apobec[apobec.snvtype=='DBS'].shape[0])
            self.apobecindellist.append(apobec[apobec.snvtype=='INDEL'].shape[0])
            self.mutvcflist.append(self.vcf)
    
    def barplot(self, group1, group2, group3, ylabel = 'APOBEC signature(%)', sig1 = 'ns', sig2 = '*p < 0.1', sig3 = '***p < 0.01'):
        '''
        Make bar plot for three groups.
        '''
        meangroup1 = np.mean(group1)
        meangroup2 = np.mean(group2)
        meangroup3 = np.mean(group3)
        
        semgroup1 = np.std(group1, ddof=1) / np.sqrt(np.size(group1))
        semgroup2 = np.std(group2, ddof=1) / np.sqrt(np.size(group2))
        semgroup3 = np.std(group3, ddof=1) / np.sqrt(np.size(group3))
        
        means = [meangroup1, meangroup2, meangroup3]
        sems = [semgroup1, semgroup2, semgroup3]
        
        # Set up the plot
        x_labels = ['WT', 'E255A', 'A3B']
        x = np.arange(len(x_labels))  # the label locations
        width = 0.25  # the width of the bars
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Plot bars with error bars
        bars = ax.bar(x, means, width, yerr=sems, capsize=10, color='white', edgecolor='black')
        
        dotx = np.concatenate((np.random.rand(len(group1)) * 0.2 - 0.1, np.random.rand(len(group2)) * 0.2 + 1 - 0.1, np.random.rand(len(group3)) * 0.2 + 2 - 0.1))
        doty = np.concatenate((group1, group2, group3))
        
        # Individual data points
        ax.scatter(dotx, doty, color='black', alpha=0.6)
        
        maxpoint1 = np.max(group1)
        maxpoint2 = np.max(group2)
        maxpoint3 = np.max(group3)
        maxpointall = np.max([maxpoint1, maxpoint2, maxpoint3])
        
        ax.plot([0, 0], [maxpoint1 * 1.1, np.max([maxpoint1, maxpoint2]) * 1.2], color='black')
        ax.plot([1, 1], [maxpoint2 * 1.1, np.max([maxpoint1, maxpoint2]) * 1.2], color='black')
        ax.plot([0, 1], [np.max([maxpoint1, maxpoint2]) * 1.2, np.max([maxpoint1, maxpoint2]) * 1.2], color='black')
        
        ax.plot([1, 1], [np.max([maxpoint1, maxpoint2]) * 1.25, maxpointall * 1.3], color='black')
        ax.plot([2, 2], [maxpoint3 * 1.1, maxpointall * 1.3], color='black')
        ax.plot([1, 2], [maxpointall * 1.3, maxpointall * 1.3], color='black')
        
        ax.plot([0, 0], [np.max([maxpoint1, maxpoint2]) * 1.25, maxpointall * 1.4], color='black')
        ax.plot([2, 2], [maxpointall * 1.35, maxpointall * 1.4], color='black')
        ax.plot([0, 2], [maxpointall * 1.4, maxpointall * 1.4], color='black')
        
        plt.text(0.5, np.max([maxpoint1, maxpoint2]) * 1.23, sig1, fontsize=20, ha='center', color='black')
        plt.text(1.5, maxpointall * 1.33, sig2, fontsize=20, ha='center', color='black')
        plt.text(1, maxpointall * 1.43, sig3, fontsize=20, ha='center', color='black')
        
        # Add labels and title
        ax.set_ylabel(ylabel, fontsize = 20)
        #ax.set_title('Mean Values with Standard Error for Each Group')
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, fontsize = 20)
        #ax.grid(axis='y')
        # Set y-axis limits
        ax.set_ylim(0, maxpointall * 1.53)  # Correct way to set y limits
        
        # Set y-ticks font size
        ax.tick_params(axis='y', labelsize=20)  # Correct way to set y-tick font size
        
        # Display the plot
        plt.tight_layout()
        plt.show()
    
    def read_sv(self, filepath = '', filenamelist = '', somatic_score_cutoff = 39):
        '''
        Read structural variation files.
        '''
        self.svvcflist = []
        self.countlist = []
        self.somaticscorelist = []
        self.dellist = []
        self.duplist = []
        self.invlist = []
        self.bndlist = []
        self.conditionlist = []
        for filename in filenamelist:
            print('Read', filename)
            self.read_vcf(filepath + filename)
            self.vcf = self.pretreat(self.vcf, posstr=True)
            self.vcf['somaticscore'] = self.vcf.apply(lambda row : int(re.search(r'SOMATICSCORE=(\d+)', row[self.vcf.columns[7]]).group(1)), axis=1)
            self.vcf['svtype'] = self.vcf.apply(lambda row : row[self.vcf.columns[7]].split(';')[0].split('=')[1] if row[self.vcf.columns[7]][:6] == 'SVTYPE' else row[self.vcf.columns[7]].split(';')[1].split('=')[1], axis=1)
            self.vcf['pr'] = self.vcf.apply(lambda row : int(row[self.vcf.columns[10]].split(':')[0].split(',')[1]), axis=1)
            self.vcf['sr'] = self.vcf.apply(lambda row : int(row[self.vcf.columns[10]].split(':')[1].split(',')[1] if len(row[self.vcf.columns[10]].split(':')) == 2 else row[self.vcf.columns[10]].split(',')[1]), axis=1)
            self.conditionlist.append(self.vcf[self.vcf.somaticscore > somatic_score_cutoff].shape[0])
            ### filter with somaticscore > somatic_score_cutoff
            self.vcf = self.vcf[self.vcf.somaticscore > somatic_score_cutoff].copy()
            self.vcf.reset_index(drop = True, inplace = True)
            self.countlist.append(self.vcf.shape[0])
            self.somaticscorelist.append(self.vcf['somaticscore'])
            self.dellist.append(self.vcf[self.vcf.svtype=='DEL'].shape[0])
            self.duplist.append(self.vcf[self.vcf.svtype=='DUP'].shape[0])
            self.invlist.append(self.vcf[self.vcf.svtype=='INV'].shape[0])
            self.bndlist.append(self.vcf[self.vcf.svtype=='BND'].shape[0])
            self.vcf = self.vcf[self.vcf.CHROM.isin(self.chrlist)].copy()
            self.vcf.reset_index(drop = True, inplace = True)
            self.svvcflist.append(self.vcf)
    
    def plot_somatic_score(self):
        '''
        Plot somatic score histogram.
        '''
        somaticscore = self.somaticscorelist
        
        x_values_list = []
        y_values_list = []
        for i in range(len(somaticscore)):
            x_values = list(i * 2 + np.random.rand(somaticscore[i].shape[0]))
            y_values = list(somaticscore[i])
            x_values_list = x_values_list + x_values
            y_values_list = y_values_list + y_values
        
        # Create a dot plot
        plt.figure(figsize=(20, 8))
        plt.scatter(x_values_list, y_values_list, color='blue', s=3)  # 's' is the size of the dots
        plt.title('Dot Plot of Somatic Score of SV', fontsize=20)
        plt.xlabel('Sample', fontsize=20)
        plt.ylabel('Somatic score', fontsize=20)
        plt.grid(False)
        
        # Set x-axis ticks and labels
        plt.xticks(ticks=np.arange(0.5, 44.5, 2), labels=[str(i) for i in range(1, 23)])
        
        # Customize tick labels
        plt.tick_params(axis='x', labelsize=20)  # Font size for x-axis tick labels
        plt.tick_params(axis='y', labelsize=20)  # Font size for y-axis tick labels
        
        # Show the plot
        plt.show()
    
    def plotsv(self, values, labels, color = 'skyblue', label = 'SV'):
        '''
        Plot counts of SV, values = vcf.countlist, labels = samplename
        '''
        # Create a bar plot
        plt.bar(range(22), values, color=color)
        
        # Add numbers on top of the bars
        for i, value in enumerate(values):
            plt.text(i, value + 0.5, str(value), ha='center', va='bottom', fontsize=16)
        
        # Labeling
        #plt.xlabel('Categories')
        plt.ylabel('# of ' + label, fontsize=20)
        #plt.title('SV called by Manta')
        
        # Set x-axis ticks and labels
        #plt.xticks(ticks=range(22), labels=[str(i) for i in range(1, 23)])
        plt.xticks(ticks=range(22), labels=labels)
        
        # Customize tick labels
        plt.tick_params(axis='x', labelsize = 20, labelrotation = 45)  # Font size for x-axis tick labels
        plt.tick_params(axis='y', labelsize = 20)  # Font size for y-axis tick labels
        
        plt.ylim(0, max(values) * 1.1)
        # Show plot
        plt.show()
    
    def plot_correlation(self, list1, list2, xlabel, ylabel):
        '''
        Plot correlation
        '''
        pearson_correlation, p_value = pearsonr(list1, list2)
        print('Pearson correlation coefficient:', pearson_correlation)
        print('p value:', p_value)
        
        allvalue = linregress(list1, list2)
        slope = allvalue[0]
        intercept = allvalue[1]
        
        plt.scatter(list1, list2, color='blue', s=20)
        plt.plot([0, max(list1)], [intercept, max(list1) * slope + intercept], color='black')
        plt.tick_params(axis='x', labelsize = 20)
        plt.tick_params(axis='y', labelsize = 20)
        plt.xlabel(xlabel, fontsize=20)
        plt.ylabel(ylabel, fontsize=20)
        plt.show()
    
    def count_enrichment(self, indexlist, mutvcflist, svvcflist, tcw = False, upstream = -100000, downstream = 100000):
        '''
        Count enrichment of mutations in SV
        '''
        samplepoints = []
        for i in indexlist:
            print('Sample No.', i)
            mutdf = mutvcflist[i].copy()
            if tcw:
                mutdf = mutdf[(mutdf.context == 'tca') | (mutdf.context == 'tct')].copy()
            mutdf['POS'] = mutdf['POS'].astype(int)
            mutdf.sort_values(by=['CHROM', 'POS'], ascending=True, inplace=True)
            mutdf.reset_index(drop = True, inplace = True)
            
            svdf = svvcflist[i].copy()
            svdf['POS'] = svdf['POS'].astype(int)
            svdf.sort_values(by=['CHROM', 'POS'], ascending=True, inplace=True)
            svdf.reset_index(drop = True, inplace = True)
            
            allpoints = []
            for svindex in svdf.index:
                targetdistance = mutdf[mutdf.CHROM==svdf.loc[svindex, 'CHROM']].POS - svdf.loc[svindex, 'POS'].copy()
                getpoints = targetdistance[(targetdistance > upstream) & (targetdistance < downstream)]
                allpoints = allpoints + list(getpoints)
            
            samplepoints = samplepoints + allpoints
        return samplepoints
    
    def plot_enrichment(self, samplepoints1, samplepoints2, samplepoints3):
        '''
        Plot enrichment
        '''
        # Create a figure and an array of subplots
        fig, axs = plt.subplots(3, 1, figsize=(8, 10))  # 3 rows, 1 column
        
        axs[0].hist(samplepoints1, bins=200, alpha=0.7, color='#AAAA00AA', edgecolor='black')
        axs[0].tick_params(axis='x', labelsize = 20)
        axs[0].tick_params(axis='y', labelsize = 20)
        #axs[0].set_xlabel('SV -100K to +100K', fontsize=20)
        #axs[0].set_ylabel('# of mutations', fontsize=20)
        
        axs[1].hist(samplepoints2, bins=200, alpha=0.7, color='#AA00AAAA', edgecolor='black')
        axs[1].tick_params(axis='x', labelsize = 20)
        axs[1].tick_params(axis='y', labelsize = 20)
        #axs[1].set_xlabel('SV -100K to +100K', fontsize=20)
        #axs[1].set_ylabel('# of mutations', fontsize=20)
        
        axs[2].hist(samplepoints3, bins=200, alpha=0.7, color='#00AAAAAA', edgecolor='black')
        axs[2].tick_params(axis='x', labelsize = 20)
        axs[2].tick_params(axis='y', labelsize = 20)
        #axs[2].set_xlabel('SV -100K to +100K', fontsize=20)
        #axs[2].set_ylabel('# of mutations', fontsize=20)
        
        # Adjust the layout
        plt.tight_layout()
        
        # Show the plot
        plt.show()
