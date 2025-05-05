
from mutlib import *

filepath = '/data/a3b4nqo/combined2/passvcf/'
savepath= '/data/a3b4nqo/combined2/caller2plus/'
filenamelist = ['E4_NQO_450_tumor_S8.muse.pass.vcf', 'E4_NQO_450_tumor_S8.mut.pass.vcf', 'E4_NQO_450_tumor_S8.strelka.pass.vcf', 'E4_NQO_450_tumor_S8.varscan.pass.vcf',
                'E4_NQO_453_tumor_S2.muse.pass.vcf', 'E4_NQO_453_tumor_S2.mut.pass.vcf', 'E4_NQO_453_tumor_S2.strelka.pass.vcf', 'E4_NQO_453_tumor_S2.varscan.pass.vcf',
                'E4_NQO_459_tumor_S4.muse.pass.vcf', 'E4_NQO_459_tumor_S4.mut.pass.vcf', 'E4_NQO_459_tumor_S4.strelka.pass.vcf', 'E4_NQO_459_tumor_S4.varscan.pass.vcf',
                'E4_NQO_469_tumor_S6.muse.pass.vcf', 'E4_NQO_469_tumor_S6.mut.pass.vcf', 'E4_NQO_469_tumor_S6.strelka.pass.vcf', 'E4_NQO_469_tumor_S6.varscan.pass.vcf',
                'CD_4NQO_175_tumor_S24.muse.pass.vcf', 'CD_4NQO_175_tumor_S24.mut.pass.vcf', 'CD_4NQO_175_tumor_S24.strelka.pass.vcf', 'CD_4NQO_175_tumor_S24.varscan.pass.vcf', 
                'CD_4NQO_178_tumor_S26.muse.pass.vcf', 'CD_4NQO_178_tumor_S26.mut.pass.vcf', 'CD_4NQO_178_tumor_S26.strelka.pass.vcf', 'CD_4NQO_178_tumor_S26.varscan.pass.vcf', 
                'CD_4NQO_247_tumor_S28.muse.pass.vcf', 'CD_4NQO_247_tumor_S28.mut.pass.vcf', 'CD_4NQO_247_tumor_S28.strelka.pass.vcf', 'CD_4NQO_247_tumor_S28.varscan.pass.vcf', 
                'CD_4NQO_271_tumor_S6.muse.pass.vcf', 'CD_4NQO_271_tumor_S6.mut.pass.vcf', 'CD_4NQO_271_tumor_S6.strelka.pass.vcf', 'CD_4NQO_271_tumor_S6.varscan.pass.vcf', 
                'CD_4NQO_276_tumor_S31.muse.pass.vcf', 'CD_4NQO_276_tumor_S31.mut.pass.vcf', 'CD_4NQO_276_tumor_S31.strelka.pass.vcf', 'CD_4NQO_276_tumor_S31.varscan.pass.vcf', 
                'CD_4NQO_286_tumor_S33.muse.pass.vcf', 'CD_4NQO_286_tumor_S33.mut.pass.vcf', 'CD_4NQO_286_tumor_S33.strelka.pass.vcf', 'CD_4NQO_286_tumor_S33.varscan.pass.vcf', 
                'CD_4NQO_333_tumor_S1.muse.pass.vcf', 'CD_4NQO_333_tumor_S1.mut.pass.vcf', 'CD_4NQO_333_tumor_S1.strelka.pass.vcf', 'CD_4NQO_333_tumor_S1.varscan.pass.vcf', 
                'CD_4NQO_339_tumor_S9.muse.pass.vcf', 'CD_4NQO_339_tumor_S9.mut.pass.vcf', 'CD_4NQO_339_tumor_S9.strelka.pass.vcf', 'CD_4NQO_339_tumor_S9.varscan.pass.vcf', 
                'CD_4NQO_347_tumor_S2.muse.pass.vcf', 'CD_4NQO_347_tumor_S2.mut.pass.vcf', 'CD_4NQO_347_tumor_S2.strelka.pass.vcf', 'CD_4NQO_347_tumor_S2.varscan.pass.vcf', 
                'CD_4NQO_401_tumor_S35.muse.pass.vcf', 'CD_4NQO_401_tumor_S35.mut.pass.vcf', 'CD_4NQO_401_tumor_S35.strelka.pass.vcf', 'CD_4NQO_401_tumor_S35.varscan.pass.vcf', 
                'CD_4NQO_456_tumor_S3.muse.pass.vcf', 'CD_4NQO_456_tumor_S3.mut.pass.vcf', 'CD_4NQO_456_tumor_S3.strelka.pass.vcf', 'CD_4NQO_456_tumor_S3.varscan.pass.vcf', 
                'CD_4NQO_457_tumor_S13.muse.pass.vcf', 'CD_4NQO_457_tumor_S13.mut.pass.vcf', 'CD_4NQO_457_tumor_S13.strelka.pass.vcf', 'CD_4NQO_457_tumor_S13.varscan.pass.vcf', 
                'CD_4NQO_463_tumor_S15.muse.pass.vcf', 'CD_4NQO_463_tumor_S15.mut.pass.vcf', 'CD_4NQO_463_tumor_S15.strelka.pass.vcf', 'CD_4NQO_463_tumor_S15.varscan.pass.vcf', 
                'CD_4NQO_468_mouth_tumor_S4.muse.pass.vcf', 'CD_4NQO_468_mouth_tumor_S4.mut.pass.vcf', 'CD_4NQO_468_mouth_tumor_S4.strelka.pass.vcf', 'CD_4NQO_468_mouth_tumor_S4.varscan.pass.vcf', 
                'CD_4NQO_468_tongue_tumor_S17.muse.pass.vcf', 'CD_4NQO_468_tongue_tumor_S17.mut.pass.vcf', 'CD_4NQO_468_tongue_tumor_S17.strelka.pass.vcf', 'CD_4NQO_468_tongue_tumor_S17.varscan.pass.vcf', 
                'CD_4NQO_474_tumor_S19.muse.pass.vcf', 'CD_4NQO_474_tumor_S19.mut.pass.vcf', 'CD_4NQO_474_tumor_S19.strelka.pass.vcf', 'CD_4NQO_474_tumor_S19.varscan.pass.vcf', 
                'CD_4NQO_475_tumor_S21.muse.pass.vcf', 'CD_4NQO_475_tumor_S21.mut.pass.vcf', 'CD_4NQO_475_tumor_S21.strelka.pass.vcf', 'CD_4NQO_475_tumor_S21.varscan.pass.vcf', 
                'CD_4NQO_476_tumor_S5.muse.pass.vcf', 'CD_4NQO_476_tumor_S5.mut.pass.vcf', 'CD_4NQO_476_tumor_S5.strelka.pass.vcf', 'CD_4NQO_476_tumor_S5.varscan.pass.vcf']

vcf = Mutlib()

###########
### integration of 4 callers

for sampleno in range(len(filenamelist)/4):
    vcf.fourcallers(sampleno = sampleno, filenamelist = filenamelist, filepath = filepath, savepath=savepath)

print('...END...')
