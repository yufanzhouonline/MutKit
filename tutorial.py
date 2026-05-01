### for integration of four-caller mutations 

from mutlib import *

filepath = '/your passed mutation vcf path/passvcf/'
savepath= '/your saved mutation path/caller2plus/'
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

#####################

### integration of 4 callers

for sampleno in range(len(filenamelist)/4):
    vcf.fourcallers(sampleno = sampleno, filenamelist = filenamelist, filepath = filepath, savepath=savepath)

##########################
### for integration of four-caller structural variations


from mutsv import *

filepath = '/your passed structure variation vcf path/passvcf/'
savepath = '/your saved structure variation path/call2plus/'

filenamelist1 = [
                'CD_4NQO_271_tumor_S6.gridss.pass.vcf',
                'CD_4NQO_333_tumor_S1.gridss.pass.vcf',
                'CD_4NQO_347_tumor_S2.gridss.pass.vcf',
                'E4_NQO_450_tumor_S8.gridss.pass.vcf',
                'CD_4NQO_463_tumor_S15.gridss.pass.vcf',
                'CD_4NQO_475_tumor_S21.gridss.pass.vcf',
                'CD_4NQO_175_tumor_S24.gridss.pass.vcf',
                'CD_4NQO_178_tumor_S26.gridss.pass.vcf',
                'CD_4NQO_339_tumor_S9.gridss.pass.vcf',
                'CD_4NQO_401_tumor_S35.gridss.pass.vcf',
                'E4_NQO_453_tumor_S2.gridss.pass.vcf',
                'CD_4NQO_456_tumor_S3.gridss.pass.vcf',
                'CD_4NQO_457_tumor_S13.gridss.pass.vcf',
                'E4_NQO_459_tumor_S4.gridss.pass.vcf',
                'CD_4NQO_468_mouth_tumor_S4.gridss.pass.vcf',
                'CD_4NQO_468_tongue_tumor_S17.gridss.pass.vcf',
                'E4_NQO_469_tumor_S6.gridss.pass.vcf',
                'CD_4NQO_474_tumor_S19.gridss.pass.vcf',
                'CD_4NQO_476_tumor_S5.gridss.pass.vcf',
                'CD_4NQO_247_tumor_S28.gridss.pass.vcf',
                'CD_4NQO_276_tumor_S31.gridss.pass.vcf',
                'CD_4NQO_286_tumor_S33.gridss.pass.vcf'
]

filenamelist2 = [
                'CD_4NQO_271_tumor_S6.manta.pass.vcf',
                'CD_4NQO_333_tumor_S1.manta.pass.vcf',
                'CD_4NQO_347_tumor_S2.manta.pass.vcf',
                'E4_NQO_450_tumor_S8.manta.pass.vcf',
                'CD_4NQO_463_tumor_S15.manta.pass.vcf',
                'CD_4NQO_475_tumor_S21.manta.pass.vcf',
                'CD_4NQO_175_tumor_S24.manta.pass.vcf',
                'CD_4NQO_178_tumor_S26.manta.pass.vcf',
                'CD_4NQO_339_tumor_S9.manta.pass.vcf',
                'CD_4NQO_401_tumor_S35.manta.pass.vcf',
                'E4_NQO_453_tumor_S2.manta.pass.vcf',
                'CD_4NQO_456_tumor_S3.manta.pass.vcf',
                'CD_4NQO_457_tumor_S13.manta.pass.vcf',
                'E4_NQO_459_tumor_S4.manta.pass.vcf',
                'CD_4NQO_468_mouth_tumor_S4.manta.pass.vcf',
                'CD_4NQO_468_tongue_tumor_S17.manta.pass.vcf',
                'E4_NQO_469_tumor_S6.manta.pass.vcf',
                'CD_4NQO_474_tumor_S19.manta.pass.vcf',
                'CD_4NQO_476_tumor_S5.manta.pass.vcf',
                'CD_4NQO_247_tumor_S28.manta.pass.vcf',
                'CD_4NQO_276_tumor_S31.manta.pass.vcf',
                'CD_4NQO_286_tumor_S33.manta.pass.vcf'
]

filenamelist3 = [
                'CD_4NQO_271_tumor_S6.svaba.pass.vcf',
                'CD_4NQO_333_tumor_S1.svaba.pass.vcf',
                'CD_4NQO_347_tumor_S2.svaba.pass.vcf',
                'E4_NQO_450_tumor_S8.svaba.pass.vcf',
                'CD_4NQO_463_tumor_S15.svaba.pass.vcf',
                'CD_4NQO_475_tumor_S21.svaba.pass.vcf',
                'CD_4NQO_175_tumor_S24.svaba.pass.vcf',
                'CD_4NQO_178_tumor_S26.svaba.pass.vcf',
                'CD_4NQO_339_tumor_S9.svaba.pass.vcf',
                'CD_4NQO_401_tumor_S35.svaba.pass.vcf',
                'E4_NQO_453_tumor_S2.svaba.pass.vcf',
                'CD_4NQO_456_tumor_S3.svaba.pass.vcf',
                'CD_4NQO_457_tumor_S13.svaba.pass.vcf',
                'E4_NQO_459_tumor_S4.svaba.pass.vcf',
                'CD_4NQO_468_mouth_tumor_S4.svaba.pass.vcf',
                'CD_4NQO_468_tongue_tumor_S17.svaba.pass.vcf',
                'E4_NQO_469_tumor_S6.svaba.pass.vcf',
                'CD_4NQO_474_tumor_S19.svaba.pass.vcf',
                'CD_4NQO_476_tumor_S5.svaba.pass.vcf',
                'CD_4NQO_247_tumor_S28.svaba.pass.vcf',
                'CD_4NQO_276_tumor_S31.svaba.pass.vcf',
                'CD_4NQO_286_tumor_S33.svaba.pass.vcf'
]

filenamelist4 = [
                'CD_4NQO_271_tumor_S6.delly.pass.vcf',
                'CD_4NQO_333_tumor_S1.delly.pass.vcf',
                'CD_4NQO_347_tumor_S2.delly.pass.vcf',
                'E4_NQO_450_tumor_S8.delly.pass.vcf',
                'CD_4NQO_463_tumor_S15.delly.pass.vcf',
                'CD_4NQO_475_tumor_S21.delly.pass.vcf',
                'CD_4NQO_175_tumor_S24.delly.pass.vcf',
                'CD_4NQO_178_tumor_S26.delly.pass.vcf',
                'CD_4NQO_339_tumor_S9.delly.pass.vcf',
                'CD_4NQO_401_tumor_S35.delly.pass.vcf',
                'E4_NQO_453_tumor_S2.delly.pass.vcf',
                'CD_4NQO_456_tumor_S3.delly.pass.vcf',
                'CD_4NQO_457_tumor_S13.delly.pass.vcf',
                'E4_NQO_459_tumor_S4.delly.pass.vcf',
                'CD_4NQO_468_mouth_tumor_S4.delly.pass.vcf',
                'CD_4NQO_468_tongue_tumor_S17.delly.pass.vcf',
                'E4_NQO_469_tumor_S6.delly.pass.vcf',
                'CD_4NQO_474_tumor_S19.delly.pass.vcf',
                'CD_4NQO_476_tumor_S5.delly.pass.vcf',
                'CD_4NQO_247_tumor_S28.delly.pass.vcf',
                'CD_4NQO_276_tumor_S31.delly.pass.vcf',
                'CD_4NQO_286_tumor_S33.delly.pass.vcf'
]

sv = Mutsv()
sampleno = 22
twocallerlist = sv.sv_four_callers(sampleno = sampleno, filepath = filepath, filenamelist1 = filenamelist1, filenamelist2 = filenamelist2, filenamelist3 = filenamelist3, filenamelist4 = filenamelist4, savepath= savepath, somatic_score_cutoff = 39,  qual_cutoff = 500)

