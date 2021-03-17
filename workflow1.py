
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gwf workflow found in source activate gwfenviron.
#Change only grouping1,grouping1,grouping1 to change round
This is a workflow used for genomic prediction
by first performing EMMAX+EMMA200 GWAS on a 5/6 training 
population of clover, leaving 1/6 of all individuals for testing,
then using the top 200 or 500 significant SNPs (MAF=0.05 filter)
to predict GEBVs of testing population individuals using RF or BayesB.
The method is also repeated with 200 or 500 random SNPs for comparison.
That workflow consists of four parts.

Part1: Creating a training population and a testing population.
This part outputs:
1) A GWAS-ready genotype with only training individuals

Part2: Running EMMAX+EMMA200 GWAS.
A .pval file is generated with p-values of all SNPs.

Part3: Making four genotype files.
Outputs:
1) a list of top 25 significant SNPs,
2) a list of top 200 significant SNPs,
3) a list of top 500 significant SNPs, 
4) a list of the 25 random SNPs chosen,
5) a list of the 200 random SNPs chosen,
6) A list of the 500 random SNPs chosen,
7) A shortened genotype file including only top 25 SNPs,
8) A shortened genotype file including only top 200 SNPs,
9) A shortened genotype file including only top 500 SNPs,
10) A shortened genotype file including only 25 random SNPs 
11) A shortened genotype file including only 200 random SNPs 
12) A shortened genotype file including only 500 random SNPs.

Part4: Genomic prediction using BayesB
Genomic prediction is performed using the R-package
BGLR and using Random Forest (Caret package in R)

Part5: All predictions are combined into one file to create one correlation coefficient
"""

from gwf import Workflow

gwf=Workflow()




def script_caller(groups, geno, outputs, script_name):
    inputs = [groups,geno]
    outputs = outputs
    options = {
        'cores': 1,
        'account': 'nchain',
        'memory': '8g',
        'walltime': '08:00:00'
    }

    spec = '''
    python {script_name} {groups} {geno}
    '''.format(script_name=script_name, groups=groups,geno=geno)

    return inputs, outputs, options, spec


def GWAS_caller(phenofile, genofile,dicname,outputs):
    inputs = [phenofile,genofile]
    outputs = outputs
    options = {
        'cores': 9,
        'account': 'nchain',
        'memory': '8g',
        'walltime': '08:00:00'
    }
    

    spec = '''
    source ~/miniconda2/etc/profile.d/conda.sh
    conda activate myproject
    mkdir {dicname}
    cd {dicname}
    python /home/cks/NChain/faststorage/GWASimplementation/atgwas/src/gwa.py -o "" -a emmax -m 6 -r {phenofile} -f {genofile} --data_format="diploid_int" 
    '''.format(dicname=dicname,phenofile=phenofile,genofile=genofile)
    
    return inputs, outputs, options, spec


def Genofile_from_topSNPs(pval, geno,top25,top200,top500,random25,random200,random500,genotop25,genotop200,genotop500,genorandom25,genorandom200,genorandom500):
    inputs = [pval,geno]
    outputs = [top25,top200,top500,random25,random200,random500,genotop25,genotop200,genotop500,genorandom25,genorandom200,genorandom500]
    options = {
        'cores': 1,
        'account': 'nchain',
        'memory': '100g',
        'walltime': '01:00:00'
    }

    spec = '''
    python Make_genotypes_based_on_specified_SNPs.py {pval} {geno} {top25} {top200} {top500} {random25} {random200} {random500} {genotop25} {genotop200} {genotop500} {genorandom25} {genorandom200} {genorandom500} 
    '''.format(pval=pval,geno=geno,top25 = top25,top200=top200,top500=top500,random25=random25, random200=random200,random500=random500,genotop25=genotop25,genotop200=genotop200,genotop500=genotop500,genorandom25=genorandom25,genorandom200=genorandom200,genorandom500=genorandom500)
    
    
    return inputs, outputs, options, spec


def GP_script_caller(pheno,groupingsystem,genogroup1,genogroup2,genogroup3,genogroup4,genogroup5,genogroup6,predictions,correlations,weights):
    inputs= [pheno,groupingsystem,genogroup1,genogroup2,genogroup3,genogroup4,genogroup5,genogroup6]
    outputs=[predictions,correlations,weights]
    options = {
        'cores': 16,
        'account': 'nchain',
        'memory': '100g',
        'walltime': '02:00:00'    
    }
    
    spec='''
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate Rprogram
    Rscript GP_ML_RandomForest.R {pheno} {groupingsystem} {genogroup1} {genogroup2} {genogroup3} {genogroup4} {genogroup5} {genogroup6} {predictions} {correlations} {weights}     
    '''.format(pheno=pheno,groupingsystem=groupingsystem,genogroup1=genogroup1,genogroup2=genogroup2,genogroup3=genogroup3,genogroup4=genogroup4,genogroup5=genogroup5,genogroup6=genogroup6,predictions=predictions,correlations=correlations, weights=weights)
    
    return inputs, outputs, options, spec




gwf.target_from_template("Training_pop_generator",
                         script_caller(groups="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                       outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining1.csv","/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining2.csv","/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining3.csv","/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining4.csv","/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining5.csv","/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining6.csv"],
                                       script_name="GWAS_creating_trainingpop_20201123.py"))

gwf.target_from_template("GWAS_1",
                         GWAS_caller(dicname="Round1_tstpop1",phenofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Phenofile_iSize.csv",
                         genofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining1.csv",
                         outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop1/results/_pid1_iSizeMean_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_2",
                         GWAS_caller(dicname="Round1_tstpop2",phenofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Phenofile_iSize.csv",
                         genofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining2.csv",
                         outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop2/results/_pid1_iSizeMean_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_3",
                         GWAS_caller(dicname="Round1_tstpop3",phenofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Phenofile_iSize.csv",
                         genofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining3.csv",
                         outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop3/results/_pid1_iSizeMean_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_4",
                         GWAS_caller(dicname="Round1_tstpop4",phenofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Phenofile_iSize.csv",
                         genofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining4.csv",
                         outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop4/results/_pid1_iSizeMean_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_5",
                         GWAS_caller(dicname="Round1_tstpop5",phenofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Phenofile_iSize.csv",
                         genofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining5.csv",
                         outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop5/results/_pid1_iSizeMean_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_6",
                         GWAS_caller(dicname="Round1_tstpop6",phenofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Phenofile_iSize.csv",
                         genofile="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/GWAStraining6.csv",
                         outputs=["/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop6/results/_pid1_iSizeMean_emmax_none_t75.pvals"]))

gwf.target_from_template("Geno_Files_Generator_group1",
                         Genofile_from_topSNPs(pval="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop1/results/_pid1_iSizeMean_emmax_none_t75.pvals",
                                           geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                           top25="Round1_tstpop1_top25snps.txt",
                                           top200="Round1_tstpop1_top200snps.txt",
                                           top500="Round1_tstpop1_top500snps.txt",
                                           random25="Round1_tstpop1_random25snps.txt", 
                                           random200="Round1_tstpop1_random200snps.txt",
                                           random500="Round1_tstpop1_random500snps.txt",
                                           genotop25="Geno_based_on_top25SNPs_Round1_tstpop1.txt",
                                           genotop200="Geno_based_on_top200SNPs_Round1_tstpop1.txt",
                                           genotop500="Geno_based_on_top500SNPs_Round1_tstpop1.txt",
                                           genorandom25 ="Geno_based_on_random25SNPs_Round1_tstpop1.txt",
                                           genorandom200="Geno_based_on_random200SNPs_Round1_tstpop1.txt",
                                           genorandom500="Geno_based_on_random500SNPs_Round1_tstpop1.txt"))

gwf.target_from_template("Geno_Files_Generator_group2",
                         Genofile_from_topSNPs(pval="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop2/results/_pid1_iSizeMean_emmax_none_t75.pvals",
                                           geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                           top25 = "Round1_tstpop2_top25snps.txt", 
                                           top200="Round1_tstpop2_top200snps.txt",
                                           top500="Round1_tstpop2_top500snps.txt",
                                           random25="Round1_tstpop2_random25snps.txt",
                                           random200="Round1_tstpop2_random200snps.txt",
                                           random500="Round1_tstpop2_random500snps.txt",
                                           genotop25="Geno_based_on_top25SNPs_Round1_tstpop2.txt",
                                           genotop200="Geno_based_on_top200SNPs_Round1_tstpop2.txt",
                                           genotop500="Geno_based_on_top500SNPs_Round1_tstpop2.txt",
                                           genorandom25 = "Geno_based_on_random25SNPs_Round1_tstpop2.txt",
                                           genorandom200="Geno_based_on_random200SNPs_Round1_tstpop2.txt",
                                           genorandom500="Geno_based_on_random500SNPs_Round1_tstpop2.txt"))

gwf.target_from_template("Geno_Files_Generator_group3",
                         Genofile_from_topSNPs(pval="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop3/results/_pid1_iSizeMean_emmax_none_t75.pvals",
                                           geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                           top25="Round1_tstpop3_top25snps.txt",
                                           top200="Round1_tstpop3_top200snps.txt",
                                           top500="Round1_tstpop3_top500snps.txt",
                                           random25="Round1_tstpop3_random25snps.txt",
                                           random200="Round1_tstpop3_random200snps.txt",
                                           random500="Round1_tstpop3_random500snps.txt",
                                           genotop25="Geno_based_on_top25SNPs_Round1_tstpop3.txt",
                                           genotop200="Geno_based_on_top200SNPs_Round1_tstpop3.txt",
                                           genotop500="Geno_based_on_top500SNPs_Round1_tstpop3.txt",
                                           genorandom25="Geno_based_on_random25SNPs_Round1_tstpop3.txt",
                                           genorandom200="Geno_based_on_random200SNPs_Round1_tstpop3.txt",
                                           genorandom500="Geno_based_on_random500SNPs_Round1_tstpop3.txt"))

gwf.target_from_template("Geno_Files_Generator_group4",
                         Genofile_from_topSNPs(pval="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop4/results/_pid1_iSizeMean_emmax_none_t75.pvals",
                                           geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                           top25="Round1_tstpop4_top25snps.txt",
                                           top200="Round1_tstpop4_top200snps.txt",
                                           top500="Round1_tstpop4_top500snps.txt",
                                           random25="Round1_tstpop4_random25snps.txt",
                                           random200="Round1_tstpop4_random200snps.txt",
                                           random500="Round1_tstpop4_random500snps.txt",
                                           genotop25="Geno_based_on_top25SNPs_Round1_tstpop4.txt",
                                           genotop200="Geno_based_on_top200SNPs_Round1_tstpop4.txt",
                                           genotop500="Geno_based_on_top500SNPs_Round1_tstpop4.txt",
                                           genorandom25="Geno_based_on_random25SNPs_Round1_tstpop4.txt",
                                           genorandom200="Geno_based_on_random200SNPs_Round1_tstpop4.txt",
                                           genorandom500="Geno_based_on_random500SNPs_Round1_tstpop4.txt"))

gwf.target_from_template("Geno_Files_Generator_group5",
                         Genofile_from_topSNPs(pval="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop5/results/_pid1_iSizeMean_emmax_none_t75.pvals",
                                           geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                           top25="Round1_tstpop5_top25snps.txt",
                                           top200="Round1_tstpop5_top200snps.txt",
                                           top500="Round1_tstpop5_top500snps.txt",
                                           random25="Round1_tstpop5_random25snps.txt",
                                           random200="Round1_tstpop5_random200snps.txt",
                                           random500="Round1_tstpop5_random500snps.txt",
                                           genotop25="Geno_based_on_top25SNPs_Round1_tstpop5.txt",
                                           genotop200="Geno_based_on_top200SNPs_Round1_tstpop5.txt",
                                           genotop500="Geno_based_on_top500SNPs_Round1_tstpop5.txt",
                                           genorandom25="Geno_based_on_random25SNPs_Round1_tstpop5.txt",
                                           genorandom200="Geno_based_on_random200SNPs_Round1_tstpop5.txt",
                                           genorandom500="Geno_based_on_random500SNPs_Round1_tstpop5.txt"))

gwf.target_from_template("Geno_Files_Generator_group6",
                         Genofile_from_topSNPs(pval="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_20201120_Averages/GWAS_iSize_followedByRF_20210201/Round1_tstpop6/results/_pid1_iSizeMean_emmax_none_t75.pvals",
                                           geno="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/Preparation_of_genotypeFile/Geno_LDfiltered0.5and1.00_20200728.csv",
                                           top25="Round1_tstpop6_top25snps.txt",
                                           top200="Round1_tstpop6_top200snps.txt",
                                           top500="Round1_tstpop6_top500snps.txt",
                                           random25="Round1_tstpop6_random25snps.txt",
                                           random200="Round1_tstpop6_random200snps.txt",
                                           random500="Round1_tstpop6_random500snps.txt",
                                           genotop25="Geno_based_on_top25SNPs_Round1_tstpop6.txt",
                                           genotop200="Geno_based_on_top200SNPs_Round1_tstpop6.txt",
                                           genotop500="Geno_based_on_top500SNPs_Round1_tstpop6.txt",
                                           genorandom25="Geno_based_on_random25SNPs_Round1_tstpop6.txt",
                                           genorandom200="Geno_based_on_random200SNPs_Round1_tstpop6.txt",
                                           genorandom500="Geno_based_on_random500SNPs_Round1_tstpop6.txt"))

gwf.target_from_template("RandomForest_iSize_top25_Round1",
                         GP_script_caller(pheno="Phenofile_iSize.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_top25SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_top25SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_top25SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_top25SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_top25SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_top25SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_iSize_RF_top25SNPs_grouping1.txt",
                                          correlations="Correlations_iSize_RF_top25SNPs_grouping1.txt",
                                          weights="Weights_iSize_RF_top25SNPs_grouping1.txt"))



gwf.target_from_template("RandomForest_iSize_top200_Round1",
                         GP_script_caller(pheno="Phenofile_iSize.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_top200SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_top200SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_top200SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_top200SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_top200SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_top200SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_iSize_RF_top200SNPs_grouping1.txt",
                                          correlations="Correlations_iSize_RF_top200SNPs_grouping1.txt",
                                          weights="Weights_iSize_RF_top200SNPs_grouping1.txt"))


gwf.target_from_template("RandomForest_iSize_top500_Round1",
                         GP_script_caller(pheno="Phenofile_iSize.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_top500SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_top500SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_top500SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_top500SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_top500SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_top500SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_iSize_RF_top500SNPs_grouping1.txt",
                                          correlations="Correlations_iSize_RF_top500SNPs_grouping1.txt",
                                          weights="Weights_iSize_RF_top500SNPs_grouping1.txt"))



gwf.target_from_template("RandomForest_iSize_random25_Round1",
                         GP_script_caller(pheno="Phenofile_iSize.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_random25SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_random25SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_random25SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_random25SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_random25SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_random25SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_iSize_RF_random25SNPs_grouping1.txt",
                                          correlations="Correlations_iSize_RF_random25SNPs_grouping1.txt",
                                          weights="Weights_iSize_RF_random25SNPs_grouping1"))



gwf.target_from_template("RandomForest_iSize_random200_Round1",
                         GP_script_caller(pheno="Phenofile_iSize.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_random200SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_random200SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_random200SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_random200SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_random200SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_random200SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_iSize_RF_random200SNPs_grouping1.txt",
                                          correlations="Correlations_iSize_RF_random200SNPs_grouping1.txt",
                                          weights="Weights_iSize_RF_random200SNPs_grouping1"))

gwf.target_from_template("RandomForest_iSize_random500_Round1",
                         GP_script_caller(pheno="Phenofile_iSize.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_random500SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_random500SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_random500SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_random500SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_random500SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_random500SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_iSize_RF_random500SNPs_grouping1.txt",
                                          correlations="Correlations_iSize_RF_random500SNPs_grouping1.txt",
                                          weights="Weights_iSize_RF_random500SNPs_grouping1.txt"))


gwf.target_from_template("RandomForest_gpd_top25_Round1",
                         GP_script_caller(pheno="Phenofile_gpd.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_top25SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_top25SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_top25SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_top25SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_top25SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_top25SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_gpd_RF_top25SNPs_grouping1.txt",
                                          correlations="Correlations_gpd_RF_top25SNPs_grouping1.txt",
                                          weights="Weights_gpd_RF_top25SNPs_grouping1.txt"))


gwf.target_from_template("RandomForest_gpd_top200_Round1",
                         GP_script_caller(pheno="Phenofile_gpd.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_top200SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_top200SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_top200SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_top200SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_top200SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_top200SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_gpd_RF_top200SNPs_grouping1.txt",
                                          correlations="Correlations_gpd_RF_top200SNPs_grouping1.txt",
                                          weights="Weights_gpd_RF_top200SNPs_grouping1.txt"))

gwf.target_from_template("RandomForest_gpd_top500_Round1",
                         GP_script_caller(pheno="Phenofile_gpd.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_top500SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_top500SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_top500SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_top500SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_top500SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_top500SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_gpd_RF_top500SNPs_grouping1.txt",
                                          correlations="Correlations_gpd_RF_top500SNPs_grouping1.txt",
                                          weights="Weights_gpd_RF_top500SNPs_grouping1.txt"))

gwf.target_from_template("RandomForest_gpd_random25_Round1",
                         GP_script_caller(pheno="Phenofile_gpd.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_random25SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_random25SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_random25SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_random25SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_random25SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_random25SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_gpd_RF_random25SNPs_grouping1.txt",
                                          correlations="Correlations_gpd_RF_random25SNPs_grouping1.txt",
                                          weights="Weights_gpd_RF_random25SNPs_grouping1.txt"))

gwf.target_from_template("RandomForest_gpd_random200_Round1",
                         GP_script_caller(pheno="Phenofile_gpd.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_random200SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_random200SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_random200SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_random200SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_random200SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_random200SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_gpd_RF_random200SNPs_grouping1.txt",
                                          correlations="Correlations_gpd_RF_random200SNPs_grouping1.txt",
                                          weights="Weights_gpd_RF_random200SNPs_grouping1.txt"))

gwf.target_from_template("RandomForest_gpd_random500_Round1",
                         GP_script_caller(pheno="Phenofile_gpd.csv",
                                          groupingsystem="/home/cks/NChain/faststorage/WHITE_CLOVER/GP_V3_July2020_CKS_LessQualityFiltering/GP_2020-0826/gpd_ResCor_laterAndBetter_iSizeCor_AllRounds/grouping1.txt",
                                          genogroup1="Geno_based_on_random500SNPs_Round1_tstpop1.txt",
                                          genogroup2="Geno_based_on_random500SNPs_Round1_tstpop2.txt",
                                          genogroup3="Geno_based_on_random500SNPs_Round1_tstpop3.txt",
                                          genogroup4="Geno_based_on_random500SNPs_Round1_tstpop4.txt",
                                          genogroup5="Geno_based_on_random500SNPs_Round1_tstpop5.txt",
                                          genogroup6="Geno_based_on_random500SNPs_Round1_tstpop6.txt",
                                          predictions="Predictions_gpd_RF_random500SNPs_grouping1.txt",
                                          correlations="Correlations_gpd_RF_random500SNPs_grouping1.txt",
                                          weights="Weights_gpd_RF_random500SNPs_grouping1.txt"))