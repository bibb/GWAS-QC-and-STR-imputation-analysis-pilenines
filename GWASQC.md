## Individual level QCs for GWAS data include:
- Sex check: sex misclasification
- Heterozygosity rate: high or low heterozygosity rates
- Individual genotype missingness: removal of individuals with more than 10% of missing genotypes
- Identity by descent: cryptic relationships
- PCA: Population substructure


### Sex check
```
plink --bfile file --check-sex --out file.sexcheck
```
```
FID	IID	PEDSEX	SNPSEX	STATUS	F
001_10	001_10	1	2	PROBLEM	-1.316
003-08	003-08	1	2	PROBLEM	0.002542
003_10	003_10	1	0	PROBLEM	0.4857
004-08	004-08	1	0	PROBLEM	0.4699
004-10	004-10	1	0	PROBLEM	0.4892
004_06	004_06	2	2	OK	-0.5289
006-10	006-10	2	2	OK	0.002542
```
Manual inspection of outputfile: males with F values < 0.5 and females > 0.5 should be removed from analysis. Then, create file fail-sexcheck.txt with
the IDs of the removed individuals.

### Heterozygosity rate
```
plink --bfile file --het --out file.het
```
```
FID	IID	O(HOM)	E(HOM)	N(NM)	F
001_10	001_10	36370	3.61E+04	46379	0.02218
003-08	003-08	37546	3.77E+04	48338	-0.01449
003_10	003_10	36519	3.68E+04	47170	-0.02498
004-08	004-08	37585	3.75E+04	48141	0.004292
004-10	004-10	37370	3.75E+04	48116	-0.01368
004_06	004_06	36776	3.68E+04	47174	-0.0006843
006-10	006-10	37570	3.76E+04	48259	-0.00583
```
Calculate the observed heterozygosity rate per individual using the formula 
```
(N(NM) - O(Hom))/N(NM)
```
Create a graph where the observed heterozygosity rate per individual is plotted on the x-axis and the proportion of missing SNPs per individuals is plotted on the y-axis. 
This step is open to interpretation. A good rule of thum is to detect outlying individuals with more and less than 3 standard deviations from the mean for the result from (N(NM) - O(Hom))/N(NM). If you remove too many individuals, let's say > 10%, use a less strict cutoff, like 5 SD. 

Save the IDs of the outlying individuals in a file > fail-hetrate.txt

### Individual genotype missingness
```
plink --bfile file --missing --out file.missing
```
Examine the file file.missing.imiss and remove individuals with more than 10% of missing genotypes and save their ids to > fail-missingness.txt

### Identity by descent

This is a LD sensitive analysis, therefore prunning is necesary.
```
plink --bfile file --indep-pairwise 50 5 0.2 --out file
```
Select the file file.prune.in with the indidependent SNPs (r2 < 0.2) and run the genome command in plink.
```
plink --bfile file --extract file.prune.in --genome --out file.prune.ibd
```
In order to detect pairs of individuals with cryptic relatedness (pi_hat > 0.187) and keep only 1, you need to provide the missinges file, because is a good rule of thumb to select the individual with the highest genotyping rate. For this, the .genome and .imiss files should have the same name.
```
cp file.missing.imiss file.ibd.imiss
cp file.prune.ibd.genome file.ibd.genome
```
Then run the script ***run-IBD-QC.pl*** located in the resources folder to detect the related individuals and to generate the file fail-ibd.txt.
```
perl run-IBD-QC.pl file.ibd
```
This step si similar to the heterozygosity one, because if you remove too many individuals with this threshold, you can be less strict and rais the pi_hat cutoff up to 0.25, i.e. removing individuals with 2 degree of relatioship ans avobe. Explore your results and make a decision.

### Population substructure

In order to detect individuals with outliying population structure, we will run a PCA analysis using the 1000 genomes project
as a reference population. 

First, download the 1000 genomes phase3 genotypes and sites database from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

```
for i in {1..22}
do
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
done

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi
```

Now get the common variants between your GWAS dataset and the 1000 genomes project using [provided scripts](GWAS_QC_scripts/).
Manually generate a bim file from the 1000 genomes project sites vcf file. Also change rsnumbers in the ID column to a more unique identifier such as chr_bp_ref_alt
```
zcat ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | grep -v "#" | awk '{print $1 "\t" $1"_"$2"_"$4"_"$5 "\t" "0" "\t" $2 "\t" $4 "\t" $5}' > 1kg3p3.bim


perl bimAnnotationUpdate.pl file.bim < 1kg3p3.bim > file.bim.info
```
Examine SNPs carefully
```
cut -f 7,8 file.bim.info | sort | uniq -c
 121236 id      1       # same ID, same position, same alleles
   6739 id      -1      # same ID, same position, different  alleles
 359469 id      2       # same ID, same position, different order of alleles
     49 id      -9      # same ID, different position and alleles
  24280 no      0       # Snps not found in the reference panel
   7099 pos     1       # same position, same alleles, different ID
   1205 pos     -1      # same position, different ID and Alleles
  21423 pos     2       # same position, different order of alleles, different ID 
    201 pos     -9      # same position, different Alleles, different ID
```
Get only good SNPs
```
awk '$8>0{print $2}' file.bim.info > file.snps 

plink --bfile file --extract file.snps --make-bed --out  file.goodsnps

make a security copy of your original .bim file

cp file.goodsnps.bim file.goodsnps.bim.orig
```
Verify results
```
perl bimAnnotationUpdate.pl file.goodsnps.bim.orig < 1kg3p3.bim > file.goodsnps.bim.orig.info
```
Generate a new bim file for your GWAS, corrected for SNP ID and allele coding/orientation according to the 1kg3 reference
```
awk '{if($8==2){a=$14;$14=$15;$15=a};print $10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' \
  file.goodsnps.bim.orig.info > file.goodsnps.bim
```
Now change the SNP ids on the 1kg3 genotypes vcfs

```
for i in {1..22}
do
zcat ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep '^#' > ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.header

zcat ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '^#' | awk 'BEGIN {OFS="\t"}{$3=$1"_"$2"_"$4"_"$5} {print}' | cat ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.header - | gzip -c >  ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.vcf.gz

rm ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.header
done
```
Create a file with only SNP IDs from your GWAS data (file.ID.txt) and use it to extract them from the edited 1kg3p3 files with Vcftools
```
for i in {1..22}
do
vcftools --gzvcf ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.vcf.gz --snps file.ID.txt --recode --recode-INFO-all --out ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered
done 
```
Concatenate all chromosomes 1 to 22 from the resulting VCFs from 1kg3p3. First create a list of the VCF files:
```
for i in {1..22}
do
echo ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered.recode.vcf >> vcf_input.txt
done
```
Concatenate with vcf-concat tool from Vcftools (I would recommend to export the VCFtools directory into your PERL PATH enviroment first, i.e. export PERL5LIB=/path/to/your/installation/perl)
```
vcf-concat `cat vcf_input.txt` | gzip -c > chr1-22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered.recode.vcf.gz
```
Convert the concatenated VCF file to plink format
```
plink --vcf chr1-22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered.recode.vcf.gz --double-id --biallelic-only strict list --vcf-half-call missing --make-bed --out chr1-22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered
```
Merge your GWAS data with the 1kg3p3 file
```
plink --bfile file.goodsnps --bmerge chr1-22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered --make-bed --out file.1kg3p3.merged
```
Prune your data similiarly to how was done for the IBD analysis
```
plink --bfile file.1kg3p3.merged --indep-pairwise 50 5 0.2 --out file.1kg3p3.merged

plink --bfile file.1kg3p3.merged --extract file.1kg3p3.merged.prune.in --make-bed --out file.1kg3p3.merged.pruned
```
Now prepare files for running PCA
- Create a file named pca-populations.txt :
```
3
4
5
6
7
```
Each number corresponds to the 5 superpopulations present in the 1000 genomes project phase 3 population, you can extract the information from this [spreadsheet](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx). In this example, we correlate individuals ID in the .fam file with their superpopulation number in this way, 3=EUR, 4=AMR, 5=AFR, 6=EAS, 7=SAS. Modify the phenotype column in the file.1kg3p3.merged.fam file (6th column) so these numbers are reflected in the population. Superpopulation information is [here](http://www.internationalgenome.org/category/population/)
```
...
NA21128 NA21128 0 0 0 7
NA21129 NA21129 0 0 0 7
NA21130 NA21130 0 0 0 7
NA21133 NA21133 0 0 0 7
NA21135 NA21135 0 0 0 7
NA21137 NA21137 0 0 0 7
NA21141 NA21141 0 0 0 7
NA21142 NA21142 0 0 0 7
NA21143 NA21143 0 0 0 7
NA21144 NA21144 0 0 0 7
...
```
Create the .pedind and .pedsnp files from the plink files
```
cp file.1kg3p3.merged.pruned.fam file.1kg3p3.merged.pruned.pedind
cp file.1kg3p3.merged.pruned.bim file.1kg3p3.merged.pruned.pedsnp
```
Run the smartpca program from the EIGENSOFT software. The following command runs the PCA requesting 10 eigenvecrtors and running with 10 cores. It is recommended to export the EIGNESOFT bin folder to your PATH enviroment first.
```
smartpca.perl -i file.1kg3p3.merged.pruned.bed -a file.1kg3p3.merged.pruned.pedsnp -b file.1kg3p3.merged.pruned.pedind -o file.1kg3p3.merged.pruned.pca -p file.1kg3p3.merged.pruned.plot -e file.1kg3p3.merged.pruned.eval -l file.1kg3p3.merged.pruned.log -k 10 -t 10 -w pca-populations.txt
```
Take the file file.1kg3p3.merged.evec and calculate the average and plus/minos 6 standard deviations for each of the 10 eigenvectors, separatedly, and identify outlier individuals. Create the file fail-PCA.txt with the outlying individual IDs.

### Remove all QC-failing individuals from the original plink file
```
cat fail-sexcheck.txt fail-hetrate.txt fail-missingness.txt fail-PCA.txt | sort | uniq | awk '{print $1 "\t" $1}' > fail-ALLQC.txt

plink --bfile file --remove fail-ALLQC.txt --make-bed --out file.AllIndivQC
```


