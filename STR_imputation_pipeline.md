## STR imputation pipeline
This pipeline takes directions from [here](http://gymreklab.com/2018/03/05/snpstr_imputation.html)
- Download STR reference panel for each chromosome (1-22) from [here](http://gymreklab.com/2018/03/05/snpstr_imputation.html)
- Harmonize your target GWAS data with the reference panel
- Perform imputation with Beagle
- QC on imputed genotypes

### Download STR reference panel for each chromosome
```
for i in {1..22}
do
wget https://s3.amazonaws.com/snp-str-imputation/1000genomes/1kg.snp.str.chr${i}.vcf.gz
wget https://s3.amazonaws.com/snp-str-imputation/1000genomes/1kg.snp.str.chr${i}.vcf.gz.tbi
done
```
### Harmonize your target GWAS data with the reference panel
**Make sure your GWAS data is in VCF format
```
for i in {1..22}
do
  java -jar -Xmx110g conform-gt.24May16.cee.jar \
  gt=GWAS.vcf.gz\
  ref=1kg.snp.str.chr${i}.vcf.gz \
  chrom=${i} \ 
  match=POS \ 
  out=snp.chr${i}.consistent
done
```
### Perform imputation with Beagle
```
for i in {1..22}
do
  java -jar -Xmx110g beagle.28Sep18.793.jar \
  gt=snp.chr${i}.consistent.vcf.gz \
  ref=1kg.snp.str.chr${i}.vcf.gz \
  nthreads=20 \
  window=10    \
  overlap=1.0 \
  out=snp.str.chr${i}.vcf.gz
done
```
To filter the resulting VCF is kind of tricky, it has a DR2 measurement for how well the imputation went for each site, but there is no program use this variable consistently for SNPs, biallelic STRs and multiallelic STRs. Here I use the default of DR2 > 0.3 as filtering criteria, using a custom method. Variants with DR2 < 0.3 are set as missing.
First, I separate SNPs from STRs in different VCF files

Extract SNPs and STRs in different files
```
for i in {1..22}
do
vcftools --gzvcf snp.str.chr${i}.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > onlySNPs.chr${i}.vcf.gz

vcftools --gzvcf snp.str.chr${i}.vcf.gz --keep-only-indels  --recode --recode-INFO-all --stdout | gzip -c > onlySTRs.chr${i}.vcf.gz
done
```
Manually change the VCF headers on the STRs files to reflect this:

old
```
##INFO=<ID=DR2,Number=1
```
new
```
##INFO=<ID=DR2,Number=A
```
```
for i in {1..22}
do
zcat onlySTRs.chr${i}.vcf.gz | sed 's/##INFO=<ID=DR2,Number=1/##INFO=<ID=DR2,Number=A/g' | gzip -c > DR2A.onlySTRs.chr${i}.vcf.gz
done
```
Perform multiallelic splitting on STRs with the VT software
```
for i in {1..22}
do
vt decompose -s DR2A.onlySTRs.chr${i}.vcf.gz -o DR2A.onlySTRs.chr${i}.split.vcf.gz
done
```
Change the STR IDs to STR_123_# where # is the variant number from first to last in the VCF file. With this we can get unique IDs for each variants, because when multipliting occurs, if 1 variant has many observed alleles (ALT), you will that many variants with the same ID and that we would like to avoid for the annotation steps later.
```
for i in {1..22}
do
zcat DR2A.onlySTRs.chr${i}.split.vcf.gz | grep -F "#" > DR2A.onlySTRs.chr${i}.split.uniqueID.vcf

zcat splitted.DR2A.STR_only.noQC.chr${i}.vcf.gz | grep -F -v "#" | awk 'BEGIN {OFS="\t"} {$3=$3"_"NR} {print}' >> DR2A.onlySTRs.chr${i}.split.uniqueID.vcf
done
```
Filter STRs by DR2 > 0.3
```
for i in {1..22}
do
grep -v "#" DR2A.onlySTRs.chr${i}.split.uniqueID.vcf | cut -f3,8 | sed 's/DR2=/\t/g' | sed 's/;/\t/g' | awk '$2 > 0.3{print $1}' > chr${i}_newsplit_Dr2_03.temp 

vcftools --vcf DR2A.onlySTRs.chr${i}.split.uniqueID.vcf --snps chr${i}_newsplit_Dr2_03.temp --recode --recode-INFO-all --stdout | gzip -c > DR2A.onlySTRs.chr${i}.split.uniqueID.DR2_03.vcf.gz
done
```

Filter SNPs by DR2 > 0.3
```
for i in {1..22}
do
grep -v "#" onlySNPs.chr${i}.vcf.gz | cut -f3,8 | sed 's/DR2=/\t/g' | sed 's/;/\t/g' | awk '$2 > 0.3{print $1}' > chr${i}_newsplit_Dr2_03.SNPS.temp 

vcftools --vcf onlySNPs.chr${i}.vcf.gz --snps chr${i}_newsplit_Dr2_03.SNPS.temp --recode --recode-INFO-all --stdout | gzip -c > onlySNPs.chr${i}.DR2_03.vcf.gz
done
```

Then, the next steps is to transform each file to Plink format and perform the desired associations with the target phenotype

### Bonus
Before splitting of the VCF files containing STRs only, we can divide the file into biallelic STRs and multiallelic STRs
```
for i in {1..22}
do
bcftools view --max-alleles 2 onlySTRs.chr${i}.vcf.gz -Oz -o onlySTRs.chr${i}.biallelic.vcf.gz
done
```
Then you can follow the filtering criteria for the biallelic ones
