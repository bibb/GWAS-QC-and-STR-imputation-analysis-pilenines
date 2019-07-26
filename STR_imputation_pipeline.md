## STR imputation pipeline
- Download STR reference panel for each chromosome (1-22) from http://gymreklab.com/2018/03/05/snpstr_imputation.html
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
Filter the resulting VCF is kind of tricky, it has a DR2 measurement for how well the imputation went for each site, but there is no program use this variable consistently for SNPs, biallelic STRs and multiallelic STRs. Here I use the default of DR2 > 0.3 as filtering criteria, using a custom method. Variants with DR2 < 0.3 are set as missing.
First, I separate SNPs from STRs in different VCF files

Extract SNPs and STRs in different files
```
for i in {1..22}
do
vcftools --gzvcf snp.str.chr${i}.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > onlySNPs.chr${i}.vcf.gz

vcftools --gzvcf snp.str.chr${i}.vcf.gz --keep-only-indels  --recode --recode-INFO-all --stdout | gzip -c > onlySTRs.chr${i}.vcf.gz
done
```
Perform multiallelic splitting
