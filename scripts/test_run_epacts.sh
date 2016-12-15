#!/bin/bash
set -e
mkdir --p out/
DIR=`dirname $0`
VCF=${DIR}/../share/EPACTS/1000G_exome_chr20_example_softFiltered.calls.vcf.gz
PED=${DIR}/../share/EPACTS/1000G_dummy_pheno.ped
OUT=out/test
GRP=out/1000G_exome_chr20_example_softFiltered.calls.sites.anno.grp 

mkdir --p out/
echo "Performing EPACTS Binary Single Variant Score Test.."
${DIR}/epacts single --vcf ${VCF} --ped ${PED} --min-maf 0.001 --chr 20 --pheno DISEASE --cov AGE --cov SEX --test b.score --out ${OUT}.single.b.score --run 5 --anno --topzoom 1

echo "Creating Empirical Kinship matrix..."
${DIR}/epacts make-kin --vcf ${VCF} --ped ${PED} --min-maf 0.01 --chr 20 --out ${OUT}.single.q.emmax.kinf --run 5

echo "Running EMMAX single variant test..."
${DIR}/epacts single --vcf ${VCF} --ped ${PED} --min-maf 0.001 --chr 20 --pheno QT --cov AGE --cov SEX --test q.emmax --out ${OUT}.single.q.emmax --kinf ${OUT}.single.q.emmax.kinf --run 5 --anno

echo "Taking the site list from the VCF file";
zcat ${VCF} | cut -f 1-8 | ${DIR}/bgzip -c > out/1000G_exome_chr20_example_softFiltered.calls.sites.vcf.gz

echo "Annotate the VCF file";
${DIR}/epacts anno --in out/1000G_exome_chr20_example_softFiltered.calls.sites.vcf.gz --out out/1000G_exome_chr20_example_softFiltered.calls.sites.anno.vcf.gz

echo "Create the input group file for burden test";
${DIR}/epacts make-group --vcf out/1000G_exome_chr20_example_softFiltered.calls.sites.anno.vcf.gz --out ${GRP} -format epacts -nonsyn

echo "Performing SKAT-O test";
${DIR}/epacts group --vcf ${VCF} --ped ${PED} --max-maf 0.05 --groupf ${GRP} --pheno QT --cov AGE --cov SEX --test skat --out ${OUT}.gene.skat --skat-o --unit 100 --run 5
