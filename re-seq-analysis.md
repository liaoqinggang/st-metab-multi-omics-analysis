Re-sequencing (re-seq) Data Analysis Pipeline

This pipeline outlines the workflow for aligning and variant calling of re-seq data. It follows standard GATK Best Practices, utilizing Snakemake for automation where applicable.

```
1.we align re-seq data with snakemake and create bam files

1.snake_BWA_mem.py

2.we call the variants with gatk using # only 用一个部分来作为例子

HaplotypeCaller -R DM8.1_genome.fasta -I C001.mk.sorted.bam -O C001.g.vcf --native-pair-hmm-threads 5 -ERC GVCF

3.combine the gvcfs

gatk CombineGVCFs \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V C001.g.vcf \
-V C002.g.vcf \
-V C003.g.vcf \
-V C006.g.vcf \
-V C007.g.vcf \
-V C009.g.vcf \
-V C010.g.vcf \
-V C012.g.vcf \
-V C014.g.vcf \
-V C015.g.vcf \
-V C016.g.vcf \
-V C017.g.vcf \
-V C018.g.vcf \
-V C019.g.vcf \
-V C020.g.vcf \
-V C021.g.vcf \
-V C022.g.vcf \
-V C024.g.vcf \
-V C025.g.vcf \
-V C026.g.vcf \
-V C027.g.vcf \
-V C028.g.vcf \
-V C029.g.vcf \
-V C031.g.vcf \
-V C032.g.vcf \
-V C033.g.vcf \
-V C034.g.vcf \
-V C035.g.vcf \
-V C036.g.vcf \
-V C037.g.vcf \
-O merge_1.vcf

4.combine the paritial vcfs

gatk CombineGVCFs \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V C001.g.vcf \
-V C002.g.vcf \
-V C003.g.vcf \
-V C006.g.vcf \
-V C007.g.vcf \
-V C009.g.vcf \
-V C010.g.vcf \
-V C012.g.vcf \
-V C014.g.vcf \
-V C015.g.vcf \
-V C016.g.vcf \
-V C017.g.vcf \
-V C018.g.vcf \
-V C019.g.vcf \
-V C020.g.vcf \
-V C021.g.vcf \
-V C022.g.vcf \
-V C024.g.vcf \
-V C025.g.vcf \
-V C026.g.vcf \
-V C027.g.vcf \
-V C028.g.vcf \
-V C029.g.vcf \
-V C031.g.vcf \
-V C032.g.vcf \
-V C033.g.vcf \
-V C034.g.vcf \
-V C035.g.vcf \
-V C036.g.vcf \
-V C037.g.vcf \
-O merge_1.vcf

5.convert gvcf and get SNP and INDEL

gatk GenotypeGVCFs \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all.vcf \
-O merge_all_raw.vcf

gatk VariantFiltration \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_raw.vcf \
-O merge_all_filter.vcf \
--filter-expression "QUAL < 30.0" \
--filter-name "LowQual" \

gatk SelectVariants \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_filter.vcf \
-O merge_all_filter_snp.vcf \
--select-type-to-include SNP

gatk SelectVariants \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_filter.vcf \
-O merge_all_filter_indel.vcf \
--select-type-to-include INDEL

6.filter SNP and INDEL

gatk VariantFiltration \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_filter_snp.vcf \
-O merge_all_filter2_snp.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 60.0" \
-filter-name "MQ_filter" -filter "MQ < 40.0" \
-filter-name "SOR_filter" -filter "SOR > 4.0" \
-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk SelectVariants \
--exclude-filtered \
-V merge_all_filter2_snp.vcf \
-O merge_all_filter3_snp.vcf

gatk VariantFiltration \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_filter_indel.vcf \
-O merge_all_filter2_indel.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 200.0" \
-filter-name "SOR_filter" -filter "SOR > 10.0"

gatk SelectVariants \
--exclude-filtered \
-V merge_all_filter2_indel.vcf \
-O merge_all_filter3_indel.vcf

6.Base Quality Score Recalibration # only 用一个bam作为代表

gatk BaseRecalibrator \
 -R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
 -I C001.mk.sorted.bam \
 --known-sites merge_all_filter3_snp.vcf \
 --known-sites merge_all_filter3_indel.vcf \
 -O recal_data.table

gatk ApplyBQSR \
 -R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
 -I C001.mk.sorted.bam \
 -bqsr recal_data.table \
 -O C001_recal_reads.bam \

7.call variants

gatk HaplotypeCaller -R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta -I C001_recal.bam -OC001.g.vcf --native-pair-hmm-threads 5 -ERC GVCF

8.combing gvcf

gatk CombineGVCFs \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V C001.g.vcf \
-V C002.g.vcf \
-V C003.g.vcf \
-V C006.g.vcf \
-V C007.g.vcf \
-V C009.g.vcf \
-V C010.g.vcf \
-V C012.g.vcf \
-V C014.g.vcf \
-V C015.g.vcf \
-V C016.g.vcf \
-V C017.g.vcf \
-V C018.g.vcf \
-V C019.g.vcf \
-V C020.g.vcf \
-V C021.g.vcf \
-V C022.g.vcf \
-V C024.g.vcf \
-V C025.g.vcf \
-V C026.g.vcf \
-V C027.g.vcf \
-V C028.g.vcf \
-V C029.g.vcf \
-V C031.g.vcf \
-V C032.g.vcf \
-V C033.g.vcf \
-V C034.g.vcf \
-V C035.g.vcf \
-V C036.g.vcf \
-V C037.g.vcf \
-O merge_1.vcf

gatk CombineGVCFs \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_1.vcf \
-V merge_2.vcf \
-V merge_3.vcf \
-V merge_4.vcf \
-V merge_5.vcf \
-V merge_6.vcf \
-O merge_all.vcf

9.convert gvcf

gatk GenotypeGVCFs \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all.vcf \
-O merge_all_raw.vcf

10.初步过滤

gatk SelectVariants \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_raw.vcf \
-O merge_all_snp.vcf \
--select-type-to-include SNP

gatk VariantFiltration \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_snp.vcf \
-O merge_filtered_snps_final1.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 60.0" \
-filter-name "MQ_filter" -filter "MQ < 40.0" \
-filter-name "SOR_filter" -filter "SOR > 4.0" \
-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk SelectVariants \
-V merge_filtered_snps_final1.vcf \
--reference /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
--restrict-alleles-to BIALLELIC \
-O merge_filtered_snps_biallelic_final2.vcf

gatk VariantFiltration \
 -V ../merge_filtered_snps_biallelic_final2.vcf \
 --filter-expression "AN < 302 || DP < 850" \
 --filter-name "LowCoverageOrHighMissing" \
 -O merge_filtered_snps_biallelic_final3.vcf
gatk SelectVariants \
 -V merge_filtered_snps_biallelic_final3.vcf \
 --exclude-filtered \
 -O merge_filtered_snps_biallelic_filtered.vcf

gatk VariantFiltration \
 -V merge_filtered_snps_biallelic_filtered.vcf \
 --filter-expression "AF < 0.01 || AF > 0.99" \
 --filter-name "LowMAF" \
 -O merge_filtered_snps_biallelic_final4.vcf
gatk SelectVariants \
 -V merge_filtered_snps_biallelic_final4.vcf \
 --exclude-filtered \
 -O merge_filtered_snps_biallelic_filtered2.vcf

gatk VariantFiltration \
-R /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/01DNA-reseq/03_bamfiles/03_vcf_8.1/01_part1/DM8.1_genome.fasta \
-V merge_all_indel.vcf \
-O merge_filtered_indels_final.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 200.0" \
-filter-name "SOR_filter" -filter "SOR > 10.0"

bcftools filter -i 'TYPE!="indel" || (ABS(STRLEN(REF)-STRLEN(ALT)) <= 5)' ../merge_filtered_indels_final.vcf -Ov -o merge_filtered_5bp_indels_final.vcf
```

xw
