
files, = glob_wildcards("analyses/trimmed/{sample}.trimmed.fastq.gz")

include: "bowtie2.smk"
include: "prepare_databases.smk"

genomeloc=""
location=""


hisatgenome= genomeloc + "hg38"
humangenome= genomeloc + "hg38.fa"


known_sites1=location + "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
known_sites2=location + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_sites3=location + "hapmap_3.3.hg38.vcf.gz"
known_sites4=location + "Homo_sapiens_assembly38.known_indels.vcf.gz"
known_sites5=location + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"

#with open("chromosomes.txt") as chr:
#	chromosomes=chr.readlines()

rule all_variant:
	input:
		expand("results/haplo_vcf/{sample}.g.vcf.gz",sample=files),
		"results/main/output.vcf.gz",
		"results/main/output_snp.vcf.gz",
		"results/main/output_indel.vcf.gz",
		"results/main/output_snp.filtered.vcf.gz"


rule hisat_mapping_reads:
	input:
		"analyses/trimmed/{sample}.trimmed.fastq.gz"
	
	output:
		"analyses/hisat2/{sample}.sorted.bam"

   
	threads: 20
	
	shell:
		"""
		
        hisat2 -k 1 -p {threads} -x {hisatgenome} -U {input} | samtools view -bS - | samtools sort - -o {output}
        samtools index {output}

				
		"""


rule mark_duplicates:
	input:
		"analyses/hisat2/{sample}.sorted.bam"

	output:
		"analyses/picard/{sample}.mark_dup.bam"

	threads: 5

	resources:
		mem_mb=5000

	shell:
		"""
		java -Xmx4g -jar picard.jar MarkDuplicates I={input} O=analyses/picard/{wildcards.sample}.mark_dup.bam M=analyses/picard/{wildcards.sample}.marked_dup_metrics.txt

		"""

	
rule split_cigar_reads:
	input:
		"analyses/picard/{sample}.mark_dup.bam"

	output:
		temp("analyses/gatk/{sample}.tmp.bam")

	
	threads: 2

	resources:
		mem_mb=5000
	shell:
		"""

		gatk SplitNCigarReads -R {humangenome} -I {input} -O {output}
		


		"""

rule addreadgroups:
	input:
		"analyses/gatk/{sample}.tmp.bam"
	output:
		"analyses/gatk/{sample}.bam"
	threads: 2

	resources:
		mem_mb=5000
	shell:
		"""
		#java -jar picard.jar AddOrReplaceReadGroups I={input}       O={output} RGID=4       RGLB=lib1       RGPL=illumina       RGPU=unit1       RGSM=20
		java -jar picard.jar AddOrReplaceReadGroups I={input}       O={output} RGID=4       RGLB=lib1       RGPL=illumina       RGPU=unit1 RGSM={wildcards.sample}
		"""


rule recalibrate:
	input:
		"analyses/gatk/{sample}.bam",
		known_sites1,
		known_sites2,
		known_sites3,
		known_sites4,
		known_sites5,
	output:
		"analyses/gatk/recal/{sample}.recal_data.table"

	threads: 2

	resources:
		mem_mb=5000
	shell:
		"""
		gatk BaseRecalibrator -I {input[0]} -R {humangenome} --known-sites {known_sites1} 	--known-sites {known_sites2} --known-sites {known_sites3} --known-sites {known_sites4} --known-sites {known_sites5} -O {output}


		"""

rule applycalibrate:
	input:
		"analyses/gatk/{sample}.bam",
		"analyses/gatk/recal/{sample}.recal_data.table"
	output:
		"analyses/filtered/{sample}.bam"

	threads: 2

	resources:
		mem_mb=5000
	shell:
		"""
		gatk ApplyBQSR \
   -R {humangenome} \
   -I {input[0]} \
   --bqsr-recal-file {input[1]} \
   -O {output}

   samtools index {output}
		"""


rule haplocaller:
	input:
		"analyses/filtered/{sample}.bam"
	output:
		"results/haplo_vcf/{sample}.g.vcf.gz"

	threads: 3

	resources:
		mem_mb=21000

	shell:
		"""
		
		gatk --java-options "-Xmx20g" HaplotypeCaller  --native-pair-hmm-threads {threads}    -R {humangenome}    -I {input}    -O {output}  -ERC GVCF

		"""

rule genomicsdbi:
	input:
		expand("results/haplo_vcf/{sample}.g.vcf.gz",sample=files)
	output:
		"analyses/dgbi/vcfheader.vcf"

	threads: 15

	resources:
		mem_mb=21000
	shell:
		"""
		rm -r analyses/dgbi/
		paramaters=$(for i in {input}; do printf " -V %s " $i; done)
		gatk --java-options "-Xmx20g -Xms20g" GenomicsDBImport --reader-threads {threads}  $paramaters --genomicsdb-workspace-path analyses/dgbi  -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr1 -L chr20 -L chr21 -L chr22 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chrM -L chrX -L chrY
		
		"""


rule mergevcfs:
	input:
		"analyses/dgbi/vcfheader.vcf"
	output:
		"results/main/output.vcf.gz"
	resources:
		mem_mb=21000
	shell:
		"""
		gatk --java-options "-Xmx20g -Xms20g" GenotypeGVCFs -R {humangenome} -V gendb://analyses/dgbi/ -O {output}

		"""

rule snpsandindels:
	input:
		"results/main/output.vcf.gz"
	output:
		snp="results/main/output_snp.vcf.gz",
		indel="results/main/output_indel.vcf.gz"
	resources:
		mem_mb=21000
	shell:
		"""
		 gatk --java-options "-Xmx20g -Xms20g" SelectVariants -R {humangenome} -V {input} --select-type-to-include SNP -O {output.snp}

	 	 gatk --java-options "-Xmx20g -Xms20g" SelectVariants -R {humangenome} -V {input} --select-type-to-include INDEL -O {output.indel}

		"""
	

		
rule hardfiltersnp:
	input:
		"results/main/output_snp.vcf.gz"
	output:
		"results/main/output_snp.filtered.vcf.gz"
	resources:
		mem_mb=21000
	shell:
		"""
		gatk --java-options "-Xmx20g -Xms20g" VariantFiltration -V {input} 		-filter "QD < 2.0" --filter-name "QD2"		-filter "QUAL < 30.0" --filter-name "QUAL30"		-filter "SOR > 3.0" --filter-name "SOR3" 		-filter "FS > 60.0" --filter-name "FS60"		-filter "MQ < 40.0" --filter-name "MQ40" 		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" 		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 		-O {output}

		"""

