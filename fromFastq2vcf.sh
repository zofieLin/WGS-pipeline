
#!/bin/bash
cd $PBS_O_WORKDIR

dir="/path/project_dir"
gatk_resources="/path/GATK/bundle/hg19"
refseq="/path/2.8/hg19/ucsc.hg19.fasta"
fastq_dir="/path/primary_seq" 

#### Step 1, BWA alignment #############################################################################################################
if [ ! -d "$dir/sh.e.o/01_bwamem" ]; then mkdir -p $dir/sh.e.o/01_bwamem; fi

for sample in `cat $dir/sample.list`
do
echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load bwa/0.7.17 samtools/1.11
temp_dir=\"$dir/$sample/tmp\"
log=\"$dir/$sample/${sample}.log\"
if [ ! -d "$dir/$sample/alignment" ]; then mkdir -p $dir/$sample/alignment; fi
if [ ! -d "$dir/$sample/log" ]; then mkdir -p $dir/$sample/log; fi

## Detection of RG, for further VerifyBamID
if [ ! -s "$dir/$sample/alignment/${sample}.RG " ]
then
      echo \"Detecting RGs in fastq file started on \`date\`\" >> \$log
      zcat $fastq_dir/${sample}_1.fastq.gz | awk -F \":\" '/^@/ { print \$3\".\"\$4 }' | sort | uniq  >$dir/$sample/alignment/${sample}.RG
      echo \"Detecting RGs in fastq file ended on \`date\`\" >> \$log
fi

rgheader=\"\"
for rg in \`cat $dir/$sample/alignment/${sample}.RG\`; do
    if [[ -z \$rgheader ]]; then sep=\"\"; else sep=\"\n\"; fi
        rgheader=\${rgheader}\${sep}\"@RG\tID:\$rg\tPL:ILLUMINA\tSM:$sample\"
done

## Alignment using bwamem
echo \"Mapping reads using bwa-mem started on \`date\`\" >> \$log
    bwa mem -M -t 10 -R \"@RG\\tID:BA\\tPL:ILLUMINA\\tSM:$sample\" $refseq $fastq_dir/${sample}_1.fastq.gz $fastq_dir/${sample}_2.fastq.gz | awk -v rgheader=\$rgheader -F \":\" '/^@RG/ { \$0=rgheader } !/^@/ { rg=\$3\".\"\$4; gsub(\"RG:Z:BA\",\"RG:Z:\"rg,\$0) }{ print }' | samtools view -bS -o $dir/$sample/alignment/${sample}.bam -  && \\
echo \"Mapping reads using bwa-mem ended on \`date\`\" >> \$log

## Sort the bam file
if [ ! -d "\$temp_dir" ]; then mkdir \$temp_dir; fi
echo \"Sorting compressed bam file started on \`date\`\" >> \$log
    /software/sambamba/bin/sambamba sort -t 10 -m 30G --tmpdir=\$temp_dir -o $dir/$sample/alignment/${sample}.sorted.bam $dir/$sample/alignment/${sample}.bam  && \\
samtools index $dir/$sample/alignment/${sample}.sorted.bam
echo \"Sorting compressed bam file ended on \`date\`\" >> \$log"    >$dir/sh.e.o/01_bwamem/step1_${sample}_bwa.sh
done
#### Step 1, BWA alignment #############################################################################################################


#### Step 2, Mark duplicated reads #####################################################################################################
if [ ! -d "$dir/sh.e.o/02_dedup" ]; then mkdir -p $dir/sh.e.o/02_dedup; fi

for sample in `cat $dir/sample.list`
do
echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load java/8.0_161 Picard/2.18.9
temp_dir=\"$dir/$sample/tmp\"
log=\"$dir/$sample/${sample}.log\"

if [ -e $dir/\$sample/alignment/\${sample}.bam ]; then rm $dir/$sample/alignment/${sample}.bam; fi
echo \"Mark duplicates started on \`date\`\" >> \$log
    java -XX:ParallelGCThreads=8 -Xmx30g -jar /software/Picard/2.18.9/picard.jar MarkDuplicates INPUT=$dir/$sample/alignment/${sample}.sorted.bam OUTPUT=$dir/$sample/alignment/${sample}.sorted.dedup.bam METRICS_FILE=$dir/$sample/alignment/lane.dedup.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=\$temp_dir
echo \"Mark duplicates ended on \`date\`\" >> \$log"   >$dir/sh.e.o/02_dedup/step2_${sample}_dedup.sh
done
#### Step 2, Mark duplicated reads #####################################################################################################


#### Step 3, Realignment ###############################################################################################################
## this step have beed discarded in GATK, can consider to remove this step as it take long time ##
if [ ! -d "$dir/sh.e.o/03_realign" ]; then mkdir -p $dir/sh.e.o/03_realign; fi

for sample in `cat $dir/sample.list`
do
echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load java/8.0_161 Picard/2.18.9 GenomeAnalysisTK/3.8.1.0
#Local realignment
#    - (a) generate regional target for realingmentlsls
#    - (b) realignment performed on target regions

temp_dir=\"$dir/$sample/tmp\"
log=\"$dir/$sample/${sample}.log\"
rm $dir/$sample/alignment/${sample}.sorted.bam $dir/$sample/alignment/${sample}.sorted.bam.bai

echo \"Realignment a started on \`date\`\" >> \$log
    java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 8 -R $refseq -I $dir/$sample/alignment/${sample}.sorted.dedup.bam -o $dir/$sample/alignment/${sample}.sorted.dedup.intervals -known $gatk_resources/1000G_phase1.indels.hg19.sites.vcf -known $gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -S LENIENT && \\
echo \"Realignment a ended on \`date\`\" >> \$log

echo \"Realignment b started on \`date\`\" >> \$log
    java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T IndelRealigner -filterNoBases -R $refseq -targetIntervals $dir/$sample/alignment/${sample}.sorted.dedup.intervals -I $dir/$sample/alignment/${sample}.sorted.dedup.bam -o $dir/$sample/alignment/${sample}.sorted.dedup.realn.bam -S LENIENT
echo \"Realignment b ended on \`date\`\" >> \$log"     >$dir/sh.e.o/03_realign/step3_${sample}_realign.sh
done
#### Step 3, Realignment ###############################################################################################################


#### Step 4, Base Quality Score Recalibration ##########################################################################################
for sample in `cat $dir/sample.list`
do
if [ ! -d "$dir/sh.e.o/04_bqsr" ]; then mkdir -p $dir/sh.e.o/04_bqsr; fi
echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load java/8.0_161 Picard/2.18.9 GenomeAnalysisTK/3.8.1.0
temp_dir=\"$dir/$sample/tmp\"
log=\"$dir/$sample/${sample}.log\"

rm $dir/$sample/alignment/${sample}.sorted.dedup.bam $dir/$sample/alignment/${sample}.sorted.dedup.bai
echo \"BQSR started on \`date\`\" >> \$log
    java -XX:ParallelGCThreads=8 -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refseq -knownSites $gatk_resources/1000G_phase1.indels.hg19.sites.vcf -knownSites $gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites $gatk_resources/dbsnp_138.hg19.vcf -I $dir/$sample/alignment/${sample}.sorted.dedup.realn.bam -o $dir/$sample/alignment/${sample}.sorted.dedup.realn.recal.grp -nct 8 -nt 1 -S LENIENT && \\
    java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T PrintReads -R $refseq -I $dir/$sample/alignment/${sample}.sorted.dedup.realn.bam -BQSR $dir/$sample/alignment/${sample}.sorted.dedup.realn.recal.grp -o $dir/$sample/alignment/${sample}.sorted.dedup.realn.recal.bam -nct 8 -nt 1 -S LENIENT && \\
echo \"BQSR ended on \`date\`\" >> \$log"     >$dir/sh.e.o/04_bqsr/step4_${sample}_bqsr.sh
done
#### Step 4, Base Quality Score Recalibration ##########################################################################################


#### Step 5, Variant calling by chromosome ##############################################################################################
if [ ! -d "$dir/sh.e.o/05_byChr" ]; then mkdir -p $dir/sh.e.o/05_byChr; fi
for chr in `seq 22` X Y
do
    if [ ! -d "$dir/sh.e.o/05_byChr/chr$chr" ]; then mkdir -p $dir/sh.e.o/05_byChr/chr$chr; fi
    for sample in `cat $dir/sample.list`
    do
echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load java/8.0_161 GenomeAnalysisTK/3.8.1.0
temp_dir=\"$dir/$sample/tmp\"
log=\"$dir/$sample/${sample}.log\"
## Generate vcf file using HaplotypeCaller and GenotypeGVCFs
if [ ! -d "$dir/$sample/variantCalling/byChr/chr$chr" ]; then mkdir -p $dir/$sample/variantCalling/byChr/chr$chr; fi
echo \"HC GVCF calling for chr$chr started on \`date\`\" >> \$log
    java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T HaplotypeCaller -R $refseq -I $dir/$sample/alignment/${sample}.sorted.dedup.realn.recal.bam -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --dbsnp $gatk_resources/dbsnp_138.hg19.vcf -A StrandOddsRatio -A Coverage -A QualByDepth -A FisherStrand -A MappingQualityRankSumTest -A ReadPosRankSumTest -A RMSMappingQuality -o $dir/$sample/variantCalling/byChr/chr$chr/chr${chr}.gvcf.vcf.gz -L chr$chr && \\   
echo \"HC GVCF calling for chr$chr ended on \`date\`\" >> \$log

echo \"Genotype GVCF calling of chr$chr started on \`date\`\" >> \$log
    java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $refseq --variant $dir/$sample/variantCalling/byChr/chr$chr/chr${chr}.gvcf.vcf.gz -A StrandOddsRatio -A Coverage -A QualByDepth -A FisherStrand -A MappingQualityRankSumTest -A ReadPosRankSumTest -A RMSMappingQuality -o $dir/$sample/variantCalling/byChr/chr$chr/chr${chr}.gatkHC.vcf.gz --dbsnp $gatk_resources/dbsnp_138.hg19.vcf -stand_call_conf 30.0 -newQual -L chr$chr && \\
echo \"Genotype GVCF calling of chr$chr ended on \`date\`\" >> \$log"     >$dir/sh.e.o/05_byChr/chr$chr/step5_${sample}_chr${chr}_vcf.sh
    done
done
#### Step 5, Variant calling by chromosome ##############################################################################################


#### Step 6, Generating individual gVCF and VCF ######################################################################################### 
for sample in `cat $dir/sample.list`
do
      if [ ! -d "$dir/sh.e.o/06_catgvcf" ]; then mkdir -p $dir/sh.e.o/06_catgvcf; fi

echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load java/8.0_161 GenomeAnalysisTK/3.8.1.0 vcftools/0.1.15

temp_dir=\"$dir/$sample/tmp\"
log=\"$dir/$sample/${sample}.log\"

gvcf_list=\"\"
vcf_list=\"\"
for chr in \`seq 22\` X Y
do
      gvcf_list=\"\$gvcf_list -V $dir/$sample/variantCalling/byChr/chr\$chr/chr\${chr}.gvcf.vcf.gz\"
      vcf_list=\"\$vcf_list $dir/$sample/variantCalling/byChr/chr\$chr/chr\${chr}.gatkHC.vcf.gz\"
done

## Combined all byCHR gvcf
echo \"CatVariants chr gvcf started on \`date\`\" >> \$log
      java -Xmx30g -cp /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants --assumeSorted -R $refseq \$gvcf_list -out $dir/$sample/variantCalling/${sample}.gvcf.vcf.gz && \\
      tabix -f -p vcf $dir/$sample/variantCalling/${sample}.gvcf.vcf.gz
echo \"CatVariants chr gvcf ended on \`date\`\" >> \$log

## Combined all byCHR gatkHC.vcf
echo \"Combined chr gatkHC.vcf started on \`date\`\" >> \$log
      perl /software/vcftools/0.1.15/bin/vcf-concat \$vcf_list | bgzip -c > $dir/$sample/variantCalling/${sample}.gatkHC.vcf.gz && \\
      tabix -f -p vcf $dir/$sample/variantCalling/${sample}.gatkHC.vcf.gz
echo \"Combined chr gatkHC.vcf ended on \`date\`\" >> \$log"      >$dir/sh.e.o/06_catgvcf/step6_${sample}_catgvcf.sh
done
#### Step 6, Generating individual gVCF and VCF #########################################################################################


#### Step 7, QC on bam and vcf ##########################################################################################################
for sample in `cat $dir/sample.list`
do
        if [ ! -d "$dir/$sample/qc_stat" ]; then mkdir -p $dir/$sample/qc_stat; fi
        if [ ! -d "$dir/sh.e.o/07_qc" ]; then mkdir -p $dir/sh.e.o/07_qc; fi

echo "#!/bin/bash
cd \$PBS_O_WORKDIR

module load java/8.0_161 GenomeAnalysisTK/3.8.1.0 Picard/2.18.9
temp_dir="$dir/$sample/tmp"
log="$dir/$sample/${sample}.log"

# QC1. Variant evaluation
echo \"Variant evaluation started on \`date\`\" >> \$log
      java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar -T VariantEval -R $refseq --eval $dir/$sample/variantCalling/${sample}.gatkHC.vcf.gz --dbsnp $gatk_resources/dbsnp_138.hg19.excluding_sites_after_129.vcf -noST -ST Sample -ST Novelty-ST Filter -noEV -EV CountVariants -noEV -EV TiTvVariantEvaluator -noEV -EV IndelSummary -noEV -EV MultiallelicSummary -o $dir/$sample/qc_stat/${sample}.varEval.grp -nt 2 && \\
echo \"Variant evaluation ended on \`date\`\" >> \$log

# QC2. VerifyBamID for detecting contamination from bam files
echo \"VerifyBamID started on \`date\`\" >> \$log
      /software/VerifyBamID/1.1.3/verifyBamID --vcf /software/verifyBamID/Omni25_genotypes_1525_samples_v2.hg19.PASS.ALL.sites.vcf.gz --bam $dir/$sample/alignment/${sample}.sorted.dedup.realn.recal.bam --maxDepth 1000 --precise --out $dir/$sample/qc_stat/${sample}.verifyBamID && \\
echo \"VerifyBamID ended on \`date\`\" >> \$log

# QC3. Whole genome metrics - including DP
if [ ! -e "$dir/$sample/qc_stat/${sample}.wgsMetric.txt" ]
then
echo \"Picard CollectWgsMetrics started on `date`\" >> \$log
      java -Xmx30g -Djava.io.tmpdir=\$temp_dir -jar /software/Picard/2.18.9/picard.jar CollectWgsMetrics INPU
T=$dir/$sample/alignment/${sample}.sorted.dedup.realn.recal.bam OUTPUT=$dir/$sample/qc_stat/${sample}.wgsMetric.txt REFERENCE_SEQUENCE=$refseq VALIDATION_STRINGENCY=LENIENT
echo \"Picard CollectWgsMetrics ended on \`date\`\" >> \$log"     >$dir/sh.e.o/07_qc/step7_${sample}_qc.sh
done
#### Step 7, QC on bam and vcf ##########################################################################################################
