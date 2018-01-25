#!/bin/bash

Sample_list=path/to/Samplelist.txt
FASTQ_dir=path/to/FASTQ_dir

FastQC=path/to/FastQC_0.11.5/fastqc
CutAdapt=/path/to/cutadapt
Bwa=/path/to/bwa-0.7.12/bwa
Bwa_index=/path/to/bwa-0.7.12/Index
ref_fa=${Bwa_index}/ref.fa
Samtools=/path/to/samtools-1.6/samtools
		
Single() {
	for SamplePath in `cat $Sample_list`;do
		SamplePath=${SamplePath/\r/}
		Sample=${SamplePath/\//_}
		echo ${Sample}
		echo ${SamplePath}
		
		FastQC_input=${FASTQ_dir}/${SamplePath}.fastq.gz
		FastQC_outdir=path/to/FastQC
		$FastQC -f fastq -o ${FastQC_outdir} -d ${FastQC_outdir} -t 10 ${FastQC_input}

		Output=/path/to/Result/CutAdapt
		FiveAdapt=TATTATAAGGCCTGCTGAAA
		ThreeAdapt=TTCAGAATCATTTTGTGGACGAATATGATCCAA # Tumor
		#ThreeAdapt=AGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAA
		
		$CutAdapt -g ${FiveAdapt} -o ${Output}/${Sample}.fastq.gz ${FASTQ_dir}/${SamplePath}
		$CutAdapt -a ThreePrimeCut=${ThreeAdapt} -o ${Output}/${Sample}_{name}.fastq.gz ${Output}/${Sample}.fastq.gz
		
		Bwa_output=/path/to/Result/BWA
		Bwa_input=/path/to/Result/CutAdapt/${Sample}_ThreePrimeCut.fastq.gz
		
		${Bwa} aln ${Bwa_index}/Result ${Bwa_input} -f ${Bwa_output}/${Sample}.sai -t 10
		
		${Bwa} samse ${Bwa_index}/Result ${Bwa_output}/${Sample}.sai ${Bwa_input} | ${Samtools} view -bS -q 1- > ${Bwa_output}/${Sample}.bam 
		${Samtools} sort -o ${Bwa_output}/Sorted_${Sample}.bam ${Bwa_output}/${Sample}.bam
		${Samtools} index ${Bwa_output}/Sorted_${Sample}.bam 
	done
}

Pair() {

	for Sample_pair in `cat $Sample_list`;do
		for i in 1 2;do
			Sample=${Sample_pair}_R${i}_001
			echo $Sample
		
			FastQC_input=${FASTQ_dir}/${Sample}.fastq.gz
			FastQC_outdir=/path/to/Result/FastQC
			$FastQC -f fastq -o ${FastQC_outdir} -d ${FastQC_outdir} -t 10 ${FastQC_input}

			Output=/path/to/Result/CutAdapt
			FiveAdapt=TATTATAAGGCCTGCTGAAA
			ThreeAdapt=TTCAGAATCATTTTGTGGACGAATATGATCCAA # Tumor
			#ThreeAdapt=AGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAA
			$CutAdapt -g ${FiveAdapt} -o ${Output}/${Sample}.fastq.gz ${FASTQ_dir}/${Sample}.fastq.gz
			$CutAdapt -a ThreePrimeCut=${ThreeAdapt} -o ${Output}/trimmed-${Sample}_{name}.fastq.gz ${Output}/${Sample}.fastq.gz

			Bwa_output=/path/to/Result/BWA
			Bwa_input=/path/to/Result/CutAdapt/trimmed-${Sample}_ThreePrimeCut.fastq.gz

			${Bwa} aln ${Bwa_index}/WJK ${Bwa_input} -f ${Bwa_output}/${Sample}.sai -t 10
		
		Sample=${Sample_pair}
		Bwa_input_1=/path/to/Result/CutAdapt/trimmed-${Sample}_R1_001_ThreePrimeCut.fastq.gz
		Bwa_input_2=/path/to/Result/CutAdapt/trimmed-${Sample}_R2_001_ThreePrimeCut.fastq.gz
		
		${Bwa} sampe ${Bwa_index}/Result ${Bwa_output}/${Sample}_R1_001.sai ${Bwa_output}/${Sample}_R2_001.sai ${Bwa_input_1} ${Bwa_input_2}| ${Samtools} view -bS -q 1- > ${Bwa_output}/${Sample}.bam 
		${Samtools} sort -o ${Bwa_output}/Sorted_${Sample}.bam ${Bwa_output}/${Sample}.bam
		${Samtools} index ${Bwa_output}/Sorted_${Sample}.bam 

		done
	done
}

## Select Single(single-end sequencing reads) or Pair(paired-end sequencing reads)
Single
#Pair
