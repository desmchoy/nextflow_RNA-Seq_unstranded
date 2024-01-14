#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.logname = "logname"
params.jobID = "jobID"
params.sample = "sample"
params.read1 = "read1.fastq"
params.read2 = "read2.fastq"
params.transcriptome = "transcriptome"
params.GTF = "GTF"
params.fiveprime = "fiveprime"
params.threeprime = "threeprime"
params.bowtie2_ref = "rRNA_ref"
params.RSEM_ref = "RSEM_ref"

params.FastQC_outdir = "$baseDir/results_FastQC/sample"
params.cutadapt_outdir = "$baseDir/inter_cutadapt/sample"
params.bowtie2_outdir = "$baseDir/inter_Bowtie2_rRNA/sample"
params.STAR_outdir = "$baseDir/results_STAR/sample"
params.RSEM_outdir = "$baseDir/results_RSEM/sample"


logname_ch = Channel.value(params.logname)
jobID_ch = Channel.value(params.jobID)
read1_ch = Channel.fromPath(params.read1, checkIfExists: true)
read2_ch = Channel.fromPath(params.read2, checkIfExists: true)
transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)
GTF_ch = Channel.fromPath(params.GTF, checkIfExists: true)
threeprime_ch = Channel.value(params.threeprime)
fiveprime_ch = Channel.value(params.fiveprime)


log.info """\


        RNA-Seq (unstranded library)    
        ===================================
        User         : ${params.logname}
        HPC Job ID   : ${params.jobID}
        transcriptome: ${params.transcriptome}
        GTF          : ${params.GTF}
        Read 1       : ${params.read1}
        Read 2       : ${params.read2}
        5' Adapter   : ${params.fiveprime}
        3' Adapter   : ${params.threeprime}
        """
        .stripIndent()



process FASTQC {

    tag "FASTQC on $sample"
    cpus 1
    debug true

    input:
    val sample

    script:
    """
    mkdir ${baseDir}/results_FastQC/${sample}
    for file in ${baseDir}/RAW/${sample}_*.fastq.gz ; do
        fastqc --noextract --outdir ${baseDir}/results_FastQC/${sample} \${file} > ${baseDir}/results_FastQC/${sample}/\$(basename \${file}).fastqc.log 2> ${baseDir}/results_FastQC/${sample}/\$(basename \${file}).fastqc.error
    done
    """

}



process SEQTK {

    tag "SEQTK on $sample"
    cpus 1
    debug true
    
    input:
    path read1
    path read2

    output:
    path "${sample}_1.temp.fastq", emit: seqtk_1 
    path "${sample}_2.temp.fastq", emit: seqtk_2
    
    script:
    """
    seqtk seq -q20 -n N $read1 > ${sample}_1.temp.fastq
    seqtk seq -q20 -n N $read2 > ${sample}_2.temp.fastq
    """

}



process CUTADAPT {

    tag "CUTADAPT on $sample"
    publishDir params.cutadapt_outdir, mode: 'move'
    cpus 1
    debug true

    input:
    val threeprime
    val fiveprime
    path seqtk_1
    path seqtk_2

    output:
    path "${sample}.trimmed.1.fastq", emit: cutadapt_1
    path "${sample}.trimmed.2.fastq", emit: cutadapt_2
    path "${sample}.cutadapt.log"
    path "${sample}.cutadapt.error"

    script:
    """
    cutadapt -a ${threeprime} -A ${fiveprime} -O 6 -e 0.1 -m 40 --max-n 10 -o ${sample}.trimmed.1.fastq -p ${sample}.trimmed.2.fastq ${seqtk_1} ${seqtk_2} > ${sample}.cutadapt.log 2> ${sample}.cutadapt.error
    """

}



process BOWTIE2 {

    tag "BOWTIE2 on $sample"
    publishDir params.bowtie2_outdir, mode: 'move'
    cpus 4
    debug true

    input:
    path cutadapt_1
    path cutadapt_2

    output:
    path "*", emit: bowtie2_out

    script:
    """
    bowtie2 -p ${task.cpus} -x /SAN/colcc/cellranger_references/rRNA/human/bowtie2-2.4.1/human_rRNA -1 ${baseDir}/inter_cutadapt/${sample}/${sample}.trimmed.1.fastq -2 ${baseDir}/inter_cutadapt/${sample}/${sample}.trimmed.2.fastq -S ${sample}.rRNA.mapped.sam --un-conc ${sample}.rRNA.filtered.%.fastq > ${sample}.bowtie2.log 2> ${sample}.bowtie2.error
    """

}



process STAR {

    tag "STAR on $sample"
    publishDir params.STAR_outdir, mode: 'move'
    cpus 6
    debug true

    input:
    path transcriptome
    path GTF
    path bowtie2_out

    output:
    path "*", emit: STAR_out

    script:
    """
    STAR --runThreadN ${task.cpus} --genomeDir ${transcriptome} --readFilesIn ${baseDir}/inter_Bowtie2_rRNA/${sample}/${sample}.rRNA.filtered.1.fastq ${baseDir}/inter_Bowtie2_rRNA/${sample}/${sample}.rRNA.filtered.2.fastq --sjdbGTFfile ${GTF} --outFileNamePrefix ${sample}. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic > ${sample}.STAR.log 2> ${sample}.STAR.error

    samtools flagstat ${sample}.Aligned.sortedByCoord.out.bam > ${sample}.Aligned.sortedByCoord.out.flagstat
    samtools index ${sample}.Aligned.sortedByCoord.out.bam
    samtools view ${sample}.Aligned.sortedByCoord.out.bam > ${sample}.Aligned.sortedByCoord.out.sam

    cut -f9 ${sample}.Aligned.sortedByCoord.out.sam | awk '\$1>0' > ${sample}_STAR.insert_size

    rm ${sample}.Aligned.sortedByCoord.out.sam
    """

}



process RSEM {

    tag "RSEM on $sample"
    publishDir params.RSEM_outdir, mode: 'move'
    cpus 6
    debug true

    input:
    val logname
    val jobID
    path STAR_out

    output:
    path "*", emit: RSEM_out

    script:
    """
    rsem-calculate-expression -p ${task.cpus} --temporary-folder /scratch0/${logname}/${jobID} --strandedness none --paired-end --bam --estimate-rspd --calc-ci ${SGE_O_WORKDIR}/results_STAR/${sample}/${sample}.Aligned.toTranscriptome.out.bam /SAN/colcc/cellranger_references/GRCh38/RSEM-1.3.3/hg38 ${sample} > ${sample}.RSEM.log 2> ${sample}.RSEM.error
    """

}


workflow {

    FASTQC(sample)
    SEQTK(read1_ch, read2_ch)
    CUTADAPT(threeprime_ch, fiveprime_ch, SEQTK.out.seqtk_1.collect(), SEQTK.out.seqtk_2.collect())
    BOWTIE2(CUTADAPT.out.cutadapt_1.collect(), CUTADAPT.out.cutadapt_2.collect())
    STAR(transcriptome_ch, GTF_ch, BOWTIE2.out.bowtie2_out.collect())
    RSEM(logname_ch, jobID_ch, STAR.out.STAR_out.collect())

}


