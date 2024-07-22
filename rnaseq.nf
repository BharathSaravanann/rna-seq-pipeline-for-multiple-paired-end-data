

process trimming {
publishDir "TRIMMED",mode:'copy'
input:
    tuple val(sample_id),path (reads)

output:
    path "*"
    path "*trimmed*.fq.gz",emit:trimmed
script:

"""
trim_galore --paired -q 20 --gzip --basename ${sample_id}_trimmed ${reads}
"""

}

process qc {
publishDir "qualityreport", mode:'copy'
input:
    path (reads)


output:
    path "*"


script:

"""
fastqc ${reads}
multiqc *fastqc*
mkdir FASTQC
mv *fastqc* FASTQC
"""



}

process STAR_INDEX{

publishDir "INDEX",mode:'copy'
input:
    path (fasta)
    path (gtf)


output:
path "*",emit:index

script:

"""
STAR --runThreadN 8 \\
--runMode genomeGenerate \\
--genomeDir index \\
--genomeFastaFiles ${fasta} \\
--sjdbGTFfile ${gtf} \\
--genomeSAindexNbases 12
"""

}

process STAR_MAPPING {
publishDir "MAPPING",mode:'copy'

cpus 4

input:
    tuple val(sample_id),path (read1),path (read2),path (index)

output:
    path "*"
    path "*.bam",emit:bams
script:

"""
STAR --runThreadN 8 \\
--genomeDir ${index} \\
--readFilesIn ${read1} ${read2} \\
--outSAMtype BAM SortedByCoordinate \\
--outFileNamePrefix ${sample_id} \\
--readFilesCommand zcat

"""

}

process FEATURE_COUNT {
publishDir "FEATURECOUNT",mode:'copy'

input:
    path (bams)
    path (gtf)
    val(strand)
    
output:
    path "*"

script:

"""

featureCounts -T 8 -s ${strand} -p --countReadPairs -t exon \\
-g gene_id -Q 10 -a ${gtf} -o gene_count ${bams}

multiqc gene_count*

"""


}

workflow {

ref_fasta=Channel.fromPath(params.ref_fasta)
ref_gtf= Channel.fromPath(params.ref_gtf)
fastq_ch=Channel.fromFilePairs(params.reads)
strand=Channel.of(params.strand)
trimming(fastq_ch).set{trimmed}

raw_fastq=fastq_ch.map{items -> items[1]}.flatten().collect()
trim_fastq=trimmed.trimmed.flatten().collect()
raw_fastq.mix(trim_fastq).collect() | qc
STAR_INDEX(ref_fasta,ref_gtf).set{star_index}
trimmed.trimmed.map{read1,read2 -> tuple ("${read1.getFileName()}".split("_trimmed")[0],read1,read2)}
|combine(star_index.index)|STAR_MAPPING
|set{bams}

bams.bams.collect().set{finalbams}

FEATURE_COUNT(finalbams,ref_gtf,strand)         

}
