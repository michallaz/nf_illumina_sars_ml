process indelQual {
    tag "indelQual:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "forvariants.bam*"

    input:
    tuple val(sampleId), path(bam), path(bai)

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    script:
    """
    lofreq indelqual --ref \${GENOME_FASTA} \
                     --out forvariants.bam \
                     --dindel ${bam}
    samtools sort -@ ${params.threads} \
                  -o forvariants.bam \
                  forvariants.bam
    samtools index forvariants.bam
    """
}
