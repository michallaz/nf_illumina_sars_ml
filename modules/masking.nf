process masking {
    tag "masking:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai)
    val(primers)
    val(pairs)

    output:
    tuple val(sampleId), path('ivar_trimmed_all.bam')

    script:
    """
    length=`echo "${params.length} - 40" | bc -l`
    ivar trim -i ${bam} \
             -b ${primers} \
             -m \${length} \
             -f ${pairs} \
             -q ${params.quality_initial} \
             -e \
             -p ivar_trimmed_all
    """
}
