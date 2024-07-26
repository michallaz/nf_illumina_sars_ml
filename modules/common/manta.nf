process manta {
    tag "manta:$sampleId"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy'
    container 'nf_illumina_sars-3.0-manta:latest'

    input:
    tuple val(sampleId), path(bam_file), path(bai_file), path(consensus_masked_fasta)
    tuple val(sampleId), path(GENOME_FASTA)

    output:
    tuple val(sampleId), path(consensus_masked_fasta)

    script:
    """
    # This is dummy module. It is not working. We only rewrite input to output
    for i in ${consensus_masked_fasta}; do
        cp \$i \${i%.fasta}_SV.fasta
    done

#    split_fasta() {
#        # Split fasta file into multiple files
#        # File names are based on the fasta header and provided suffix <header><sufix>.fasta
#
#        local input_file="\$1"
#        awk '
#        /^>/ {
#            if (seq)
#                print seq > file
#            file = substr(\$0, 2) "_manta.fasta"
#            print \$0 > file
#            seq = ""
#            next
#        }
#        {
#            seq = seq \$0 "\\n"
#        }
#        END {
#            if (seq)
#              print seq > file
#        }
#        ' "\$input_file"
#    }
#
#    samtools faidx ${GENOME_FASTA}
#    ILE_ODCZYTOW=`samtools view ${bam_file} | wc -l`
#    if [  \${ILE_ODCZYTOW} -lt 1000 ]; then
#        # pusty bam, nie puszczamy manty, po prostu tworzymy kopie plikow z poprawionymi nazwami
#        HEADER=`head -1 ${consensus_masked_fasta}`
#        NEW_HEADER=`echo -e "\${HEADER}_SV"`
#        cat ${consensus_masked_fasta} | sed s"/\${HEADER}/\${NEW_HEADER}/"g > output_consensus_masked_SV.fa
#    else
#        python /opt/docker/manta/bin/configManta.py --bam ${bam_file} --reference ${GENOME_FASTA} --runDir Manta_results
#        python Manta_results/runWorkflow.py -j ${params.threads} --quiet
#
#        if [ -e Manta_results/results/variants/diploidSV.vcf.gz ]; then
#            # Wywalamy skomplikowane SV jak translokacje itd typun BND
#            bcftools view -O z -o manta_results.vcf.gz -i '(FILTER="PASS" | FILTER="MaxDepth" | FILTER="NoPairSupport") && SVTYPE != "BND"' Manta_results/results/variants/diploidSV.vcf.gz
#                tabix manta_results.vcf.gz
#
#            ILE_SV=`zcat manta_results.vcf.gz  | grep -v '#' | grep -v cont | wc -l`
#
#            if [ \${ILE_SV} -gt 0 ]; then
#                cat ${GENOME_FASTA} | bcftools consensus -s - manta_results.vcf.gz  > output_manta.fa
#                split_fasta output_manta.fa
#                for fasta_file in `ls *_manta.fasta`; do
#                    insert_SV_python2.py ${consensus_masked_fasta} \${fasta_file} "\${fasta_file%.fasta}_SV.fasta"
#                done
#            else
#                HEADER=`head -1 ${consensus_masked_fasta}`
#                NEW_HEADER=`echo -e "\${HEADER}_SV"`
#                cat ${consensus_masked_fasta} | sed s"/\${HEADER}/\${NEW_HEADER}/"g > output_consensus_masked_SV.fa
#            fi
#        else
#            echo "Error brak pliku Manta_results/results/variants/diploidSV.vcf.gz"
#            exit 1
#        fi
#    fi
    """
}
