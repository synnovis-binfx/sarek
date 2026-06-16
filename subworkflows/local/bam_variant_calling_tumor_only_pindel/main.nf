include { PINDEL_PINDEL as PINDEL              } from '../../../modules/nf-core/pindel/pindel/main'
include { GATK4_MERGEVCFS  as MERGE_PINDEL     } from '../../../modules/nf-core/gatk4/mergevcfs/main.nf'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_PINDEL {
    take:
    input          // channel: [mandatory] [ meta, tumor_bam, tumor_bai ]
    fasta          // channel: [mandatory] [ fasta ]
    fai            // channel: [mandatory] [ fasta_fai ]
    pindel_targets 
    dict           // channel: /path/to/reference/fasta/dictionary

    main:
    versions = Channel.empty()

    // Pindel wont work with split intervals

    PINDEL(input, fasta, fai, pindel_targets) // Call SV variants in Pindel 

    // only one  vcf from the pindel with no intervals /task split upstream ; so no collec/merge of vcf
    
    //just run merge with single vcf
    //MERGE_PINDEL(PINDEL.out.vcf, dict)

    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata
    //vcf = Channel.empty().mix(MERGE_PINDEL.out.vcf).map{ meta, vcf -> [ meta - meta.submap + [ variantcaller:'pindel' ], vcf ] }
    //tbi = Channel.empty().mix(MERGE_PINDEL.out.tbi).map{ meta, tbi -> [ meta - meta.subMap + [ variantcaller:'pindel' ], tbi ] }

    versions = versions.mix(PINDEL.out.versions)

    emit:
    //vcf
    //tbi
    versions
}
