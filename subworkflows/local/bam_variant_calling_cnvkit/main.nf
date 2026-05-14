//
// CNVKIT calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH                                     } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_CALL                                      } from '../../../modules/nf-core/cnvkit/call/main'
include { CNVKIT_CUSTOM_CALL_AND_PLOT                      } from '../../../modules/local/cnvkit_custom/custom_call_and_plot/main'
include { CNVKIT_EXPORT                                    } from '../../../modules/nf-core/cnvkit/export/main'
include { CNVKIT_GENEMETRICS                               } from '../../../modules/nf-core/cnvkit/genemetrics/main'

workflow BAM_VARIANT_CALLING_CNVKIT {
    take:
    cram                // channel: [mandatory] meta, bam, vcf, normal_  // @asmith
    fasta               // channel: [mandatory] meta, fasta
    fasta_fai           // channel: [optional]  meta, fasta_fai
    targets             // channel: [mandatory] meta, bed
    reference           // channel: [optional]  meta, cnn
    cnvkit_plot_targets

    main:
    versions = Channel.empty()
    generate_pon = false

    // vcf file now included from input list and tuble element in module @asmith
    CNVKIT_BATCH(cram, fasta, fasta_fai, targets, reference, generate_pon)

    // right now option to use input VCF to improve the calling of B alleles under specific module and using mutect outputs from tumor only workflow only  @asmith
    // based on SNV frequencies from the VCF file
    //CNVKIT_CALL(CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2], []]})

    // Call and plot module for custom calling copy number and plotting (using output from call) and input VCF for LOH/BAF @asmith
    if (params.cnv_custom_call_and_plot) { 
        cns_in = CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2]]}
        vcf_bam_collect = cns_in.join(cram).join(CNVKIT_BATCH.out.cnr).map{ meta, cns, bam_, vcf_, normal_, cnr -> [meta, cns, cnr, vcf_]}
        vcf_bam_collect.view()
        CNVKIT_CUSTOM_CALL_AND_PLOT(vcf_bam_collect, cnvkit_plot_targets) // now using cnvkit LOI for targeted plotting
        versions = versions.mix(CNVKIT_CUSTOM_CALL_AND_PLOT.out.versions)
        CNVKIT_EXPORT(CNVKIT_CUSTOM_CALL_AND_PLOT.out.cns)
        cnv_calls_out = CNVKIT_CUSTOM_CALL_AND_PLOT.out.cns
        }
    else {
        CNVKIT_CALL(CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2], []]})
        versions = versions.mix(CNVKIT_CALL.out.versions)
        CNVKIT_EXPORT(CNVKIT_CALL.out.cns)
        cnv_calls_out = CNVKIT_CALL.out.cns
    }

        

    ch_genemetrics = CNVKIT_BATCH.out.cnr.join(CNVKIT_BATCH.out.cns).map{ meta, cnr, cns -> [meta, cnr, cns[2]]}
    CNVKIT_GENEMETRICS(ch_genemetrics)

    versions = versions.mix(CNVKIT_BATCH.out.versions)
    versions = versions.mix(CNVKIT_GENEMETRICS.out.versions)
    versions = versions.mix(CNVKIT_EXPORT.out.versions)

    emit:
    cnv_calls               = cnv_calls_out                   // channel: [ meta, cns ]
    cnv_calls_export        = CNVKIT_EXPORT.out.output         // channel: [ meta, export_format ]
    versions                                                   // channel: [ versions.yml ]
}
