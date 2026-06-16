include { VARDICTJAVA as VARDICT               } from '../../../modules/nf-core/vardictjava/main.nf'
include { GATK4_MERGEVCFS  as MERGE_VARDICT    } from '../../../modules/nf-core/gatk4/mergevcfs/main.nf'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_VARDICTJAVA {
    take:
    input     // channel: [mandatory] [ meta, tumor_cram, tumor_crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ]
    dict      // channel: /path/to/reference/fasta/dictionary

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map
        .map {meta, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [ num_intervals:num_intervals ], tumor_cram, tumor_crai, intervals]}

    VARDICT(input_intervals, fasta, fai) // Call variants with  vardict

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = VARDICT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = VARDICT.out.tbi.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()

    MERGE_VARDICT(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata
    vcf = Channel.empty().mix(MERGE_VARDICT.out.vcf, vcf_branch.no_intervals).map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'vardict' ], vcf ] }
    tbi = Channel.empty().mix(MERGE_VARDICT.out.tbi, tbi_branch.no_intervals).map{ meta, tbi -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'vardict' ], tbi ] }

    versions = versions.mix(MERGE_VARDICT.out.versions)

    emit:
    vcf
    tbi
    versions
}
