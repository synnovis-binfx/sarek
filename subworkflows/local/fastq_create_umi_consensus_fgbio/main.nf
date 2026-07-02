//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLUMICONSENSUS } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FASTQTOBAM                  as FASTQTOBAM       } from '../../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI  } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FASTQ_ALIGN                       as ALIGN_UMI        } from '../fastq_align/main'
include { SAMTOOLS_MERGE                    as MERGE_CONSENSUS  } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_BAM2FQ                   as BAM2FASTQ        } from '../../../modules/nf-core/samtools/bam2fq/main.nf'
include { FASTP as FASTP_UMI                                    } from '../../../modules/nf-core/fastp/main'
include { CUTADAPT as CUTADAPT_UMI                              } from '../../../modules/nf-core/cutadapt/main'

workflow FASTQ_CREATE_UMI_CONSENSUS_FGBIO {
    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
    fai                       // channel: [optional] /path/to/reference/fasta_fai, needed for Sentieon
    map_index                 // channel: [mandatory] Pre-computed mapping index
    groupreadsbyumi_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"
    adapter_fasta_r1          // channel: [optional] collected file path for r1 adapter sequences 
    adapter_fasta_r2          // channel: [optional] collected file path for r2 adapter sequences

    main:
    ch_versions = Channel.empty()

    // params.umi_read_structure is passed out as ext.args
    // FASTQ reads are converted into a tagged unmapped BAM file (uBAM)
    FASTQTOBAM(reads)

    // in order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    // TODO check if DRAGMAP works well with BAM inputs
    // but keep the appropriate UMI tags in the FASTQ comment field and produce
    // an interleaved FASQT file (hence, split = false)
    split = false
    BAM2FASTQ(FASTQTOBAM.out.bam, split)

    // Trimming prior to consensus calling @asmith
    if (params.trim_fastq_umi) {

         save_trimmed_fail = false
         save_merged = false
         FASTP_UMI(
             BAM2FASTQ.out.reads,
             [], // we are not using any adapter fastas at the moment
             false, // we don't use discard_trimmed_pass at the moment
             save_trimmed_fail,
             save_merged
         )
         reads_for_alignment = FASTP_UMI.out.reads
         //reports = reports.mix(FASTP_UMI.out.json.collect{ _meta, json -> json })
         //reports = reports.mix(FASTP_UMI.out.html.collect{ _meta, html -> html })


         ch_versions = ch_versions.mix(FASTP_UMI.out.versions)
    }

    // Trimming prior to consensus calling with cutadapt @asmith
    if (params.trim_fastq_umi_cutadapt) {

         CUTADAPT_UMI(
             BAM2FASTQ.out.reads,  // currently interleaved fastq
             adapter_fasta_r1, // place holder for adapter fasta for r1; rather than stating in ext.args
             adapter_fasta_r2, // place holder for adapter fasta for r2; rather than stating in ext.args
         )
         reads_for_alignment = CUTADAPT_UMI.out.reads
         //reports = reports.mix(CUTADAPT_UMI.out.json.collect{ _meta, json -> json })
         //reports = reports.mix(CUTADAPT_UMI.out.html.collect{ _meta, html -> html })


         ch_versions = ch_versions.mix(CUTADAPT_UMI.out.versions)


    } else {
         reads_for_alignment = BAM2FASTQ.out.reads
    }
    

    // appropriately tagged interleaved FASTQ reads are mapped to the reference
    // bams will not be sorted (hence, sort = false)
    sort = false
    ALIGN_UMI(reads_for_alignment, map_index, sort, fasta, fai)

    bams_to_merge = ALIGN_UMI.out.bam
    // id currently includes the lane, so swap to just id=sample and groupKey to avoid blocking
        .map {meta, bam ->
            tuple( groupKey(meta + [id:meta.sample], meta.num_lanes), bam)
            }
        .groupTuple()
        // undo the groupKey, else the meta map is not a normal map.
        .map{meta, bam -> tuple(meta.target, bam)}
        .branch { meta, bam ->
            single: meta.num_lanes <= 1
            return [meta, bam[0]]
            multiple: meta.num_lanes > 1
        }

    // Merge across runs/lanes for the same sample
    MERGE_CONSENSUS(bams_to_merge.multiple, [[], []], [[], []])

    bams_all = MERGE_CONSENSUS.out.bam.mix(bams_to_merge.single)

    // appropriately tagged reads are now grouped by UMI information
    GROUPREADSBYUMI(bams_all, groupreadsbyumi_strategy)

    // Using newly created groups
    // To call a consensus across reads in the same group
    // And emit a consensus BAM file
    // TODO: add params for call_min_reads and call_min_baseq
    call_min_reads = 1
    call_min_baseq = 10
    CALLUMICONSENSUS(GROUPREADSBYUMI.out.bam, call_min_reads, call_min_baseq)

    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)
    ch_versions = ch_versions.mix(ALIGN_UMI.out.versions)
    ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions)
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)
    ch_versions = ch_versions.mix(MERGE_CONSENSUS.out.versions)

    emit:
    umibam         = FASTQTOBAM.out.bam             // channel: [ val(meta), [ bam ] ]
    groupbam       = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam   = CALLUMICONSENSUS.out.bam       // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
    umigrouphist   = GROUPREADSBYUMI.out.histogram  // channel: [ val(meta), [ histogram ]  ]
}
