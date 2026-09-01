//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'
include { SENTIEON_BWAMEM        } from '../../../modules/nf-core/sentieon/bwamem/main'
include { FGBIO_ZIPPERBAMS        } from '../../../modules/nf-core/fgbio/zipperbams/main'

// for umi changes would need to have reads + bam prior to umi coming in here from updated tuple in fastq_gatk workflow - meta, reads, bam (or []) if umi ///reads = this tuple ; else reads = current tuple; then tuple for bwa and new module fgbio zipperbams (same for index)
workflow FASTQ_REALIGN_UMI {
    take:
    reads_bams // channel: [mandatory] meta, reads, ubams
    index // channel: [mandatory] index
    sort  // boolean: [mandatory] true -> sort, false -> don't sort
    fasta
    fasta_fai

    main:

    versions = Channel.empty()
    reports = Channel.empty()
    
    reads = reads_bams.map{ meta, reads, ubam -> [ meta, reads ] }
    ubams = reads_bams.map{ meta, reads, ubam -> [ meta, ubam ] }

    // Only one of the following should be run
    BWAMEM1_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(reads, index, [[id:'no_fasta'], []], sort) // If aligner is dragmap
    // The sentieon-bwamem-module does sorting as part of the conversion from sam to bam.
    SENTIEON_BWAMEM(reads, index, fasta, fasta_fai) // If aligner is sentieon-bwamem

    // Get the bam files from the aligner
    // Only one aligner is run
    bam = Channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)
    bam = bam.mix(DRAGMAP_ALIGN.out.bam)
    bam = bam.mix(SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bam ] })

    bai = SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bai ] }

    //fgbio bamzipper
    ubam_bam = bam.join(ubams).map{ meta, bam, ubam -> [ meta, bam, ubam ] }
    refs_ch =fasta.combine(fasta_fai)
        .combine(index)
        .map { fa_meta, fa_file, fai_meta, fai_file, bwa_meta, bwa_files ->
            tuple(
                [meta_id: 'id'],
                fa_file,
                fai_file,
                bwa_files
            )
        }
    refs_ch.view()
    FGBIO_ZIPPERBAMS(ubam_bam, refs_ch)
    bam = FGBIO_ZIPPERBAMS.out.bam

    // Gather reports of all tools used
    reports = reports.mix(DRAGMAP_ALIGN.out.log)

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(DRAGMAP_ALIGN.out.versions)
    versions = versions.mix(SENTIEON_BWAMEM.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    reports
    versions // channel: [ versions.yml ]
}
