# Background to Chosen Custom Configuration Settings - Trinity/Aviti Data 

Configuration settings based around data generated with Qiagen-trinity-Element library prep and sequencing.

## 📑 Other Background Elements of Note Include:
*   Data generated to date is showing low duplication levels - in-line with moderate sequencing levels, 40-100ng input DNA and lower PCR cycles (as capture done on Aviti FC, no post capture PCR / PCR subsampling effects). 
*   **NOTE:** High off-target levels (approx 60% seen) also mean lower average dup levels due to those off-target being spread genome-wide and with singleton/low depth & thus unlikely duplication.
*   Low dup level combined with higher innate base quality scores for Element vs Illumina have influenced custom configurations here.

---

## ⚙️ Config List with Background as per panHaem_params.yaml Supplied at Run-Time

### `input`
```yaml
input: './resources/synnovis_configs/samplesheet_test_panhaem.csv'
```
*   SampleSheet generated as per Sarek documentation - with new meta input for chemistry type to allow optional switches between eg Illumina and Element tool config as required.

### `outdir`
```yaml
outdir: './test_results/'
```
*   N/A

### `split_fastq`
```yaml
split_fastq: 0
```
*   We want to keep fastq / bams as one, in-line with UMI actions - splitting is turned off inside the umi sub-workflow anyway. Currently, this size panel will not benefit from splitting regardless, especially with additional overhead of spinning up new instances when run in cloud.

### `rp_intervals`
```yaml
rp_intervals: '/assay_targets/panHaem/QIAseq_xHYB.CXHS-11459Z-00.roi-covered-QGEN-ONLY_sorted_with_calr_csf3r_merged.bed'
```
*   No comments. "rp" prefix denotes relative path to be converted to full path via params.resources - for flexibility of external files/path provision at run time (cloud running provision).

### `wes`
```yaml
wes: true
```
*   As per operation instructions for panels.

### `aligner`
```yaml
aligner: 'bwa-mem'
```
*   Standard selection. `bwa-mem2` not used here - requires ~40GB to 50GB of RAM to run and double for index generation. TBC.

### `tools`
```yaml
tools: 'lofreq,mutect2,cnvkit,vardict,pindel'
```
*   Selected callers for validation. Lofreq unknown quantity and is more at the dev stage - may not get to validation stage.

### `trim_fastq` & `trim_fastq_umi`
```yaml
trim_fastq: false
trim_fastq_umi: false
```
*   For in-situ pipeline fastp trimming after umi pipeline run. Not used here as replaced with our custom Cutadapt addition - see below.

### `trim_fastq_umi_cutadapt`
```yaml
trim_fastq_umi_cutadapt: true
```
*   Custom Cutadapt addition: for greater control of trimming after UMI extracted to read name, allowing removal of 5' and 3' adapter sequences and phased adapter. Better practice to remove all traces of adapter prior to umi grouping/consensus calling.

### `split_fastq_umi`
```yaml
split_fastq_umi: 0
```
*   Actually set to false currently in umi workflow. Not good to subgroup data when doing UMI grouping.

### `interleaved_umi`
```yaml
interleaved_umi: true
```
*   Output of umi workflow is interleaved.

### `auto_trim_adapter_umi`
```yaml
auto_trim_adapter_umi: true
```
*   For fastp operation; not used for Trinity/Aviti workflow with Cutadapt trimming.

### `rp_adapter_fasta_r1` & `rp_adapter_fasta_r2`
```yaml
rp_adapter_fasta_r1: '/assay_targets/trinity_adapters/trinity_phased_adapters.fa'
rp_adapter_fasta_r2: '/assay_targets/trinity_adapters/trinity_phased_adapters_r2.fa'
```
*   Adapters files to supply Cutadapt can now be supplied, which is the case for trinity chemistry. "rp" prefix denotes relative path to be converted to full path via params.resources - for flexibility of external files/path provision at run time (cloud running provision). Adapters themselves, are linked adapters (eg...) for all different phased adapter variations. - 5' required 3' not required. 3' kept short currently otherwise combinations would be 5x more. - TBC.

### `trim_fileR1_operation`, `trim_fileR2_operation`, & `discard_untrimmed_umi`
```yaml
trim_fileR1_operation: "-a"
trim_fileR2_operation: "-A"
discard_untrimmed_umi: true
```
*   Action type associated with the Cutadapt trimming (for read 1 and read 2 respectively). For our phased adapter combinations, this means the 5' part of the linked adapter pattern is required but the 3' is not required. As reads expect at least the 5' adapter (anchored to 5' end), after umi sequence removal to read name by Fgbio before Cutadapt action, we filter out sequences that don't have this. In practice, approx 98% reads retained as expected.

### `ignore_soft_clipped_bases`
```yaml
ignore_soft_clipped_bases: false
```
*   Mutect2 to use soft clips. - best practices for Mutect2.

### `normalize_vcfs`
```yaml
normalize_vcfs: true
```
*   Left align, and remove duplicates via bcftools only. Also uses `--multi-allelics -both` to split multi-allelics to separate records.

### `filter_vcfs`
```yaml
filter_vcfs: true
```
*   Filtering of vcfs via bcftools. Now using variant caller type from meta to pick VC specific filtering action. default = bcftools_filter_criteria below.

### `bcftools_filter_mutect2`, `bcftools_filter_criteria`, `bcftools_filter_lofreq`, & `bcftools_filter_vardict`
```yaml
bcftools_filter_mutect2: "-e "(FILTER!='PASS' && FILTER!~'haplotype' && FILTER!~'germline' && FILTER!~'strand_bias' && FILTER!~'panel_of_normals') ||  FILTER~'base_qual' ||  FILTER~'clustered_events' || FILTER~'multiallelic' ||  FILTER~'weak_evidence' || FILTER~'orientation' || FILTER~'slippage' || MAX(AF[*]) < 0.02""
bcftools_filter_criteria: "-e "AF < 0.02""
bcftools_filter_lofreq: "-e "HRUN > 10 || AF < 0.02""
bcftools_filter_vardict: "-e "(FILTER!='PASS' && FILTER!~'Bias' && FILTER!~'MSI15') || FILTER~'Q10' || FILTER~'q25' || FILTER~'NM5.25' || (FORMAT/AF < 0.02) || (FORMAT/VD < 6) || ( (FORMAT/VD < 10) && ((MQ < 55.0 && NM > 1.0) || (MQ < 60.0 && NM > 2.0) || (FORMAT/DP < 100) || (QUAL < 40.0)) )""
```
*   Options selected for variant filtering via bcftools. Based on current multi-modal pipeline set-up and initial VC options (see below), with changes to reflect seq chemistry including. **TO BE FURTHER ASSESSED IN VALIDATIONS.**
*   Base filtering of < 'q25', to reflect higher Aviti base qualities with lower duplication levels in mind (score partially downgraded by FGbio - although this is restricted with FGbio options below).
*   Currently not filtering strand-bias (bias) - will be keeping tag in DSS for assessment (should be limited numbers due to random fragmentation type chemistry).
*   `NM5.25` vardict filter (avg number mismatches in reads) to be assessed in validation (with aviti in mind).
*   Vardict secondary set of filters above for lower quality reads - to be assessed for aviti chemistry.
*   Variant support as per qiaseq/multimodal settings. Needs to be reduced for small validation samples (to 2-3).

### `skip_tools`
```yaml
skip_tools: 'baserecalibrator,baserecalibrator_report,markduplicates,vcftools'
```
*   No comment.

### CNVkit Options
```yaml
#cnvkit reference below is a default for current mgp4 panel; requires updating
rp_cnvkit_reference: '/assay_targets/panHaem/referencePON_filtered_43_liftover_to_38_reannotate.cnn'
rp_cnvkit_plot_targets_bed: '/assay_targets/panHaem/PanHaem_cnvkit_targets.bed'
cnv_custom_call_and_plot: true
```
*   PON for Cnvkit to be generated as part of validation. CNVkit options then further considered and added here (currently in modules config ext.args).

### `umi_read_structure`
```yaml
umi_read_structure: '7M+T 7M+T'
```
*   As per trinity library structure & UMI then template. Phased adapter actually comes after the UMI but this is trimmed sequentially by Cutadapt in our custom UMI workflow as described above (beyond capability of fgbio).

### `seq_platform`
```yaml
seq_platform: "element"
```
*   Added meta value for sequencing type. This can be used to switch between illumina and other technologies such as element here, for specific tool configs. Eg used to lofreq_comp custom module to switch between realignment or not. Lofreq realignment being assessed and currently used.

### `vardict_args` & `vardictfilter_args`
```yaml
vardict_args: "-f 0.02 -z -c 1 -S 2 -E 3 -g 4 -I 200 -X 1 -Q 20 --mfreq 0.15 "
vardictfilter_args: "-A -f 0.02 -p 6 -I 15 -q 25"
```
*   Args for Vardict. Based on Qiaseq/multimodal pipeline. `-X 2` (default) will treat MNPs as separates if separated by more than 2 base but use in realignment if only 2 or less base - conflict with HGVS but likely safer than splitting if residing in different codons. Other disadvantage, is may miss database entry as these likely to be split (as per GATK best practices).
*   To be discussed with lab at validation time. May be best for Vardict and Mutect to extend out further to always capture all forms of adjacent calls at expense of being less HGVS-compliant - although VEP etc should give you the correct call but in combined form (one call). `-X 1` Guarantees a Same-Codon Merge however according to Vardict logic and in line with Mutect `--max-mnp-distance 2` (default 1).
*   `--pcr-snv-qual 45` and model overrides act as the mathematical filter, preventing singletons (as per this Aviti low dup data) from calling false positive errors. See below also. To be assessed in validation.

### `mutect2_args` & `mutect_filters`
```yaml
mutect2_args: "--max-reads-per-alignment-start 0 --max-mnp-distance 2 --genotype-germline-sites true --genotype-pon-sites true --min-base-quality-score 25 --pcr-indel-model CONSERVATIVE --pcr-snv-qual 45 --pcr-indel-qual 42"
mutect_filters: "--min-median-base-quality 25 --min-median-mapping-quality 10 --max-alt-allele-count 5 --min-allele-fraction 0.02 --max-events-in-region 8 --min-median-read-position 1 --pcr-slippage-rate 0.075"
```
*   MNP discussed above for Vardict. 
*   `q25` in line with VarDict and higher qualities for Aviti data (taking into account some downgrade to scores for low duplication, fgbio consensus corrections (see below)).

### `lofreq_args`, `lofreq_indelqual_args`, & `lofreq_filter_args`
```yaml
lofreq_args: "--call-indels --min-alt-bq 20 --min-bq 20 --sig 1"
lofreq_indelqual_args: "--uniform 35,35"
lofreq_filter_args: "--verbose -v 10 -a 0.02 --sb-incl-indels --sb-alpha 0.001 --sb-mtc fdr --print-all"
```
*   Lofreq in dev/unknown quality to be assessed in validations. Settings currently to reduce filtering in first instance. Uniform base qualities picked to fit in with aviti data (not using pre-calculated system dindel which is based on Illumina data). Blunt standardization to 35 for indels should capture most calls and be representative of good aviti data.

### `rp_pindel_targets_bed` & `pindel_args`
```yaml
rp_pindel_targets_bed: '/assay_targets/panHaem/Panhaem_pindel_targets.bed'
pindel_args: '-M 6 -k'
```
*   Using restricted targets for Pindel currently - eg flt3, KMT2A. TBC with lab. Permissive settings as base for validation. 6 supporting reads and calling breakpoints. TBC.

### `mosdepth_thresholds` & `mosdepth_quantize`
```yaml
mosdepth_thresholds: "20,100,200,300,400,500,600,800,1000"
mosdepth_quantize: "0:5:100:200:300:400:500:600:700:800"
```
*   Thresholds to produce custom outputs for amount of supplied ROI covered and what to output to Multiqc for display.

### `samtools_depth`
```yaml
samtools_depth: "-q 15 -J -d 0"
```
*   As per Qiaseq/multimodal pipelines. Indels included into depth calc currently. And full read compliment used (no subsample capping at 8000 - default).
*   `-q15` currently permissive. May raise to `q25` in line with base callers. To be assessed in validations against variant caller depth calculations.

### `snv_consensus_calling` & `consensus_min_count`
```yaml
snv_consensus_calling: true
consensus_min_count: 2
```
*   Testing consensus of the variant callers (Vardict, Mutect, and Lofreq). Likely would require more work due to different VC call nuances, especially for indels.

### `fgbio_consensus`, `fgBioFilter_minReads`, `fgBioFilter_minBaseq`, & `fgBioFilter_maxBaseErrorRate`
```yaml
fgbio_consensus: '--error-rate-pre-umi 45 --error-rate-post-umi 40'
fgBioFilter_minReads: '1'
fgBioFilter_minBaseq: '25'
fgBioFilter_maxBaseErrorRate: '0.15'
json_sqvd: true
```
*   Best estimate of appropriate settings for low duplication level aviti data and how downstream variant callers use data. To be assessed in validation work.
*   `min reads = 1` per consensus as this is low dupl data with very little duplication. Still duplication may occur heavily at specific target regions warranting use. Off target also bringing this avg global figure down. TBC with lab on increasing duplication level (eg use of less input, less number per seq batch to increase seq depth).
*   `maxBaseErrorRate`: default for this parameter is 0.1 (a strict 10% maximum error rate). By raising it to 0.15 (15%), makes consensus filtering more permissive to address low dup singleton data with more chance of low quality scores by fgbio downgrading.
*   Singletons/low duplication and Aviti high quality: pre-UMI error prior is set to 45, no singleton read will exit Fgbio with a quality score higher than Q45, keeping its trust realistically capped. Post-umi 40 also with Aviti high quality data in mind (unlikely even with singleton that error happened post umi/adapter ligation eg in any pcr or sequencing reactions). Actually this is default when using min-reads = 1. To be validated against.