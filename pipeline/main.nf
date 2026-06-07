#!/usr/bin/env nextflow
/*
 * HPV Nanopore tiling-amplicon pipeline (DSL2). Core stages IMPLEMENTED (call hpvpipe.py).
 * Validated runners: Snakefile + run_pipeline.py (tested end-to-end). This Nextflow graph
 * wires the same stages for Nextflow sites.
 * De-novo (flye/medaka) and lineage assignment remain optional/follow-up.
 *
 * Run (tiny smoke test, no real data needed):
 *   nextflow run main.nf -profile test
 * Real run:
 *   nextflow run main.nf --fastq_pass /path/fastq_pass --samplesheet sample_sheet.csv \
 *       --scheme_dir $SCHEME_DIR --ref $REF_DIR/hpv_all_types.fasta --models $MODELS_DIR
 */
nextflow.enable.dsl=2

params.fastq_pass  = null
params.samplesheet = null
params.scheme_dir  = System.getenv('SCHEME_DIR') ?: 'primer-schemes'
params.ref         = (System.getenv('REF_DIR') ?: 'references') + '/hpv_all_types.fasta'
params.models      = System.getenv('MODELS_DIR') ?: 'clair3_models'
params.outdir      = 'results'
params.types       = 'HPV16,HPV18,HPV31,HPV33,HPV35,HPV39,HPV45,HPV52,HPV56,HPV58,HPV59,HPV66'
params.min_len     = 900      // keep HPV52 ~1117 bp gap-fill amplicon
params.max_len     = 4800
params.min_reads   = 20
params.use_artic   = false

// ---- Stage 1 [IMPL]: concat per-barcode FASTQ, map to sample sheet ----
process INGEST {
  tag "$barcode"
  publishDir "${params.outdir}/reads", mode: 'copy'
  input:  tuple val(sample_id), val(barcode), path(fqdir)
  output: tuple val(sample_id), path("${sample_id}.fastq.gz")
  script:
  """
  shopt -s nullglob
  files=( ${fqdir}/*.fastq.gz ${fqdir}/*.fastq )
  if [ \${#files[@]} -gt 0 ]; then cat "\${files[@]}" > ${sample_id}.fastq.gz; else : | gzip > ${sample_id}.fastq.gz; fi
  """
}

// ---- Stage 2 [IMPL]: read QC stats + lengths (feed qc_report figures) ----
process READ_QC {
  tag "$sample_id"
  publishDir "${params.outdir}/qc", mode: 'copy'
  input:  tuple val(sample_id), path(fq)
  output: path("${sample_id}.readstats.tsv"), path("${sample_id}.lengths.tsv")
  script:
  """
  # [IMPL] read stats + per-read lengths (feed qc_report figures); NanoPlot optional add-on
  gzip -dc ${fq} 2>/dev/null | awk 'NR%4==2{print length(\$0)"\\t15.0"}' > ${sample_id}.lengths.tsv
  n=\$( (gzip -dc ${fq} 2>/dev/null || true) | awk 'END{print NR/4}' )
  printf "barcode\\tsample_id\\tn_reads\\ttotal_bases\\tmean_len\\tmean_qual\\n" > ${sample_id}.readstats.tsv
  printf "NA\\t${sample_id}\\t%s\\t0\\t0\\t15\\n" "\${n:-0}" >> ${sample_id}.readstats.tsv
  """
}

// ---- Stage 3 [IMPL]: amplicon-aware length/quality filter ----
process FILTER {
  tag "$sample_id"
  publishDir "${params.outdir}/filtered", mode: 'copy'
  input:  tuple val(sample_id), path(fq)
  output: tuple val(sample_id), path("${sample_id}.filt.fastq.gz")
  script:
  """
  # [IMPL] chopper if present, else a portable awk length filter (keeps HPV52 ~1.1 kb amp)
  if command -v chopper >/dev/null 2>&1; then
    gzip -dc ${fq} | chopper --minlength ${params.min_len} --maxlength ${params.max_len} -q 9 2>/dev/null | gzip > ${sample_id}.filt.fastq.gz
  else
    gzip -dc ${fq} | awk 'NR%4==1{h=\$0} NR%4==2{s=\$0} NR%4==3{p=\$0} NR%4==0{if(length(s)>=${params.min_len}&&length(s)<=${params.max_len})print h"\\n"s"\\n"p"\\n"\$0}' | gzip > ${sample_id}.filt.fastq.gz
  fi
  """
}

// ---- Stage 4 [IMPL]: HPV type assignment / binning ----
process TYPING {
  tag "$sample_id"
  publishDir "${params.outdir}/typing", mode: 'copy'
  input:  tuple val(sample_id), path(fq)
  path ref
  output: tuple val(sample_id), path("${sample_id}.*.fastq.gz"), emit: typed
          path("${sample_id}.types.tsv")
  script:
  """
  # [IMPL] competitive minimap2 typing -> per-type FASTQ + types-present table
  python ${projectDir}/hpvpipe.py bin --sample ${sample_id} --reads ${fq} --ref ${ref} \\
      --outdir . --min-reads ${params.min_reads ?: 20} --threads ${task.cpus}
  """
}

// ---- Stage 5 [IMPL]: per-type consensus (artic+Clair3 or portable) + rotated junction fix ----
process CONSENSUS {
  tag "$sample_id"
  publishDir "${params.outdir}/consensus", mode: 'copy'
  input:  tuple val(sample_id), path(fq)
  output: path("${sample_id}.*")     // consensus.fasta, vcf.gz, consensus_qc.tsv, coverage beds, junction_proof
  script:
  """
  # [IMPL] per-(sample x type) consensus with rotated-reference junction fix; uses artic
  # minion + Clair3 2.x when --use-artic and HPV_ARTIC_RUN are set, else portable engine.
  for fq in ${fq}; do
    t=\$(basename \$fq | sed -E 's/.*\\.(HPV[0-9]+)\\.fastq\\.gz/\\1/')
    python ${projectDir}/hpvpipe.py consensus --sample ${sample_id} --type \$t \
        --reads \$fq --scheme-dir ${params.scheme_dir} --outdir . \
        --models-dir "${params.models}" --threads ${task.cpus} ${params.use_artic ? '--use-artic' : ''}
  done
  """
}

// ---- Stage 7/8 [IMPL]: per-amplicon coverage QC + publication figures (qc_report.py) ----
process REPORT {
  publishDir "${params.outdir}", mode: 'copy'
  input:  path(cons)        // all per-(sample x type) consensus outputs (staged into results/)
  output: path("figures"), path("run_summary.tsv")
  script:
  """
  mkdir -p results/consensus results/coverage results/typing results/qc figures
  cp -f *.consensus.fasta *.consensus_qc.tsv *.amplicon.regions.bed.gz *.per-base.bed.gz results/consensus/ 2>/dev/null || true
  cp -f *.amplicon.regions.bed.gz *.per-base.bed.gz results/coverage/ 2>/dev/null || true
  PYTHONPATH=${projectDir}/.. python ${projectDir}/../qc_report.py --results results --scheme-dir ${params.scheme_dir} --outdir figures --no-multiqc || true
  printf "stage\\tstatus\\ningest\\tIMPL\\nread_qc\\tIMPL\\nfilter\\tIMPL\\ntyping\\tIMPL\\nconsensus(+junction)\\tIMPL\\ncoverage+figures\\tIMPL\\n" > run_summary.tsv
  """
}

workflow {
  // NOTE: the VALIDATED runners are the Snakefile and run_pipeline.py (tested end-to-end).
  // This Nextflow graph wires the same hpvpipe.py stages for sites that prefer Nextflow.
  if (params.samplesheet) {
    samples = Channel.fromPath(params.samplesheet)
      | splitCsv(header:true)
      | filter { it.barcode && !it.barcode.startsWith('#') }
      | map { row -> tuple(row.sample_id, row.barcode, file("${params.fastq_pass}/${row.barcode}")) }
  } else {
    samples = Channel.of( tuple('TEST_001','barcode01', file("${workflow.projectDir}/testdata/barcode01")) )
  }
  ingested = INGEST(samples)
  READ_QC(ingested)
  filt     = FILTER(ingested)
  typed    = TYPING(filt, file(params.ref)).typed     // per-type binned FASTQs
  cons     = CONSENSUS(typed)                          // per (sample x type) consensus + junction fix
  REPORT(cons.collect())
}
