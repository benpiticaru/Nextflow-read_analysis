#!/usr/bin/env nextflow

// Define pipeline parameters
params.sample_sheet = "/{path_to}/sample_sheet.csv" // Path to the sample sheet CSV file
params.outdir = "exp_4.out" // Output directory
params.genome_reference = "/{path_to}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" // Path to genome reference
params.bed_reference = "/{path_to}/Homo_sapiens_GRCh38_110.bed" // Path to BED reference file

// Process: readLengthandCountAnalysis
// This process calculates read lengths and counts from processed FASTQ files and combines them with sequencing metadata.
process readLengthandCountAnalysis {
    label "smallJob"

    input:
        tuple val(sample_id), val(replicate), val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count)
        file(sequence_info)

    output:
        path("final_info.csv")

    script:
    """
    python3 scripts/read_analysis.py read_length_and_count_analysis \
        --processed_fastq ${processed_fastq} \
        --sequence_info ${sequence_info} \
        --replicate ${replicate} \
        --condition ${condition} \
        --output_file final_info.csv
    """
}

// Process: splitReads
// This process combines FASTQ files, removes adaptors, and splits the reads into smaller chunks for downstream processing.
process splitReads {
    label "split_reads"

    input:
        tuple val(sample_id), val(replicate), val (condition), path(fast5), path(fastq_dir)
    output:
        path ("split_locations.tsv")

    script:
    """
    ml seqkit
    ml scipy-stack

    cat ${fastq_dir}/pass/*.fastq.gz > combined.fastq.gz
    fastqwiper --fastq_in combined.fastq.gz --fastq_out ${sample_id}.fastq.gz
    seqkit split -s 500000 -j ${task.cpus} ${sample_id}.fastq.gz
    new_files=\$(ls ${sample_id}.fastq.gz.split/${sample_id}.part_*.fastq.gz)
    for file in \$new_files; do
        part_number="\${file##*_}"
        part_number="\${part_number%%.*}"
        part_number=\$(echo "\$part_number" | sed 's/^0*//')
        echo -e "${sample_id}\t${replicate}\t${condition}\t\$PWD/${fast5}\t\$PWD/${fastq_dir}\t\$PWD/\${file}\t\${part_number}" >> split_locations.tsv
    done
    """
}

// Process: removeAdaptors
// This process removes adaptors from raw FASTQ files using Porechop.
process removeAdaptors {
    label "remove_adaptor"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), val(output_file_count)

    output:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path("${sample_id}.${output_file_count}.trimmed.fastq"), val(output_file_count), emit: data

    script:
    """
    porechop -i ${raw_fastq} -o "${sample_id}.${output_file_count}.trimmed.fastq" --threads ${task.cpus} 
    """
}

// Process: alignAndSort
// This process aligns reads to the genome reference, sorts the alignments, and generates read metrics.
process alignAndSort {
    label "aligning"

    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: "*.bam*"

    input:
        tuple val(sample_id), val(replicate), val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count)
        file genome_reference
        file reference_bed

    output:
        tuple val(sample_id), val(replicate), val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count), emit: read_info
        path("${sample_id}.${output_file_count}.read_metrics.tsv"), emit: read_metrics
        path("sequence_info.csv"), emit: sequencing_metrics
        tuple path("${sample_id}.${output_file_count}.dna.bam.bai"), path("${sample_id}.${output_file_count}.dna.bam"), emit: dna_alignments

    script:
    """
    ml minimap2 samtools bedtools scipy-stack seqtk
    seqtk seq -a ${processed_fastq} > reads.fasta
    minimap2 -t ${task.cpus} -ax splice ${genome_reference} reads.fasta > aligned.sam
    samtools sort -@ ${task.cpus} aligned.sam | samtools view -@ ${task.cpus} -F 260 -q 60 -O BAM -o ${sample_id}.${output_file_count}.dna.bam -
    samtools index -@ ${task.cpus} ${sample_id}.${output_file_count}.dna.bam

    grep -w "gene" ${reference_bed} | sort > gene_only.bed
    bedtools bamtobed -i ${sample_id}.${output_file_count}.dna.bam > reads.bed
    samtools view ${sample_id}.${output_file_count}.dna.bam | awk '{print length(\$10)}' | paste reads.bed - > read_lengths.bed
    bedtools intersect -f 0.90 -r -a read_lengths.bed -b gene_only.bed -wb | \
        awk -v OFS='\t' '{print \$4, \$11, \$10-\$9, \$7}' | \
        sed 's/"//; s/";//' -  | sort -t \$'\t' -k 1 > gene_intersects.bed

    grep "gene_biotype" gene_only.bed | \
        sed -E 's/.*gene_id "([^"]*)".*gene_biotype "([^"]*)".*/\\1\t\\2/' | uniq > gene_id_biotype.tsv

    python3 scripts/read_analysis.py process_align_and_sort \
        --gene_id_biotype gene_id_biotype.tsv \
        --gene_intersects gene_intersects.bed \
        --sample_id ${sample_id} \
        --replicate ${replicate} \
        --condition ${condition} \
        --output_metrics ${sample_id}.${output_file_count}.read_metrics.tsv \
        --output_sequence_info sequence_info.csv
    """
}

// Process: polyATailAnalysis
// This process analyzes polyA tail lengths using Nanopolish.
process polyATailAnalysis {
    label "polyA"

    publishDir "${params.outdir}/polya", mode: 'copy', pattern: "*polya_length.tsv"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count)
        tuple file(bam_bai), file(bam)
        file(genome_reference)

    output:
        path("${sample_id}.${output_file_count}.polya_length.tsv"), emit: polya_length_tsv

    script:
    """
    ml StdEnv/2020 gcc/9.3.0 nanopolish/0.13.2
    nanopolish index --directory=${fast5} --sequencing-summary=${fastq_dir}/sequencing_summary.txt ${raw_fastq}
    nanopolish polya --threads=1 --reads=${raw_fastq} --bam=${bam} --genome=${genome_reference} > polya.output
    grep -w "PASS" polya.output | \
        awk -v OFS='\t' -v sample_id=${sample_id} -v replicate=${replicate} -v condition=${condition} '{print sample_id, replicate, condition, \$1, \$9}' \
            > ${sample_id}.${output_file_count}.polya_length.tsv
    """
}

// Process: polyaGraphing
// This process generates boxplots and statistical analysis for polyA tail lengths.
process polyaGraphing {
    label "smallJob"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "PolyA_length_boxplot.png"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: "anova_tukey_output.txt"

    input:
        file polya_tsv

    output:
        file("PolyA_length_boxplot.png")
        file("anova_tukey_output.txt")

    script:
    """
    python3 scripts/read_analysis.py generate_polya_graphing \
        --polya_tsv ${polya_tsv} \
        --output_boxplot PolyA_length_boxplot.png \
        --output_anova anova_tukey_output.txt
    """
}

// Process: combiningTsvFiles
// This process combines sequencing metrics from multiple samples into a summary TSV file.
process combiningTsvFiles {
    label "smallJob"

    publishDir "${params.outdir}", mode: 'copy', pattern: "experiments_summary.tsv"

    input:
        file sequence_info

    output:
        path("exp_summary.tsv"), emit: experiments_summary
        path("experiments_summary.tsv")

    script:
    """
    python3 scripts/read_analysis.py combine_tsv_files \
        --sequence_info ${sequence_info} \
        --output_summary exp_summary.tsv \
        --output_experiments experiments_summary.tsv
    """
}

// Process: readLengthHistogram
// This process generates histograms and statistical summaries for read lengths.
process readLengthHistogram {
    label "smallJob"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.tsv"

    input:
        file experiments_summary

    output:
        file("Read_length_histogram_by_condition.png")
        file("summary_statistics.tsv")
        file("mann_whitney_u_test_result.tsv")

    script:
    """
    python3 scripts/read_analysis.py generate_read_length_histogram \
        --experiments_summary ${experiments_summary} \
        --output_histogram Read_length_histogram_by_condition.png \
        --output_stats summary_statistics.tsv \
        --output_test mann_whitney_u_test_result.tsv
    """
}

// Process: readLengthScatterplot
// This process generates scatterplots for read lengths and TPM values.
process readLengthScatterplot {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"

    input:
        path(read_lengths_tsv)

    output:
        file("read_length_scatterplot.png")
        file("TPM_scatterplot.png")

    script:
    """
    python3 scripts/read_analysis.py generate_read_length_scatterplot \
        --read_lengths_tsv ${read_lengths_tsv} \
        --output_scatterplot read_length_scatterplot.png \
        --output_tpm_scatterplot TPM_scatterplot.png
    """
}

// Process: geneMakeupGraph
// This process generates bar plots showing the gene makeup of the library.
process geneMakeupGraph {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "gene_makeup_graph.png"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.tsv"

    input:
        file read_gene_ids

    output:
        file("gene_makeup_graph.png")
        file("gene_makeup_mrna.tsv")
        file("gene_makeup_other.tsv")

    script:
    """
    python3 scripts/read_analysis.py generate_gene_makeup_graph \
        --read_gene_ids ${read_gene_ids} \
        --output_graph gene_makeup_graph.png \
        --output_mrna gene_makeup_mrna.tsv \
        --output_other gene_makeup_other.tsv
    """
}

// Process: upsetPlot
// This process generates an UpSet plot to visualize gene overlaps between samples.
process upsetPlot {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"

    input:
        file read_gene_ids

    output:
        file("upset_plot.png")

    script:
    """
    python3 scripts/read_analysis.py generate_upset_plot \
        --read_gene_ids ${read_gene_ids} \
        --output_plot upset_plot.png
    """
}

// Workflow definition
workflow {
    genome_reference = file(params.genome_reference)
    bed_reference = file(params.bed_reference)
    sample_sheet = file(params.sample_sheet)

    Channel
        .fromPath(params.sample_sheet)
        .splitCsv( header: true, sep: ',' )
        .map { row -> tuple( row.sample_id,row.replicate,row.condition,row.fast5,row.fastq_dir) }
        .set { data_ch }

    splitReads(data_ch).collectFile(name: 'split_reads.tsv',newLine: true)
        .splitCsv( header: false, sep: '\t' )
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6]) }
        .set { split_data_ch }
    
    removeAdaptors(split_data_ch)
    alignAndSort(removeAdaptors.out.data,genome_reference,bed_reference)
    read_info_tsv = alignAndSort.out.read_metrics.collectFile(name: 'all_read_info.tsv', newLine: true)
    readLengthandCountAnalysis(alignAndSort.out.read_info,alignAndSort.out.sequencing_metrics)
    polyATailAnalysis(alignAndSort.out.read_info,alignAndSort.out.dna_alignments,genome_reference)
    combiningTsvFiles(readLengthandCountAnalysis.out.collectFile(name: 'sequence_info.csv', newLine: true))
    readLengthHistogram(combiningTsvFiles.out.experiments_summary)
    polyaGraphing(polyATailAnalysis.out.polya_length_tsv.collectFile(name: 'all_polya_lengths.tsv', newLine: true))
    readLengthScatterplot(read_info_tsv)
    geneMakeupGraph(read_info_tsv)
    upsetPlot(read_info_tsv)
}