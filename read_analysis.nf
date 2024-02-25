#!/usr/bin/env nextflow
params.sample_sheet = "/home/bdpitica/scratch/RNA_vs_mRNA/scripts/final_sample_sheet.csv"
params.outdir = "final_exp_4.out"
params.genome_reference = "/home/bdpitica/scratch/RNA_vs_mRNA/reference/Homo_sapiens.GRCh38.cdna.all.fa"
params.cds_reference = "/home/bdpitica/scratch/RNA_vs_mRNA/reference/Homo_sapiens.GRCh38.cds.all.fa"
params.bed_reference = "/home/bdpitica/scratch/RNA_vs_mRNA/reference/Homo_sapiens_GRCh38_110.bed"

process readLengthandCountAnalysis {
    label "smallJob"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(processed_fastq), val(output_file_count)
        tuple file(read_lengths), file(sequence_info)


    output:
        path("final_info.csv")

    script:
    $/
    python3 - <<EOF

    import os
    from Bio import SeqIO
    import pandas as pd

    def calculate_read_lengths(file_name):
        read_lengths = []
        with open(file_name, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_lengths.append(len(record.seq))
        return read_lengths

    read_lengths = calculate_read_lengths("${processed_fastq}")
    read_counts = len(read_lengths)

    df = pd.read_csv("${sequence_info}",sep=",",header=None)

    df["read_counts"] = read_counts
    df["read_lengths"] = [read_lengths]
    df["replicate"] = "${replicate}"
    df["condition"] = "${condition}"

    df.to_csv("final_info.csv",sep=",",header=None,index=False)
    EOF
    /$
}

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
    fastqwiper --fastq_in combined.fastq.gz --fastq_out ${sample_id}.cleaned.fastq.gz
    seqkit split -s 1000000 -j ${task.cpus} ${sample_id}.cleaned.fastq.gz
    new_files=\$(ls ${sample_id}.cleaned.fastq.gz.split/${sample_id}.cleaned.part_*.fastq.gz)
    for file in \$new_files; do
        part_number="\${file##*_}"
        part_number="\${part_number%%.*}"
        part_number=\$(echo "\$part_number" | sed 's/^0*//')
        echo -e "${sample_id}\t${replicate}\t${condition}\t\$PWD/${fast5}\t\$PWD/${fastq_dir}\t\$PWD/\${file}\t\${part_number}" >> split_locations.tsv
    done
    """

}

process removeAdaptors {
    label "remove_adaptor"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), val(output_file_count)

    output:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path("${sample_id}.${output_file_count}.trimmed.fastq"), val(output_file_count)

    script:
    """
    porechop -i ${raw_fastq} -o "${sample_id}.${output_file_count}.trimmed.fastq" --threads ${task.cpus}
    """
}

process alignAndSort {
    label "aligning"

    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: "*.bam*"
    publishDir "${params.outdir}", mode: 'copy', pattern: "read_lengths.tsv"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(processed_fastq), val(output_file_count)
        file genome_reference
        file cds_reference

    output:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(processed_fastq), val(output_file_count), emit: read_info
        tuple path("read_lengths.tsv"),path("sequence_info.csv"), emit: read_metrics
        tuple path("${sample_id}.${output_file_count}.cds.bam.bai"), path("${sample_id}.${output_file_count}.cds.bam"), emit: cds_alignments
        tuple path("${sample_id}.${output_file_count}.dna.bam.bai"), path("${sample_id}.${output_file_count}.dna.bam"), emit: dna_alignments

    script:
    """
    minimap2 -t ${task.cpus} -ax splice ${genome_reference} ${processed_fastq} | samtools sort -@ ${task.cpus} | samtools view -F 2308 -o ${sample_id}.${output_file_count}.dna.bam -
    samtools index -@ ${task.cpus} ${sample_id}.${output_file_count}.dna.bam
    echo -e "sample_id\treplicate\tcondition\tcontig\tread_length\tmapping_score" > read_lengths.tsv
    samtools view ${sample_id}.${output_file_count}.dna.bam | awk -v sample_id="${sample_id}" -v replicate="${replicate}" -v condition="${condition}" \
        '{if (length(\$10) > 1) print sample_id"\t"replicate"\t"condition"\t"\$3"\t"length(\$10)"\t"\$5}' >> read_lengths.tsv

    minimap2 -t ${task.cpus} -ax splice ${cds_reference} combined.fastq.gz | samtools sort -@ ${task.cpus} | samtools view -F 2308 -o ${sample_id}.${output_file_count}.cds.bam -
    samtools index -@ ${task.cpus} ${sample_id}.${output_file_count}.cds.bam

    primary_alignments=\$(samtools view -c ${sample_id}.${output_file_count}.dna.bam)
    cds_alignments=\$(samtools view -c ${sample_id}.${output_file_count}.cds.bam)
    bases_mapped=\$(samtools view -F 4 ${sample_id}.${output_file_count}.cds.bam | awk '{sum += length(\$10)} END {print sum}')

    echo -ne "${sample_id},\$primary_alignments,\$cds_alignments,\$bases_mapped" > sequence_info.csv
    """
}

process polyATailAnalysis {
    label "polyA"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(processed_fastq), val(output_file_count)
        tuple file(bam_bai), file(bam)
        file(genome_reference)

    output:
        file("${sample_id}.${output_file_count}.polya_length.tsv")

    script:
    """
    nanopolish index -d ${fast5} -s ${fastq_dir}/sequencing_summary.txt ${processed_fastq}
    nanopolish polya --threads ${task.cpus} --reads ${processed_fastq} --bam ${bam} --genome ${genome_reference} > ${sample_id}.${output_file_count}.polya_length.tsv
    """
}

process addColumnstoPolya {
    label "smallJob"

    publishDir "${params.outdir}/polya", mode: 'copy', pattern: "*polya_length.tsv"
    publishDir "${params.outdir}/read_lengths", mode: 'copy', pattern: "*read_length.tsv"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(processed_fastq), val(output_file_count)
        file polya_tsv
        file read_info

    output:
        path("${sample_id}.${output_file_count}.polya_length.tsv"), emit: polya_length_tsv

    script:
        $/
        python3 - <<EOF
        import os
        import pandas as pd
        columns = ['readname','transcript_id','read_length','transcript_length','gene_type','sample_id','replicate','condition']
        df = pd.read_csv("${polya_tsv}",sep='\t',header=0)
        df = df[df["qc_tag"] == "PASS"]

        read_gene_df = pd.read_csv("${read_info}",sep='\t',header=0,names=columns)

        merged_df = pd.merge(df, read_gene_df, left_on='readname', right_on='readname', how='left')
        merged_df.drop(columns=['readname'], inplace=True)

        merged_df["sample_id"] = "${sample_id}"
        merged_df["replicate"] = "${replicate}"
        merged_df["condition"] = "${condition}"
        merged_df = merged_df[['sample_id','replicate','condition','transcript_id','polya_length']]
        merged_df.to_csv("${sample_id}.${output_file_count}.polya_length.tsv", index=False, sep='\t',header=False)
        EOF
        /$
}

process polyaGraphing {
    label "smallJob"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"

    input:
        file polya_tsv

    output:
        file("*.png")
    
    script:
        $/
        python3 - <<EOF
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        columns = ['sample_id','replicate','condition','reference_name','polya_length']
        df = pd.read_csv("${polya_tsv}",sep='\t',header=None,names=columns)

        selected_df = df[df["condition"] == "selected"]
        unselected_df = df[df["condition"] == "unselected"]
        selected_df = selected_df.groupby(['reference_name', 'sample_id', 'replicate', 'condition'])['polya_length'].mean().reset_index(name='average_polya_length')
        unselected_df = unselected_df.groupby(['reference_name', 'sample_id', 'replicate', 'condition'])['polya_length'].mean().reset_index(name='average_polya_length')
        merged_df = pd.merge(selected_df, unselected_df, on='reference_name', how='inner')

        plt.figure(figsize=(8, 6))

        sns.regplot(x=merged_df['average_polya_length_x'], y=merged_df['average_polya_length_y'])
        spearman_corr = merged_df['average_polya_length_x'].corr(merged_df['average_polya_length_y'], method='spearman')
        plt.text(0.1, 0.9, f'Spearman R: {spearman_corr:.2f}', transform=plt.gca().transAxes)

        plt.xlabel('Selected PolyA Length')
        plt.ylabel('Unselected PolyA Length')
        plt.title('')

        plt.savefig("polyA_scatterplot.png",format="png")
        EOF
        /$
}

process combiningTsvFiles {
    label "smallJob"

    publishDir "${params.outdir}", mode: 'copy', pattern: "experiments_summary.tsv"

    input:
        file sequence_info

    output:
        file sequence_info
        path("experiments_summary.tsv"), emit: experiments_summary
    
    script:
    $/
    python3 - <<EOF
    import os
    import pandas as pd

    columns = ['sample_id' , 'primary_alignments' , 'cds_alignments' , 'bases_mapped', 'read_counts' , 'read_lengths', 'replicate','condition']

    result_dataframe = pd.read_csv("${sequence_info}",sep=',',header=None)

    result_dataframe.columns = columns

    def combine_lists(series):
        combined_list = []
        for lst in series:
            combined_list.extend(lst)
        return combined_list

    aggregated_df = result_dataframe.groupby(['sample_id', 'replicate', 'condition'], as_index=False).agg({
        'read_lengths': combine_lists,
        'primary_alignments': 'sum',
        'cds_alignments': 'sum',
        'bases_mapped': 'sum',
        'read_counts': 'sum',
    })

    result_dataframe.to_csv("experiments_summary.tsv", index=False, sep='\t')



    EOF
    /$

}

process readLengthHistogram {
    label "smallJob"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"

    input:
        file experiments_summary

    output:
        file("*.png")
    
    script:
    $/
    python3 - <<EOF
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import ast

    df = pd.read_csv("${experiments_summary}",sep='\t',header=0)
    selected_read_lengths = ast.literal_eval(df.loc[df["condition"] == "selected", "read_lengths"].values[0])
    unselected_read_lengths = ast.literal_eval(df.loc[df["condition"] == "unselected", "read_lengths"].values[0])
    selected_read_lengths = [x for x in selected_read_lengths if x <= 7000]
    unselected_read_lengths = [x for x in unselected_read_lengths if x <= 7000]
    plt.hist(unselected_read_lengths, bins=50, alpha=0.5, label='Unselected Reads')
    plt.hist(selected_read_lengths, bins=50, alpha=0.5, label='Selected Reads')
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.title('Histogram of Read lengths Between Conditons')
    plt.legend()
    plt.savefig("Read_length_histogram.png",format="png")
    plt.close()
    EOF
    /$

}

process geneDiversityAnalysis {
    label "smallJob"
    
    publishDir "${params.outdir}", mode: 'copy', pattern: "gene_biotype_counts.tsv"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(processed_fastq), val(output_file_count)
        tuple path(bam_bai), path(bam)
        file(bed_reference)
        
    output:
        path("read_info.tsv")
    script:
        $/
        python3 - <<EOF
        import pysam
        import pandas as pd

        def extract_read_info_and_get_gene_names(bam_file):

            df = pd.DataFrame({"readname":[],"transcript_id":[],"read_length":[],"transcript_length":[]})

            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam:
                    read_id = read.query_name
                    reference_name = read.reference_name
                    reference_name = reference_name.split('.')[0]
                    read_length = read.query_length
                    gene_length = read.reference_length
                    new_row = pd.DataFrame({"readname":[read_id],"transcript_id":[reference_name], \
                                            "read_length":[read_length],"transcript_length":[gene_length]})
                    df = pd.concat([df,new_row], ignore_index=True)
            return df

        genes_df = extract_read_info_and_get_gene_names("${bam}")
        gene_types = pd.read_csv("${bed_reference}",sep='\t',header=None)

        s = gene_types.iloc[:, 9]
        gene_info = []
        for item in s:
            attributes = item.split('; ')
            gene_biotype = None
            transcript_id = None
            for attribute in attributes:
                attribute = attribute.rstrip(';')
                if 'gene_biotype' in attribute:
                    gene_biotype = attribute.split(' ')[1].replace('"', '')
                if 'transcript_id' in attribute:
                    transcript_id = attribute.split(' ')[1].replace('"', '')
            if gene_biotype and transcript_id:
                gene_info.append((transcript_id, gene_biotype))

        gene_info_df = pd.DataFrame(gene_info, columns=['transcript_id', 'gene_type']).drop_duplicates()

        merged_df = pd.merge(genes_df, gene_info_df, on='transcript_id', how='inner')
        merged_df.fillna('Unknown', inplace=True)
        merged_df["sample_id"] = "${sample_id}"
        merged_df["replicate"] = "${replicate}"
        merged_df["condition"] = "${condition}"
        merged_df.to_csv('read_info.tsv', sep='\t',header=False,index=False)
        EOF
        /$
}

process readLengthScatterplot {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${read_lengths_tsv}"

    input:
        path(read_lengths_tsv)
    output:
        file("*.png")
        path(read_lengths_tsv)

    script:
    $/
    python3 - <<EOF
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    columns = ['readname','transcript_id','read_length','transcript_length','gene_type','sample_id','replicate','condition']
    df = pd.read_csv("${read_lengths_tsv}", sep='\t', header=None,names=columns)
    selected_df = df[df["condition"] == "selected"]
    unselected_df = df[df["condition"] == "unselected"]
    selected_df = selected_df.groupby(['transcript_id', 'sample_id','condition'])['read_length'].mean().reset_index(name='average_read_length')
    unselected_df = unselected_df.groupby(['transcript_id', 'sample_id','condition'])['read_length'].mean().reset_index(name='average_read_length')
    merged_df = pd.merge(selected_df, unselected_df, on='transcript_id', how='inner')

    plt.figure(figsize=(8, 6))
    sns.regplot(x=merged_df['average_read_length_x'], y=merged_df['average_read_length_y'])
    spearman_corr = merged_df['average_read_length_x'].corr(merged_df['average_read_length_y'], method='spearman')
    plt.text(0.1, 0.9, f'Spearman R: {spearman_corr:.2f}', transform=plt.gca().transAxes)
    plt.xlabel('Selected Average Read Length')
    plt.ylabel('Unselected Average Read Length')
    plt.title('')
    plt.savefig("read_length_scatterplot.png",format="png")
    plt.close()

    counts_df = df.value_counts(subset=['condition','replicate','transcript_length','transcript_id']).reset_index(name='counts')
    for col in counts_df.select_dtypes(include=['object']).columns:
        counts_df[col] = counts_df[col].str.strip()
    counts_df["counts"] /= (counts_df["transcript_length"] / 1_000)
    read_num = counts_df.groupby(["condition","replicate"])['counts'].sum().div(1_000_000).reset_index(name='rpm_count')
    counts_df = pd.merge(counts_df, read_num, on=['condition','replicate'], how='left')
    counts_df["counts"] /= (counts_df["rpm_count"])
    counts_df = counts_df.groupby(['condition','transcript_id']).mean().reset_index()
    selected_df = counts_df[counts_df["condition"] == "selected"]
    unselected_df = counts_df[counts_df["condition"] == "unselected"]
    merged_df = pd.merge(selected_df, unselected_df, on='transcript_id', how='inner')

    plt.figure(figsize=(8, 6))
    sns.regplot(x=merged_df['counts_x'], y=merged_df['counts_y'])
    spearman_corr = merged_df['counts_x'].corr(merged_df['counts_y'], method='spearman')
    plt.text(0.1, 0.9, f'Spearman R: {spearman_corr:.2f}', transform=plt.gca().transAxes)
    plt.xlabel('Selected TPM')
    plt.ylabel('Unselected TPM')
    plt.title('')
    plt.savefig("TPM_scatterplot.png",format="png")
    plt.close()
    EOF
    /$
}

process geneMakeupGraph {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "gene_makeup_graph.png"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.tsv"

    input:
        file gene_type_coutns

    output:
        file("gene_makeup_graph.png")
        file("*.tsv")
    
    script:
        $/
        python3 - <<EOF
        import numpy as np
        import pandas as pd
        from scipy import stats
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from matplotlib import pyplot , patches
        import seaborn as sb
        import statistics

        pallete = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
            '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#ff9896', '#98df8a', 
            '#f7b6d2', '#c5b0d5', '#c49c94', '#dbdb8d', '#9edae5', '#ad494a']

        columns = ['readname','transcript_id','read_length','transcript_length','gene_type','sample_id','replicate','condition']
        merged_df = pd.read_csv("${gene_type_coutns}",sep='\t',header=None,names=columns)

        counts_df = merged_df.value_counts(subset=['gene_type','condition','replicate','transcript_length','transcript_id']).reset_index(name='counts')
        for col in counts_df.select_dtypes(include=['object']).columns:
            counts_df[col] = counts_df[col].str.strip()
        counts_df["counts"] /= (counts_df["transcript_length"] / 1_000)
        read_num = counts_df.groupby(["condition","replicate"])['counts'].sum(numeric_only=True).div(1_000_000).reset_index(name='rpm_count')
        counts_df = pd.merge(counts_df, read_num, on=['condition','replicate'], how='left')
        counts_df["counts"] /= (counts_df["rpm_count"])
        counts_df = counts_df.groupby(['gene_type','condition','replicate']).sum(numeric_only=True).reset_index()

        read_num = counts_df.groupby(["condition","replicate"])['counts'].sum().reset_index(name='total_count')
        counts_df = pd.merge(counts_df, read_num, on=['condition','replicate'], how='left')
        counts_df['counts'] = (counts_df['counts'] / counts_df['total_count']) * 100  

        counts_df.drop(columns=["total_count"],inplace=True)
        others_df = counts_df[counts_df["counts"] < 0.5]
        others_df = others_df.groupby(["condition","replicate"])['counts'].sum(numeric_only=True).reset_index(name='counts')
        others_df['gene_type'] = "Other"
        result_df = pd.concat([counts_df[counts_df['counts'] >= 0.5], others_df], ignore_index=True)
        pivot_df = result_df.pivot_table(index=["condition","replicate"], columns='gene_type', values='counts', aggfunc='sum')
        mrna_subset = pivot_df[['protein_coding']].copy()
        mrna_subset.loc[:,"Other"] = mrna_subset.loc[:,'protein_coding'].apply(lambda x: 100 - x)
        mrna_subset.reset_index(inplace=True)

        pivot_df.drop(columns="protein_coding",inplace=True)
        sorted_df = pivot_df.sort_index(axis=1,ascending=False)
        sorted_df.reset_index(inplace=True)
        sorted_df.drop(columns=['replicate'], inplace=True)
        mrna_subset.drop(columns=['replicate'], inplace=True)
        types_num = (len(sorted_df.columns) - 1) * 1.2
        fig, axes = plt.subplots(nrows=1, ncols=2,sharex=True,figsize=(6,5))
        plt.subplots_adjust(left=.1375, bottom=.2, right=0.95, top=0.9, wspace=None, hspace=None)
        mrna_subset.plot(x="condition",kind='bar',stacked=True,ax=axes[0],linewidth=0,color=["#195190FF","#A2A2A1FF"]).legend(loc="right",bbox_to_anchor=(2.975, 0.95),frameon=False)
        sorted_df.plot(x="condition",kind="bar",stacked=True,ax=axes[1],linewidth=0,color=pallete).legend(loc="right",bbox_to_anchor=(2.42, 0.625),frameon=False)
        for ax in axes:
            ax.set_xlabel("")
        fig.suptitle("")
        fig.supylabel("Percentage of Library (%)")
        plt.text(-1.5,-8,"Sample",fontsize=14.5)
        rectangle = patches.Rectangle((4.885, (18 - types_num)), .51, types_num,edgecolor='#A2A2A1FF',facecolor="#A2A2A1FF", linewidth=0,clip_on=False)
        ax.add_patch(rectangle)
        fig.savefig("gene_makeup_graph.png",format='png',bbox_inches='tight')
        mrna_subset.to_csv("gene_makeup_mrna.tsv",sep='\t')
        sorted_df.to_csv("gene_makeup_other.tsv",sep='\t')
        EOF
        /$
}

process upsetPlot {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}", mode: 'copy', pattern: "test.tsv"

    input:
        file read_gene_ids

    output:
        file("*.png")
        file "test.tsv"
    
    script:
        $/
        python3 - <<EOF
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from upsetplot import UpSet,plot,from_contents
        columns = ['readname','transcript_id','read_length','transcript_length','gene_type','sample_id','replicate','condition']

        df = pd.read_csv("${read_gene_ids}",sep='\t',header=0,names=columns)

        df.to_csv("test.tsv",sep='\t')

        contig_dict = df[['transcript_id', 'sample_id']].drop_duplicates().groupby('sample_id')['transcript_id'].agg(list).to_dict()
        binary_matrix = from_contents(contig_dict)
        ax_dict = UpSet(binary_matrix, subset_size="count", sort_by='-cardinality',show_counts=True,element_size=None).plot()
        plt.savefig("upset_plot.png", format="png")
        plt.close()
        EOF
        /$
}


workflow {
    genome_reference = file(params.genome_reference)
    cds_reference = file(params.cds_reference)
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
    alignAndSort(removeAdaptors.out,genome_reference,cds_reference)
    geneDiversityAnalysis(alignAndSort.out.read_info,alignAndSort.out.dna_alignments,bed_reference)
    readLengthandCountAnalysis(alignAndSort.out.read_info,alignAndSort.out.read_metrics)
    polyATailAnalysis(alignAndSort.out.read_info,alignAndSort.out.cds_alignments,genome_reference)
    addColumnstoPolya(alignAndSort.out.read_info,polyATailAnalysis.out,geneDiversityAnalysis.out)
    combiningTsvFiles(readLengthandCountAnalysis.out.collectFile(name: 'sequence_info.csv', newLine: true))
    readLengthHistogram(combiningTsvFiles.out.experiments_summary)
    polyaGraphing(addColumnstoPolya.out.polya_length_tsv.collectFile(name: 'all_polya_lengths.tsv', newLine: true))
    readLengthScatterplot(geneDiversityAnalysis.out.collectFile(name: 'all_read_lengths.tsv', newLine: true))
    geneMakeupGraph(geneDiversityAnalysis.out.collectFile(name: 'all_read_info.tsv', newLine: true))
    upsetPlot(geneDiversityAnalysis.out.collectFile(name: 'all_gene_names.tsv', newLine: true))
}