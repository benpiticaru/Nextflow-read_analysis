#!/usr/bin/env nextflow
params.sample_sheet = "/{path_to}/sample_sheet.csv"
params.outdir = "exp_4.out"
params.genome_reference = "/{path_to}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.bed_reference = "/{path_to}/Homo_sapiens_GRCh38_110.bed"

process readLengthandCountAnalysis {
    label "smallJob"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count)
        file(sequence_info)


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

process alignAndSort {
    label "aligning"

    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: "*.bam*"

    input:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count)
        file genome_reference
        file reference_bed

    output:
        tuple val(sample_id) , val(replicate) , val(condition), path(fast5), path(fastq_dir), path(raw_fastq), path(processed_fastq), val(output_file_count), emit: read_info
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

    python3 - <<EOF
    import pandas as pd
    df1 = pd.read_csv("gene_id_biotype.tsv",sep='\t',header=None,names=["gene_id","gene_type"])
    df2 = pd.read_csv("gene_intersects.bed",sep='\t',header=None,names=["readname","gene_id","gene_length","read_length"])
    merged_df = pd.merge(df1, df2, on="gene_id")
    merged_df["sample_id"] = "${sample_id}"
    merged_df["replicate"] = "${replicate}"
    merged_df["condition"] = "${condition}"
    merged_df.to_csv("${sample_id}.${output_file_count}.read_metrics.tsv",sep="\t",header=False,index=False)
    EOF

    primary_alignments=\$(samtools view -c ${sample_id}.${output_file_count}.dna.bam)
    cds_alignments=\$(grep -o "protein_coding" ${sample_id}.${output_file_count}.read_metrics.tsv | wc -l)
    bases_mapped=\$(grep -w "protein_coding" ${sample_id}.${output_file_count}.read_metrics.tsv | awk '{ sum += \$5 } END { print sum }')

    echo -ne "${sample_id},\$primary_alignments,\$cds_alignments,\$bases_mapped" > sequence_info.csv
    """
}

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
        $/
        python3 - <<EOF
        import pandas as pd
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        from ast import literal_eval
        from scipy.stats import f_oneway
        from statsmodels.stats.multicomp import pairwise_tukeyhsd

        def find_limits(data):
                q1, q3 = np.percentile(data, [25, 75])
                iqr = q3 - q1
                lwr = q1 - 1.5 * iqr
                upr = q3 + 1.5 * iqr
                return lwr , upr

        df = pd.read_csv("${polya_tsv}",sep='\t',header=None,names=['sample_id','replicate','condition',"readname","tail_length"])
        pivoted_df = df.pivot_table(index='readname',columns='sample_id',values='tail_length')

        for sample in pivoted_df.columns:
            lwr, upr = find_limits(pivoted_df[sample].dropna())
            pivoted_df = pivoted_df[(pivoted_df[sample] > lwr) & (pivoted_df[sample] <= upr) | pivoted_df[sample].isna()]

        data = []
        groups = []

        for column in pivoted_df.columns:
            column_data = pivoted_df[column].dropna().tolist()
            data.extend(column_data)
            groups.extend([column] * len(column_data))

        anova_result = f_oneway(*[pivoted_df[column].dropna().tolist() for column in pivoted_df.columns])
        anova_pvalue = anova_result.pvalue

        if anova_pvalue < 0.05:
            posthoc = pairwise_tukeyhsd(data, groups, alpha=0.05)
            posthoc_output = posthoc.summary()
        else:
            posthoc_output = "No significant difference between groups."

        with open("anova_tukey_output.txt", "w") as file:
            file.write(f"ANOVA p-value: {anova_pvalue}\n\n")
            file.write("Post-hoc Tukey HSD Test:\n")
            file.write(str(posthoc_output))

        plt.figure(figsize=(10, 6))
        sns.boxplot(data=pivoted_df, fill=None,color="black")
        plt.xlabel('Sample')
        plt.ylabel('Tail Length')
        plt.title('')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig("PolyA_length_boxplot.png", format="png")
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
        path("exp_summary.tsv"), emit: experiments_summary
        path("experiments_summary.tsv")
    
    script:
    $/
    python3 - <<EOF
    import os
    import pandas as pd
    import ast

    columns = ['sample_id' , 'primary_alignments' , 'cds_alignments' , 'bases_mapped', 'read_counts' , 'read_lengths', 'replicate','condition']
    result_dataframe = pd.read_csv("${sequence_info}",sep=',',header=None,names=columns)
    
    def combine_lists(series):
        combined_list = []
        for lst in series:
            combined_list.extend(ast.literal_eval(lst))
        return combined_list

    aggregated_df = result_dataframe.groupby(['sample_id', 'replicate', 'condition'], as_index=False).agg({
        'read_lengths': combine_lists,
        'primary_alignments': 'sum',
        'cds_alignments': 'sum',
        'bases_mapped': 'sum',
        'read_counts': 'sum',
    })

    final_df = aggregated_df[['sample_id' , 'primary_alignments' , 'cds_alignments' , 'bases_mapped', 'read_counts', 'replicate','condition']]
    aggregated_df.to_csv("exp_summary.tsv", index=False, sep='\t')
    final_df.to_csv("experiments_summary.tsv", index=False, sep='\t')
    EOF
    /$

}

process readLengthHistogram {
    label "smallJob"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "*.png"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.tsv"

    input:
        file experiments_summary

    output:
        file("*.png")
        file("mann_whitney_u_test_result.tsv")
        file("summary_statistics.tsv")
    
    script:
    $/
    python3 - <<EOF
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from ast import literal_eval
    from scipy.stats import mannwhitneyu

    df = pd.read_csv("${experiments_summary}", sep='\t', header=0)
    df['read_lengths'] = df['read_lengths'].apply(literal_eval)

    combined_read_lengths = []

    for index, row in df.iterrows():
        condition = row['condition']
        read_lengths = row['read_lengths']
        for length in read_lengths:
            combined_read_lengths.append({'condition': condition, 'read_length': length})

    combined_read_lengths_df = pd.DataFrame(combined_read_lengths)

    combined_read_lengths_df = combined_read_lengths_df[combined_read_lengths_df['read_length'] <= 7000]

    summary_stats = combined_read_lengths_df.groupby('condition')['read_length'].describe()
    summary_stats.to_csv("summary_statistics.tsv", sep='\t',index=True)
    print("Summary Statistics:")
    print(summary_stats)


    plt.figure(figsize=(10, 6))
    sns.histplot(data=combined_read_lengths_df, x='read_length', hue='condition', bins=50)
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.title('')

    selected_lengths = combined_read_lengths_df[combined_read_lengths_df['condition'] == 'selected']['read_length']
    unselected_lengths = combined_read_lengths_df[combined_read_lengths_df['condition'] == 'unselected']['read_length']
    statistic, p_value = mannwhitneyu(unselected_lengths, selected_lengths,alternative='greater')

    print("Mann-Whitney U Test:")
    print(f"Statistic: {statistic}")
    print(f"P-value: {p_value}")

    test_result_df = pd.DataFrame({"Statistic": [statistic], "P-value": [p_value]})
    test_result_df.to_csv("mann_whitney_u_test_result.tsv", sep='\t', index=False)
    plt.savefig("Read_length_histogram_by_condition.png", format="png")
    plt.close()
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
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    def log_mean_read_length(df):
        length_df = df.groupby(['gene_id','condition'])['read_length'].mean().reset_index(name='average_read_length')
        length_df['average_read_length'] = np.log(length_df['average_read_length'])
        pivoted_df = length_df.pivot(index='gene_id', columns='condition', values='average_read_length').dropna()
        return pivoted_df

    def plot_scatter(df, x_col, y_col, xlabel, ylabel, title, filename):
        plt.figure(figsize=(8, 6))
        plt.scatter(df[x_col], df[y_col], color="black", s=2, label='_nolegend_')  
        slope, intercept = np.polyfit(df[x_col], df[y_col], 1) 
        plt.plot(df[x_col], slope * df[x_col] + intercept, color='#1f77b4', label='Regression Line')  
        spearman_corr = df[x_col].corr(df[y_col], method='spearman')
        plt.text(0.03, 0.875, f'Spearman R: {spearman_corr:.2f}', transform=plt.gca().transAxes, va='top', ha='left') 
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title('')
        plt.gca().set_aspect('equal', adjustable='box')
        max_val = max(df[x_col].max(), df[y_col].max())
        min_val = min(df[x_col].min(), df[y_col].min())
        plt.plot([min_val, max_val], [min_val, max_val], color='red', label='x = y')
        plt.legend(loc='upper left') 
        plt.savefig(filename, format="png", bbox_inches='tight')
        plt.close()


    def log_tpm(df):
        counts_df = df.value_counts(subset=['condition', 'gene_length', 'gene_id']).reset_index(name='counts')
        counts_df = counts_df[counts_df["counts"] >= 5]
        counts_df["counts_cpk"] = counts_df["counts"] / (counts_df["gene_length"] / 1000)
        rpm_count = counts_df.groupby(["condition"])['counts_cpk'].sum().div(1_000_000).reset_index(name='rpm_count')
        counts_df = pd.merge(counts_df, rpm_count, on=['condition'], how='left')
        counts_df["tpm"] = (counts_df["counts_cpk"] / counts_df["rpm_count"])
        counts_df["log_tpm"] = np.log(counts_df["tpm"])
        counts_df.replace(-np.inf, np.nan, inplace=True)
        pivoted_df = counts_df.pivot(index='gene_id', columns='condition', values='log_tpm').dropna()
        return pivoted_df

    def main():
        columns = ['gene_id','gene_biotype','readname','gene_length','read_length','sample_id','replicate','condition']
        df = pd.read_csv("${read_lengths_tsv}", sep='\t', header=None,names=columns)

        length_df = log_mean_read_length(df)
        plot_scatter(length_df, 'unselected', 'selected', 
                    'Unselected Average Read Length (log)', 'Selected Average Read Length (log)', 
                    'Read Length Scatterplot', 'read_length_scatterplot.png')

        tpm_df = log_tpm(df)
        plot_scatter(tpm_df, 'unselected', 'selected', 
                    'Unselected TPM (log)', 'Selected TPM (log)', 
                    'TPM Scatterplot', 'TPM_scatterplot.png')

    if __name__ == "__main__":
        main()
    EOF
    /$
}

process geneMakeupGraph {
    label "python"

    publishDir "${params.outdir}/graphs", mode: 'copy', pattern: "gene_makeup_graph.png"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.tsv"

    input:
        file read_gene_ids

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

        columns = ['gene_id','gene_biotype','readname','gene_length','read_length','sample_id','replicate','condition']
        merged_df = pd.read_csv("${read_gene_ids}",sep='\t',header=None,names=columns)

        counts_df = merged_df.value_counts(subset=['gene_biotype','condition','replicate','gene_length','gene_id']).reset_index(name='counts')
        for col in counts_df.select_dtypes(include=['object']).columns:
            counts_df[col] = counts_df[col].str.strip()
        read_num = counts_df.groupby(["condition","replicate"])['counts'].sum(numeric_only=True).div(1_000_000).reset_index(name='rpm_count')
        counts_df = pd.merge(counts_df, read_num, on=['condition','replicate'], how='left')
        counts_df["counts"] /= (counts_df["rpm_count"])
        counts_df = counts_df.groupby(['gene_biotype','condition','replicate']).sum(numeric_only=True).reset_index()

        read_num = counts_df.groupby(["condition","replicate"])['counts'].sum().reset_index(name='total_count')
        counts_df = pd.merge(counts_df, read_num, on=['condition','replicate'], how='left')
        counts_df['counts'] = (counts_df['counts'] / counts_df['total_count']) * 100  

        counts_df.drop(columns=["total_count"],inplace=True)
        others_df = counts_df[counts_df["counts"] < 0.5]
        others_df = others_df.groupby(["condition","replicate"])['counts'].sum(numeric_only=True).reset_index(name='counts')
        others_df['gene_biotype'] = "Other"
        result_df = pd.concat([counts_df[counts_df['counts'] >= 0.5], others_df], ignore_index=True)
        pivot_df = result_df.pivot_table(index=["condition","replicate"], columns='gene_biotype', values='counts', aggfunc='sum')
        mrna_subset = pivot_df[['protein_coding']].copy()
        mrna_subset.loc[:,"Other:"] = mrna_subset.loc[:,'protein_coding'].apply(lambda x: 100 - x)
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
        sorted_df.plot(x="condition",kind="bar",stacked=True,ax=axes[1],linewidth=0,color=pallete).legend(loc="right",bbox_to_anchor=(1.75, 0.8),frameon=False)
        for ax in axes:
            ax.set_xlabel("")
        fig.suptitle("")
        fig.supylabel("Percentage of Library (%)")
        plt.text(-1.5,-5,"Sample",fontsize=14.5)
        rectangle = patches.Rectangle((3.8, (16 - types_num)), .51, types_num,edgecolor='#A2A2A1FF',facecolor="#A2A2A1FF", linewidth=0,clip_on=False)
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

    input:
        file read_gene_ids

    output:
        file("*.png")
    
    script:
        $/
        python3 - <<EOF
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from upsetplot import UpSet,plot,from_contents
        columns = ['gene_id','gene_biotype','readname','gene_length','read_length','sample_id','replicate','condition']
        df = pd.read_csv("${read_gene_ids}",sep='\t',header=0,names=columns)

        contig_dict = df[['gene_id', 'sample_id']].drop_duplicates().groupby('sample_id')['gene_id'].agg(list).to_dict()
        binary_matrix = from_contents(contig_dict)
        ax_dict = UpSet(binary_matrix, subset_size="count", sort_by='-cardinality',show_counts=True,element_size=None).plot()
        plt.savefig("upset_plot.png", format="png")
        plt.close()
        EOF
        /$
}


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