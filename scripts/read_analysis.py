import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.stats import f_oneway, mannwhitneyu
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import seaborn as sns
import matplotlib.pyplot as plt
from ast import literal_eval


def calculate_read_lengths(file_name):
    read_lengths = []
    with open(file_name, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_lengths.append(len(record.seq))
    return read_lengths


def read_length_and_count_analysis(processed_fastq, sequence_info, replicate, condition, output_file):
    read_lengths = calculate_read_lengths(processed_fastq)
    read_counts = len(read_lengths)

    df = pd.read_csv(sequence_info, sep=",", header=None)
    df["read_counts"] = read_counts
    df["read_lengths"] = [read_lengths]
    df["replicate"] = replicate
    df["condition"] = condition

    df.to_csv(output_file, sep=",", header=None, index=False)


def combine_tsv_files(sequence_info, output_summary, output_experiments):
    columns = ['sample_id', 'primary_alignments', 'cds_alignments', 'bases_mapped', 'read_counts', 'read_lengths', 'replicate', 'condition']
    result_dataframe = pd.read_csv(sequence_info, sep=',', header=None, names=columns)

    def combine_lists(series):
        combined_list = []
        for lst in series:
            combined_list.extend(literal_eval(lst))
        return combined_list

    aggregated_df = result_dataframe.groupby(['sample_id', 'replicate', 'condition'], as_index=False).agg({
        'read_lengths': combine_lists,
        'primary_alignments': 'sum',
        'cds_alignments': 'sum',
        'bases_mapped': 'sum',
        'read_counts': 'sum',
    })

    final_df = aggregated_df[['sample_id', 'primary_alignments', 'cds_alignments', 'bases_mapped', 'read_counts', 'replicate', 'condition']]
    aggregated_df.to_csv(output_summary, index=False, sep='\t')
    final_df.to_csv(output_experiments, index=False, sep='\t')


def generate_read_length_histogram(experiments_summary, output_histogram, output_stats, output_test):
    df = pd.read_csv(experiments_summary, sep='\t', header=0)
    df['read_lengths'] = df['read_lengths'].apply(literal_eval)

    combined_read_lengths = []
    for _, row in df.iterrows():
        condition = row['condition']
        read_lengths = row['read_lengths']
        for length in read_lengths:
            combined_read_lengths.append({'condition': condition, 'read_length': length})

    combined_read_lengths_df = pd.DataFrame(combined_read_lengths)
    combined_read_lengths_df = combined_read_lengths_df[combined_read_lengths_df['read_length'] <= 7000]

    summary_stats = combined_read_lengths_df.groupby('condition')['read_length'].describe()
    summary_stats.to_csv(output_stats, sep='\t', index=True)

    plt.figure(figsize=(10, 6))
    sns.histplot(data=combined_read_lengths_df, x='read_length', hue='condition', bins=50)
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.title('')
    plt.savefig(output_histogram, format="png")
    plt.close()

    selected_lengths = combined_read_lengths_df[combined_read_lengths_df['condition'] == 'selected']['read_length']
    unselected_lengths = combined_read_lengths_df[combined_read_lengths_df['condition'] == 'unselected']['read_length']
    statistic, p_value = mannwhitneyu(unselected_lengths, selected_lengths, alternative='greater')

    test_result_df = pd.DataFrame({"Statistic": [statistic], "P-value": [p_value]})
    test_result_df.to_csv(output_test, sep='\t', index=False)


def generate_read_length_scatterplot(read_lengths_tsv, output_scatterplot, output_tpm_scatterplot):
    def log_mean_read_length(df):
        length_df = df.groupby(['gene_id', 'condition'])['read_length'].mean().reset_index(name='average_read_length')
        length_df['average_read_length'] = np.log(length_df['average_read_length'])
        pivoted_df = length_df.pivot(index='gene_id', columns='condition', values='average_read_length').dropna()
        return pivoted_df

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

    def plot_scatter(df, x_col, y_col, xlabel, ylabel, filename):
        plt.figure(figsize=(8, 6))
        plt.scatter(df[x_col], df[y_col], color="black", s=2, label='_nolegend_')
        slope, intercept = np.polyfit(df[x_col], df[y_col], 1)
        plt.plot(df[x_col], slope * df[x_col] + intercept, color='#1f77b4', label='Regression Line')
        spearman_corr = df[x_col].corr(df[y_col], method='spearman')
        plt.text(0.03, 0.875, f'Spearman R: {spearman_corr:.2f}', transform=plt.gca().transAxes, va='top', ha='left')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.gca().set_aspect('equal', adjustable='box')
        max_val = max(df[x_col].max(), df[y_col].max())
        min_val = min(df[x_col].min(), df[y_col].min())
        plt.plot([min_val, max_val], [min_val, max_val], color='red', label='x = y')
        plt.legend(loc='upper left')
        plt.savefig(filename, format="png", bbox_inches='tight')
        plt.close()

    columns = ['gene_id', 'gene_biotype', 'readname', 'gene_length', 'read_length', 'sample_id', 'replicate', 'condition']
    df = pd.read_csv(read_lengths_tsv, sep='\t', header=None, names=columns)

    length_df = log_mean_read_length(df)
    plot_scatter(length_df, 'unselected', 'selected',
                 'Unselected Average Read Length (log)', 'Selected Average Read Length (log)',
                 output_scatterplot)

    tpm_df = log_tpm(df)
    plot_scatter(tpm_df, 'unselected', 'selected',
                 'Unselected TPM (log)', 'Selected TPM (log)',
                 output_tpm_scatterplot)


def generate_polya_graphing(polya_tsv, output_boxplot, output_anova):
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats import f_oneway
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    def find_limits(data):
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        lwr = q1 - 1.5 * iqr
        upr = q3 + 1.5 * iqr
        return lwr, upr

    df = pd.read_csv(polya_tsv, sep='\t', header=None, names=['sample_id', 'replicate', 'condition', "readname", "tail_length"])
    pivoted_df = df.pivot_table(index='readname', columns='sample_id', values='tail_length')

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

    with open(output_anova, "w") as file:
        file.write(f"ANOVA p-value: {anova_pvalue}\n\n")
        file.write("Post-hoc Tukey HSD Test:\n")
        file.write(str(posthoc_output))

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=pivoted_df, fill=None, color="black")
    plt.xlabel('Sample')
    plt.ylabel('Tail Length')
    plt.tight_layout()
    plt.savefig(output_boxplot, format="png")


def generate_gene_makeup_graph(read_gene_ids, output_graph, output_mrna, output_other):
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib import patches

    pallete = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
               '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#ff9896', '#98df8a',
               '#f7b6d2', '#c5b0d5', '#c49c94', '#dbdb8d', '#9edae5', '#ad494a']

    columns = ['gene_id', 'gene_biotype', 'readname', 'gene_length', 'read_length', 'sample_id', 'replicate', 'condition']
    merged_df = pd.read_csv(read_gene_ids, sep='\t', header=None, names=columns)

    counts_df = merged_df.value_counts(subset=['gene_biotype', 'condition', 'replicate', 'gene_length', 'gene_id']).reset_index(name='counts')
    counts_df["counts"] /= counts_df.groupby(["condition", "replicate"])['counts'].transform('sum')
    counts_df = counts_df.groupby(['gene_biotype', 'condition', 'replicate']).sum(numeric_only=True).reset_index()

    others_df = counts_df[counts_df["counts"] < 0.005]
    others_df = others_df.groupby(["condition", "replicate"])['counts'].sum().reset_index(name='counts')
    others_df['gene_biotype'] = "Other"
    result_df = pd.concat([counts_df[counts_df['counts'] >= 0.005], others_df], ignore_index=True)

    pivot_df = result_df.pivot_table(index=["condition", "replicate"], columns='gene_biotype', values='counts', aggfunc='sum')
    mrna_subset = pivot_df[['protein_coding']].copy()
    mrna_subset.loc[:, "Other"] = 1 - mrna_subset['protein_coding']
    mrna_subset.reset_index(inplace=True)

    pivot_df.drop(columns="protein_coding", inplace=True)
    sorted_df = pivot_df.sort_index(axis=1, ascending=False)
    sorted_df.reset_index(inplace=True)

    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(10, 5))
    mrna_subset.plot(x="condition", kind='bar', stacked=True, ax=axes[0], linewidth=0, color=["#195190FF", "#A2A2A1FF"])
    sorted_df.plot(x="condition", kind="bar", stacked=True, ax=axes[1], linewidth=0, color=pallete)

    for ax in axes:
        ax.set_xlabel("")
    fig.supylabel("Percentage of Library (%)")
    plt.tight_layout()
    plt.savefig(output_graph, format='png')

    mrna_subset.to_csv(output_mrna, sep='\t', index=False)
    sorted_df.to_csv(output_other, sep='\t', index=False)


def generate_upset_plot(read_gene_ids, output_plot):
    import pandas as pd
    import matplotlib.pyplot as plt
    from upsetplot import from_contents, UpSet

    columns = ['gene_id', 'gene_biotype', 'readname', 'gene_length', 'read_length', 'sample_id', 'replicate', 'condition']
    df = pd.read_csv(read_gene_ids, sep='\t', header=0, names=columns)

    contig_dict = df[['gene_id', 'sample_id']].drop_duplicates().groupby('sample_id')['gene_id'].agg(list).to_dict()
    binary_matrix = from_contents(contig_dict)

    ax_dict = UpSet(binary_matrix, subset_size="count", sort_by='-cardinality', show_counts=True, element_size=None).plot()
    plt.savefig(output_plot, format="png")
    plt.close()


def process_align_and_sort(gene_id_biotype, gene_intersects, sample_id, replicate, condition, output_metrics, output_sequence_info):
    import pandas as pd

    # Read and merge data
    df1 = pd.read_csv(gene_id_biotype, sep='\t', header=None, names=["gene_id", "gene_type"])
    df2 = pd.read_csv(gene_intersects, sep='\t', header=None, names=["readname", "gene_id", "gene_length", "read_length"])
    merged_df = pd.merge(df1, df2, on="gene_id")

    # Add metadata columns
    merged_df["sample_id"] = sample_id
    merged_df["replicate"] = replicate
    merged_df["condition"] = condition

    # Save read metrics
    merged_df.to_csv(output_metrics, sep="\t", header=False, index=False)

    # Calculate summary statistics
    primary_alignments = len(merged_df)
    cds_alignments = len(merged_df[merged_df["gene_type"] == "protein_coding"])
    bases_mapped = merged_df[merged_df["gene_type"] == "protein_coding"]["read_length"].sum()

    # Save sequence info
    with open(output_sequence_info, "w") as f:
        f.write(f"{sample_id},{primary_alignments},{cds_alignments},{bases_mapped}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Read Analysis Script")
    subparsers = parser.add_subparsers(dest="command")

    # Subcommand: read_length_and_count_analysis
    parser_read_length = subparsers.add_parser("read_length_and_count_analysis")
    parser_read_length.add_argument("--processed_fastq", required=True)
    parser_read_length.add_argument("--sequence_info", required=True)
    parser_read_length.add_argument("--replicate", required=True)
    parser_read_length.add_argument("--condition", required=True)
    parser_read_length.add_argument("--output_file", required=True)

    # Subcommand: combine_tsv_files
    parser_combine = subparsers.add_parser("combine_tsv_files")
    parser_combine.add_argument("--sequence_info", required=True)
    parser_combine.add_argument("--output_summary", required=True)
    parser_combine.add_argument("--output_experiments", required=True)

    # Subcommand: generate_read_length_histogram
    parser_histogram = subparsers.add_parser("generate_read_length_histogram")
    parser_histogram.add_argument("--experiments_summary", required=True)
    parser_histogram.add_argument("--output_histogram", required=True)
    parser_histogram.add_argument("--output_stats", required=True)
    parser_histogram.add_argument("--output_test", required=True)

    # Subcommand: generate_read_length_scatterplot
    parser_scatterplot = subparsers.add_parser("generate_read_length_scatterplot")
    parser_scatterplot.add_argument("--read_lengths_tsv", required=True)
    parser_scatterplot.add_argument("--output_scatterplot", required=True)
    parser_scatterplot.add_argument("--output_tpm_scatterplot", required=True)

    # Subcommand: generate_polya_graphing
    parser_polya = subparsers.add_parser("generate_polya_graphing")
    parser_polya.add_argument("--polya_tsv", required=True)
    parser_polya.add_argument("--output_boxplot", required=True)
    parser_polya.add_argument("--output_anova", required=True)

    # Subcommand: generate_gene_makeup_graph
    parser_gene_makeup = subparsers.add_parser("generate_gene_makeup_graph")
    parser_gene_makeup.add_argument("--read_gene_ids", required=True)
    parser_gene_makeup.add_argument("--output_graph", required=True)
    parser_gene_makeup.add_argument("--output_mrna", required=True)
    parser_gene_makeup.add_argument("--output_other", required=True)

    # Subcommand: generate_upset_plot
    parser_upset = subparsers.add_parser("generate_upset_plot")
    parser_upset.add_argument("--read_gene_ids", required=True)
    parser_upset.add_argument("--output_plot", required=True)

    # Subcommand: process_align_and_sort
    parser_align_sort = subparsers.add_parser("process_align_and_sort")
    parser_align_sort.add_argument("--gene_id_biotype", required=True)
    parser_align_sort.add_argument("--gene_intersects", required=True)
    parser_align_sort.add_argument("--sample_id", required=True)
    parser_align_sort.add_argument("--replicate", required=True)
    parser_align_sort.add_argument("--condition", required=True)
    parser_align_sort.add_argument("--output_metrics", required=True)
    parser_align_sort.add_argument("--output_sequence_info", required=True)

    args = parser.parse_args()

    if args.command == "read_length_and_count_analysis":
        read_length_and_count_analysis(args.processed_fastq, args.sequence_info, args.replicate, args.condition, args.output_file)
    elif args.command == "combine_tsv_files":
        combine_tsv_files(args.sequence_info, args.output_summary, args.output_experiments)
    elif args.command == "generate_read_length_histogram":
        generate_read_length_histogram(args.experiments_summary, args.output_histogram, args.output_stats, args.output_test)
    elif args.command == "generate_read_length_scatterplot":
        generate_read_length_scatterplot(args.read_lengths_tsv, args.output_scatterplot, args.output_tpm_scatterplot)
    elif args.command == "generate_polya_graphing":
        generate_polya_graphing(args.polya_tsv, args.output_boxplot, args.output_anova)
    elif args.command == "generate_gene_makeup_graph":
        generate_gene_makeup_graph(args.read_gene_ids, args.output_graph, args.output_mrna, args.output_other)
    elif args.command == "generate_upset_plot":
        generate_upset_plot(args.read_gene_ids, args.output_plot)
    elif args.command == "process_align_and_sort":
        process_align_and_sort(args.gene_id_biotype, args.gene_intersects, args.sample_id, args.replicate, args.condition, args.output_metrics, args.output_sequence_info)


if __name__ == "__main__":
    main()
