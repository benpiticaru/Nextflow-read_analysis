# Nextflow Read Analysis Pipeline

This pipeline processes sequencing data to analyze read lengths, polyA tail lengths, and gene makeup. It generates various outputs, including summary statistics, visualizations, and alignment files.

---

## Features

- **Read Length Analysis**: Calculates read lengths and generates histograms.
- **PolyA Tail Analysis**: Analyzes polyA tail lengths using Nanopolish.
- **Gene Makeup Analysis**: Visualizes the gene composition of the library.
- **Alignment and Sorting**: Aligns reads to a reference genome and generates BAM files.
- **Statistical Analysis**: Performs ANOVA, Tukey HSD, and Mann-Whitney U tests.
- **Scatterplots**: Generates scatterplots for read lengths and TPM values.
- **UpSet Plot**: Visualizes gene overlaps between samples.

---

## Requirements

### Software
- [Nextflow](https://www.nextflow.io/)
- Python 3.8 or higher
- Required Python libraries:
  - `pandas`
  - `numpy`
  - `seaborn`
  - `matplotlib`
  - `scipy`
  - `statsmodels`
  - `upsetplot`
- Required tools:
  - `minimap2`
  - `samtools`
  - `bedtools`
  - `seqtk`
  - `nanopolish`
  - `porechop`

### Input Files
1. **Sample Sheet**: A CSV file with the following columns:
   - `sample_id`: Unique identifier for the sample.
   - `replicate`: Replicate number.
   - `condition`: Experimental condition.
   - `fast5`: Path to the FAST5 directory.
   - `fastq_dir`: Path to the FASTQ directory.
2. **Genome Reference**: FASTA file of the reference genome.
3. **BED Reference**: BED file with gene annotations.

---

## Usage

### Running the Pipeline
1. Edit the `params` section in `read_analysis.nf` to specify the paths to your input files.
2. Run the pipeline:
   ```bash
   nextflow run read_analysis.nf
   ```

### Output
The pipeline generates the following outputs in the specified output directory:
- **Alignments**: BAM files and their indices.
- **Read Metrics**: TSV files with read length and alignment statistics.
- **PolyA Analysis**: Boxplots and statistical summaries for polyA tail lengths.
- **Gene Makeup**: Bar plots showing the gene composition of the library.
- **Histograms**: Read length histograms by condition.
- **Scatterplots**: Scatterplots for read lengths and TPM values.
- **UpSet Plot**: Visualizes gene overlaps between samples.

---

## Example Input and Output

### Example Sample Sheet (CSV)
```csv
sample_id,replicate,condition,fast5,fastq_dir
sample1,1,control,/path/to/fast5,/path/to/fastq
sample2,1,treatment,/path/to/fast5,/path/to/fastq
```

### Example Output
- `alignments/`: Contains BAM files and indices.
- `graphs/`: Contains visualizations (e.g., histograms, scatterplots, boxplots, UpSet plots).
- `tables/`: Contains statistical summaries (e.g., ANOVA, Mann-Whitney U test results, gene makeup tables).

---

## Python Script Integration

The pipeline uses a Python script (`scripts/read_analysis.py`) to handle various data processing and visualization tasks. The script includes the following subcommands:
- `read_length_and_count_analysis`: Calculates read lengths and counts.
- `combine_tsv_files`: Combines sequencing metrics into summary TSV files.
- `generate_read_length_histogram`: Creates histograms and statistical summaries for read lengths.
- `generate_read_length_scatterplot`: Generates scatterplots for read lengths and TPM values.
- `generate_polya_graphing`: Produces boxplots and ANOVA results for polyA tail lengths.
- `generate_gene_makeup_graph`: Creates bar plots for gene composition.
- `generate_upset_plot`: Generates UpSet plots for gene overlaps.
- `process_align_and_sort`: Processes alignment data and generates read metrics.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Contact

For questions or support, please contact the project maintainer.