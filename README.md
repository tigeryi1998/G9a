# G9a

## EHMT2 (G9a) Breast Cancer Isoforms Differential Gene Expression Analysis


## Structure of G9a project folder

```
G9a/
├── scripts/               # All bash / qsub / pipeline scripts
│   ├── run_star_alignment.txt
│   ├── run_featureCounts.txt
│   ├── run_build_gene_map.txt
│   ├── build_star_index.txt
│   └── trimmomatic.txt
│
├── data/                 # Generated intermediate data (NOT raw FASTQ)
│   ├── gene_counts.txt
│   ├── gene_counts_cleaned.txt
│   ├── gene_counts.txt.summary
│   ├── ensembl_gene_map.csv
│   ├── metadata.csv
│   ├── samples.txt
│   ├── splicing_isoform_signals.txt
│   ├── signatures.csv
│   └── ...
│
├── results/              # DESeq2 and downstream analysis results
│   ├── rescue.csv
│   ├── all_AFEPS.csv
│   ├── DE_CSV/           # Subfolder for multiple DE output CSVs
│   └── ...
│
├── figures/              # Output plots
│   ├── Heatmap_log2FC.png
│   ├── EHMT2_log2FC.png
│   ├── volcano_plots/
│   └── ...
│
├── notes/                # Project documentation, description, etc.
│   ├── notes.txt
│   ├── description.txt
│   └── qsub.txt          
│
├── README.md
├── LICENSE.md
├── .gitignore
├── DE.R                 # R script for differential gene expression analysis
```

