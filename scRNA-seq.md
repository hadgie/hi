# Single Cell RNA 101

## welcome to my workflow
I will go over the tools & commands I used to pre-process, cluster single cell RNA data.
Alignment was done by `cellranger count` from the `cellranger` pipeline. You can follow the basic steps for installing cellranger from their website. [Click here to visit.](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in)


## 1. Data align using `cellranger count`

Before running `cellranger count`, ***you need to make sure the file name is fixed like the following:***

```
name_S0_L000_R0_000.fastq.gz

# example:
pbmc_1k_v3_S1_L001_R2_001.fastq.gz
pbmc_1k_v3_S1_L002_I1_001.fastq.gz
pbmc_1k_v3_S1_L001_R1_001.fastq.gz
pbmc_1k_v3_S1_L002_R1_001.fastq.gz
pbmc_1k_v3_S1_L002_R2_001.fastq.gz
pbmc_1k_v3_S1_L001_I1_001.fastq.gz
```



`cellranger count` can be run by the following commands:
```
cellranger count   --id=SRR19882110 \ 
                   --transcriptome=/home/harold/yard/reference/refdata-gex-mm10-2020-A.tar.gz \
                   --fastqs=/home/harold/sras \
                   --sample=mysample \
                   --localcores=16 \
                   --localmem=64
```



`--id` will be the name of the sample you will be running. `--transcriptome` is the reference genome(can download it from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)), `--fastqs` is the directory where the raw data is located. `--sample` is the name of the project, `--localcores` and `--localmem` are computational resource assigned to the tool.

The output directory of `cellranger count` will be located where the input file is located with following data:
|File/Folder |	Description|
|------------|-------------|
|analysis/ |	Folder containing the results of graph-based clusters and K-means clustering 2-10; differential gene expression analysis between clusters; and PCA, t-SNE, and UMAP dimensionality reduction. [Learn more](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/analysis)|
|cloupe.cloupe |	A [Loupe Browser](https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/tutorial) visualization and analysis file.|
|filtered_feature_bc_matrix/ | 	Contains only detected cell-associated barcodes in MEX format. Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column), as described in the [feature-barcode matrix page](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices). This file can be input into third-party packages and allows users to wrangle the barcode-feature matrix (e.g. to filter outlier cells, run dimensionality reduction, normalize gene expression).|
|filtered_feature_bc_matrix.h5 |	Same information as `filtered_feature_bc_matrix` in [HDF5 format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices).|
|metrics_summary.csv |	Run summary metrics file in CSV format, described in the [Gene Expression metrics page](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/gex-metrics).|
|molecule_info.h5 |	Contains per-molecule information for all molecules that contain a valid barcode, valid UMI, and were assigned with high confidence to a gene or Feature Barcode. This file is a required input to run `cellranger aggr`. [Learn more](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info) |
|possorted_genome_bam.bam |	Indexed BAM file containing position-sorted reads aligned to the genome and transcriptome, as well as unaligned reads, annotated with barcode information. [Learn more](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam)|
|possorted_genome_bam.bam.bai |	Index file for the `possorted_genome_bam.bam`. In cases where the reference transcriptome is generated from a genome with very long chromosomes (>512 Mbp), Cell Ranger v7.0+ generates a `possorted_genome_bam.bam.csi` index file instead.|
|raw_feature_bc_matrix/ |	Contains all detected barcodes in MEX format. Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column), as described in the [feature-barcode matrix page](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices).|
|raw_feature_bc_matrix_h5 |	Same information as `raw_feature_bc_matrix` in [HDF5 format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices). |
|web_summary.html |	Run summary metrics and charts in HTML format, described in the `count` [web summary page](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/summary).|

The processed matrix data will be in `outs` directory.

## 2. Matrix data processing using `Seurat`

After completing the alignment process, 3 files are needed for the next step: `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`.
**I strongly recommend to run `seurat` from local, so that you can see plots from each step.**
The following steps are based on the vignette of Satija Lab. [Click here to visit.](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

You can see my script [here](https://github.com/hoonbiolab/HL-SCRNA/blob/main/sandbox/HS/Seurat.R).
