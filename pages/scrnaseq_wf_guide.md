---
title: scRNAseq Processing Workflows
type: guide
contributors: 
description: How-to guide for scRNAseq workflows on Galaxy Australia
affiliations: QCIF
toc: false
---

- [Analysis overview](#analysis-overview)
- [Example output](#example-output)
- [User guide](#user-guide)
    + [Running a single sample workflow](#running-a-single-sample-workflow)
    + [Running a multi sample experiment](#running-a-multi-sample-experiment)
    + [Next steps](#next-steps)
- [Background and Tutorials](#background-and-tutorials)

This document describes how to use some scanpy-based scRNAseq workflows on galaxy Australia. 

The aim of these workflows is to handle the routine ‘boring’ part of single cell RNAseq data processing. It will produces an ‘AnnData’ object, which can then be used as a base for downstream analysis – either within galaxy or outside of it. AnnData is a standard format used by the ‘scanpy’ python package. 

These workflows represent just one way of processing data for a ‘typical’ scRNAseq experiment – there are many other options!  



This document describes 3 sub-workflows for processing single cell RNAseq data with scanpy 

* **Load counts matrix:** [link](https://usegalaxy.org.au/u/s.williams/w/scrnaseq-load-counts-matrix-subworkflow)   This workflow adds a sample name, which enables multi-sample analyses 
* **Single cell QC:** [link](https://usegalaxy.org.au/u/s.williams/w/scrnaseqcellqc)  This workflow generates some basic QC plots and applies filtering 
* **Single cell QC to basic processing:** [link](https://usegalaxy.org.au/u/s.williams/w/copy-of-scanpyqcplusdraft ) This generates a UMAP, does clustering and calculates cluster marker genes. 

 
For single sample experiments, there is a streamlined workflow that runs all 3 sub-workflows all at once 
* **Single sample workflow:** [link](https://usegalaxy.org.au/u/s.williams/w/copy-of-scrnaseqcountsmatrixtoqc) This workflow loads counts matrix, does some basic processing, suitable for a single sample.


These workflows are all available on galaxy australia.

# Analysis overview

![Processing flowchart](./images/workflow_diagram_simple.png)

1. Start with fastq files  
2. Cellranger|starSOLO will align the reads to the genome, and make a table of the number of times each gene is counted per cell. 
3. Perform some basic QC on the counts matrix, and filter out ‘cells’ that have too little RNA counts or too much mitochondrial gene content 
4. Run some basic single cell analyses: Normalisation, PCA, UMAP, clustering and identify cluster markesr 
5. The resulting AnnData object can be analysed further. See ‘Next steps’ 


# Example output
 
When run in full, these workflows produce the following main outputs

* A processed [AnnData](https://anndata.readthedocs.io/en/latest/) file, which contains gene expression and annotation information ready for downstream analysis.
* Tables of 'marker' gene information - to aid detemination of cell types present in the experiment
* Two report summaries. 
  + Single cell QC report [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-ede7b160ea86b66e) : This shows QC metrics at the cell level, to evaluate filtering thresholds and data quality. 
  + Single cell basic processing report [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-21aa7559fbcd167e) : Some basic UMAP, clustering and cluster marker results to begin an analysis with. 


# User guide

## Running a single sample workflow

When there is only a single biological sample in a study, there is a streamlined workflow. 

* [Single sample workflow](https://usegalaxy.org.au/u/s.williams/w/copy-of-scrnaseqcountsmatrixtoqc)

1. Input *will be* fastq files.
2. Search for the *scRNAseq Single sample workflow* under 'workflow' and run it. You'll be prompted to customise any filtering parameters, and choose a sensible name for the biological sample. 

![Single Sample Launch prompt](./images/screen_single_sample_launch.png)

3. This pipeline should take a [few minutes/few hours] to run. 

4. Return to galaxy, and look up the run, or invocation, of this workflow: 

_User Menu > Workflow Invocations_

This brings up the history of workflow invocations. This particular workflow runs a number of sub-workflows. The **Cell QC** and **QC to Basic Processing** subworkflows produce reports, which can be viewed shared and saved.

* Cell QC : Generates the cell level QC plots  
* QC to Basic Processing: Shows the subsetquent UMAP, clustering and marker information. Includes links to download the processed AnnData object for downstream work 

![workflow invocations screen](./images/screen_invocations.png)


## Running a multi sample experiment

With multi-sample experiments, each sample is loaded independently and then combined. The overall method is the same, but the QC and processing steps are run separately. 

![subworkflows](./images/workflow_diagram_large.jpg)

1. Upload the raw data for one sample. Then run The ‘scRNAseq: Load counts matrix’ workflow – this will prompt you for a sample name that will be used throughout. 

![Load counts matrix launch](./images/screen_load_counts_matrix_launch.png)

```
**IF USING the 3 file input, describe how to change barcodes.tsv to tabular, and note that features == genes, and that .tsv.gz will work fine. ** 
```

This part of the workflow will load the counts matrix into an anndata object, and then adds an extra column in the metadata called ‘sample’. This means the sample information can be tracked when multiple samples are combined. 

The AnnData object that it produces in your history will probably be named something like ‘Manipulate AnnData (add_annotation) on data 20 and data 17’. You may choose to rename this object via the ‘edit attributes’ option in the history panel, so its easier to find later. 

![rename](./images/screen_rename.png)


2. In the same history, repeat for all other samples. 

3. Next, join all samples with the ‘Manipulate AnnData object’ Tool (search on the tools pane on the left). 

This tool can do several different operations – listed under ‘Function to manipulate the object’, but we want the default; “Concatenate along the observation axis”. This combines cells (observations) from multiple sample runs. 
 
Choose the AnnData object of one of your samples in the  ‘Annotated data matrix’ dropdown. Then, choose the rest of your samples under ‘Annotated matrix to add’. Use ctrl-select / option-select to highlight multiple samples. 

_Note:_ Be careful not to select the sample in the ‘Annotated data matrix’ dropdown again – else it will be joined to itself! In this example there is only two samples, so only pbmc8k is selected to be added to pbmc1k. 

![join anndatas](./images/screen_multi_merge.png)

A combined anndata object will created in your history. 

4. Next, run the **scRNAseq Cell QC** workflow on your combined anndata object.  

This workflow plots some basic cell-level QC thresholds, and applies the QC thresholds to produce a filtered Anndata object.  

![](./images/screen_cellQC_launch.png)

5. Once it finishes running, view the report (Go to User menu > Invocations to find it).  

You’ll notice you can see each sample plotted separately in the QC plots. You may elect to rerun with tweaked thresholds (e.g. higher minimum counts threshold) once you've seen this output.

![](./images/cell_qc_plot.png)


6. If you are happy with the filtering thresholds, you can  launch the next workflow, **scRNAseq QC to Basic Processing** to do some routine single cell calculations.  It only asks for the filtered Anndata Object (typically the last AnnData in your history, with a name ending in ‘Filtered Cells Anndata’.) 

![launch basic processing](./images/screen_qc_to_basic_processing_launch.png)


7. This takes a few minutes to run. Once finished, return to the invocations page to see the QC to Basic Processing report, as per the single sample workflow. 

The first umap now shows the different samples that make up the data. 

![umap samples](./images/umap_by_sample.png)


## Next steps

The AnnData object generated is ready for analysis! Options include

* CellXgene : CellXgene is a tool for browsing and exploring single cell data. It use the AnnData object, and all of the annotation stored within it (expression, clusters, sample names). See 
[CellXgene documentation](https://cellxgene.cziscience.com/docs/01__CellxGene), and a [user guide](https://icbi-lab.github.io/cellxgene-user-guide/). This can be launched within galaxy.
* Scanpy : [Scanpy](https://scanpy.readthedocs.io/en/stable/) is a python based suite of tool for working with single cell data. 
  + Many scanpy functions are available within galaxy (and are used within this workflow). Explore the [galaxy scRNAseq tutorials](https://training.galaxyproject.org/training-material/topics/single-cell/) for more information
  + Or, you download your AnnData object and work with the scanpy toolkit directly, using python on your computer or other server. This can provide more fine-grained control of your analysis. There are [many tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html) to work from.

Note that there are toolkits other than scanpy (e.g. Seurat, SingleCellExperiment objects) which may not be directly compatible without conversions.



## STAR Solo vs Cell Ranger

If you are unable to use cellranger, an alternate version of the workflow use STAR solo (which uses a MIT liscence and can be configured to support different sequencing technologies)

* Your fastq files
* Whitelist files: The list of expected barcodes in the kit - which varys by technology and chemistry. For 10X chromium data see here; https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-
* Genome and matching annotation reference (see below).

StarSOLO will require a genome sequence file (fasta format), and a .gtf or .gff file of the gene positions. Take care to ensure these are from the same genome version. A good source of suitable refernece/annotation pairs for a wide range of species is the [ensembl download index](https://asia.ensembl.org/info/data/ftp/index.html). For example, human reference data at ensembl v109:

* GRCH38 primary assembly: [Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz)
* Gene Annotation : [Homo_sapiens.GRCh38.109.gtf.gz](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz)

You can supply links to that data directly to a galaxy history via `upload data > Paste/Fetch data`, so there's no need to download/upload large files from your computer. And since you'll likely want to reuse the same refernece in new analyses, [its possible to copy to new histories as needed](https://training.galaxyproject.org/training-material/faqs/galaxy/histories_copy_dataset.html)

Note that STARsolo will not produce a .cloupe object for the cell loupe browser.


# Background and Tutorials 

For more general information about **single cell RNAseq processing on galaxy**; there are some excellent tutorials to be found here on the [galaxy training website scRNA section](https://training.galaxyproject.org/training-material/topics/single-cell/). The workflow implemented here is heavily influenced by the [Clustering 3kPBMCs with Scanpy tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html )

More general information on **using galaxy** can be found on the [galaxy training website](https://training.galaxyproject.org/training-material/)
 
There are many general resources online about the princials of single cell analysis. The [Scanpy preprocessing and clustering tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) may be of particular use becuase it describes the scanpy methods used in this workflow. Even if you don't use the python code, it works through and explains many of the plots these workflows generate.  



# License(s)


Note that Cellranger is subject to a [custom licence](https://github.com/10XGenomics/cellranger/blob/master/LICENSE). Most notably, it can only be used with 10X technology. Any use of these workflows that use these cellRanger must adhere to that licence. 

An alternative option is the STARsolo workflow, distributed under an MIT licence. 

Otherwise, useage of these workflows is depenant on the (generally permissive) licences of the underlying tools and platforms; including;

* Cell Ranger : https://github.com/10XGenomics/cellranger/blob/master/LICENSE
* STARSolo : https://github.com/alexdobin/STAR/blob/master/LICENSE
* Galaxy : https://galaxyproject.org/admin/license/
* Galaxy australia terms of service: https://site.usegalaxy.org.au/about#terms-of-service
* Scanpy : https://github.com/scverse/scanpy/blob/master/LICENSE\
* Scanpy Scripts: https://github.com/ebi-gene-expression-group/scanpy-scripts/blob/develop/LICENSE




# Acknowledgements/citations/credits

The workflow implemented here is heavily influenced by the [Clustering 3kPBMCs with Scanpy tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html )
