---
title: scRNAseq Processing Workflows
type: guide
description: How-to guide for scRNAseq workflows on Galaxy Australia
affiliations: QCIF
toc: false
---

This document describes how to use some scanpy-based scRNAseq workflows on galaxy Australia. 

The aim of these workflows is to handle the routine ‘boring’ part of single cell RNAseq data processing. It will produces an ‘AnnData’ object, which can then be used as a base for downstream analysis – either within galaxy or outside of it. AnnData is a standard format used by the ‘scanpy’ python package. 

These workflows represent just one way of processing data for a ‘typical’ scRNAseq experiment – there are many other options!  


This document describes 3 sub-workflows for processing single cell RNAseq data with scanpy 

* **Load counts matrix (CellRanger):**    This workflow adds a sample name, which enables multi-sample analyses. NB: Does not yet run cellranger. 
* **Count and load (starSOLO):** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=466&trs_version=2)   This workflow takes fastq files, uses starSOLO to generate a counts matrix, loads the data into a standard AnnData format and adds a sample name. 
* **Single cell QC:** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=467&trs_version=2)  This workflow generates some basic QC plots and applies filtering 
* **Single cell QC to basic processing:** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=468&trs_version=2)  This generates a UMAP, does clustering and calculates cluster marker genes. 

 
For single sample experiments, there are streamlined workflows that runs all 3 sub-workflows all at once 
* **Single sample workflow (CellRanger):** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=464&trs_version=2)  This workflow loads counts matrix, does some basic processing, suitable for a single sample.
* **Single sample workflow (StarSOLO):** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=465&trs_version=2)  This workflow loads counts matrix, does some basic processing, suitable for a single sample.


These workflows are all available on galaxy australia.

# Background

There are many general resources online about the principals of single cell analysis. The [Scanpy preprocessing and clustering tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) may be of particular use because it describes the scanpy methods used in this workflow. Even if you don't use the python code, it works through and explains many of the plots these workflows generate.  

General information on **using galaxy** can be found on the [galaxy training website](https://training.galaxyproject.org/training-material/)
 
For more general information about **single cell RNAseq processing on galaxy**; there are some excellent tutorials to be found here on the [galaxy training website scRNA section](https://training.galaxyproject.org/training-material/topics/single-cell/). The workflow implemented here is heavily influenced by the [Clustering 3kPBMCs with Scanpy tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html )


{% include callout.html type="note" content="These workflows will not currently work with the fastq data from 10X's probe based 'Gene expression Flex' kits, only sequence based kits." %}

# Analysis overview

![Processing flowchart](./images/workflow_diagram_simple.png)

1. Start with fastq files. (Or skip to step 3 with a counts matrix)
2. CellRanger or starSOLO will align the reads to the genome, and make a table of the number of times each gene is counted per cell. 
3. Perform some basic QC on the counts matrix, and filter out ‘cells’ that have too little RNA counts or too much mitochondrial gene content 
4. Run some basic single cell analyses: Normalisation, PCA, UMAP, clustering and identify cluster markers 
5. The resulting AnnData object can be analysed further. See ‘Next steps’ 

When run in full, these workflows produce the following main outputs

* A processed [AnnData](https://anndata.readthedocs.io/en/latest/) file, which contains gene expression and annotation information ready for downstream analysis.
* Tables of 'marker' gene information - to aid determination of cell types present in the experiment
* Two report summaries. 
  + Single cell QC report [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-ede7b160ea86b66e) : This shows QC metrics at the cell level, to evaluate filtering thresholds and data quality. 
  + Single cell basic processing report [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-21aa7559fbcd167e) : Some basic UMAP, clustering and cluster marker results to begin an analysis with. 



# User guide

This diagram shows how the different workflows work together - the choice depends on whether you have a single or multiple samples, or are using CellRanger or StarSOLO. 

{% include image.html file="/workflow_options.png" max-width="800px" %}



## Input options? STAR Solo vs Cell Ranger vs counts matrix

Within these workflows, there are two options for generating a counts matrix from fastq sequence data (counting copies per gene per cell); 
This changes how the first part of analysis is run, but subsequent steps are the same. 

* CellRanger: Developed by 10X genomics, this is a widely used package that works for 10X genomics data. Usage is subject to licence conditons (see Licenses section) and available only to registered users. Currenlty only the inbuilt human and mouse references are avaliable.  CellRanger will also output a '.cloupe' file, which may be used with the 'cell loupe' desktop program.
* StarSOLO: An open source alternative, which can be configured to support different sequencing technologies. The current implementation works multiple species (subject to reference annotations).
* Counts matrix: Sometimes you already have a counts matrix - these are available for many public datasets already. Or, if you are using data from probe-based kits (e.g. 10X fixed rna kit), you'll need to get your counts matrix as described by the manufactuer.




## Prepare your fastq inputs

Single cell sequencing data is typically generated as paired-end sequencing data.

In galaxy, we can store paired sequencing data in a 'paired collection' - that way the R1 and matched R2 will always be together. For more information what a galaxy _collection_ is see [guide to collections in galaxy](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html)  

These workflows expect a paired collection of data for one sample. Its possible to have multiple fastq pairs for one biological sample - but this workflow will usually be run on a paired collection of one pair.

{% include callout.html type="note" content="In galaxy, collections (paired or not) are usually used group multiple samples. In this workflow however, we will have each sample in its own collection" %}

You may also have I1 and I2 fastq datasets, but these are indicies not used in these pipelines. If necessary, data should be demultiplexed before starting - each R1/R2 file pair should only contain data from one sample. 

1. Load the _R1 and _R2 fastq files into your galaxy history. See [guide to loading data](https://galaxyproject.org/support/loading-data/) for details.

2. Select the fastq files that make up one biological sample. Typically this may be 2 files (R1 and R2), in this case, there are 4. 

    a. Click the 'tick in a box' select files button at the top of the history (dataset) panel.
    b. Checkboxes appear for all history items - check the specific files needed.
    c. There will be a dropdown button at the top of the history panel - 'selected' - hit this, and select 'Build list of dataset pairs' 

{% include image.html file="/pair_collection_fastq_select.png" alt="Selecting fastqs for collection" max-width="500px" %}

3. This brings up a new window to indicate the forward and reverse reads. It may initially say that it cannot create any pairs. In this case, we need to update the 'Forward' box to '_R1' and the reverse to '_R2'. These identifiers are then found in the filenames.

{% include image.html file="/pair_collection_pairing.png" alt="Pairing the collection" max-width="800px" %}

4. Hit 'auto-pair' and now they are highlighed green as below.

{% include image.html file="/pair_collection_paired.png" alt="Paired collection" max-width="800px" %}



5. Give it a sensible name, then hit 'create collection', and now that paired collection will be visible in your history as below;

{% include image.html file="/fastqcollection.png" alt="collection in history" max-width="500px" %}


## Running a single sample workflow

When there is only a single biological sample in a study, there is a streamlined workflow, one for CellRanger and one for StarSolo. 

### Using CellRanger



### With a Counts matrix

If you have a counts matrix. 


{% include image.html file="/screen_single_sample_launch.png" alt="Single Sample Launch prompt.png" max-width="800px" %}



### Using StarSOLO

The input for the star solo workflow is:

* Your fastq files, prepared as above
* Whitelist files: The list of expected barcodes in the kit - which varies by technology and chemistry. For 10X chromium data see here; https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-
* Genome and matching annotation reference (see below).

StarSOLO will require a genome sequence file (fasta format), and a .gtf or .gff file of the gene positions. Take care to ensure these are from the same genome version. A good source of suitable refernece/annotation pairs for a wide range of species is the [ensembl download index](https://asia.ensembl.org/info/data/ftp/index.html). For example, human reference data at ensembl v109:

* GRCH38 primary assembly: [Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz)
* Gene Annotation : [Homo_sapiens.GRCh38.109.gtf.gz](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz)

You can supply links to that data directly to a galaxy history via `upload data > Paste/Fetch data`, so there's no need to download/upload large files from your computer. And since you'll likely want to reuse the same reference in new analyses, [its possible to copy to new histories as needed](https://training.galaxyproject.org/training-material/faqs/galaxy/histories_copy_dataset.html). Note that the current configuration of StarSOLO will build a reference every time it is run, which can take some time. 

1. Load your input files into galaxy.

2. Import the **Single sample workflow (StarSOLO):** (listed above) This workflow will include a couple of sub workflows.

3. Hit run to bring up the following launch form. The first option _'paired fastqs for one sample'_ should be the paired collection you created and named above. You'll also be prompted to customise any filtering parameters, and choose a sensible name for the biological sample.

{% include image.html file="/screen_single_sample_starsolo_launch.png" alt="Single Sample Launch prompt" max-width="800px" %}




Once ready hit the 'run workflow' button. This pipeline should take an hour or several to run. 

4. Return to galaxy, and look up the run, or invocation, of this workflow: 

_User Menu > Workflow Invocations_

This brings up the history of workflow invocations. This particular workflow runs a number of sub-workflows. The **Cell QC** and **QC to Basic Processing** subworkflows produce reports, which can be viewed shared and saved.

* Cell QC : Generates the cell level QC plots. [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-ede7b160ea86b66e) 
* QC to Basic Processing:  Shows the subsequent UMAP, clustering and marker information. Includes links to download the processed AnnData object for downstream work . [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-21aa7559fbcd167e)

{% include image.html file="/screen_invocations.png" alt="workflow invocations screen" max-width="500px" %}




## Running a multi sample experiment

With multi-sample experiments, each sample is loaded independently and then combined. The overall method is the same, but the QC and processing steps are run separately. 


1. Upload the raw data for one sample. There are 3 options for the first step, but after that the process is the same.

-------

**Option 1: Cell Ranger** 

TBA.

**Option 2: StarSOLO** If you have fastq files and want to use starSOLO. 

Run the 'scrnaseq: Count and Load with starSOLO workflow.' The options are a subset of those listed for the single-sample workflow, with the same inputs. 

{% include image.html file="/screen_counts_and_load_starsolo_launch.png" alt="Load counts matrix launch" max-width="800px" %}


**Option 3: Counts matrix** If you already have a counts matrix then run The ‘scRNAseq: Load counts matrix’ workflow – this will prompt you for a sample name that will be used throughout. 

![Load counts matrix launch](./images/screen_load_counts_matrix_launch.png)

{% include callout.html type="note" content="Older datasets will have genes.tsv(.gz), where as newer dataset will have features.tsv(.gz). Either works" %}

{% include callout.html type="note" content="Sometimes, if you have an uncompressed barcodes.tsv file, galaxy will determine that it is a text file, rather than the desired '.tabular', which will prevent it from being selected as input (will be missing from the dropdown). This is fixed by manually telling galaxy to change the datatype to .tabular. [Instructions](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_change_datatype.html)" %}

This part of the workflow will load the counts matrix into an AnnData object, and then adds an extra column in the metadata called ‘sample’. This means the sample information can be tracked when multiple samples are combined. 

------
Any option outputs an AnnData object in you history - containing your counts and annotated with the sample name internally.

That object in your history will probably be named something like ‘Manipulate AnnData (add_annotation) on data 20 and data 17’. You may like to rename this object via the ‘edit attributes’ option in the history panel, so its easier to find later. 

{% include image.html file="/screen_rename.png" alt="Rename" max-width="500px" %}


2. In the same history, repeat for all other samples. 

3. Next, join all samples with the ‘Manipulate AnnData object’ Tool (search on the tools pane on the left). 

This tool can do several different operations – listed under ‘Function to manipulate the object’, but we want the default; “Concatenate along the observation axis”. This combines cells (observations) from multiple sample runs. 
 
Choose the AnnData object of one of your samples in the  ‘Annotated data matrix’ dropdown. Then, choose the rest of your samples under ‘Annotated matrix to add’. Use ctrl-select / option-select to highlight multiple samples. 


{% include callout.html type="important" content="Be careful not to select the sample in the ‘Annotated data matrix’ dropdown again – else it will be joined to itself! In this example there is only two samples, so only pbmc8k is selected to be added to pbmc1k." %}

{% include image.html file="/screen_multi_merge.png" alt="join anndatas" max-width="500px" %}

A combined AnnData object will created in your history. 

4. Next, run the **scRNAseq Cell QC** workflow on your combined AnnData object.  

This workflow plots some basic cell-level QC thresholds, and applies the QC thresholds to produce a filtered AnnData object. Configure thresholds as appropriate.

{% include image.html file="/screen_cellQC_launch.png" alt="cell QC_launch" max-width="800px" %}


5. Once it finishes running, view the report (Go to User menu > Invocations to find it).  

You’ll notice you can see each sample plotted separately in the QC plots. You may elect to rerun with tweaked thresholds (e.g. higher minimum counts threshold) once you've seen this output.

{% include image.html file="/cell_qc_plot.png" alt="Cell qc quality plot" max-width="500px" %}

6. If you are happy with the filtering thresholds, you can launch the next workflow, **scRNAseq QC to Basic Processing** to do some routine single cell calculations. It only asks for the filtered AnnData Object (typically the last AnnData in you history, which may be named something like 'scanpy scrublet on data x') 

{% include image.html file="/screen_qc_to_basic_processing_launch.png" alt="launch basic processing" max-width="800px" %}

7. This takes a few minutes to run. Once finished, return to the invocations page to see the QC to Basic Processing report, as per the single sample workflow. 

The first umap now shows the different samples that make up the data. 

{% include image.html file="/umap_by_sample.png" alt="UMAP coloured by sample" max-width="500px" %}


## Next steps

The AnnData object generated is ready for analysis! Options include

* CellXene : CellXgene is a tool for browsing and exploring single cell data. It use the AnnData object, and all of the annotation stored within it (expression, clusters, sample names). See 
[CellXgene documentation](https://cellxgene.cziscience.com/docs/01__CellxGene), and a [user guide](https://icbi-lab.github.io/cellxgene-user-guide/). This can be launched within galaxy.
* Scanpy : [Scanpy](https://scanpy.readthedocs.io/en/stable/) is a python based suite of tool for working with single cell data. 
  + Many scanpy functions are available within galaxy (and are used within this workflow). Explore the [galaxy scRNAseq tutorials](https://training.galaxyproject.org/training-material/topics/single-cell/) for more information
  + Or, you download your AnnData object and work with the scanpy toolkit directly, using python on your computer or other server. This can provide more fine-grained control of your analysis. There are [many tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html) to work from.

Note that there are toolkits other than scanpy (e.g. Seurat, SingleCellExperiment objects) which are not directly compatible without conversions.




# Licenses


Note that CellRanger is subject to a [custom license](https://github.com/10XGenomics/cellranger/blob/master/LICENSE). Most notably, it can only be used with 10X technology. Any use of these workflows that use these cellRanger must adhere to that license. 

An alternative option is the STARsolo, distributed under an MIT license. 

Otherwise, usage of these workflows is dependant on the (generally permissive) licenses of the underlying tools and platforms; including;

* [Cell Ranger](https://github.com/10XGenomics/cellranger/blob/master/LICENSE)
* [STARSolo](https://github.com/alexdobin/STAR/blob/master/LICENSE)
* [Galaxy](https://galaxyproject.org/admin/license/)
* [Galaxy Australia terms of service](https://site.usegalaxy.org.au/about#terms-of-service)
* [Scanpy](https://github.com/scverse/scanpy/blob/master/LICENSE\)
* [Scanpy Scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/blob/develop/LICENSE)



