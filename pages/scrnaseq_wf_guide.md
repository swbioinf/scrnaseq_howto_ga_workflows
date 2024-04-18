---
title: scRNAseq Processing Workflows
type: guide
description: How-to guide for scRNAseq workflows on Galaxy Australia
toc: true
sidebar: false
---


This document describes how to use some scanpy-based scRNAseq workflows on galaxy Australia. 

The aim of these workflows is to handle the routine ‘boring’ part of single cell RNAseq data processing. It will produces an ‘AnnData’ object, which can then be used as a base for downstream analysis – either within galaxy or outside of it. AnnData is a standard format used by the ‘scanpy’ python package. 

These workflows represent just one way of processing data for a ‘typical’ scRNAseq experiment – there are many other options!  




# Workflows

This document describes a set of workflows for processing single cell RNAseq data with scanpy 

These workflows are all available on galaxy australia.

{% include callout.html type="note" content="To import these workflows, you must be logged in to galaxy australia with an account, before clicking the import links." %}



## Standard workflows

With 3 ways to load the data:

* **Load counts matrix**: [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=512&trs_version=1)   This workflow takes a counts matrix (of form counts.mtx + genes.csv + features.csv), loads the data into a standard AnnData format and adds a sample name.
* **Count and load (starSOLO)**: [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=513&trs_version=1)  This workflow takes fastq files, uses starSOLO to generate a counts matrix, loads the data into a standard AnnData format and adds a sample name. 
* **Count and Load (Cell Ranger)**: [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=646&trs_version=1)  This workflow takes fastq files, uses Cell Ranger to generate a counts matrix, loads the data into a standard AnnData format and adds a sample name.

<!-- CHecking galaxy to repo Nov 23: ok, ok, current. -->


Followed by two steps of routine processing: 

* **Single cell QC:** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=467&trs_version=3)  This workflow generates some basic QC plots and applies filtering 
* **Single cell QC to basic processing:** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=468&trs_version=2)  This generates a UMAP, does clustering and calculates cluster marker genes. 

<!-- CHecking galaxy to repo Nov 23: ok ok-->

## Single sample streamlined workflows
 
For single sample experiments, there are streamlined workflows that the steps all at once.

* **Single sample workflow (Cell Ranger):** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=647&trs_version=1)  This workflow starts from fastq files, generates counts with Cell Ranger, QCs and does some basic processing. It suitable for a single sample.
* **Single sample workflow (StarSOLO):** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=465&trs_version=3)  This workflow starts from fastq files, generates counts with starSOLO, QCs and does some basic processing. It suitable for a single sample.
* **Single sample workflow (Counts matrix):** [import](https://usegalaxy.org.au/workflows/trs_import?trs_server=workflowhub.eu&trs_id=514&trs_version=2) This workflow starts from a counts matrix, QCs and does some basic processing. It suitable for a single sample.

<!-- CHecking galaxy to repo Nov 23: new. DIFFERENT,  DIFFERENT-->


# Background

There are many general resources online about the principals of single cell analysis. The [Scanpy preprocessing and clustering tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) may be of particular use because it describes the *scanpy methods* used in this workflow. Even if you don't use the python code, it works through and explains many of the plots these workflows generate.  

General information on *using galaxy* can be found on the [galaxy training website](https://training.galaxyproject.org/training-material/)
 
For more general information about *single cell RNAseq processing on galaxy*; there are some excellent tutorials to be found here on the [galaxy training website scRNA section](https://training.galaxyproject.org/training-material/topics/single-cell/). The workflow implemented here is heavily influenced by the [Clustering 3kPBMCs with Scanpy tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html )


{% include callout.html type="note" content="These workflows will not currently work with the fastq data from 10X's probe-based 'Gene expression Flex' kits, only sequence based kits." %}

# Analysis overview

![Processing flowchart](./images/workflow_diagram_simple.png)

1. Start with fastq files.
2. Cell Ranger or starSOLO will align the reads to the genome, and make a table of the number of times each gene is counted per cell. 
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

This diagram shows how the different workflows work together - initial part of the process depends on whether you have a single or multiple samples, or are using Cell Ranger, StarSOLO or a counts matrix. The subsequent process is the same.

What are the 3 different options?

* **Cell Ranger**: Developed by 10X genomics, this is a widely used package that works for 10X genomics data. Usage is subject to licence conditions (see Licenses section) and available only to registered users. Currently only the inbuilt human and mouse references are available (it is possible to build and use a custom reference on galaxy, but not within this workflow). 
* **StarSOLO**: An open source alternative, which can be configured to support different sequencing technologies or custom references. 
* **Counts matrix**: Sometimes you already have a counts matrix - these are available for many public datasets already. Or, if you are using data from probe-based kits (e.g. 10X fixed rna kit), you'll need to get your counts matrix as described by the manufacturer.


{% include image.html file="/workflow_options.png" max-width="800px" %}


## Configuration options

There are some QC thresholds that may need to be adjusted. Note that there is no ideal threshold for minimum counts and genes per cell - thresholds vary wildly between technologies and experiment. 

* **Mitochondrial Prefix**: Prefix on mitochondrial gene names. Default MT- typical for human, mt- for mouse. Using the wrong form will results in plotting errors - check your gene names if this occurs.
If mitochondrial genes are not prefixed, mitochondrial flagging, plotting and filtering steps will need to be removed - else they'll generated errors.
* **Min Count Per Cell**: Cells with fewer than this number of total counts will be removed. Check your QC plots to make sure your thresholds are reasonable!
* **Min Genes Per Cell**: Cells with fewer than this many distinct genes observed will be removed. Again, check your QC plots to make sure your thresholds are resonable!
* **MaxMTpc** : Cells with more than this percentage of counts coming from mitochondrial genes will be removed. Generally high MT-content cells are uninteresting, but again, the threshold can differ between experiments.

## Prepare your fastq inputs

Single cell sequencing data is typically generated as paired-end sequencing data. This is handled by 'collections' in galaxy. For more information what a galaxy _collection_ is see [guide to collections in galaxy](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html)  

Cell Ranger and starSOLO expect fastq files to be arranged in slightly differnet collections, as described below.


### For Cell Ranger


Cell Ranger expects fastq files in a simple collection.

1. Load the _R1 and _R2 fastq files into your galaxy history. See [guide to loading data](https://galaxyproject.org/support/loading-data/) for details.

{% include callout.html type="important" content="The Cell Ranger tool expects fastq files to be named in a bcl2fastq naming convention *[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz* .  This is typical for illumina sequencers. If your data is of a different naming convention (e.g. SRA downloads), you may need to rename your files. More detail [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input)." %}

2. Select the fastq files that make up one biological sample. Typically this may be 2 files (R1 and R2), but there can be more (e.g. often multiple lanes).

    a. Click the 'tick in a box' select files button at the top of the history (dataset) panel.
    b. Checkboxes appear for all history items - check the specific files needed.
    c. There will be a dropdown button at the top of the history panel - 'selected' - hit this, and select 'Build dataset list' 

{% include image.html file="/screen_flat_fastq_collection.png" alt="Screenshot of selecting fastqs for cell ranger collection" max-width="500px" %}

However it expects the collection to be named in a particular way. Fastq files will typically have names like: *SomeSample_S2_L001_R1_001.fastq.gz*, and the collection should be named with everything before *_S2_L001_R1_001.fastq.gz*. Here that is *'500_PBMC_3p_LT_Chromium_X'*.

{% include callout.html type="important" content="The name of your collection matters!" %}

{% include callout.html type="note" content="If your run fails with the error: `ERROR: Requested sample(s) not found in fastq directory. Available samples: ...`; this indicates a mismatch between the collection name and filenames. Consider renaming your collection to what is listed as an available sample." %}




### For starSOLO 


In galaxy, we can store paired sequencing data in a 'paired collection' - that way the R1 and matched R2 will always be together. 

These workflows expect a paired collection of data for one sample. Its possible to have multiple fastq pairs for one biological sample - but this workflow will usually be run on a paired collection of one pair.

{% include callout.html type="note" content="In galaxy, collections (paired or not) are usually used group multiple samples. In this workflow however, we will have each sample in its own collection" %}

You may also have I1 and I2 fastq datasets, but these are indices not used in these pipelines. If necessary, data should be demultiplexed before starting - each R1/R2 file pair should only contain data from one sample. 

1. Load the _R1 and _R2 fastq files into your galaxy history. See [guide to loading data](https://galaxyproject.org/support/loading-data/) for details.

2. Select the fastq files that make up one biological sample. Typically this may be 2 files (R1 and R2), in this case, there are 4. 

    a. Click the 'tick in a box' select files button at the top of the history (dataset) panel.
    b. Checkboxes appear for all history items - check the specific files needed.
    c. There will be a dropdown button at the top of the history panel - 'selected' - hit this, and select 'Build list of dataset pairs' 

{% include image.html file="/pair_collection_fastq_select.png" alt="Selecting fastqs for collection" max-width="500px" %}

3. This brings up a new window to indicate the forward and reverse reads. It may initially say that it cannot create any pairs. In this case, we need to update the 'Forward' box to '_R1' and the reverse to '_R2'. These identifiers are then found in the filenames.

{% include image.html file="/pair_collection_pairing.png" alt="Pairing the collection" max-width="800px" %}

4. Hit 'auto-pair' and now they are highlighted green as below.

{% include image.html file="/pair_collection_paired.png" alt="Paired collection" max-width="800px" %}


5. Give it a sensible name, then hit 'create collection', and now that paired collection will be visible in your history as below;

{% include image.html file="/fastqcollection.png" alt="collection in history" max-width="500px" %}


## Running a single sample workflow

When there is only a single biological sample in a study, there is a streamlined workflow.

{% include callout.html type="important" content="Please note: with a recent update we no longer see the two sub-reports when running the single sample workflows. This documentation is in the process of being updated." %}



### Using Cell Ranger


The input for the star solo workflow is:

* **Fastq files**: prepared as above
* **Reference**: Only if not working with human or mouse. 


{% include callout.html type="note" content=" If you're using a different species, you can make your reference with the mkgtf and mkref Cell Ranger subtools. Refer to [10X custom reference documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr). For use of these workflows consider editing the workflow, or running the Cell Ranger tool separately to use a counts matrix input." %}


1. If you haven't already, apply for access to Cell Ranger - see [Getting access to cell ranger](#getting-access-to-cell-ranger). Otherwise Cell Ranger won't be visible in your tools list!

2. Import the **Single sample workflow (Cell Ranger):** (listed above). This workflow will include a couple of sub workflows.

3. Open the workflow menu, find the workflow and hit run to bring up the following launch form. The first option _'Fastqs for one sample'_ should be the paired collection you created and named above. You'll also be prompted to customise any filtering parameters, and choose a sensible name for the biological sample.

{% include image.html file="/screen_cellranger_ss_launch.png" alt="Single Sample Launch prompt" max-width="800px" %}

Once ready hit the 'run workflow' button. This pipeline should take an hour or several to run. 


### Using StarSOLO

The input for the star solo workflow is:

* **Fastq files**: prepared as above
* **Whitelist files**: The list of expected barcodes in the kit - which varies by technology and chemistry. For 10X chromium data see here; https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-  . The file you need for each kit is listed, and the files themselves can be obtained by downloading the Cell Ranger software package.
* **Reference**: Genome and matching annotation reference (see below).

StarSOLO will require a genome sequence file (fasta format), and a .gtf or .gff file of the gene positions. Take care to ensure these are from the same genome version. A good source of suitable refernece/annotation pairs for a wide range of species is the [ensembl download index](https://asia.ensembl.org/info/data/ftp/index.html). For example, human reference data at ensembl v109:

* GRCH38 primary assembly: [Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz)
* Gene Annotation : [Homo_sapiens.GRCh38.109.gtf.gz](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz)

You can supply links to that data directly to a galaxy history via `upload data > Paste/Fetch data`, so there's no need to download/upload large files from your computer. And since you'll likely want to reuse the same reference in new analyses, [its possible to copy to new histories as needed](https://training.galaxyproject.org/training-material/faqs/galaxy/histories_copy_dataset.html). Note that the current configuration of StarSOLO will build a reference every time it is run, which can take some time. 

1. Load your input files into galaxy.

2. Import the **Single sample workflow (StarSOLO):** (listed above) This workflow will include a couple of sub workflows.

3. Open the workflow menu, find the workflow and hit run to bring up the following launch form. The first option _'paired fastqs for one sample'_ should be the paired collection you created and named above. You'll also be prompted to customise any filtering parameters, and choose a sensible name for the biological sample.

{% include image.html file="/screen_single_sample_starsolo_launch.png" alt="Single Sample Launch prompt" max-width="800px" %}

Once ready hit the 'run workflow' button. This pipeline should take an hour or several to run. 


### With a Counts matrix

If you have a counts matrix you can instead run the **scRNAseq Single Sample Processing Counts Matrix** workflow.

**Step 1** Load your input files into galaxy. The counts matrix is expected in a [matrix market](https://networkrepository.com/mtx-matrix-market-format.html) format, which consists of 3 files. It is a common format for single cell data, generated by tools like Cell Ranger or StarSOLO, and is often available for public datasets. Files may be gzipped (suffixed .gz), there is no need to uncompress them. 

 * **matrix.mtx (.gz)** : A file with the number of counts per cell
 * **genes.tsv (.gz)** or **features.tsv (.gz)** : A file listing all the genes or features in the experiment
 * **barcodes.tsv (.gz)** : A file listing the cellular barcodes in the experiment. 

You may have a 'filtered' and a 'raw' set of counts matrices - generally, the 'filtered' set should be used. The raw set may include large numbers of barcodes with very low counts - real cells are expected to be in the 'filtered' set.

**Step 2.** Import the **scRNAseq Single Sample Processing (Counts Matrix)** workflow (listed above)

**Step 3.** Hit run to bring up the following launch form.  There is no need for the reference genome, since genes have already been counted.

Once ready hit the 'run workflow' button. This pipeline should usually run in less than an hour. 

{% include image.html file="/screen_single_sample_launch.png" alt="Single Sample Launch prompt.png" max-width="800px" %}


## Check the results

Once your workflow has had time to run, return to galaxy and look up its 'invocation'.

_User Menu > Workflow Invocations_

This brings up the history of workflow invocations. This particular workflow runs a number of sub-workflows. The **Cell QC** and **QC to Basic Processing** subworkflows produce reports, which can be viewed shared and saved.

* Cell QC : Generates the cell level QC plots. [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-ede7b160ea86b66e) 
* QC to Basic Processing:  Shows the subsequent UMAP, clustering and marker information. Includes links to download the processed AnnData object for downstream work . [(example)](https://usegalaxy.org.au/u/s.williams/p/invocation-report-21aa7559fbcd167e)

{% include image.html file="/screen_invocations.png" alt="workflow invocations screen" max-width="500px" %}



## Running a multi sample experiment

With multi-sample experiments, each sample is loaded independently and then combined before running the downstream processing. 
This means that the overall process is the same, but the steps are launched manually. This is also useful when something goes wrong.


**Step 1.** Upload the raw data for one sample. 

**Step 2.** There are 3 options for this first step, but after that the process is the same.

  *Option 1: Cell Ranger* 

  Run the 'scrnaseq: Count and Load with Cell Ranger workflow.' The options are a subset of those listed for the single-sample workflow, with the same inputs. 

  {% include image.html file="/screen_count_and_load_cellranger.png" alt="Cellranger load and count launch prompt.png" max-width="800px" %}

  *Option 2: StarSOLO* If you have fastq files and want to use starSOLO. 

  Run the 'scrnaseq: Count and Load with starSOLO workflow.' The options are a subset of those listed for the single-sample workflow, with the same inputs. 

  {% include image.html file="/screen_counts_and_load_starsolo_launch.png" alt="Load and count with starSOLO launch" max-width="800px" %}

  *Option 3: Counts matrix* If you already have a counts matrix then run The ‘scRNAseq: Load counts matrix’ workflow – this will prompt you for a sample name that will be used throughout. 

  <!--![Load counts matrix launch](./images/screen_load_counts_matrix_launch.png)-->
  {% include image.html file="/screen_load_counts_matrix_launch.png" alt="Load counts matrix launch" max-width="800px" %}



  {% include callout.html type="note" content="Older datasets will have genes.tsv(.gz), where as newer dataset will have features.tsv(.gz). Either works" %}

  {% include callout.html type="note" content="Sometimes, if you have an uncompressed barcodes.tsv file, galaxy will determine that it is a text file, rather than the desired '.tabular', which will prevent it from being selected as input (will be missing from the dropdown). This is fixed by manually telling galaxy to change the datatype to .tabular. [Instructions](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_change_datatype.html)" %}

  This part of the workflow will load the counts matrix into an AnnData object, and then adds an extra column in the metadata called ‘sample’. This means the sample information can be tracked when multiple samples are combined. 


**Step 3.** Any of these 3 options outputs an AnnData object in you history - containing your counts and annotated with the sample name internally.

  That object in your history will probably be named something like ‘Manipulate AnnData (add_annotation) on data 20 and data 17’. You may like to rename this object via the ‘edit attributes’ option in the history panel, so its easier to find later. 

  {% include image.html file="/screen_rename.png" alt="Rename" max-width="500px" %}


**Step 4.** In the same history, repeat for all other samples. 

**Step 5.**  Next, join all samples with the ‘Manipulate AnnData object’ Tool (search on the tools pane on the left). 

  {% include callout.html type="important" content="There's currently a versioning conflict in the tools installed on galaxy. When using this concatenate tool, we currently need to use the older version '0.7.5+galaxy1' for compatability with downstream tools [How to select tool version](https://training.galaxyproject.org/training-material/faqs/galaxy/tools_change_version.html) [Issue details](https://github.com/bioconda/bioconda-recipes/issues/45164). You may find that AnnData objects saved with current versions of the AnnData library will fail." %}

  This tool can do several different operations – listed under ‘Function to manipulate the object’, but we want the default; “Concatenate along the observation axis”. This combines cells (observations) from multiple sample runs. 
    
  Choose the AnnData object of one of your samples in the  ‘Annotated data matrix’ dropdown. Then, choose the rest of your samples under ‘Annotated matrix to add’. Use ctrl-select / option-select to highlight multiple samples. 


  {% include callout.html type="important" content="Be careful not to select the sample in the ‘Annotated data matrix’ dropdown again – else it will be joined to itself! In this example there is only two samples, so only pbmc8k is selected to be added to pbmc1k." %}

  {% include image.html file="/screen_multi_merge.png" alt="join anndatas" max-width="500px" %}

  A combined AnnData object will created in your history. 

**Step 6.** Next, run the **scRNAseq Cell QC** workflow on your combined AnnData object.  

  This workflow plots some basic cell-level QC thresholds, and applies the QC thresholds to produce a filtered AnnData object. Configure thresholds as appropriate.

  {% include image.html file="/screen_cellQC_launch.png" alt="cell QC_launch" max-width="800px" %}

**Step 7.** Once it finishes running, view the report (Go to User menu > Invocations to find it).  

  You’ll notice you can see each sample plotted separately in the QC plots. You may elect to rerun with tweaked thresholds (e.g. higher minimum counts threshold) once you've seen this output.

  {% include image.html file="/cell_qc_plot.png" alt="Cell qc quality plot" max-width="500px" %}

**Step 8.** If you are happy with the filtering thresholds, you can launch the next workflow, **scRNAseq QC to Basic Processing** to do some routine single cell calculations. It only asks for the filtered AnnData Object (typically the last AnnData in you history, which may be named something like 'scanpy scrublet on data x') 

  {% include image.html file="/screen_qc_to_basic_processing_launch.png" alt="launch basic processing" max-width="800px" %}

**Step 9.** This takes a few minutes to run. Once finished, return to the invocations page to see the QC to Basic Processing report, as per the single sample workflow. 

  The first umap now shows the different samples that make up the data. 

  {% include image.html file="/umap_by_sample.png" alt="UMAP coloured by sample" max-width="500px" %}


# Next steps

The AnnData object generated by these workflows is intended for further analysis. Options include

* CellXene : CellXgene is a tool for browsing and exploring single cell data. It use the AnnData object, and all of the annotation stored within it (expression, clusters, sample names). See 
[CellXgene documentation](https://cellxgene.cziscience.com/docs/01__CellxGene), and a [user guide](https://icbi-lab.github.io/cellxgene-user-guide/). This can be launched within galaxy.
* Scanpy : [Scanpy](https://scanpy.readthedocs.io/en/stable/) is a python based suite of tool for working with single cell data. 
  + Many scanpy functions are available within galaxy (and are used within this workflow). Explore the [galaxy scRNAseq tutorials](https://training.galaxyproject.org/training-material/topics/single-cell/) for more information
  + Or, you download your AnnData object and work with the scanpy toolkit directly, using python on your computer or other server. This can provide more fine-grained control of your analysis. There are [many tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html) to work from.

Note that there are toolkits other than scanpy (e.g. Seurat, SingleCellExperiment objects) which are not directly compatible without conversions.


Note that the Cell Ranger tool will also output a '.cloupe' file, which may be used with the 'cell loupe' desktop program. If you download the .cloupe file, you'll need to change its file ending from .binary to .cloupe.  It is hidden by default within these workflows. 




# Licenses

Cell Ranger is proprietry software form 10X, and is only useable on galaxy australia once you have applied for access - see [Getting access to cell ranger](#getting-access-to-cell-ranger)

An alternative option is the STARsolo, distributed under an MIT license. 

Otherwise, usage of these workflows is dependant on the (generally permissive) licenses of the underlying tools and platforms; including;

* [Cell Ranger](https://github.com/10XGenomics/cellranger/blob/master/LICENSE)
* [STARSolo](https://github.com/alexdobin/STAR/blob/master/LICENSE)
* [Galaxy](https://galaxyproject.org/admin/license/)
* [Galaxy Australia terms of service](https://site.usegalaxy.org.au/about#terms-of-service)
* [Scanpy](https://github.com/scverse/scanpy/blob/master/LICENSE\)
* [Scanpy Scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/blob/develop/LICENSE)


## Getting access to Cell Ranger

When you first log in to Galaxy, the Cell Ranger tool with _not_ be visible. 

Cell Ranger is subject to a [custom license](https://github.com/10XGenomics/cellranger/blob/master/LICENSE). Most notably, it can only be used with 10X technology. 

To get access on your account, apply here (simply requires agreeing to the licence): https://site.usegalaxy.org.au/request/access/cellranger

Once successful, Cell Ranger will appear in your tool list under _MISCELLANEOUS TOOLS > Local Tools_ Once you receive a confirmation email, you may need to log out and log back in for it to take effect.


# Finding help

For general galaxy training, see the [Galaxy Training Network](https://training.galaxyproject.org/)

For error problems with galaxy tools on galaxy australia - please check out the [error debugging guide](https://galaxyproject.org/support/tool-error/) and use the built-in bug reporting functionality. Or, for more general or involved galaxy help try the [galaxy help forum](https://help.galaxyproject.org/).

To get in contact with problems you find specifically with this workflow; please report a ticket on [this github repository](https://github.com/swbioinf/scrnaseq_howto_ga_workflows/issues)
