---
title: scRNAseq Processing Workflow landing page
type: guide
contributors: 
description: How-to guide for scRNAseq workflows on Galaxy Australia.
affiliations: 
---


## About 

This document describes how to use some scanpy-based scRNAseq workflows on galaxy Australia. 

The aim of these workflows is to handle the routine ‘boring’ part of single cell RNAseq data processing. It will produces an ‘AnnData’ object, which can then be used as a base for downstream analysis – either within galaxy or outside of it. AnnData is a standard format used by the ‘scanpy’ python package. 

These workflows represent just one way of processing data for a ‘typical’ scRNAseq experiment – there are many other options!  

The how-to guide is available [here](scrnaseq_wf_guide.html)



## Acknowledgements

The workflows are based on the excellent [Clustering 3kPBMCs with Scanpy tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html)

This guide makes use of the ELIXIR toolkit theme: [![theme badge](https://img.shields.io/badge/ELIXIR%20toolkit%20theme-jekyll-blue?color=0d6efd)](https://github.com/ELIXIR-Belgium/elixir-toolkit-theme)

## References

These workflows depend on the following tools and resources;

* Scanpy : (https://scanpy.readthedocs.io/en/stable/) [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0)
* Cell Ranger (10X): (https://support.10xgenomics.com/docs/citations)
* STARSolo : [publication](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full)
* Galaxy : (https://galaxyproject.org/) 
* Galaxy Australia: (https://usegalaxy.org.au/)
* Scanpy Scripts: (https://github.com/ebi-gene-expression-group/scanpy-scripts)


