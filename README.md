### Introduction
A work-in-progress pipeline for QC of genotyping-by-sequencing (GBS) generated by National Genomics Infrastructure in Stockholm. GBS is a method closely related to RAD-seq (Restriction-site Associated DNA sequencing) that allows for deep yet sparsely sampled sequencing of many individuals in a highly multiplexed manner. Typical applications of this methods includes QTL mapping, GWAS studies, high resolution population differentiation/phylogeny, pedigree reconstruction and SNP discovery for other more high throughput assays. The Stockholm protocol is yet to be published, but this pipeline runs several steps that might be applicable to other RAD-seq flavours. These steps are mainly designed to only characterize the data and attempt to correct defects (*e.g.* adapter contamination and restriction-site sequencing errors) for further downstream analysis, but not to draw any biologically relevant conclusions:

* FastQC
* Adapter trimming using trimmomatic
* Flash to calculate % read overlap
* Jellyfish (testing for now)
* Stacks de novo pipeline with the most lenient parameters

### Implementation
The pipeline is implemented in [Nextflow](https://www.nextflow.io/). Implementation and usage details to follow.


### To-do

 Some ideas for continued development

* Module integrations for [MultiQC](http://multiqc.info/)
* Second pipeline for working with organisms where a reference genome is available
* Some form of support for DDRad
* More tweakable input files -- ie. sample-names and populations structures
* More Stacks output options (*e.g.* VCF and plink files) and general tweakability.
