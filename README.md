# CircCode2

CircCode2 is a Python3-base pipeline for the prediction and visualization of circRNA and ORF.

# Requirement
## Data:

- Transcript sequences (fasta format)
- Candidate circRNA sequences (fasta format)
- rRNA sequences (fasta format)
- Ribosome profiling data (sra format) 
- circRNA annotation file (txt format): This gtf file should be like this type(make by yourself), samples and gene_id are optional, you can input NaN if you don't need.

```

circRNA_id      strand      circRNA_lenth      gene_id      samples

TC-hsa-RERE_0160      -      185      ENSG00000142599.17      NaN

TC-hsa-UBE2J2_0010      -      607      ENSG00000160087.20      NaN

```

- riboseq_adapters fasta file

## Software:

- cutadapt (v 1.18)
- bowtie (v.1.2.3)
- STAR (v.2.7.8)
- bedtools (v.2.27.1)
- ViennaRNA (v 2.4.18)

## python3 package:

- Biopython
- Pandas
- numpy
- yaml
- argparse
- time
- os
- math
- PIL
- itertools
- torch
- sklearn
- tqdm
- matplotlib

## R package:

- UpSetR
- edgeR
- ggplot2

# Quick Start
Usually you can download the package from github simply,and then:
```
git clone git://github.com/Lilab-SNNU/CircCode2.git
cd CircCode2
tar -zxvf requiredSoft.tar.gz
tar -zxvf test_data.tar.gz
```


# Usage

Attention: Before you begin to use this package, you need to make sure that you have install the required software and add them to the environment variables. Besides, we provided the test data, you can verify your environment by test it.


1. Fill the config file (config.yaml), input absolute path of each required file.


2. Run CircCode2 step by step:


  - Deal raw ribo-seq data

  ```
   python deal_raw_reads.py -y config.yaml
  ```
  
  - Filter and obtain the mapped reads and circRNA

  ```
   python filter_coding_reads_circ.py -y config.yaml
  ```
  
  - Identify the ORF from mapped circRNAs by some features

  ```
  python filter_coding_orf.py -y config.yaml
  ```
  
  - Visualize the IRES region secondary structure, the circ annoated graph, the word annoated graph, the express_analysis of the circ and the distribution of the predicted ORFs
  
  ```
  python visual_circ_orf.py -y config.yaml
  ```


### How to fill in the config.yaml file?

When opening the config file in text format, there are some lines that need to be filled in, they are:

 - transcript_fasta: Fill in the absolute path of the transcript sequences related to circRNA(not the relative path!).

 - raw_reads: Fill in the absolute path of the Ribo-Seq data(One or two tissue samples data in the same species) related to your interest species(not the relative path!).
   
 - ribosome_fasta: Fill in the absolute path of the rRNA data related to your interest species(not the relative path!).
   
 - circrnas: Fill in the absolute path of the candidate circRNA(not the relative path!).You need to ensure that the id format of your circRNA is the same as the circrna_gtf file.
   
 - circrna_gtf: Fill in the absolute path of the candidate circRNA annoated file(not the relative path!). This gtf file should be like this type(make by yourself), samples are optional, you can input NaN if you don't need.

```

circRNA_id      strand      circRNA_lenth      gene_id      samples

TC-hsa-RERE_0160      -      185      ENSG00000142599.17      NaN

TC-hsa-UBE2J2_0010      -      607      ENSG00000160087.20      NaN

```
 
 - result_file_location: Fill in the absolute path of a folder to save the final results(not the relative path!).
 
 - tmp_file_location: Fill in the absolute path of a folder to save the temporary files(not the relative path!).
 
 - thread: Fill the number of threads running.
 
 - ribotype: The type of Ribo-seq data(sra or fastq).
 
  
## Citation

CircCode2: A useful workflow for Identification and Visualization of circRNAs with coding potential

## Contact us

If you encounter any problems while using CircCode2, please send an email to mr1997@snnu.edu.cn or glli@snnu.edu.cn, and we will resolve it as soon as possible.
