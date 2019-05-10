# HIPS: Hiding In Plain Site
## Unannotated structural variants in public genomic data sets

## Contributors:
### Rebecca Torene, Alexandrea Stylianou, Jitong Cai, Ariel Gershman, Alexandra Weber, and Nancy Hansen

With thanks to former UMUC Bioinformatics master program students:
* Mary Schramko
* Upasana Pandey
* Eric Keller
* Mike Fox
* Al Koroma
* David Jalali
---
*Premise:*
---
Public data, such as gnomAD, makes VCFs of short nucleotide variants (SNVs) available. There is anecdotal evidence that some apparent “insertions” are actually part of larger structural variants

*Goal:*
---
Use VCFs of SNVs to call structural variants. In doing so, can we identify SVs that have otherwise escaped detection.

*Description:*
---
Large, aggregate, genomics data sets of short nucleotide variants (SNVs) such as [gnomAD/ExAC](https://gnomad.broadinstitute.org/) are an invaluable resource to the research and clinical genetics communities. Such data sets are often used as a control to determine the frequencies of variants in unaffected populations and gnomAD has [recently](https://www.biorxiv.org/content/10.1101/578674v1) released a structural variant call set on a subset of their samples. There is evidence, however, that additional SVs exist in the data, but haven't yet been annotated. For example, [this insertion](https://gnomad.broadinstitute.org/variant/21-18612332-A-ACCCAGGCAAACAGCGTCTGGAGTGGACCTCCAGGAAACAGGGTCTGGAGTGGACCTCCAGCAGACCTGCAGCAGAGGCACCTGTT) is actually indicative of a [known deletion](http://dgv.tcag.ca/dgv/app/variant?id=esv3646472&ref=hg19). Such anecdotal findings of deletions, processed pseudogenes, inversions, tandem duplications, and mobile element insertions exist in publicly available data, hiding there in plain sight.

Without individual alignment files, the broader genomics community cannot evaluate structural variants in the complete data set. We, therefore, developed **HIPS**, a structural variant caller that uses VCF of SNVs as input and then outputs structural variant calls in bed format. HIPS starts from VCF, parses insertion sequences and converts to fastA, aligns by BLAST, and filters BLAST results for possible structural variants. Running the snakefile then merges similar calls and adds annotation. The [Shiny app](https://hidinginplainsight.shinyapps.io/HidingInPlaneSight/) takes gnomAD data processed by HIPS and provides a web interface. 

*Suggested Future Improvements:*
---
This project may be on GitHub and it may have a great name (HIPS), but we cannot deny the fact that this is the product of a 3-day hackathon and that it is imperfect. We welcome the GitHub community to make improvements and suggest the following as future improvements:

* fastaToBed.py can be modified to split VCFs by position (already splits by chromosome) so that parallelization is more efficient
* fastaToBed.py arguments could be refined (e.g., meiOnly could be boolean flag instead of string input)
* stitch together fastaToBed.py and snakefile to have a seamless pipeline from VCF to annotated SV calls
* Add more flexibility for file names, species, etc.
* Option to output SV calls as VCF (currently outputs BED)
* Add a low quality flag for SVs based on bitscore and based on overlap of breakpoints with simple repeats

---
![flowchart](https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/blob/master/resources/prelim_flowchart.png)

# Web Interface
hosted at https://hidinginplainsight.shinyapps.io/HidingInPlaneSight/

## Graph View
![graph_view](https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/blob/master/resources/graph_view_screenshot.png)

## Table View
![graph_view](https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/blob/master/resources/table_view_screenshot.png)

---

## [Results](https://hidinginplainsight.shinyapps.io/HidingInPlaneSight/) from applying to gnomAD SNV VCFs:
![counts](https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/blob/master/resources/Num_SVs_by_chrom.jpeg)


---
## Getting Started
### For front end development -  download R and Rstudio - more info [here](https://www.ics.uci.edu/~sternh/courses/210/InstallingRandRStudio.pdf)

Once Rstudio is installed, import:
```
shiny
dplyr
shinyWidgets
devtools
```
To host the shiny app we will use:
```
https://www.shinyapps.io
```

### For back end development - follow the steps below to create your environment and install an editor.
I prefer the free version of Pycharm, found [here](https://www.jetbrains.com/pycharm/)

1. [Download miniconda](https://docs.conda.io/en/latest/miniconda.html) and choose Python 3.7 64-bit installer.
```
After downloading it make sure to open a new terminal and test with `which conda`.
```
For example:
```
MDMBASTYLIANOU:~ astylianou900045$ which conda
/Users/astylianou900045/miniconda3/bin/conda
```
2. Add conda channels - do this to get all the packages needed for this project:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Then check that everything installed correctly with: conda config --show channels

For example:
```
MDMBASTYLIANOU:~ astylianou900045$ conda config --show channels
channels:
  - conda-forge
  - bioconda
  - defaults
```
Now that conda is installed, clone the repo
```
3. git clone https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets.git
```
Cloning will make a new directory. After it is cloned, enter that directory with 'cd' 

For example:
```
4. cd Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets
```
Time to create and activate the environment for this project:

5. ```conda env create -f environment.yml```

6. `conda activate sv_env`

You will know that it was activated because sv_env will show in your command prompt, for example:
```
(sv_env) MDMBASTYLIANOU:~ astylianou900045$
```
7. Run test data in MEI-only mode. Should create a bed file with MEI calls.
```
python /path/to/repo/vcfToBed.py \
-inFile /path/to/repo/test/test.vcf.gz \
-outFile out.fa -chr 21 \
-dir /path/to/repo/resources -meiOnly True
```
Where:
* `-inFile` is a gzip'ed VCF with SNV calls in hg19.
* `-outFile` is the name of the output for the fastA
* `-chr` is the chromosome to work on (without leading 'chr')
* `-dir` is the directory where BLAST dictionaries are located. The MEI BLAST dictionary is located in the `resources` directory of this repo
* `-meiOnly` is whether to detect only MEI structural variants. 

## Running on other SV types (deletions, inversions, and tandem duplications)
1. Download human reference sequences for each chromosome for hg19
ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/README.txt

2. Make BLAST database for each chromosome
```
for f in *.fa.gz; do gunzip -k $f; done
for f in *.fa; do makeblastdb -in $f -dbtype nucl; done
```
3. Make sure your MEI BLAST database and chromosome BLAST databases are in the same directory

4. Run test data
```
python /path/to/repo/vcfToBed.py \
-inFile /path/to/repo/test/test.vcf.gz \
-outFile out.fa -chr 21 \
-dir /path/to/repo/resources -meiOnly ""
```
