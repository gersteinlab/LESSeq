# LESSeq
Local event-based analysis of alternative splicing using RNA-Seq data

## Overview
The Local Event-based analysis of alternative Splicing using RNA-Seq (or LESSeq) is a Linux-based processing pipeline for analyzing alternative splicing events from RNA-Seq data.

## Software requirements
1. Linux with g++ compiler
2. C++ Boost library: http://archive.gersteinlab.org/boost/
3. R: [http://cran.us.r-project.org/](http://cran.us.r-project.org/)

## Installation guide

Either download the package by clicking the "Clone or download" button, unzipping file in desired location, and renaming the directory "LESSeq" OR use the command line ``` git clone https://github.com/gersteinlab/LESSeq ```.

Add the following lines to your ~/.bashrc (without quotations):

1. "export LD_LIBRARY_PATH=PATH_TO_PACKAGE/gsl/lib/:PATH_TO_PACKAGE/cppunit/lib/"
2. "export PATH=PATH_TO_PACKAGE/bin/:$PATH"
     
Where PATH_TO_PACKAGE is the absolute path to the LESSeq/ folder.

### From pre-compiled executables
All C++ executables and R scripts required to run the LESSeq pipeline are found within LESSeq/bin/ folder.

### Compile from source
From the command line, enter the following commands:

1. ``` cd LESSeq/classify/ ```
2. ``` make ```

Repeat for LESSeq/count/ (steps 3-4) and LESSeq/solve/ (5-6) folders. Then move all executables to the LESSeq/bin folder:

7. ``` cd LESSeq/ ```
8. ``` mv classify/bin/* bin/ ``` 

Repeat for LESSeq/count/bin/ (steps 9-10) and LESSeq/solve/bin/ (11-12) folders.

### C++ Boost library

Download and unpackaged the library in the LESSeq/ folder:

1. ``` cd LESSeq/ ```
2. ``` wget http://archive.gersteinlab.org/boost/boost_1_34_1.zip ```
3. ``` unzip boost_1_34_1.zip ```

Finally in your ~/.bashrc file, add the following lines (without quotations):

1. "export LD_LIBRARY_PATH= PATH_TO_PACKAGE/gsl/lib/:PATH_TO_PACKAGE/cppunit/lib/:PATH_TO_PACKAGE/boost_1_34_1/boost/"
2. "export PATH=PATH_TO_PACKAGE/bin/:$PATH"

Where PATH_TO_PACKAGE is the absolute path to the LESSeq/ folder.

Then re-login or source your ~./bashrc file: ``` source ~/.bashrc ```

## Quick start

Below are instructions to run the four major steps described in the LESSeq manuscript.

### Step 1 --- Refining gene models using RNA-Seq data

The computational pipeline `Cufflinks' is used in this step of LESSeq. The Cufflinks executable is included within the LESSeq/ folder, but can be replaced with other tools should the user prefer.

   Cufflinks GitHub repo: [https://github.com/cole-trapnell-lab/cufflinks](https://github.com/cole-trapnell-lab/cufflinks)
   
Cufflinks manual: [http://cufflinks.cbcb.umd.edu/](http://cufflinks.cbcb.umd.edu/)

To run Cufflinks, type the following in the terminal: 
``` cufflinks```
	 
### Step 2 --- Identifying local events

To identify local alternative splicing events, type the following in the terminal:
```
classify log_level proj_name out_prefix isoform_format g2i_format g2i_path gene_begin_idx gene_end_idx
positional arguments:
		log_level			determines how much information to send to stdout while running
		proj_name			the name of the project given by the user
		out_prefix			directory name for output files
		isoform_format		file format for each form of the local events
		isoforms_path		path to the isoform_format file
		g2i_format			file format for grouping local event forms
		g2i_path			path to the g2i_format file
		gene_begin_idx		index of the first local event to be quantified
		gene_end_idx 		index of the last local event to be analyzed
```
Notes on positional arguments:
1. ``` log_level ``` should be an integer (e.g., 0, 1, 2)
2. ```out_prefix``` should already exist and string must include '/' at the end
, and the current choice is LH_GENE_TXT (which is the same as 
3. ``` isoform_format ``` is used to specify the coordinates of different forms of local events (default value--- LH_GENE_TXT). LH_GENE_TXT is equivalent to the `interval' format define [here](http://info.gersteinlab.org/RSEQtools#Interval).
4. ``` g2i_format ``` indicates which local event forms belong to the same local event (default value --- UCSC_GENE2ISOFORM). UCSC_GENE2ISOFORM is equivalent to the `knownIsoforms.txt' files [here](http://info.gersteinlab.org/RSEQtools#mergeTranscripts).
 read_format is the format of aligned reads, the current choice is MRF_SINGLE (single-end reads in MRF format [http://info.gersteinlab.org/RSEQtools#Mapped_Read_Format_.28MRF.29](http://info.gersteinlab.org/RSEQtools#Mapped_Read_Format_.28MRF.29))

``` classify ``` generates splicing graphs for each gene.  "Events.r" can then be used to generate local events from the splicing graphs.

The output from this step are eight sets of annotation files for the eight local event types specified in the LESSeq manuscript:
1. Skipped Exon (SE)
2. Retained Intron (RI)
3. Alternative 5' Splice Site (A5SS)
4. Alternative 3' Splice Site (A3SS)
5. Mutually eXclusive Exon (MXE) 
6. Alternative First Exon (AFE)
7. Alternative Last Exon (ALE)
8. Tandem 3' UTRs (T3)

For each local event type, two files are generated:
1. LH_GENE_TXT (or `interval') formatted file containing the annotation information of local event forms
2. UCSC_GENE2ISOFORM formatted file containing the grouping information of local event forms

Two executables are provided in the LESSeq/bin/ folder for generating the above two file formats from GTF/GFF files.
1. ``` parseGencode``` converts a GTF/GFF file generated by ``` cufflinks ``` into the LH_GENE_TXT (or 'interval') annotation file format. Its usage is as follows:
``` cat GTF/GFF_INPUT_FILE_NAME | parseGencode > OUTPUT_FILE_NAME.interval ``` 
2. ``` gencodeIsoformMap ``` converts the LH_GENE_TXT file to UCSC_GENE2ISOFORM grouping file format. Its usage is as follows:
``` cut -f1 OUTPUT_FILE_NAME.interval | gencodeIsoformMap > OUTPUT_FILE_NAME.map``` 
 
### Step 3 --- Counting reads compatible with alternative forms of local events and estimating their relative expression levels

The ``` count ``` command will associate raw read counts of reads compatible with alternative forms of local event.
``` 
count
positional arguments:
	log_level				determines how much information to send to stdout while running
	proj_name				the name of the project given by the user
	out_prefix				directory name for output files
	isoform_format			file format for each form of the local events
	isoforms_path			path to the isoform_format file
	g2i_format				file format for grouping local event forms
	g2i_path				path to the g2i_format file
	gene_begin_idx			index of the first local event to be quantified
	gene_end_idx 			index of the last local event to be analyzed
	read_formats			format of aligned reads
	read_type				the type of reads
	reads_paths				path to the read alignment file
	expected_read_lengths	the average read length
```
Notes on positional arguments:
1. ``` read_format ``` default value --- MRF_SINGLE, which is the single-end reads in MRF format described [here](http://info.gersteinlab.org/RSEQtools#Mapped_Read_Format_.28MRF.29)
2. ``` read_type ``` default value --- SHORT_READ

The output of ``` count ``` is a four column table, where columns are as follows:
1. grouping ID of local events
2. total number of reads mapped to a local event column
3. ID of a specific form of a local event column
4. number of reads compatible with the specific form of local event

Using the ``` solve ``` command will estimate the relative expression levels.
```
solve
positional arguments:
	log_level				determines how much information to send to stdout while running
	proj_name				the name of the project given by the user
	out_prefix				directory name for output files
	isoform_format			file format for each form of the local events
	isoforms_path			path to the isoform_format file
	g2i_format				file format for grouping local event forms
	g2i_path				path to the g2i_format file
	gene_begin_idx			index of the first local event to be quantified
	gene_end_idx 			index of the last local event to be analyzed
	read_formats			format of aligned reads
	read_type				the type of reads
	reads_paths				path to the read alignment file
	expected_read_lengths	the average read length
	total_read_bases		the total number of bases in the alignment file
```

The output of ```solve``` is a six column table, where columns are as follows:
1. grouping ID of local events
2. total number of reads mapped to a local event
3. ID for a specific form of a local event
4. relative expression level for the specific form of a local event
5. Reads Per Kilobase Million (RPKM) value for the specific form of a local event
6. log likelihood statistic of estimation method


For both ```count``` and ```solve``` commands, output will be printed to standard output (stdout), but can be directed to a file using '>':
```count arguments > count_output.txt```

### Step 4 --- Test differential alternative splicing events
Multiple R packages are required to run the log-linear model method found in the LESSeq/bin/Test_AS.r R script:

1. [epicalc](https://cran.r-project.org/src/contrib/Archive/epicalc/)
2. [lmtest](https://cran.r-project.org/web/packages/lmtest/index.html)
3. [MASS](https://cran.r-project.org/web/packages/MASS/MASS.pdf)
4. [multtest](https://www.bioconductor.org/packages/release/bioc/html/multtest.html)
5. [xtable](https://cran.r-project.org/web/packages/xtable/xtable.pdf)

To run LESSeq/bin/Test_AS.r:

``` Rscript LESSeq/bin/Test_AS.r ```

## License:

The MIT License

Copyright (c) 2019 Gerstein Lab, Mark B. Gerstein

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
