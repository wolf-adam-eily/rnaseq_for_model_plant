# RNA-Seq for Model Plant (Arabidopsis thaliana)

This repository is a usable, publicly available tutorial for analyzing differential expression data and creating topological gene networks. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools, such as nano, vim, or emacs.  If you are new to Linux, please use <a href="https://bioinformatics.uconn.edu/unix-basics/">this</a> handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as <a href="https://en.wikipedia.org/wiki/FASTA_format">FASTA</a>, <a href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a>, <a href="https://en.wikipedia.org/wiki/SAM_(file_format)">SAM/BAM</a>, and <a href="https://en.wikipedia.org/wiki/General_feature_format">GFF3/GTF</a>. You can learn even more about each file format <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>. If you do notp have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one <a href="https://bioinformatics.uconn.edu/contact-us/">here</a>.
	
<div id="toc_container">
<p class="toc_title">Contents</p>
<ul class="toc_list">
<li><a href="#First_Point_Header">1 Introduction and programs</>
<li><a href="#Second_Point_Header">2 Accessing the data using sra-toolkit</a></li>
<li><a href="#Third_Point_Header">3 Quality control using sickle</a></li>
<li><a href="#Fourth_Point_Header">4 Aligning reads to a genome using hisat2</a></li>
<li><a href="#Fifth_Point_Header">5 Transcript assembly and quantification with StringTie</a></li>
<li><a href="#Sixth_Point_Header">6 Differential expression analysis using ballgown</a></li>
<li><a href="#Seventh_Point_Header">Topological networking using cytoscape</a></li>
 <li><a href="#Integration">8 Integrating the DE Results with the Annotation Results</a></li>
<li><a href="#Citation">Citations</a></li>
</ul>
</div>

<h2 id="First_Point_Header">Introduction and Programs</h2>

In this tutorial, we will be analyzing thale cress (Arabidopsis thaliana) RNA-Seq data from various parts of the plant (roots, stems). Perhaps one of the most common organisms for genetic study, the aggregrate wealth of genetic information of the thale cress makes it ideal for new-comers to learn. Organisms such as this we call "model organisms". You may think of model organisms as a subset of living things which, under the normal conventions of analysis, behave nicely. The data we will be analyzing comes from an experiment in which various cellular RNA was collected from the roots and shoots of a single thale cress. The RNA profiles are archived in the SRA, and meta-information on each may be viewed through the SRA ID: <a href="https://www.ncbi.nlm.nih.gov/sra?term=SRX1756762">SRR3498212</a>, <a href="https://www.ncbi.nlm.nih.gov/sra/?term=SRR3498213">SRR3498213</a>, <a href="https://www.ncbi.nlm.nih.gov/sra?term=SRX1756765">SRR3498215</a>, <a href="https://www.ncbi.nlm.nih.gov/sra?term=SRX1756766">SRR3498216</a>.

The Single Read Archive, or SRA, is a publicly available database containing read sequences from a variety of experiments. Scientists who would like their read sequences present on the SRA submit a report containing the read sequences, experimental details, and any other accessory meta-data.

Our data, SRR3498212, SRR3498213, SRR3498215, SRR3498216 come from root 1, root 2, shoot 1, and shoot 2 of a single thale cress, respectively. Our objective in this analysis is to determine which genes are expressed in all samples, quantify the expression of each common gene in each sample, identify genes which are lowly expressed in roots 1 and 2 but highly expressed in shoots 1 and 2, or vice versa, quantify the relative expression of such genes, and lastly to create a visual topological network of genes with similar expression profiles.

Before beginning, we need to understand a few aspects of the Xanadu server. When first logging into Xanadu from your local terminal, you will be connected to the submit node. The submit node is the interface with which users on Xanadu may <i>submit</i> their processes to the desired compute nodes, which will run the process. Never, under any circumstance, run processes directly in the submit node. Your process will be killed and all of your work lost! This tutorial will not teach you shell script configuration to submit your tasks on Xanadu. Therefore, before moving on, read and master the topics covered in the <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/xanadu/">Xanadu tutorial</a>.

Now that we have covered the introduction and objective of our analysis, we may begin!

<h2 id="Second_Point_Header">Accessing the data using sra-toolkit </h2>

We know that the SRA contain the read sequences and accessory meta-information from experiments. Rather than downloading experimental data through a browser, we may use the <a href="https://www.ncbi.nlm.nih.gov/books/NBK158900/">sratoolkit</a>'s "fastq-dump" function to directly dump raw read data into the current terminal directory. Let's have a look at this function (it is expected that you have read the Xanadu tutorial, and are familiar with loading modules):

<pre style="color: silver; background: black;">-bash-4.2$ module load sratoolkit
  fastq-dump [options] <path> [<path>...]
  fastq-dump [options] <accession>

Use option --help for more information

fastq-dump : 2.8.2 </pre>

For our needs, we will simply be using the accession numbers to dump our experimental data into our directory. We know our accession numbers, so let's write a shell script to retrieve our raw reads. There are a variety of text editors available on Xanadu. My preferred text editor is "nano". Therefore, we will be using nano to write our shell script.

<pre style="color: silver; background: black;">-bash-4.2$ nano data_dump.sh

  GNU nano 2.3.1                                                      File: data_dump.sh                                                                                                                    

#!/bin/bash
#SBATCH --job-name=data_dump
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o data_dump_%j.out
#SBATCH -e data_dump_%j.err
#SBATCH --partition=general

mkdir /home/CAM/$USER/tmp/
export TMPDIR=/home/CAM/$USER/tmp/

module load sratoolkit

fastq-dump SRR3498212
fastq-dump SRR3498213
fastq-dump SRR3498215
fastq-dump SRR3498216
                                                                                             [ Read 20 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>

As a precautionary measure, always include your temporary directory in the environment. While not all programs require a temporary directory to work, it takes far less time including ours in the environment than it is waiting for an error! After typing our script, we press CTRL + X to exit, 'y', and then enter to save.

Now that we have our script saved, we submit it to the compute nodes with the following command:

<pre style="color: silver; background: black;">-bash-4.2$ sbatch data_dump.sh</pre>

Now we wait until we receive an email that our process has finished.

Let's take a look at one of our files:

<pre style="color: silver; background: black;">-bash-4.2$ head SRR3498212.fastq
@SRR3498212.1 SN638:767:HC555BCXX:1:1108:2396:1996 length=50
NTCAATCGGTCAGAGCACCGCCCTGTCAAGGCGGAAGCAGATCGGAAGAG
+SRR3498212.1 SN638:767:HC555BCXX:1:1108:2396:1996 length=50
#&#60;DDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR3498212.2 SN638:767:HC555BCXX:1:1108:2934:1998 length=50
NCAGTTTTGACCAGATAGGTCTCGCTAAGAAGATTGAGAAGATCGGAAGA
+SRR3498212.2 SN638:767:HC555BCXX:1:1108:2934:1998 length=50
#&#60;&#60;&#60;BDFFEEHIIIHIHEHEHIIEHHHH?&#60;GCCCHCCGICHH1&#60;GHHHIC
@SRR3498212.3 SN638:767:HC555BCXX:1:1108:3860:2000 length=50
NTCGCTTCGTAAGCGAAAGGCCGCGAGTTCGAAGATCGGAAGAGCACACG
</pre>

We see that for our first three runs we have information about the sampled read including its length followed by the nucleotide read and then a "+" sign. The "+" sign marks the beginning of the corresponding scores for each nucleotide read for the nucleotide sequence preceding the "+" sign. 

<h2 id="Third_Point_Header">Quality control using sickle</h2>

Sickle performs quality control on illumina paired-end and single-end short read data using a sliding window. As the window slides along the fastq file, the average score of all the reads contained in the window is calculated. Should the average window score fall beneath a set threshold, <a href="https://github.com/najoshi/sickle/blob/master/README.md">sickle</a> determines the reads responsible and removes them from the run. After visiting the SRA pages for our data, we see that our data are single end reads. Let's find out what sickle can do with these:

<pre style="color: silver; background: black;">-bash-4.2$ module load sickle

-bash-4.2$ sickle

Usage: sickle <command> [options]

Command:
pe	paired-end sequence trimming
se	single-end sequence trimming

--help, display this help and exit
--version, output version information and exit</pre>

We have single-end sequences. 

<pre style="color: silver; background: black;">-bash-4.2$ sickle se

Usage: sickle se [options] -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>

Options:
-f, --fastq-file, Input fastq file (required)
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)
-o, --output-file, Output trimmed fastq file (required)
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
-x, --no-fiveprime, Don't do five prime trimming.
-n, --trunc-n, Truncate sequences at position of first N.
-g, --gzip-output, Output gzipped files.
--quiet, Don't print out any trimming information
--help, display this help and exit
--version, output version information and exit</pre>

The quality may be any score from 0 to 40. The default of 20 is much too low for a robust analysis. We want to select only reads with a quality of 35 or better. Additionally, the desired length of each read is 50bp. Again, we see that a default of 20 is much too for analysis confidence. We want to select only reads whose length exceeds 45bp. Lastly, we must know the scoring type. While the quality type is not listed on the SRA pages, most SRA reads use the "sanger" quality type. Unless explicitly stated, try running sickle using the sanger qualities. If an error is returned, try illumina. If another error is returned, lastly try solexa.

Let's put all of this together for our sickle script using our downloaded fastq files:

<pre style="color: silver; background: black;">-bash-4.2$ nano sickle_run.sh

  GNU nano 2.3.1                                                      File: sickle_run.sh                                                                                                                   

#!/bin/bash
#SBATCH --job-name=sickle_run
#SBATCH --mail-user=wolf.adam.eily@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o sickle_run_%j.out
#SBATCH -e sickle_run_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/aeily/tmp/

module load sickle

sickle se -f SRR3498212.fastq -t sanger -o trimmed_SRR3498212.fastq -q 35 -l 45

sickle se -f SRR3498213.fastq -t sanger -o trimmed_SRR3498213.fastq -q 35 -l 45

sickle se -f SRR3498215.fastq -t sanger -o trimmed_SRR3498215.fastq -q 35 -l 45

sickle se -f SRR3498216.fastq -t sanger -o trimmed_SRR3498216.fastq -q 35 -l 45

                                                                                             [ Read 23 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>
<br>
<pre style="color: silver; background: black;">-bash-4.2$ sbatch sickle_run.sh</pre>

It is helpful to see how the quality of the data has changed after using sickle. To do this, we will be using the commandline versions of <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt">fastqc</a> and <a href="http://multiqc.info/docs/">MultiQC</a>. These two programs simply create reports of the average quality of our trimmed reads, with some graphs. There is no way to view a --help menu for these programs in the command-line. However, their use is quite simple, we simply run "fastqc <trimmed_fastq>" or "multiqc -f -n trimmed trimmed*". Do not worry too much about the options for MultiQC! Let's write our script:

<pre style="color: silver; background: black;">-bash-4.2$ nano quality_control.sh

  GNU nano 2.3.1                                                    File: quality_control.sh                                                                                                                

#!/bin/bash
#SBATCH --job-name=quality_control
#SBATCH --mail-user=wolf.adam.eily@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o quality_control_%j.out
#SBATCH -e quality_control_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/aeily/tmp/

module load fastqc
module load MultiQC

fastqc trimmed_SRR3498212.fastq

fastqc trimmed_SRR3498213.fastq

fastqc trimmed_SRR3498215.fastq

fastqc trimmed_SRR3498216.fastq

multiqc -f -n trimmed trimmed*


                                                                                             [ Read 26 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>
<br>
<pre style="color: silver; background: black;">-bash-4.2$ sbatch quality_control.sh</pre>

fastqc will create the files "trimmed_file_fastqc.html". To have a look at one, we need to move all of our "trimmed_file_fastqc.html" files moved into a single directory, and then <a href="https://www.techrepublic.com/article/how-to-use-secure-copy-for-file-transfer/">secure copy</a> that folder to our local directory. Then, we may open our files! If that seems like too much work for you, you may open the files directly through this github. Simply click on any "html" file and you may view it in your browser immediately. Because of this, the steps mentioned above will not be placed in this tutorial.

This script will create a directory "trimmed_data". Let's look inside of that directory:

<pre style="color: silver; background: black;">-bash-4.2$ cd trimmed_data</pre>


<h2 id="Fourth_Point_Header">Aligning reads to a genome using hisat2</h2>
<a href="https://ccb.jhu.edu/software/hisat2/manual.shtml">HISAT2</a> is a fast and sensitive aligner for mapping next generation sequencing reads against a reference genome. HISAT2 requires two arguments: the reads file being mapped and the indexed genome to which those reads are mapped. Typically, the hisat2-build command is used to make a HISAT index file for the genome. It will create a set of files with the suffix .ht2, these files together build the index. What is an index and why is it helpful? Genome indexing is the same as indexing a tome, like an encyclopedia. It is much easier to locate information in the vastness of an encyclopedia when you consult the index, which is ordered in an easily navigatable way with pointers to the location of the information you seek within the encylopedia. Genome indexing is thus the structuring of a genome such that it is ordered in an easily navigatable way with pointers to where we can find whichever gene is being aligned. Let's have a look at how the hisat2-build command works:

<pre style="color: silver; background: black;">-bash-4.2$ module load hisat2
bash-4.2$ hisat2-build

No input sequence or sequence file specified!
HISAT2 version 2.1.0 by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo)
Usage: hisat2-build [options]* <reference_in> <ht2_index_base>
    reference_in            comma-separated list of files with ref sequences
    hisat2_index_base       write ht2 data to files with this dir/basename</pre>

As you can see, we simply enter our reference genome files and the desired prefix for our .ht2 files. Now, fortunately for us, Xanadu has many indexed genomes which we may use. To see if there is a hisat2 <i>Arabidopsis thaliana</i> indexed genome we need to look at the <a href="https://bioinformatics.uconn.edu/databases/">Xanadu databases</a> page. We see that our desired indexed genome is in the location /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/. Now we are ready to align our reads using hisat2 (for hisat2, the script is going to be written first with an explanation of the options after).

<pre style="color: silver; background: black;">nano hisat2_run.sh

  GNU nano 2.3.1                                                      File: hisat2_run.sh                                                                                                                   

#!/bin/bash
#SBATCH --job-name=hisat2_run
#SBATCH --mail-user=wolf.adam.eily@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o hisat2_run_%j.out
#SBATCH -e hisat2_run_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/aeily/tmp/

module load hisat2

hisat2 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q trimmed_SRR3498212.fastq -S athaliana_root_1.sam

hisat2 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q trimmed_SRR3498213.fastq -S athaliana_root_2.sam

hisat2 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q trimmed_SRR3498215.fastq -S athaliana_shoot_1.sam

hisat2 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q trimmed_SRR3498216.fastq -S athaliana_shoot_2.sam

                                                                                             [ Read 23 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>
<br>
<pre style="color: silver; background: black;">-p : number of processors been used
--dta: report alignments tailored for transcript assemblers
-x: path to index generated from previous step
-q: query input files in fastq format
-S: output SAM file</pre>
<br>
<pre style="color: silver; background: black;">bash-4.2$ sbatch hisat2_run.sh</pre>

Once the mapping have been completed, the file structure is as follows:
<pre style="color: silver; background: black;">bash-4.2$ ls
|-- athaliana_root_1.sam
|-- athaliana_root_2.sam
|-- athaliana_shoot_1.sam
|-- athaliana_shoot_2.sam</pre>

When HISAT2 completes its run, it will summarize each of it’s alignments, and it is written to the standard error file, which can be found in the same folder once the run is completed.

<pre style="color: silver; background: black;">bash-4.2$ nano hisat2_run&#42;err

  GNU nano 2.3.1                                                   File: hisat2_run_507187.err                                                                                                              

34475799 reads; of these:
  34475799 (100.00%) were unpaired; of these:
    33017550 (95.77%) aligned 0 times
    1065637 (3.09%) aligned exactly 1 time
    392612 (1.14%) aligned >1 times
4.23% overall alignment rate
42033973 reads; of these:
  42033973 (100.00%) were unpaired; of these:
    40774230 (97.00%) aligned 0 times
    931377 (2.22%) aligned exactly 1 time
    328366 (0.78%) aligned >1 times
3.00% overall alignment rate
31671127 reads; of these:
  31671127 (100.00%) were unpaired; of these:
    31103167 (98.21%) aligned 0 times
    465131 (1.47%) aligned exactly 1 time
    102829 (0.32%) aligned >1 times
1.79% overall alignment rate
49890217 reads; of these:
  49890217 (100.00%) were unpaired; of these:
    48622480 (97.46%) aligned 0 times
    1029943 (2.06%) aligned exactly 1 time
    237794 (0.48%) aligned >1 times
2.54% overall alignment rate
</pre>
<br>
Let's have a look at a SAM file:

<pre style="color: silver; background: black;">-bash-4.2$ head -n 20 rnaseq_athaliana_root_1.sam
@HD	VN:1.0	SO:unsorted
@SQ	SN:Chr1	LN:30427671
@SQ	SN:Chr2	LN:19698289
@SQ	SN:Chr3	LN:23459830
@SQ	SN:Chr4	LN:18585056
@SQ	SN:Chr5	LN:26975502
@SQ	SN:ChrM	LN:366924
@SQ	SN:ChrC	LN:154478
@PG	ID:hisat2	PN:hisat2	VN:2.1.0	CL:"/isg/shared/apps/hisat2/2.1.0/hisat2-align-s --wrapper basic-0 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q trimmed_SRR3498212.fastq -S rnaseq_athaliana_root_1.sam"
SRR3498212.6	4	*	0	0	*	*	0	0	TTTCCAAGCCCTTTCTAGTCTGCGCTTGAGTTTGATTGCAGAGATCGGAA	DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.1	4	*	0	0	*	*	0	0	CAATCGGTCAGAGCACCGCCCTGTCAAGGCGGAAGCAGATCGGAAGAG	DDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.4	4	*	0	0	*	*	0	0	AAAGGGCGTGGGTTCAAATCCCACAGATGTCACCAGATCGGAAGAGC	DDHIIIIIIIIIEHHHIHIIIIHIIIIIIIIIIIIIIIIIIIIIIHH	YT:Z:UU
SRR3498212.8	4	*	0	0	*	*	0	0	TTAAGATTGCTGATTTTGGCCTGGCACGTGAGGTTAAGATCGGAAGAGCA	DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.19	4	*	0	0	*	*	0	0	TGGATGATGGAAAAACCAGCAAGCCCCTCTTCTTTCAAGATCGGAAGAGC	DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.23	4	*	0	0	*	*	0	0	TTTGCCTTCCAAGCAATAGACCCGGGTAGATCGGAAGAGCACACGTCTGA	DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.24	4	*	0	0	*	*	0	0	TGAAACTTCTTGGTTTTAAAGTGTGAATATAGCTGACAAAAGATTGGAAG	DDDDDIIIIIIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.12	4	*	0	0	*	*	0	0	AAGGGTGTTCTCTGCTACGGACCTCCAGATCGGAAGAGCACACGTCTGAA	DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.27	4	*	0	0	*	*	0	0	ATTGTTCCGGGCTGCCCAGTCCAAGCTGAGAGTGAAGATCGGAAGAGCAC	DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.29	4	*	0	0	*	*	0	0	TATGTCTACGCTGGTTCAAATCCAGCTCGGCCCACCAAGATCGGAAGAGC	DDDDDIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	YT:Z:UU
SRR3498212.18	4	*	0	0	*	*	0	0	CGTGGGTTCGACTCCCACTGTGGTCGCCAAGATCGGAAGAGCACACGTC	DDDCHCCHHHEIHIGIIIEGHHIIIIGHHHIIIIIIIIIIIIIIIIIII	YT:Z:UU
</pre>

All of the lines starting with an "@" symbol tell us something about the chromosomes or our input. For instance "@SQ SN:Chr1 LN:30427671" tells us that we have a sequence (@SQ) whose sequence name is Chr1 (SN:Chr1), lastly the sequence has a length of 30427671bp (LN:30427671). You may be wondering what the first line means. It is quite straightfoward! The first line is simply the header (@HD) statig that the file is unsorted (SO:unsorted). The second column in the first line is somewhat of a dummy variable, but stands for "version number". Lastly we have the "@PG" line, which, in order, keeps track of the software used to write the file (ID:hisat2), the program name used to align the reads (PN:hisat2), the version of the program used (VN:2.1.0), and lastly the user input which started the process (written in the form that the program reads, not in which we wrote it).

The alignment portion of the SAM file is much more straight-forward and may be understood by reading the SAM output formatting guide linked in the beginning of this tutorial.

Because of the density of the sam file, it is compressed to binary to create a more easily tractable file for manipulation by future programs. We convert the sam file to b<sub>inary</sub>am with the following command and sort it such that the alignments are listed in the order the genes appear in the genome. To do this we use the software <a href="https://en.wikipedia.org/wiki/SAMtools">samtools</a>

<pre style="color: silver; background: black;">-bash-4.2$ module load samtools
bash-4.2$ samtools
Usage:   samtools <command> [options]

Commands:
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate information
     reheader       replace BAM header
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags
     markdup        mark duplicates

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA

  -- Statistics
     bedcov         read depth per BED region
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)

  -- Viewing
     flags          explain BAM flags
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM
</pre>

We are truly only interested in sorting our SAM files.

<pre style="color: silver; background: black;">-bash-4.2$ samtools sort

Usage: samtools sort [options...] [in.bam]
Options:
  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
  -n         Sort by read name
  -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
  -o FILE    Write final output to FILE rather than standard output
  -T PREFIX  Write temporary files to PREFIX.nnnn.bam
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional threads to use [0]
</pre>

The sort function converts SAM files to BAM automatically. Therefore, we can cut through most of these options and do a simple "samtools sort -o <output.bam> <inupt.sam>. Let's write our script:

<pre style="color: silver; background: black;">bash-4.2$ nano sam_sort_bam.sh

  GNU nano 2.3.1                                                     File: sam_sort_bam.sh                                                                                                                  

#!/bin/bash
#SBATCH --job-name=sam_sort_bam
#SBATCH --mail-user=wolf.adam.eily@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o sam_sort_bam_%j.out
#SBATCH -e sam_sort_bam_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/aeily/tmp/

module load samtools

samtools sort -@ 16 -o athaliana_root_1.bam athaliana_root_1.sam

samtools sort -@ 16 -o athaliana_root_2.bam athaliana_root_2.sam

samtools sort -@ 16 -o athaliana_shoot_1.bam athaliana_shoot_1.sam

samtools sort -@ 16 -o athaliana_shoot_2.bam athaliana_shoot_2.sam
                                                                                             [ Read 23 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>

<pre style="color: silver; background: black;">bash-4.2$ sbatch sam_sort_bam.sh</pre>

<h2 id="Fifth_Point_Header">Transcript assembly and quantification with StringTie</h2>

This is the most intricate part of our analysis. Because of the many steps involved, this tutorial is going to walk through the steps before writing a final batch script to be submitted. 

Analysis of RNA-seq experiments requires accurate reconstructions of all the <a href="https://en.wikipedia.org/wiki/Gene_isoform">isoforms</a> expressed from each gene, as well as estimates of the relative abundance of those isoforms. Accurate quantification benefits from knowledge of precisely which reads originated from each isoform, which cannot be computed perfectly because reads are much shorter than transcripts. StringTie assembles transcripts from RNA-seq reads that have been aligned to the genome, first grouping the reads into distinct gene loci and then assembling each locus into as many isoforms as are needed to explain the data. To begin, it would be nice to have a source describing all genome features (simply a list of genes, exons, transcripts, etc., which states on which chromosome and over which base-pairs those genes are). The file format for which we're looking is GFF (which is exactly as described the sentence prior). We can download the GFF file for the thale cress with the following code:
  
<pre style="color: silver; background: black;">bash-4.2$ wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff

bash-4.2$ head TAIR_GFF3_genes.gff
Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1
Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1	TAIR10	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
Chr1	TAIR10	protein	3760	5630	.	+	.	ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	five_prime_UTR	3631	3759	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	3760	3913	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	3996	4276	.	+	2	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1
</pre>

The GFF file is quite self-explanatory. However, it'd be nice if could combine all of the pieces of information from the GFF into something better. For instance, if there are multiple overlapping, but distinct exons from a single gene, we could use that information to determine the isoforms of that gene. Then, we could make a file which gives each isoform its own track (there are other extrapolations to be made, but this is our most relevant example). Luckily for us, we can use the program "gffread" to transform our GFF file into the more useful form just stated, The output of <a href="https://github.com/gpertea/gffread">gffread --help</a> is much too dense for us to go into here, but the necessary options will be explained. Do not run this code! We are compiling this code with various other chunks into one script, be patient!

<pre style="color: silver; background: black;">bash-4.2$ module load gffread
gffread TAIR10_GFF3_genes.gff -T -o athaliana_TAIR10_genes.gtf</pre>

The option -T tells gffread to convert our input into the gtf format, and the option -o simply is how we call the output. The GTF format is simply the transcript assembly file, and is composed of exons and coding sequences. Let's have a look at the GTF file:

<pre style="color: silver; background: black;">-bash-4.2$ head athaliana_TAIR10_genes.gtf 
Chr1	TAIR10	exon	3631	3913	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	exon	3996	4276	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	exon	4486	4605	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	exon	4706	5095	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	exon	5174	5326	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	exon	5439	5899	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	CDS	3760	3913	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	CDS	3996	4276	.	+	2	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	CDS	4486	4605	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1	TAIR10	CDS	4706	5095	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";

-bash-4.2$ tail athaliana_TAIR10_genes.gtf 
ChrM	TAIR10	exon	349830	351413	.	-	.	transcript_id "ATMG01360.1"; gene_id "ATMG01360"; gene_name "ATMG01360";
ChrM	TAIR10	CDS	349830	351413	.	-	0	transcript_id "ATMG01360.1"; gene_id "ATMG01360"; gene_name "ATMG01360";
ChrM	TAIR10	exon	360717	361052	.	-	.	transcript_id "ATMG01370.1"; gene_id "ATMG01370"; gene_name "ATMG01370";
ChrM	TAIR10	CDS	360717	361052	.	-	0	transcript_id "ATMG01370.1"; gene_id "ATMG01370"; gene_name "ATMG01370";
ChrM	TAIR10	exon	361062	361179	.	-	.	transcript_id "ATMG01380.1"; gene_id "ATMG01380"; gene_name "ATMG01380";
ChrM	TAIR10	exon	361350	363284	.	-	.	transcript_id "ATMG01390.1"; gene_id "ATMG01390"; gene_name "ATMG01390";
ChrM	TAIR10	exon	363725	364042	.	+	.	transcript_id "ATMG01400.1"; gene_id "ATMG01400"; gene_name "ATMG01400";
ChrM	TAIR10	CDS	363725	364042	.	+	0	transcript_id "ATMG01400.1"; gene_id "ATMG01400"; gene_name "ATMG01400";
ChrM	TAIR10	exon	366086	366700	.	-	.	transcript_id "ATMG01410.1"; gene_id "ATMG01410"; gene_name "ATMG01410";
ChrM	TAIR10	CDS	366086	366700	.	-	0	transcript_id "ATMG01410.1"; gene_id "ATMG01410"; gene_name "ATMG01410";
</pre>

We see that whereas in our GFF file we have various untranslated regions included, as well as annotations, the GTF format contains information only on various transcripts for each gene. The "transcript_id" denoter in the last column tells us the gene and its isoform, and everything else about the GTF file is quite apparent!

Just as was stated for our conversion from gff to gtf, it would be helpful for us to perform the same operation on our aligned reads. That is, if there are multiple, overlapping but distinct reads from a single gene, we could combine these reads into one transcript isoform. Because we have the gene isoforms in the gtf file, we can re-map each assembled transcript to a gene isoform and then count how many mappings there are per isoform. This, in effect, allows us to quantify the expression rates of each isoform. We will be using the program <a href="http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual">StringTie</a> to assemble the transcripts for each sample. StringTie requires three input arguments: the BAM alignment file, the genomic GTF file, and the desired output GTF filename. Thus, our code will look like (do not run this!):

<pre style="color: silver; background: black;">bash-4.2$ module load stringtie

stringtie -p 16 athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_1.gtf

stringtie -p 16 athaliana_root_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_2.gtf

stringtie -p 16 athaliana_shoot_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_1.gtf

stringtie -p 16 athaliana_shoot_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_2.gtf</pre>

Let's have a look at the stringtie output:

<pre style="color: silver; background: black;">-bash-4.2$ athaliana_root_1.gtf
# stringtie -p 16 rnaseq_athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o rnaseq_athaliana_root_1.gtf
# StringTie version 1.3.3b
Chr1	StringTie	transcript	28500	28706	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; reference_id "AT1G01046.1"; ref_gene_id "AT1G01046"; ref_gene_name "AT1G01046"; cov "0.241546"; FPKM "3.727008"; TPM "0.747930";
Chr1	StringTie	exon	28500	28706	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; reference_id "AT1G01046.1"; ref_gene_id "AT1G01046"; ref_gene_name "AT1G01046"; cov "0.241546";
Chr1	StringTie	transcript	47494	48839	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "3.928230"; FPKM "60.611832"; TPM "12.163484";
Chr1	StringTie	exon	47494	47982	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "4.529652";
Chr1	StringTie	exon	48075	48839	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "2"; cov "3.543791";
Chr1	StringTie	transcript	50075	51199	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; reference_id "AT1G01100.2"; ref_gene_id "AT1G01100"; ref_gene_name "AT1G01100"; cov "8.437494"; FPKM "130.188904"; TPM "26.126097";
Chr1	StringTie	exon	50075	50337	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; reference_id "AT1G01100.2"; ref_gene_id "AT1G01100"; ref_gene_name "AT1G01100"; cov "6.228601";
Chr1	StringTie	exon	50419	50631	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; reference_id "AT1G01100.2"; ref_gene_id "AT1G01100"; ref_gene_name "AT1G01100"; cov "9.487524";
</pre>

While this is certainly confusing, we can still understand it. To start each row we have the chromosome for each sequence. The second column is the software used to assemble the transcript, in our case StringTie. The third column is the sequence type, either a transcript or an exon. The next two columns are the start and end bp of the feature (assuming the chromosome starts at bp 1), followed by another column which is the score. The column after the score is either "+" or "-" for forward strand and reverse strand, respectively. Our last two columns are the frame ("." means frame not determined) and the feature meta-information. 

You may be wondering about the last three columns which are not standard in GTF. The column "cov" is simply the covariance of the gene across samples (if the gene is highly or lowly expressed in both samples the covariance will be high, if it is highly expressed in one sample and lowly expressed in another sample, the covariance will be low). The FPKM column is the fragment per kilobase million value. Simply put, this is the number of times the specific feature was counted divided by how many counts there were for all features combined, in millions. That number is then divided by the length of the feature in kilobases. The reason for this being that longer features should have higher counts. For instance, when we split our mRNA in sequences of 50 or less for reading, one 5000 bp transcript will appear as 100 counts, and one 500 bp transcript will appear as 10 counts. Now, let's divide each count by the transcript length in kilobases: we have 20 as the value for the 5000 bp sequence (100/(5000/1000)) and 20 as the value for the 500 bp sequence (10/(500/1000)! Now our quantification matches the actual expression profile -- that is, both features were transcribed the same amount.

The last column, TPM, is the transcripts per feature per  million transcripts counted across all features combined. As evident from the example above, without some sort of scaling factor, this value is highly misleading.

The assembled isoforms in each sample are most likely different from those of other samples. However, we may repeat the process of determining isoforms but this time using the gtf files for all four samples. That is, if there are isoforms across the sample gtfs which overlap but are distinct, we may merge those into a single isoform. To do this we will be using the --merge option of stringtie. The --merge option of stringtie requires three input arguments: the genomic GTF file, the desired output filename, and a plain-text file containing each file to be merged separated by a return. To begin, let's make our plain-text file:

<pre style="color: silver; background: black;">bash-4.2$ nano mergelist.txt

  GNU nano 2.3.1                                                      File: mergelist.txt                                                                                                                   

athaliana_root_1.gtf
athaliana_root_2.gtf
athaliana_shoot_1.gtf
athaliana_shoot_2.gtf
                                                                                              [ Read 4 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell

</pre>

After saving our plain-text file, we can merge our samples using the following code (do not run this!):

<pre style="color: silver; background: black;">bash-4.2$ module load stringtie
stringtie --merge -p 16 -G -o merged stringtie_merged.gtf</pre>

While the options are quite self-explanatory, one thing to note is that the "-o" option is simply the output <i>prefix</i>. After running, you should see the following files in your directory (but hopefully you listened and did not run it yet!):

<pre style="color: silver; background: black;">bash-4.2$ ls
merged.annotated.gtf
merged.loci
merged.stats
merged.stringtie_merged.gtf.refmap
merged.stringtie_merged.gtf.tmap
merged.tracking</pre>

Although you have not run the code yet, let's have a look at some of the files we've generated (we will not be looking at the merged.annotated.gtf as we are already quite familiar with the gtf format!)

<pre style="color: silver; background: black;">bash-4.2$ -bash-4.2$ head merged.loci
XLOC_000001	Chr1[+]3631-5899	AT1G01010|AT1G01010.1	AT1G01010.1
XLOC_000002	Chr1[+]23146-31227	AT1G01040|AT1G01040.1,AT1G01040|AT1G01040.2	AT1G01040.1,AT1G01040.2
XLOC_000003	Chr1[+]28500-28706	AT1G01046|AT1G01046.1	AT1G01046.1
XLOC_000004	Chr1[+]44677-44787	AT1G01073|AT1G01073.1	AT1G01073.1
XLOC_000005	Chr1[+]52239-54692	AT1G01110|AT1G01110.2,AT1G01110|AT1G01110.1	AT1G01110.2,AT1G01110.1
XLOC_000006	Chr1[+]56624-56740	AT1G01115|AT1G01115.1	AT1G01115.1
XLOC_000007	Chr1[+]72339-74096	AT1G01160|AT1G01160.1,AT1G01160|AT1G01160.2	AT1G01160.1,AT1G01160.2
XLOC_000008	Chr1[+]75583-76758	AT1G01180|AT1G01180.1	AT1G01180.1
XLOC_000009	Chr1[+]88898-89745	AT1G01210|AT1G01210.1	AT1G01210.1
XLOC_000010	Chr1[+]91376-95651	AT1G01220|AT1G01220.1	AT1G01220.1</pre>

We see we have a condensed form of our various exons. The exon loci name is the first column, followed by the chromosome, strand, and bp location in the second column. The final columns are the gene ID, the transcript IDs to which the loci belong, and the isoforms to which the transcripts belong.

<pre style="color: silver; background: black;">-bash-4.2$ nano merged.stats

  GNU nano 2.3.1                                                      File: merged.stats                                                                                                                    

# gffcompare v0.10.4 | Command line was:
#gffcompare -r athaliana_TAIR10_genes.gtf -G -o merged stringtie_merged.gtf
#

#= Summary for dataset: stringtie_merged.gtf
#     Query mRNAs :   41854 in   33403 loci  (30272 multi-exon transcripts)
#            (6013 multi-transcript loci, ~1.3 transcripts per locus)
# Reference mRNAs :   41607 in   33350 loci  (30127 multi-exon)
# Super-loci w/ reference transcripts:    33184
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    99.9    |
        Exon level:   100.0     |    99.9    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |    99.5    |
  Transcript level:   100.0     |    99.4    |
       Locus level:   100.0     |    99.8    |

     Matching intron chains:   30127
       Matching transcripts:   41607
              Matching loci:   33350

          Missed exons:       0/169264  (  0.0%)
           Novel exons:      90/169578  (  0.1%)
        Missed introns:       0/127896  (  0.0%)
         Novel introns:      13/127951  (  0.0%)
           Missed loci:       0/33350   (  0.0%)
            Novel loci:      75/33403   (  0.2%)

 Total union super-loci across all input datasets: 33403
41854 out of 41854 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)

                                                                                             [ Read 30 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>

The information here is quite apparent.

<pre style="color: silver; background: black;">-bash-4.2$ head merged.stringtie_merged.gtf.refmap
ref_gene_id	ref_id	class_code	qry_id_list
AT1G01010	AT1G01010.1	=	AT1G01010|AT1G01010.1
AT1G01040	AT1G01040.1	=	AT1G01040|AT1G01040.1
AT1G01040	AT1G01040.2	=	AT1G01040|AT1G01040.2
AT1G01046	AT1G01046.1	=	AT1G01046|AT1G01046.1
AT1G01073	AT1G01073.1	=	AT1G01073|AT1G01073.1
AT1G01110	AT1G01110.2	=	AT1G01110|AT1G01110.2
AT1G01110	AT1G01110.1	=	AT1G01110|AT1G01110.1
AT1G01115	AT1G01115.1	=	AT1G01115|AT1G01115.1
AT1G01160	AT1G01160.1	=	AT1G01160|AT1G01160.1
</pre>

Here we have the gene IDs from the reference GFF, followed by the isoform IDs from the reference GTF, the "class_code" which simply tells you that the last column was matched fully ("=") or partially ("c"). Taking our first column, we see that all of isoform ATG01010.1 was matched to the gene ATG01010.

<pre style="color: silver; background: black;">-bash-4.2$ head merged.stringtie_merged.gtf.tmap
ref_gene_id	ref_id	class_code	qry_gene_id	qry_id	num_exons	FPKM	TPM		cov	len	major_iso_id	ref_match_len
AT1G01090	AT1G01090.1	=	AT1G01090	AT1G01090.1	3	0.000000	0.000000	0.000000	1627	AT1G01090.1	1627
AT1G01100	AT1G01100.2	=	AT1G01100	AT1G01100.2	4	0.000000	0.000000	0.000000	631	AT1G01100.2	631
AT1G01100	AT1G01100.1	=	AT1G01100	AT1G01100.1	4	0.000000	0.000000	0.000000	587	AT1G01100.2	587
AT1G01100	AT1G01100.4	=	AT1G01100	AT1G01100.4	4	0.000000	0.000000	0.000000	607	AT1G01100.2	607
AT1G01100	AT1G01100.3	=	AT1G01100	AT1G01100.3	5	0.000000	0.000000	0.000000	566	AT1G01100.2	566
AT1G01120	AT1G01120.1	=	AT1G01120	AT1G01120.1	1	0.000000	0.000000	0.000000	1899	AT1G01120.1	1899
AT1G01180	AT1G01180.1	=	AT1G01180	AT1G01180.1	1	0.000000	0.000000	0.000000	1176	AT1G01180.1	1176
AT1G01183	AT1G01183.1	=	AT1G01183	AT1G01183.1	1	0.000000	0.000000	0.000000	101	AT1G01183.1	101
AT1G01225	AT1G01225.1	=	AT1G01225	AT1G01225.1	2	0.000000	0.000000	0.000000	1025	AT1G01225.1	1025
</pre>

All of the information in the .tmap file may be readily understood now, knowing that "len" and "ref_match_len" are the sequence lengths and reference lengths, respectively.

Lastly,
<pre style="color: silver; background: black;">-bash-4.2$ head merged.tracking
TCONS_00000001	XLOC_000001	AT1G01010|AT1G01010.1	=	q1:AT1G01010|AT1G01010.1|6|0.000000|0.000000|0.000000|1688
TCONS_00000002	XLOC_000002	AT1G01040|AT1G01040.1	=	q1:AT1G01040|AT1G01040.1|20|0.000000|0.000000|0.000000|6251
TCONS_00000003	XLOC_000002	AT1G01040|AT1G01040.2	=	q1:AT1G01040|AT1G01040.2|20|0.000000|0.000000|0.000000|5877
TCONS_00000004	XLOC_000003	AT1G01046|AT1G01046.1	=	q1:AT1G01046|AT1G01046.1|1|0.000000|0.000000|0.000000|207
TCONS_00000005	XLOC_000004	AT1G01073|AT1G01073.1	=	q1:AT1G01073|AT1G01073.1|1|0.000000|0.000000|0.000000|111
TCONS_00000006	XLOC_000005	AT1G01110|AT1G01110.2	=	q1:AT1G01110|AT1G01110.2|5|0.000000|0.000000|0.000000|1782
TCONS_00000007	XLOC_000005	AT1G01110|AT1G01110.1	=	q1:AT1G01110|AT1G01110.1|3|0.000000|0.000000|0.000000|1439
TCONS_00000008	XLOC_000006	AT1G01115|AT1G01115.1	=	q1:AT1G01115|AT1G01115.1|1|0.000000|0.000000|0.000000|117
TCONS_00000009	XLOC_000007	AT1G01160|AT1G01160.1	=	q1:AT1G01160|AT1G01160.1|5|0.000000|0.000000|0.000000|1045
TCONS_00000010	XLOC_000007	AT1G01160|AT1G01160.2	=	q1:AT1G01160|AT1G01160.2|6|0.000000|0.000000|0.000000|1129
</pre>

the merged.tracking is the compact form of tmap file combined with the loci file. The last column has the gene, gene isoform, number of hits, FPKM, TPM, cov, and length all in that order.

Our last step is to create a count-table for the differential expression software "ballgown". A word of note about ballgown, ballgown requires that all of the files it will analyze be in their own parent directory, "ballgown", and furthermore each file is in its own subdirectory "ballgown/file_sub/file". Knowing we have four files, let's create the required structure for ballgown (run this):

<pre style="color: silver; background: black;">-bash-4.2$ mkdir ballgown
-bash-4.2$ cd ballgown
-bash-4.2$ mkdir athaliana_root_1
-bash-4.2$ mkdir athaliana_root_2
-bash-4.2$ mkdir athaliana_shoot_1
-bash-4.2$ mkdir athaliana_shoot_2
-bash-4.2$ cd -</pre>

The reason for this will become obvious soon. Now we will use StringTie to make our count tables (do not run this!):

<pre style="color: silver; background: black;">-bash-4.2$ module load stringtie

stringtie -e -B -p 16 athaliana_root_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_1/athaliana_root_1.count

stringtie -e -B -p 16 athaliana_root_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_2/athaliana_root_2.count

stringtie -e -B -p 16 athaliana_shoot_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_1/athaliana_shoot_1.count

stringtie -e -B -p 16 athaliana_shoot_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_2/athaliana_shoot_2.count
</pre>

Now we are ready to compile all of our code into a single script:

<pre style="color: silver; background: black;">-bash-4.2$ nano transcript_assembly.sh

 GNU nano 2.3.1                                                   File: transcript_assembly.sh                                                                                                             

#!/bin/bash
#SBATCH --job-name=sam_sort_bam
#SBATCH --mail-user=wolf.adam.eily@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o sam_sort_bam_%j.out
#SBATCH -e sam_sort_bam_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/aeily/tmp/

module load gffread
module load stringtie

gffread TAIR10_GFF3_genes.gff -T -o athaliana_TAIR10_genes.gtf

stringtie -p 16 athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_1.gtf

stringtie -p 16 athaliana_root_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_2.gtf

stringtie -p 16 athaliana_shoot_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_1.gtf

stringtie -p 16 athaliana_shoot_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_2.gtf

stringtie --merge -p 16 -G -o merged stringtie_merged.gtf

stringtie -e -B -p 16 athaliana_root_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_1/athaliana_root_1.count

stringtie -e -B -p 16 athaliana_root_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_2/athaliana_root_2.count

stringtie -e -B -p 16 athaliana_shoot_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_1/athaliana_shoot_1.count

stringtie -e -B -p 16 athaliana_shoot_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_2/athaliana_shoot_2.count

                                                                                             [ Read 36 lines ]
^G Get Help                       ^O WriteOut                       ^R Read File                      ^Y Prev Page                      ^K Cut Text                       ^C Cur Pos
^X Exit                           ^J Justify                        ^W Where Is                       ^V Next Page                      ^U UnCut Text                     ^T To Spell
</pre>
<br>
<pre style="color: silver; background: black;">-bash-4.2$ sbatch transcript_assembly.sh</pre>

<h2 id="Sixth_Point_Header">Differential expression analysis using ballgown</h2>
For many organisms, many of the same genes are expressed in separate cell types, with a variety of phenotype differences a result of the specific isoforms a cell will use. Therefore, when performing a differential expression analysis from different parts of one organism (not one species, but a singular organism), it is wise to perform an isoform expression analysis alongside a standard differential expression analysis and combine the results. We will only be performing the isoform expresion analysis. <a href="https://bioconductor.org/packages/release/bioc/html/ballgown.html">Ballgown</a> is a differential expression package for R via Bioconductor ideal for isoform expression analyses. Before beginning, you need to secure copy our ballgown directory from Xanadu to your local machine:

<pre style="color: silver; background: black;">-bash-4.2$ exit
logout
Connection to xanadu-submit-ext.cam.uchc.edu closed.
user:~$ scp -r YOUR.USER.NAME@xanadu-submit-ext.cam.uchc.edu:/home/CAM/aeily/rnaseq_for_model_plant/ballgown .</pre>

Now we load <a href="https://www.rstudio.com/products/rstudio/download/">RStudio</a> with administrator privileges (otherwise you cannot install packages!).

To begin we must download and load the proper packages:

<pre style="color: silver; background: black;">install.packages("devtools")
install.packages("RFLPtools")
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown", "genefilter", "dplyr", "devtools"))

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(ggplot2)
library(gplots)
library(devtools)
library(RFLPtools)</pre>

Now we need to set our working directory to the directory which contains our "ballgown" folder. For me, this is:

<pre style="color: silver; background: black;">dir <- "C://Users/Wolf/Desktop/"
setwd(dir)
list.files()</pre>

You should see the "ballgown" folder after the list.files() command.

Let's have a look at the ballgown function:

<pre style="color: silver; background: black;">help("ballgown")

<strong>constructor function for ballgown objects</strong>

<strong>Description</strong>

<em>constructor function for ballgown objects</em>

<strong>Usage</strong>

ballgown(samples = NULL, dataDir = NULL, samplePattern = NULL,
  bamfiles = NULL, pData = NULL, verbose = TRUE, meas = "all")
Arguments

samples			vector of file paths to folders containing sample-specific ballgown data (generated by tablemaker). If samples 
			is provided, dataDir and samplePattern are not used.
dataDir			file path to top-level directory containing sample-specific folders with ballgown data in them. Only used if 
			samples is NULL.
samplePattern		regular expression identifying the subdirectories of\ dataDir containing data to be loaded into the ballgown 
			object (and only those subdirectories). Only used if samples is NULL.
bamfiles		optional vector of file paths to read alignment files for each sample. If provided, make sure to sort properly
			(e.g., in the same order as samples). Default NULL.
pData			optional data.frame with rows corresponding to samples and columns corresponding to phenotypic variables.
verbose			if TRUE, print status messages and timing information as the object is constructed.
meas			character vector containing either "all" or one or more of: "rcount", "ucount", "mrcount", "cov", "cov_sd", 
			"mcov", "mcov_sd", or "FPKM". The resulting ballgown object will only contain the specified expression 	
			measurements, for the appropriate features. See vignette for which expression measurements are available for 
			which features. "all" creates the full object.</pre>

Because of the structure of our ballgown directory, we may use dataDir = "ballgown", samplePattern = "athaliana", measure = "FPKM", and pData = some_type_of_phenotype_matrix.

We want all of the objects in our arguments to be in the same order as they are present in the ballgown directory. Therefore, we want our pData matrix to have two columns -- the first column being the samples as they appear in the ballgown directory, and the second being the phenotype of each sample in the column before it (root or shoot). Let's see the order of our sample files:

<pre style="color: silver; background: black;">list.files("ballgown/")
&#35;&#35; [1] "athaliana_root_1"  "athaliana_root_2"  "athaliana_shoot_1" "athaliana_shoot_2"</pre>

Now we construct a 4x2 phenotype matrix with the first column being our samples in order and the second each sample's phenotype:

<pre style="color: silver; background: black;">pheno_data = c("athaliana_root_1", "athaliana_root_2", "athaliana_shoot_1",  "athaliana_shoot_2","root","root","shoot","shoot")

&#35;&#35;R fills the rows of each column first. Therefore, our character vector contains all of column 1 in order, followed by all of
&#35;&#35;column 2 in order. The end result will be our matrix.

pheno_matrix = matrix(pheno_data, ncol=2)

pheno_matrix = as.data.frame(pheno_matrix)

colnames(pheno_matrix) <- c("sample", "part")

rownames(pheno_matrix) <- pheno_matrix[,1]

pheno_matrix

&#35;&#35;                             sample  part
&#35;&#35;athaliana_root_1   athaliana_root_1  root
&#35;&#35;athaliana_root_2   athaliana_root_2  root
&#35;&#35;athaliana_shoot_1 athaliana_shoot_1 shoot
&#35;&#35;athaliana_shoot_2 athaliana_shoot_2 shoot</pre>

We may now create our ballgown object:

<pre style="color: silver; background: black;">
bg <- ballgown(dataDir = "ballgown", pData=pheno_matrix, samplePattern = "athaliana")
&#35;&#35;Wed Jun 13 11:12:33 2018
&#35;&#35;Wed Jun 13 11:12:33 2018: Reading linking tables
&#35;&#35;Wed Jun 13 11:12:33 2018: Reading intron data files
&#35;&#35;Wed Jun 13 11:12:35 2018: Merging intron data
&#35;&#35;Wed Jun 13 11:12:36 2018: Reading exon data files
&#35;&#35;Wed Jun 13 11:12:39 2018: Merging exon data
&#35;&#35;Wed Jun 13 11:12:40 2018: Reading transcript data files
&#35;&#35;Wed Jun 13 11:12:41 2018: Merging transcript data
&#35;&#35;Wrapping up the results
&#35;&#35;Wed Jun 13 11:12:42 2018

&#35;&#35;we filter our ballgown object to take only genes with <a href="https://en.wikipedia.org/wiki/Variance">variances</a> above 1 
&#35;&#35;using <a href="https://www.rdocumentation.org/packages/metaMA/versions/3.1.2/topics/rowVars">rowVars()</a>

??ballgown::subset</pre>

<pre>
<strong style="color: blue;">subset ballgown objects to specific samples or genomic locations</strong>

<strong style="color: grey;">"Description</strong>

<em style="color: green;">subset ballgown objects to specific samples or genomic locations</em>

<strong style="color: grey;">Usage</strong>

subset(x, ...)

## S4 method for signature 'ballgown'
subset(x, cond, genomesubset = TRUE)

<strong style="color: grey;">Arguments</strong>

x	
a 		ballgown object
...		further arguments to generic subset
cond		Condition on which to subset. See details.
genomesubset	if TRUE, subset x to a specific part of the genome. Otherwise, subset x to only include specific samples. TRUE by 
		default.

<strong style="color: grey;">Details</strong>

To use subset, you must provide the cond argument as a string representing a logical expression specifying your desired subset. The subset expression can either involve column names of texpr(x, "all") (if genomesubset is TRUE) or of pData(x) (if genomesubset is FALSE). For example, if you wanted a ballgown object for only chromosome 22, you might call subset(x, "chr == 'chr22'"). (Be sure to handle quotes within character strings appropriately).</pre>

<br>

<pre style="color: silver; background: black;">bg_filt = subset(bg, "rowVars(texpr(bg))>1", genomesubset=TRUE)

&#35;&#35;we follow the guide and subset our ballgown object under the condition that the row-variances of the expression data are
&#35;&#35;greater than one, keeping the gene names.</pre>

To perform the isoform differential expression analysis we use ballgown's "stattest" function. Let's have a look at it:
<pre style="color: silver; background: black;">??ballgown::stattest</pre>

<pre><strong style="color: blue;">statistical tests for differential expression in ballgown</strong>

<strong style="color: grey;">Description</strong>

<em style="color: green;">Test each transcript, gene, exon, or intron in a ballgown object for differential expression, using
&#35;&#35;comparisons of linear models.</em>

<strong style="color: grey;">Usage</strong>

stattest(gown = NULL, gowntable = NULL, pData = NULL, mod = NULL,
  mod0 = NULL, feature = c("gene", "exon", "intron", "transcript"),
  meas = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"),
  timecourse = FALSE, covariate = NULL, adjustvars = NULL, gexpr = NULL,
  df = 4, getFC = FALSE, libadjust = NULL, log = TRUE)

<strong style="color: grey;">Arguments</strong>

gown		name of an object of class ballgown
gowntable	matrix or matrix-like object with rownames representing feature IDs and columns representing samples, with expression 
		estimates in the cells. Provide the feature name with feature. You must provide exactly one of gown or gowntable. NB:
		gowntable is log-transformed within stattest if log is TRUE, so provide un-logged expression values in gowntable.
pData		Required if gowntable is provided: data frame giving phenotype data for the samples in the columns of gowntable. (Rows
		of pData correspond to columns of gowntable). If gown is used instead, it must have a non-null, valid pData slot (and 
		the pData argument to stattest should be left NULL).
mod		object of class model.matrix representing the design matrix for the linear regression model including covariates of
		interest
mod0		object of class model.matrix representing the design matrix for the linear regression model without the covariates of 
		interest.
feature		the type of genomic feature to be tested for differential expression. If gown is used, must be one of "gene", 
		"transcript", "exon", or "intron". If gowntable is used, this is just used for labeling and can be whatever the rows of
		gowntable represent.
meas		the expression measurement to use for statistical tests. Must be one of "cov", "FPKM", "rcount", "ucount", "mrcount",
		or "mcov". Not all expression measurements are available for all features. Leave as default if gowntable is provided.
timecourse	if TRUE, tests whether or not the expression profiles of genomic features vary over time (or another continuous
		covariate) in the study. Default FALSE. Natural splines are used to fit time profiles, so you must have more timepoints
		than degrees of freedom used to fit the splines. The default df is 4.
covariate	string representing the name of the covariate of interest for the differential expression tests. Must correspond to the
		name of a column of pData(gown). If timecourse=TRUE, this should be the study's time variable.
adjustvars	optional vector of strings representing the names of potential confounders. Must correspond to names of columns of
		pData(gown).
gexpr		optional data frame that is the result of calling gexpr(gown)). (You can speed this function up by pre-creating
		gexpr(gown).)
df		degrees of freedom used for modeling expression over time with natural cubic splines. Default 4. Only used if 
		timecourse=TRUE.
getFC		if TRUE, also return estimated fold changes (adjusted for library size and confounders) between populations. Only
		available for 2-group comparisons at the moment. Default FALSE.
libadjust	library-size adjustment to use in linear models. By default, the adjustment is defined as the sum of the sample's log
		expression measurements below the 75th percentile of those measurements. To use a different library-size adjustment,
		provide a numeric vector of each sample's adjustment value. Entries of this vector correspond to samples in in rows of
		pData. If no library size adjustment is desired, set to FALSE.
log		if TRUE, outcome variable in linear models is log(expression+1), otherwise it's expression. Default TRUE.
</pre>

We see we can determine which transcripts and genes are differentially expressed in the roots or shoots, alongside the fold changes of
&#35;&#35;each differentially expressed gene as measured in FPKM with the following code:

<pre style="color: silver; background: black;">results_transcripts = stattest(bg_filt, feature="transcript" , covariate = "part" , 
getFC = TRUE, meas = "FPKM")

results_genes = stattest(bg_filt, feature="gene" , covariate = "part" , getFC = TRUE, meas = "FPKM")</pre>

Now we want to order our results according to their pvalue, and then subset to only take results with pvalues below 0.01, writing our findings to a csv:

<pre style="color: silver; background: black;">
results_genes = arrange(results_genes,pval)
results_genes = subset(results_genes, pval < 0.01)
results_transcripts = arrange(results_transcripts, pval)
results_transcripts = subset(results_transcripts, pval < 0.01)
write.csv(results_transcripts, "transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "gene_results.csv", row.names=FALSE)
&#35;&#35;we use row.names=FALSE because currently the row names are just the numbers 1, 2, 3. . .
</pre>

Now we want to visualize our data:

<pre style="color: silver; background: black;">&#35;&#35;we want pretty colors for our visualization
tropical = c('red', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
&#35;&#35;we want to compare our genes based on their FPKM values. We know from reading ballgown's vignette that we can extract the 
&#35;&#35;expression data using texpr() and specifying a measure. 

fpkm = texpr(bg, meas = "FPKM")
&#35;&#35;let's look at the distribution
plot(density(fpkm),main="Density Plot of \nUntransformed FPKM")
&#35;&#35;We can see virtually nothing except that there are many, many genes that are lowly expressed. The reason for the sharp peak 
&#35;&#35;is that the density plot automatically scales its x-axis from the lowest expressed to the highest expressed. Let's see what 
&#35;&#35;those values are:

min(fpkm)
&#35;&#35;0
max(fpkm)
&#35;&#35;1179157

&#35;&#35;because of this crazy scaling, we cannot truly see the distribution. However, what we <i>can</i> do is to transform the data 
&#35;&#35;such that the variance is not so staggering, allowing us to see better. There are a few rules for this, all of the data must 
&#35;&#35;be transformed in a consistent and reversible manner, after transformation no data may have a negative value, and all data 
&#35;&#35;with a value of 0 must also be 0 after transformation. The reason for the second and third rules is more epistemological. For 
&#35;&#35;us, if a gene has an FPKM of 0, then for that sample the gene is unexpressed. Should we transform the data and that 
&#35;&#35;particular gene's FPKM is now above 0, we are fundamentally changing the nature of that sample -- i.e., we are now saying it 
&#35;&#35;is expresesing a gene it actually is not! Additionally, there is no such thing as negative expression, so there is no 
&#35;&#35;physical reality where we will have an FPKM beneath 0. With these three rules, we see that taking the log of all our data 
&#35;&#35;will prevent negative values, be consistent and reversible, and scale down our variance. However, log(0) = -inf! We have 
&#35;&#35;broken a cardinal rule (oddly enough, the fact that it is infinity is not a rule-breaker, but rather that it is 
&#35;&#35;<b>negative</b> infinity! Seeing this, we can simply add 1 to our data before log transforming, log(0+1) = 0. Now we have 
&#35;&#35;fulfilled all three rules.

fpkm = log2(fpkm + 1)
plot(density(fpkm),main="Density Plot of \nTransformed Data")
&#35;&#35;now we see an actual distribution. Let's see the difference in distribution between each individual part. To do this we are 
&#35;&#35;going to plot the density for each part, one by one, and watch for great changes.

dim(fpkm)
&#35;&#35;41854     4
colnames(fpkm)
&#35;&#35;"FPKM.athaliana_root_1"  "FPKM.athaliana_root_2"  "FPKM.athaliana_shoot_1" "FPKM.athaliana_shoot_2"

plot(density(fpkm[,1]),main="Density Comparison")
lines(density(fpkm[,2])
lines(density(fpkm[,3])
lines(density(fpkm[,4])

&#35;&#35;we see that overall the distributions are pretty similar, but not identical -- shoot 2 has the highest expression levels while root 2 has the lowest.

&#35;&#35;now we will generate a <a href="http://setosa.io/ev/principal-component-analysis/">PCA</a> plot. I strongly advise you read 
&#35;&#35;the PCA link before continuing if you are not familiar with Principal Component Analysis. It will not be explained in this tutorial.

&#35;&#35;let's create a vector with our PCA point names

short_names = c("r1","r2","s1","s2")
&#35;&#35;we are going to be using the <a href="https://en.wikipedia.org/wiki/Pearson_correlation_coefficient">Pearson coefficient</a> 
&#35;&#35;for our PCA plot. You may think of the Pearson coefficient simply as a measure of similarity. If two datasets are very 
&#35;&#35;similar, they will have a Pearson coefficient approaching 1 (every data compared to itself has a Pearson coefficient of 1). 
&#35;&#35;If two datasets are very dissimilar, they will have a Pearson coefficient approaching 0 Let's calculate a vector containing 
&#35;&#35;the correlation coefficient
	
r = cor(fpkm, use="pairwise.complete.obs", method="pearson")
&#35;&#35;the cor() function takes our expression matrix, as well as two other arguments. the "pairwise.complete.obs" argument tells R 
&#35;&#35;that we want to do a pairwise analysis (so compare column 1 to 2, 3, and 4, then column 2 to 3 and 4, and lastly column 3 to 
&#35;&#35;4) only on observations (genes) without missing values. Lastly, our method "pearson" tells us we're using the Pearson method 
&#35;&#35;as a similarity metric. Let's have a look at r

r

<strong style="color:blue;">				FPKM.athaliana_root_1 FPKM.athaliana_root_2 FPKM.athaliana_shoot_1 FPKM.athaliana_shoot_2
FPKM.athaliana_root_1              1.0000000             0.9127575              0.7254296		0.7639201
FPKM.athaliana_root_2              0.9127575             1.0000000              0.7247305		0.7661130
FPKM.athaliana_shoot_1             0.7254296             0.7247305              1.0000000		0.8877891
FPKM.athaliana_shoot_2             0.7639201             0.7661130              0.8877891		1.0000000</strong>

&#35;&#35;here we see each member of the diagonal is 1.000000. Of course we knew this already, as each member is 100% similar to 
&#35;&#35;itself! Then we have the similarity measures of each sample to each other sample.

&#35;&#35;rather than calculate the similarity, it would be nicer to calculate the dissimilarity or distance between each sample. We 
&#35;&#35;know that if two samples are the same, their similarity measure is 1.000000. We also know that then their dissimilarity is 
&#35;&#35;0%, or 0.000000. Here we see that if we subtract each element from 1, we get the dissimilarity matrix! Let's do it:

d = 1 - r
d
<strong style="color:blue;">
                       FPKM.athaliana_root_1 FPKM.athaliana_root_2 FPKM.athaliana_shoot_1 FPKM.athaliana_shoot_2
FPKM.athaliana_root_1             0.00000000            0.08724246              0.2745704	0.2360799
FPKM.athaliana_root_2             0.08724246            0.00000000              0.2752695	0.2338870
FPKM.athaliana_shoot_1            0.27457045            0.27526946              0.0000000	0.1122109
FPKM.athaliana_shoot_2            0.23607994            0.23388696              0.1122109	0.0000000
</strong>
&#35;&#35;R has a function which will perform the principal component analysis for us when provided with a dissimilarity matrix, 
&#35;&#35;cmdscale. Let's have a look at it:

help(cmdscale)</pre>

<pre style="color: silver; background: black;"><strong style="color: blue;">Classical (Metric) Multidimensional Scaling</strong>

<strong style="color: grey;">e="Description</strong>

<em style="color: green;">Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis (Gower, 1966).</em>

<strong style="color: grey;">Usage</strong>

cmdscale(d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE,
         list. = eig || add || x.ret)
<strong style=color: green;">Arguments</strong>

d	a distance structure such as that returned by dist or a full symmetric matrix containing the dissimilarities.
k	the maximum dimension of the space which the data are to be represented in; must be in {1, 2, …, n-1}.
eig	indicates whether eigenvalues should be returned.
add	logical indicating if an additive constant c* should be computed, and added to the non-diagonal dissimilarities such that the 
	modified dissimilarities are Euclidean.
x.ret	indicates whether the doubly centred symmetric distance matrix should be returned.
list.	logical indicating if a list should be returned or just the n * k matrix, see ‘Value:’.

<strong style="color:grey;">Details</strong>

Multidimensional scaling takes a set of dissimilarities and returns a set of points such that the distances between the points are approximately equal to the dissimilarities. (It is a major part of what ecologists call ‘ordination’.)

A set of Euclidean distances on n points can be represented exactly in at most n - 1 dimensions. cmdscale follows the analysis of Mardia (1978), and returns the best-fitting k-dimensional representation, where k may be less than the argument k.

The representation is only determined up to location (cmdscale takes the column means of the configuration to be at the origin), rotations and reflections. The configuration returned is given in principal-component axes, so the reflection chosen may differ between R platforms (see prcomp).

When add = TRUE, a minimal additive constant c* is computed such that the dissimilarities d[i,j] + c* are Euclidean and hence can be represented in n - 1 dimensions. Whereas S (Becker et al, 1988) computes this constant using an approximation suggested by Torgerson, R uses the analytical solution of Cailliez (1983), see also Cox and Cox (2001). Note that because of numerical errors the computed eigenvalues need not all be non-negative, and even theoretically the representation could be in fewer than n - 1 dimensions.</pre>

Let's perform our principal component analysis:

<pre style="color: silver; background: black;">pca = cmdscale(d, k=2)
&#35;&#35;so we expect pca to have four rows, each row corresponding to a sample, and two columns, the first column representing our 
&#35;&#35;first coordinate axis and the second dimension representing our second coordinate axis. If we plot the first column against 
&#35;&#35;the second column, the distances between points is the dissimilarity between points.

pca

<strong style="color:blue;">			[,1]         [,2]
FPKM.athaliana_root_1  -0.1229438 -0.016069664
FPKM.athaliana_root_2  -0.1225334  0.006751752
FPKM.athaliana_shoot_1  0.1454597 -0.045779068
FPKM.athaliana_shoot_2  0.1000175  0.055096980</strong>

&#35;&#35;for this next step it is assumed that you are familiar with plotting in R. If not you may look <a 
&#35;&#35;href="https://bioinformatics.uconn.edu/introduction-to-r/">here</a>

plot.new()
par(mfrow=c(1,1))
plot(pca$points, type='n', xlab="", ylab="", main="PCA plot for all libraries")
points(pca[,1], pca[,2], col="grey", cex=2, pch=16)
text(pca[,1], pca[,1], short_names, col=c("red", "red", "hotpink", "hotpink")</pre>


<h2 id="Seventh_Point_Header">Topological networking using cytoscape</h2>

<a href="https://github.com/miriamposner/cytoscape_tutorials">Cytoscape</a> is a desktop program which creates visual topological networks of data. The only requirement to create a topological network in Cytospace is to have an "edge list". An edge list is a two-column list which comprises two separate types of data. Each row of the edge list is a single piece of information (so both columns combine to form a piece of infromation for the row) that Cytospace uses. Now, let's go over the two columns of an edge list:

The two columns of an edge list represent the "sources" and "nodes". You may instead think of the two columns representing an object and its group. For instance, the edge list:

<pre style="color: silver; background: black;">A	1
	B	1
	C	4
	D	3
	A	4
	C	3</pre>
Means that object A is in both groups 1 and 4, object B is in only group 1, object C is in groups 3 and 4, and so forth. We could then say that our sources are the objects, and our nodes are the groups to which each object belongs. To make our topological network we simply put one point for each object and one point for each group on a piece of paper and draw arrows going from each object-point to each group-point. If two objects are pointing to the same group, we place those objects closer to the group. If an object is part of two or more groups, that object will have two or more arrows pointing to each distinct group to which it belong. The cardinal rule is that no object may have two arrows pointing to the same group (that is, the set for that group is pairwise disjoint). Now we need to begin thinking about our differential expression data.
