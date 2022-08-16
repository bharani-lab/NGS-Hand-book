Exome/Targeted Data Analysis
============================

Step1: Data Retrieval

**Key Learning Outcomes**

After completing this practical, the trainee should be able to retrieve
the raw data (FastQ files) from SRA database

**Database - Introduction**

The Sequence Read Archive (SRA) stores raw sequence data from
"next-generation" sequencing technologies including Illumina, 454,
IonTorrent, Complete Genomics, PacBio and OxfordNanopores. In addition
to raw sequence data, SRA now stores alignment information in the form
of read placements on a reference sequence. The SRA is NIH's primary
archive of high-throughput sequencing data and is part of the
International Nucleotide Sequence Database Collaboration (INSDC) that
includes at the NCBI Sequence Read Archive (SRA), the European
Bioinformatics Institute (EBI), and the DNA Database of Japan (DDBJ).
Data submitted to any of the three organizations are shared among them.

**Resources**

**SRA**

https://www.ncbi.nlm.nih.gov/sra

**SRA_Toolkit**

https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/

**Retrieval of raw data from SRA**

To retrieve a sequence from SRA, we will demonstrate tools called
SRA-toolkit for downloading the raw reads.

Open your terminal

1. Set your working directory as Desktop

   **cd** **Desktop/NGS_Data_Analysis**

2. At any time, help can be displayed for fastq-dump with the following
command:

   **fastq-dump -h**

3. Execute the following commands for downloading the data from SRA

   **fastq-dump SRR098401 --split-files**

(SRR098401 is SRA number and the command “--split-files” separates the
read into two files R1 and R2)

Step 2: Quality check and Read Preprocessing

**Key Learning Outcomes**

After completing this practical, the trainee should be able to:

-  Assess the overall quality of NGS sequence reads

-  Visualize the quality and other associated matrices of reads

-  to filter and trim adapter and poor quality bases in the reads based
   on the quality scores

**Resources**

**FastQC**

http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

**Fastp**

https://github.com/OpenGene/fastp

**Introduction**

All the methods in this tutorial are better applicable to Illumina data
sets only. Although Illumina high throughput sequencing provides highly
accurate sequence data, several sequence artifacts, including
base-calling errors and small insertions/deletions, poor quality reads
and a primer/adapter contamination, are quite common in the high
throughput sequencing data, including substitution errors. The error
rates can vary from 0.5-2.0% with errors mainly rising in frequency at
the 3’ ends of reads.

One way to investigate sequence data quality is to visualize the quality
scores and other metrics in a compact manner to get an idea about the
quality of a read data set. Read data sets can be improved by post
processing in different ways like trimming off low-quality bases,
cleaning up any sequencing adapters, and removing PCR duplicates. We can
also look at other statistics such as, sequence length distribution,
base composition, sequence complexity, presence of ambiguous bases etc.
to assess the overall quality of the data set.

Highly redundant coverage (>15X) of the genome can be used to correct
sequencing errors in the reads before assembly and errors. Various k-mer
based error correction methods exist but are beyond the scope of this
tutorial.

**Prepare the Environment**

To investigate sequence data quality and proproccessing, we will
demonstrate tools called FastQC and fastp. FastQC will process and
present the reports in a visual manner. Based on the results, the
sequence data can be processed using the fastp. We will use one data set
in this practical, which can be found in the QC directory on your
desktop.

1. Open the terminal and go to your directory where the data are stored:

   **cd Desktop/NGS_Data_Analysis**

2. At any time, help can be displayed for FastQC with the following
command:

   **fastqc -h**

**Quality Check and Visualization**

1. Execute the following commands to check the quality of the raw reads

   **fastqc -f fastq SRR098401.fastq**

(FastQC generates the results in the form of zipped and unzipped
directory for each input file)

2. To view the FastQC report file of your reads you need to use any web
browser such as firefox etc.

**firefox SRR098401_fastq.html**

The report file will have a Basic Statistics table and various graphs
and tables for different quality statistics. E.g.:

**Table 1:** FastQC Basic Statistics table

================= =======================
**File name**     **input.fastq**
File type         Conventional base calls
Encoding          Sanger / Illumina 1.9
Total Sequence    40000
Filtered Sequence 0
Sequence length   100
%GC               48
================= =======================

**Figure 1:** Per base sequence quality plot for example .fastq

.. image:: vertopal_4d620e5c4e1d49c18117bd01b67b88b5/media/image1.png
   :width: 5.00139in
   :height: 3.59722in

A Phred quality score (or Q-score) expresses an error probability. In
particular, it serves as a convenient and compact way to communicate
very small error probabilities. The probability that base A is wrong (P
(∼ A)) is expressed by a quality score, Q(A), according to the
relationship:

Q(A) = −10log10(P (∼ A))

**Read Trimming**

Read trimming can be done in a variety of different ways. Choose a
method that best suits your data. Here we are giving examples of
quality-based trimming.

**Quality Based Trimming**

Base call quality scores can also be used to dynamically determine the
trim points for each read. A quality score threshold and minimum read
length following trimming can be used to remove low quality data.

Run the following command to trim your data using Phred Score (q)

1. **cd Desktop/NGS_Data_Analysis**

2. **fastp -h**

3. **fastp -i SRR098401_R1.fastq.gz -I SRR098401_R2.fastq.gz -o
   SRR098401_R1_trimmed.fq.gz -O SRR098401_R2_trimed.fq.gz -h
   SRR098401_fastp.html -w 16**

4. **fastqc -f fastq SRR098401_R1_trimmed.fq.gz**

5. **fastqc -f fastq SRR098401_R2_trimmed.fq.gz**

(-q 33 indicates the input quality score are phred +33 encoded and -o
Output file name)

Run FastQC on the quality trimmed file and visualise the quality score

1. **fastqc -f SRR098401_R1_trimmed.fq.gz**

2. **firefox SRR098401_R2_trimmed.fq.gz**

The output should be like:

**Table 2:** FastQC Basic Statistics table

================= =========================
**File name**     **output_trimmed_fastqc**
File type         Conventional base calls
Encoding          Sanger / Illumina 1.9
Total Sequence    38976
Filtered Sequence 0
Sequence length   50-100
%GC               48
================= =========================

**Figure 2:** Per base sequence quality plot for the quality-trimmed
reads

.. image:: vertopal_4d620e5c4e1d49c18117bd01b67b88b5/media/image2.png
   :width: 5.57083in
   :height: 3.88819in

Step3: READ Alignment and preprocessing

**Data**

In this tutorial, we are going to use the in-house data set as below.

**Sample 1**

-  L001_R1.fastq (Forward Pair)

-  L001_R2.fastq (Reverse Pair)

**Sample 2**

-  S14_R1.fastq (Forward Pair)

-  S14_R2.fastq (Reverse Pair)

The samples were sequenced using Miseq platform. The samples were
retrieved from patients with Retinoblastoma. The user can try the step1
for both samples before Step3.

**Key Learning Outcomes**

After completing this practical, the trainee should be able to:

-  Perform the read alignment task on any one sample data set against
   human reference genome

-  Interpret and manipulate the mapping output using SAMtools

**Resources**

**BWA**

http://bio-bwa.sourceforge.net/

**Samtools**

http://samtools.sourceforge.net/

**Introduction**

The goal of this hands-on session is to perform an unspliced alignment
for a small subset of raw reads. We will align raw sequencing data
(after preprocessing) to the human genome using BWA and then we will
manipulate the SAM output in order to visualize the alignment on the IGV
browser.

**Prepare the Environment**

Open the Terminal.

1.First, go to the right folder, where the data are stored

   **cd** **Desktop/NGS_Data_Analysis**

The trimmed raw reads can be used for alignments and further steps

**Alignment**

You already know that there are a number of competing tools for short
read alignment, each with its own set of strengths, weaknesses, and
caveats. Here we will try BWA, a widely used ultrafast, memory efficient
short read aligner.

1.BWA has a number of parameters in order to perform the alignment. To
view them all type

   **BWA –help**

(BWA uses indexed genome for the alignment in order to keep its memory
footprint small. For this we need the whole human genome in FASTA
format. This can be retrieved from ncbi
`ftp <ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta>`__
site)

ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta

2. The indexed fasta file is generated using the command

   **bwa index** **Homo_sapiens_assembly19.fasta**

Now the genome is indexed, we can move on to the actual alignment. The
first argument for BWA is the base name of the index for the genome to
be searched; in our case this is Homo_sapiens_assembly. We also want to
make sure that the output is in SAM format using .sam at the end of the
output file name parameter. The last argument is the name of the FASTQ
file (SRR098401_trimmed).

   **bwa mem Homo_sapiens_assembly19.fasta L001_R1.fastq L001_R2.fq >**
   **bwa_alignment.sam**

(The above command results in the alignment in SAM format and stores
them in the file bwa_alignment.sam)

**Manipulate SAM output**

SAM files are rather big and when dealing with a high volume of NGS
data, storage space can become an issue. For that, we can convert SAM to
BAM files (their binary equivalent that are not human readable) that
occupy much less space.

Convert SAM to BAM using samtools view and store the output in the file
bwa_alignment.bam. You have to instruct samtools view that the input is
in SAM format (-S), the output should be in BAM format (-b) and the
output to be stored in the file specified by the -o option:

   **samtools view -bSo** **bwa_alignment.bam bwa_alignment.sam**

Step4: Variant calling, annotation, prioritization and Visualization

**Key Learning Outcomes**

After completing this practical, the trainee should be able to:

-  Perform the NGS data variant calling task using aligned data

-  Filter and Prioritize the variants that are associated with the
   disease

-  Visualize the alignment and variants via a standard genome browser,
   e.g. IGV browser

**Resources**

**Gatk**

https://github.com/broadinstitute/gatk

**Wannovar**

http://wannovar.wglab.org/

**IGV**

http://software.broadinstitute.org/software/igv/

**Introduction**

Once you have aligned file against the human reference genome, you
detected nucleotidie level changes in the raw reads by comparing the
reference genome using variant caller tools. There are several best
performing tools exist, such as DeepVariant, GATK, samtools and the
Strelka etc.

In this tutorial, we use GATK, one the best performing variant caller
and memory efficient. It produces a very well-annotated VCF output that
is suitable for immediate downstream filtering.

**Variant calling with GATK**

GATK is a fast and accurate variant caller optimized for germline and
somatic variants detection. In this tutorial, we used a germline method
to detect all variants from the retinoblastoma samples.

**gatk HaplotypeCaller --java-options "-Xmx100g" -R
/Reference/hg38.fasta -I alignment_RG.bam -O GATK.vcf.gz**

(--bam indicates the bam file ; --reference indicates the Reference
file; --runDir Indicated the run directory -j indicates the total number
of CPU cores -m indicated the cluster node)

**Variant Annotation:**

ANNOVAR is a rapid, efficient tool to annotate the functional
consequences of genetic variation from high-throughput sequencing data.
wANNOVAR provides easy and intuitive web-based access to the most
popular functionalities of the ANNOVAR software.

Here, we use wANNOVAR for the annotation of variants that are generated
from **GATK** tool. First, upload your final vcf file (variants.vcf.gz)
from your local computer in wANNOVAR website. The analysis will take a
while and the output file will be returned in the form of .csv. Please
provide the email ID. You can mention ‘retinoblastoma’ in the phenotype
column.

**Variant filtering and prioritization**

The output file (.csv) from the wANNOVAR can be best viewed using
windows excel. However, here we used linux excel that are from
open-source forum. Please enable the Data>filter option for the excel
data. Once the filter option is enabled, you can set the following
filter for the variant prioritization.

1. Variants that are present in the functional site are alone kept
   example Exonic region splice site are kept.

2. Non-synonymous, frameshift mutation, stop gain and stop loss mutation
   are alone further filtered

3. Allelic frequency should be lesser than 0.01 are filtered in the 1000
   genome and ExAC database

4. Among the 3 (polyphen, sift and mutpred) any two program should
   predict the amino acid substation should be deleterious.

5. CADD score should be higher than 10 for all the variants

6. GERP score should be greater than 2.5

After the above filtering step you may be seeing many variants (Figure
3), which may potentially disease-causing variants. However, we need
prioritize the variants that are associated with the phenotype. In our
case the phenotype is Retinoblastoma. The known gene panel for
retinoblastoma includes RB1 gene mutation and MYCN copy number changes,
both contribute the 80% frequency in the retinoblastoma samples.

**Figure 3:** Potential Disease-causing variants

.. image:: vertopal_4d620e5c4e1d49c18117bd01b67b88b5/media/image3.png
   :width: 6.49722in
   :height: 0.62917in

**Phenotype-Type Based Prioritization**: You can several online tools
such as phenolyzer, Varlect, exomizer for gene-prioritization based on
user-provided phenotype. In this tutorial, you can try out phenolyzer
online tool for variant prioritization.

Phenolyzer: http://phenolyzer.wglab.org/

1. Upload your gene set with potential disease-causing variants in
      phenolyzer site by enabling the gene selection option ‘yes’

2. Enter “Retinoblastoma” in the disease/phenotype box

3. keep rest of the option default

4. Submit with your email ID

**Final Manual Confirmation using IGV**

IGV is a stand-alone genome browser. Please check their website
(http://www.broadinstitute. org/igv/) for all the formats that IGV can
display. For our visualization purposes we will use the BAM formats.

When uploading a BAM file into the genome browser, the browser will look
for the index of the BAM file in the same folder where the BAM files is.
The index file should have the same name as the BAM file and the suffix
.bai. Finally, to create the index of a BAM file you need to make sure
that the file is sorted according to chromosomal coordinates.

Sort alignments according to chromosomal position and store the result
in the file with the prefix bwa_alignment.sorted:

**samtools sort bwa_alignment.bam bwa_alignment.sorted**

Index the sorted file.

**samtools index** **bwa_alignment.sorted.bam**

The indexing will create a file called bwa_alignment.sorted.bam.bai.
Note that you don’t have to specify the name of the index file when
running samtools index, it simply appends a .bai suffix to the input BAM
file.

Now we will load the data into the IGV browser for visualization. In
order to launch IGV double click on the IGV icon on your Desktop. Ignore
any warnings and when it opens you have to load the genome of interest.
on the top left of your screen choose from the drop-down menu human
hg19. Then in order to load the desire files go to:

**File > Load from File**

On the pop-up window navigate to **Desktop > Desktop/NGS_Data_Analysis**
folder and select the file **bwa_alignment.sorted.bam**

In order to see the aligned reads and the detected variants from your
BAM file, you need to zoom in to a specific region. For example, to look
at the gene RB1 in chromosome 13

**Figure 4:** Visualization of the variants in IGV

.. image:: vertopal_4d620e5c4e1d49c18117bd01b67b88b5/media/image4.png
   :alt: Rb1_workshop
   :width: 6.49583in
   :height: 2.61667in
