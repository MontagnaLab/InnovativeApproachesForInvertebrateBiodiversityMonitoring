# Bioinformatic analysis of metabarcoding data for soil invertebrates

## 2. Bioinformatic analyses

Primers used for the amplification

- **Euk575F**: `ASCYGYGGTAAYWCCAGC`
- **Euk895R**: `TCHNHGNATTTCACCNCT`

### 2.1. Preliminary settings

We start defining some variables that we will use later.
First we define the variable `RAWDIR` for the directory containing raw sequences (i.e. the fastq files) and the variable `WORKDIR` for the working directory where we want to save the results of the analyses.

```bash
RAWDIR=<path/to/rawdata/dir>
#! RAWDIR='/home/vesuvio/Desktop/StorageDisk/RawData/LUCAS_18S_RawData/Campioni_per_corso_MM'

WORKDIR=<path/to/outputs>
#! WORKDIR='/home/vesuvio/Desktop/AnalysisDisk/_prova_analisi_corso_MM'
```

We set another variable for the maximum number of processes that can be run simultaneously, usually corresponding to the maximum number of available cores.

```bash
JOBS=<number_of_cores>
#! JOBS=32
```

Finally we define a variable containing the name of the conda environment where we installed **QIIME2**.
```bash
ENV=qiime2-amplicon-2026.1
#! ENV=qiime2-amplicon-2025.10
```


### 2.2. Reads quality check
First of all we check the quality of the raw data using the softwares **FastQC** and **MultiQC**. 
If you installed these softwares you can run the commands, but since this can be a quite long step 
it may be better to look directly at the outputs [here..missing.link..](missing.link).
> [!NOTE]
> STILL HAVE TO ADD THE LINK

```bash
#! conda activate Reads_QC

# move to the directory where we want to save the outputs #
cd "$WORKDIR"
# create a directory for FastQC outputs #
mkdir -p FastQC_output

# run fastQC #
for file in "$RAWDIR"/*.fq.gz; do
    SAMPLE=$(basename "$file")
    fastqc -t $JOBS -o FastQC_output "$file"
    echo "Processed $SAMPLE"
done

# run multiqc #
cd FastQC_output
multiqc .
```

**FastQC** produces one report for each fastq file, **MultiQC** combines those reports in a single one for all samples.
Let's have a look at the reports.

⚠️
There are still a few illumina universal adapters in some reads, we will remove them later.

### 2.3 Create reference database
To taxonomically classify our sequences we need to construct a reference database.
In this case we will use [SILVA v.138.1](https://www.arb-silva.de/), a quality checked and regularly updated 
database of aligned small (16S/18S, SSU) and large subunit (23S/28S, LSU) ribosomal RNA (rRNA) sequences 
for all three domains of life (Bacteria, Archaea and Eukarya). 
This can be a very long step since the database is quite huge, so it may be better if you directly use the 
already prepared database available [here..missing.link..](missing.link).
> [!NOTE]
> STILL HAVE TO ADD THE LINK

#### 2.3.1 Obtain and clean the reference sequences

First of all we activate the conda environment containing **QIIME2**.

```bash
conda activate $ENV
```

We will use the [q2-RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) plugin that has several built-in functions
for managing and curating reference sequence databases.

**Download** the SILVA RNA sequences and the associated taxonomic labels.
```bash
qiime rescript get-silva-data \
  --p-version '138.1' \ # SILVA database version
  --p-target 'SSURef_NR99' \ # target reference sequences, in this case 99% similarity non-redundant SSU
  --p-include-species-labels \ # include species-level labels
  --p-no-rank-propagation \ # do not fill missing ranks with the taxonomy from upper-level ranks 
  --parallel \ # run in parallel
  --o-silva-sequences silva-138.1-ssu-nr99-rna-seqs.qza \ # QIIME2 output for the sequences
  --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza \ # QIIME2 output for the taxonomc labels
  --verbose
```

**Reverse-transcribe** the RNA sequences to obtain the corresponding DNA sequences.
```bash
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \ # input RNA sequences
  --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza # output DNA sequences
```

**Remove low quality sequences**, in this case those with 5 or more degenerated bases 
and/or containing homopolymers with 8 or more bases)
```bash
qiime rescript cull-seqs \
  --i-sequences silva-138.1-ssu-nr99-seqs.qza \ # input DNA sequences
  --p-num-degenerates 5 \ # number of degenerated bases
  --p-homopolymer-length 8 \ # homopolymer length
  --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \ # output clean DNA sequences
  --p-n-jobs $JOBS # number of concurrent processes
```

#### 2.3.2. Clean the reference taxonomy
Even if SILVA is a curated database the taxonomic labels can be a bit messy, we will use **R** to clean them a bit.

First we need to **export the taxonomy** in the TSV format (Tab Separated Values), so that **R** can read it.
```bash
# export taxonomy to tsv file
qiime tools export \
  --input-path silva-138.1-ssu-nr99-tax-cleaned.qza \
  --output-path exp
```

Now **open Rstudio or an R terminal** and set the working directory to the same of `$WORKDIR`.

```r
setwd("<path/to/outputs>")
```

**Load the required libraries**, install them if they are not already installed. 
```r
#install.packages("dplyr")
library(dplyr)
#install.packages("stringr")
library(stringr)
```

**Read the file** containing the SILVA taxonomy.
```r
data <- read.table("exp/taxonomy.tsv", header=T, sep="\t")
```

**Remove bad species-level ID** that do not correspond to taxon names (i.e. those starting with lowercase) 
and fill empty genus rank using species-level identifications.
```r
CleanData <- data %>%
  mutate(Taxon = case_when(
      # Case 1: First letter after "s__" is lowercase
      str_detect(Taxon, "s__[a-z]") ~ str_replace(Taxon, "s__[a-z].*", "s__"),
      # Case 2: First letter after "s__" is uppercase
      str_detect(Taxon, "s__[A-Z]") ~ str_replace(Taxon, 
                                                  "g__; s__([^_]+)",  # Capture "g__;" and genus portion after "s__"
                                                  "g__\\1; s__\\1"),  # Append genus after "g__" and keep it after "s__"
      # Default: Leave other rows unchanged
      TRUE ~ Taxon
    )
  )
```

**Remove bad genus-level ID** that do not correspond to taxon names (i.e. those starting with lowercase).
```r
CleanData2 <- CleanData %>%
  mutate(
    Taxon = case_when(
      # If the first letter after "g__" is lowercase, remove everything after "g__" and keep "s__"
      str_detect(Taxon, "g__[a-z]") ~ str_replace(Taxon, "g__[a-z].*s__", "g__; s__"),
      
      # Otherwise, leave the string as it is
      TRUE ~ Taxon
    )
  )
```

**Save** the cleaned taxonomy to a TSV file.
```r
write.table(CleanData2, "silva-138.1-ssu-nr99-tax-cleaned.tsv", sep="\t", quote = F, row.names = F)
```

Close Rstudio or the R terminal and **go back to the bash terminal** with the QIIME2 environment
to import the cleaned SILVA taxonomy as a QIIME2 artifact.
```bash
qiime tools import \
  --type FeatureData[Taxonomy] \ # type of artifact to be created
  --input-path silva-138.1-ssu-nr99-tax-cleaned.tsv \ # path to the file to be imported
  --input-format HeaderlessTSVTaxonomyFormat \ # format of the file to be imported
  --output-path silva-138.1-ssu-nr99-tax-cleaned.qza # output file
```


#### 2.3.3. Filter and dereplicate the reference database

Next we **remove** database entries with **too short SSU sequences** depending on the taxon.
In this case we use a minimum length threshold of 900 bp for Archaea, 1200 bp for Bacteria and 1400 bp for Eukaryota.
```bash
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \ # input sequences
  --i-taxonomy silva-138.1-ssu-nr99-tax-cleaned.qza \ # input taxonomy
  --p-labels Archaea Bacteria Eukaryota \ # taxonomic labels for conditional filtering
  --p-min-lens 900 1200 1400 \ # minimum length thresholds
  --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \ # output sequences
  --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza # output taxonomy
```

Next, we **dereplicate entries** with the same taxonomic ID to remove redundance.
```bash
qiime rescript dereplicate \
  --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza \ # input sequences
  --i-taxa silva-138.1-ssu-nr99-tax-cleaned.qza \ # input taxonomy
  --p-mode 'uniq' \ # retain all sequences with unique taxonomic affiliations
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \ # output sequence
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \ # output taxonomy
  --p-threads $JOBS # number of concurrent processes
```


#### 2.3.4. Extract the amplified region from the reference database

To optimize taxonomic classification and reduce database complexity we can **trim the reference sequences**
to contain only the region actually amplified by the primers used in this study.
```bash
qiime feature-classifier extract-reads \
  --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \ # input sequences
  --p-f-primer ASCYGYGGTAAYWCCAGC \ # forward primer
  --p-r-primer TCHNHGNATTTCACCNCT \ # reverse primer
  --p-identity 0.8 \ # minimum combined primer match identity threshold
  --p-min-length 150 \ # minimum amplicon length
  --p-max-length 450 \ # maximum amplicon length
  --p-n-jobs $JOBS \ # number of concurrent processes
  --p-read-orientation forward \ # orientation of primers relative to the sequences
  --o-reads silva-138.1-ssu-nr99-seqs_Euk575-895.qza # output sequences
```

**Dereplicate** again (since some of these shorter sequences may be identical now) to remove redundance.
```bash
qiime rescript dereplicate \
  --i-sequences silva-138.1-ssu-nr99-seqs_Euk575-895.qza  \ # input sequences
  --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \ #input taxonomy
  --p-mode 'uniq' \ # retain all sequences with unique taxonomic affiliations
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \ # output sequences
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \ # output taxonomy
  --p-threads $JOBS # number of concurrent processes
```


### 2.4. Import sequences in QIIME2
We are now ready to import the LUCAS sequences in QIIME2.

#### 2.4.1. Make a manifest file containing paths to raw reads
We need to create a manifest file telling QIIME2 where to find forward and reverse reads for each sample.
This is a TSV file with three columns: sample ID, path to forward reads, path to reverse reads.
First move to the directory containing the fastq files.
```bash
cd "$RAWDIR"
```

Get the **absolute paths to forward and reverse reads** and save it in separate files.
```bash
ls -1 "$PWD/"*1.fq.gz > R1.txt # paths to forward reads
ls -1 "$PWD/"*2.fq.gz > R2.txt # paths to reverse reads
```

Next we are going to extract **sample names** from the names of the fastq files.
The file are named following this scheme: "Lucas\<sampleID\>.\<sequencing-direction\>.fq.gz" 
(e.g., Lucas0001.1.fq.gz, Lucas0001.2.fq.gz, Lucas0002.1.fq.gz, Lucas0002.2.fq.gz, ...).
So the first 9 characters of the file names correspond to the sample names. 
Let's use some bash commands to extract sample names. 
```bash
# from each fastq file name extract the first 9 characters, sort them, and keep only unique values #
ls *fq.gz | cut -c 1-9 | sort | uniq > ids.txt
```

Combine sample names and absolute paths to make a manifest file.
```bash
# add column names to the manifest file #
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv

# add sample names and paths to the manifest file #
paste ids.txt R1.txt R2.txt >> manifest.tsv

# remove temporary files #
rm *.txt

# move the manifest file to the working directory #
mv manifest.tsv "$WORKDIR"/manifest.tsv
```

Let's have a look at the content of the manifest file to check that everything is ok.

#### 2.4.2. Import sequences
Now we can use the manifest file to import sequences in QIIME2. First move back to the working directory.
```bash
cd $WORKDIR
```

Import the sequences in a QZA file.
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \ # type of artifact to be created
  --input-path manifest.tsv \ # manifest file
  --output-path seqs.qza \ # output file
  --input-format PairedEndFastqManifestPhred33V2 # format of the file to be imported
```

Create a visualization (QZV) of the imported sequences.
```bash
qiime demux summarize \
  --i-data seqs.qza \ # input file
  --o-visualization seqs.qzv # output file
```
ℹ️ QIIME2 QZV files can be visualized using the command `qiime tools view <filename.qzv>` or using the [online visualizer](https://view.qiime2.org/).


If you remember the quality check reports thare are still some Illumina adapter in the sequences, it is better to remove them using [q2-cutadapt](https://github.com/qiime2/q2-cutadapt)
```bash
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences seqs.qza \ # input file
  --p-cores $JOBS \ # number of concurrent processes
  --p-adapter-f AGATCGGAAGAG \ # adapter sequence
  --p-adapter-r AGATCGGAAGAG \ # adapter sequence
  --o-trimmed-sequences seqs_trimmed.qza \ # output file
  --verbose
```

### 2.5 Denoising and sequence re-orientation

We are now ready for removing non-biological variation from our data. We use the [DADA2 algorithm](https://benjjneb.github.io/dada2/) implemented in [q2-dada2](https://github.com/qiime2/q2-dada2) that models and corrects sequencing errors to infer exact biological sequences (amplicon sequence variants, ASVs).

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs seqs_trimmed.qza \ # input file
  --p-trim-left-f 18 \ # forward primer length
  --p-trim-left-r 18 \ # reverse primer length
  --p-trunc-len-f 220 \ # position for forward reads truncation
  --p-trunc-len-r 200 \ # position for reverse reads truncation
  --p-max-ee-f 2 \ # max expected error for forward reads
  --p-max-ee-r 2 \ # max expected error for reverse reads
  --p-trunc-q 2 \ # quality threshold for reads truncation
  --p-pooling-method 'pseudo' \ # pooling method for denoising
  --p-chimera-method 'consensus' \ # method for chimera removal
  --p-n-reads-learn 1000000 \ # number of reads used for training the error model
  --p-n-threads $JOBS \ # number of concurrent processes
  --o-table table_MixedOrientation.qza \ # output file with the ASV table
  --o-representative-sequences rep-seqs_MixedOrientation.qza \ # output file with the ASV sequences
  --o-denoising-stats denoising-stats.qza \ output file with the denoising stats
  --verbose
```
❗In real life scenarios you should experiment with `--p-trunc-len-f` and `--p-trunc-len-r` parameters and compare the results (in terms of number of retained sequences per sample and sequences length) to choose the best values.

Now let's create visualizations for the stats of the DADA2 algorithm.
```bash
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```
Let's have a look at this visualization.

Since in this study barcodes and adapters were added after PCR amplification each fastq file contained both forward and reverse reads. So sequences needs to be re-orientered using the reference database as guide. We can use again the [q2-RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) plugin for doing it.
```bash
qiime rescript orient-seqs \
  --i-sequences rep-seqs_MixedOrientation.qza \ # input file with sequences in mixed orientation
  --i-reference-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \ # reference database
  --o-oriented-seqs rep-seqs.qza \ # output file with reoriented sequences 
  --o-unmatched-seqs orientation_unmatched_sequences.qza \ # output file with unmatched sequences
  --p-threads $JOBS # # number of concurrent processes
```

Exclude unmatched sequences from the ASV table.
```bash
qiime feature-table filter-features \
  --i-table table_MixedOrientation.qza \ # input table with all ASVs
  --m-metadata-file orientation_unmatched_sequences.qza \ # file with unmatched sequences
  --p-exclude-ids \ # exclude ASVs present in the metadata
  --o-filtered-table table.qza # output table without unmatched sequences
```

Create visualizations for the ASV table and ASV sequences files.
```bash
# table visualization
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

# sequences visualization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```


---

---

---

**end of formatted text**

---

---

---







---





---

## **6. Taxonomic classification**

```bash
# fit a classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --o-classifier NBclassifier.qza

# classify without a confidence threshold
qiime feature-classifier classify-sklearn \
  --i-classifier NBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence=disable \
  --p-n-jobs $JOBS \
  --o-classification taxonomy_no_confidence.qza

# generate weights
qiime clawback generate-class-weights \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --i-reference-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-samples table.qza \
  --i-taxonomy-classification taxonomy_no_confidence.qza \
  --o-class-weight class-weights.qza

# fit new classifier with weights
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --i-class-weight Euk575-895_LUCAS_class-weights.qza \
  --o-classifier weightedNBclassifier.qza

# classify using a confidence threshold of 0.9
qiime feature-classifier classify-sklearn \
  --i-classifier weightedNBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.90 \
  --p-n-jobs $JOBS \
  --o-classification taxonomy.qza

```

## **7. Final filtering and exporting the results**

```bash
# remove singletons
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 2 \
  --o-filtered-table table.qza

# remove non-metazoan ASVs
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include 'd__Eukaryota(Metazoa)' \
  --o-filtered-table table.qza

# export table (ASVs abundance per sample)
qiime tools export \
  --input-path table.qza \
  --output-path exp
biom convert -i exp/feature-table.biom -o table.tsv --to-tsv

# export taxonomy (ASVs taxonomic assignments)
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exp
mv exp/taxonomy.tsv taxonomy.tsv

```

---

































