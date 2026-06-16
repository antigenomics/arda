[![Build Status](https://travis-ci.org/mikessh/migec.svg?branch=master)](https://travis-ci.org/mikessh/migmap)

#  MiGMAP: mapper for full-length T- and B-cell repertoire sequencing

In a nutshell, this software is a smart wrapper for [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi) V-(D)-J mapping tool designed to facilitate analysis immune receptor libraries profiled using high-throughput sequencing. This package includes additional experimental modules for contig assembly, error correction and immunoglobulin lineage tree construction.

The software is distributed as an executable JAR file and a data bundle.

**NOTE** Last IgBlastWrp version is available [here](https://github.com/mikessh/migmap/releases/tag/v0.6) (source and readme are available [here](https://github.com/mikessh/migmap/tree/v0.6)), this is a completely re-written version of original software.


## Motivation

IgBlast is an excellent  of V-(D)-J mapping tool able to correctly map even severely hypermutated antibody variants. While being a gold standard, the following limitations of IgBlast v1.4.0 have driven MIGMAP development:

- It doesn't extract sequence of CDR3 region directly, neither provide coordinates for CDR3 region in reads. It reports reference Cys residue of Variable segment and Variable segment end in CDR3, but not Phe/Trp residue of J segment that marks the end of CDR3

- Output is not straightforward to parse and summarize to a readable clonotype abundance table containing CDR3 sequences, segment assignments and list of somatic hypermutations

- It doesn't account for sequence quality

- It is somewhat hard to make it running with a custom segment reference and species other than human and mouse


## Features

Present wrapper adds the following capabilities to IgBlast:

- Run on FASTQ data.

- Use a comprehensive V/D/J segment database for human, mouse, rat, rabbit and rhesus monkey.

- Speed-up by piping reads to IgBlast and parsing the output in parallel as the built-in ``--num-threads`` argument doesn't offer much optimization.

- Assemble clonotypes, apply various filtering options such as quality filtering for CDR3 N-regions and mutations, non-coding sequence filtering, etc.

- Reporting mutations (including indels) in V, D and J segments, grouped by CDR/FW region, both on nucleotide and amino-acid level.

- Frequency and parsimony-based error correction.

- Includes a post-analysis module for quantification of somatic hypermutations and building clonotype trees, output compatible with [VDJtools](https://github.com/mikessh/vdjtools) post-analysis software.


## Pre-requisites

[Java v1.8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) or higher is required to run MIGMAP. Users should then install [IgBlast v1.4.0](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) binaries that are appropriate for their system and make sure that ``igblastn`` and ``makeblastdb`` are added to ``$PATH`` or the directory that contains binaries is specified using ``--blast-dir /path/to/bin/`` argument during MiGMAP execution. IgBlast v1.4.0 binaries can also be downloaded from [here](https://github.com/mikessh/igblast-bin).

> Note that MIGMAP also works with IgBlast v1.6.1, although this was not tested extensively

A data folder named ``data/`` containing binary databases required for IgBlast to work is provided in the release bundle. It can also explicitly specify its path with ``--blast-dir /path/to/bin/`` for troubleshooting purposes.


## Installation

See [latest release](https://github.com/mikessh/migmap/releases/latest) section for MiGMAP package. For Windows you need to both install IgBlast and download the latest release. For MacOS and Linux, MIGMAP can be easily installed using [Homebrew/Linuxbrew](http://brew.sh) or [bioconda](https://bioconda.github.io/) (no need to download anything/manually install IgBlast):

```bash
brew tap mikessh/repseq
brew install migmap-macos # or migmap-linux
```

Another option is to intall MIGMAP using BIOCONDA, see corresponding [recipe](http://bioconda.github.io/recipes/migmap/README.html).

MiGMAP can be compiled from sources using [Gradle](https://gradle.org/) with ``gradle build``. Note that in order for tests to pass IgBlast binaries should be in ``$PATH`` variable, you may need to modify following part of ``build.gradle`` 

```gradle
test {
    environment "PATH", "$System.env.PATH:/usr/local/bin/:/usr/local/ncbi/igblast/bin/"
}
```

## Execution

### General

To see the full list of MiGMAP options run 

```bash
java -jar migmap.jar -h ...
```

The memory limit can be extended by using ``-Xmx`` argument (``-Xmx8G`` will be sufficient in most cases). In case installed using **Homebrew** the command to execute MIGMAP is simply ``migmap ...``.

The following command will process ``sample.fastq.gz`` file containing human Immunoglobulin Heavy Chain reads, assemble clonotypes and store them in ``out.txt``:

```bash
java -jar migmap.jar -R IGH -S human sample.fastq.gz out.txt
```

MIGMAP accepts both FASTQ and FASTA input files, raw and GZIP-compressed. MiGMAP can be also run in per-read mode and allows piping results, e.g.:

```bash
java -jar migmap.jar --by-read -R IGH -S human sample.fastq.gz - | grep "IGHV1-8" > out.txt
```

Several receptor chains can be specified, e.g. ``-R IGH,IGK,IGL``. It is always recommended to map to complete set of TCR or IG genes and filter contaminations (e.g. TRA<>TRB) later.

The list of possible options is the following:


Option              | Definition
--------------------|------------------------------------------------------------------------
``--blast-dir`` | Path to folder that contains ``igblastn`` and ``makeblastdb`` binaries, by default assume they are added to ``$PATH`` and execute them directly.
``--data-dir`` | Path to folder that contains the data bundle (``internal_data/`` and ``optional_file/`` directories). By default it the package is provided with MIGMAP binaries, that is ``install_dir/data/``.
``--custom-database`` | Path to a custom segments database. By default will use built-in database. See ``segments.txt`` and **Using your own references** for details.
``-n`` | Number of reads to analyze before stopping. Will analyze all reads by default
``-p`` | Number of threads to use. By default will use all available threads.
``--report`` | Path to the file that is going to store MIGMAP report (extraction and filtering efficiency for current input, etc). Will append report line if file exists.
``-R`` | **REQUIRED** Receptor gene and chain. Several chains can be specified, separated with commas. Allowed values are: ``IGH,IGL,IGK,TRA,TRB,TRG,TRD``.
``-S`` | **REQUIRED** Species, allowed values: ``human,mouse,rat,rabbit,rhesus_monkey``
``--all-alleles`` | Will use all alleles provided in the antigen receptor segment database. By default uses only major allele (``*01`` according to IMGT).
``--use-kabat`` | Will use KABAT nomenclature for FR/CDR region markup. Uses IMGT nomenclature by default.
``--allow-incomplete`` | Report clonotypes with partial CDR3 mapping (lacking J germline part, etc). By default those are no included into the output.
``--allow-no-cdr3`` | Report clonotypes with no CDR3 mapping. By default those are no included into the output.
``--allow-noncoding`` | Report clonotypes that have either stop codon or frameshift in their receptor sequence. By default those are no included into the output.
``--allow-noncanonical`` | Report clonotypes that have non-canonical CDR3 (do not start with C or end with F/W residues). By default those are no included into the output.
``-q`` | Quality threshold, ``2-40`` defaults to ``25``. Filter out reads that have at least one error or N-nucleotide with a quality value below the corresponding threshold.
``--details`` | Specifies the nucleotide and amino acid sequences for certain FR/CDR regions that will be added to the output. Allowed values: ``fr1X,cdr1X,fr2X,cdr2X,fr3X,cdr3X,fr4X,contigX`` where ``X`` stands for ``nt`` or ``aa``. By default all those fileds are included
``--by-read`` | Will output mapping results for each read (see **Output format** below, excluding frequency and count fields) and read headers.
``--unmapped`` | Specifies a file to store unmapped reads.


### Additional routines

There are several built-in routines implementing common result processing and analysis tasks. Note that one should use ``-cp`` instead of ``-jar`` when executing the module and specify full package name, as shown below. When using ``-cp`` (classpath) for execution always make sure that the path to executable JAR file is set correctly, otherwise JVM will throw some uninformative error message.

#### Merging contigs

``MergeContigs`` routine of MIGMAP merges clonotypes that are represented by embedded contigs. MIGMAP uses CDR3 nucleotide sequence, V and J segment names and mutation list to define a unique clonotype signature. Therefore in case of variable IG/TCR sequence coverage (can be a result of intrinsic library properties and/or read trimming) there is a chance of ambigous clonotypes: e.g. in case two reads coming from a clonotype with a mutation in FR1 region and one of the reads doesn't cover the FR1 two clonotypes will be reported by MIGMAP. To resolve these ambiguities run ``MergeContigs`` utility as follows:

```bash
java -cp migmap.jar com.antigenomics.migmap.MergeContigs out.txt out_merged.txt
```

This routine will generate a list of clonotypes represented by the set of longest completely overlapping contigs. Note: make sure that you haven't manually excluded ``contignt`` feature from output as in such case the routine will fail.

#### Error correction

PCR and sequencing errors, as well as hot-spot PCR errors in case of UMI correct data can generate a great deal of artificial (erroneous) clonotypes, especially in case of full-length IG sequence analysis. To filter erroneous sub-variants and append their read count to corresponding parent clonotypes, execute 

```bash
java -cp migmap.jar com.antigenomics.migmap.Correct out.txt corrected.txt
```

There are three error correction options that can be set by user:
- Minimal parent-to-child read ratio allowed per one mismatch can be tweaked using the ``-r`` option
- Number of mismatches between parent variant and its error to consider when looking for erroneous subvariants (aka search depth). By default this is set to ``1`` assuming uniform mutation process, i.e. if there is a pair of clonotypes that differ by 2 mismatches due to PCR/sequencing errors, intermediate variants should be present as well.
- ``--contig`` mode that will work with entire contigs and ignore region context (germline/CDR3 mutations). This option is recommended for full-length IG data analysis and should not be used for short reads.

The algorithm can deal with errors in CDR3 region as well as correctly resolve networks of erroneous sub-variants.

Note that error correction only works for clonotype tables, not by-read output. This step should be ideally executed after ``MergeContigs`` routine.

#### Post-analysis

> **IMPORTANT** It is highly recommended that you merge contigs and perform error correction prior to proceeding to post-analysis.

To summarize somatic hypermutations and generate clonotype trees run

```bash
java -cp migmap.jar com.antigenomics.migmap.Analyze -S species -R gene out.txt out
```

Species and receptor parameters are required and should be the same as the ones used in initial MIGMAP analysis. This module will generate several text files with ``out`` prefix:

- ``out.shm.txt`` a flat file with all mutations present in sample that can be processed with ``post/analyze_shm.Rmd`` [R markdown](http://rmarkdown.rstudio.com/) template.
- ``out.net.txt``, ``out.node.txt`` and ``out.edge.txt`` (network, node and edge properties) that can be imported to [Cytoscape](http://www.cytoscape.org/) using ``Import>Table>`` menu.

For more details and an example analysis of hypermutating Raji cell repertoire go to ``post/`` folder in this repository (or click [here](https://github.com/mikessh/migmap/tree/master/post)).

### Using your own references

It is possible to use your own references with MiGMAP, given they are in the same format as [internal reference file](https://github.com/mikessh/migmap/blob/master/src/main/resources/segments.txt). Note that you do **not** need to set the reference points, just put ``-1`` in corresponding column and run 

```bash
java -cp migmap.jar com.antigenomics.migmap.AnnotateSegments -S mouse -R TRB my_segments.txt my_segments_with_refs.txt
```

Here you are required to set the receptor(s) (``-R``) and species (``-S``), run the command with ``-h`` to see the list of available options.
Choose the closest receptor(s) and species, however don't worry as the Variable segment amino acid sequences are quite homologous between species and the markup should run fine in most cases.
Next, run MiGMAP with ``--custom-database my_segments_with_refs.txt`` selecting same receptor(s) and species as before.

To convert references in IMGT format into MIGEC/MiGMAP reference format use [imgtparser](https://github.com/antigenomics/imgtparser).

## Output format

The output is provided in a tab-delimited format. Note that no header column is present in piped output. Mutations are grouped by their FW/CDR region in several columns, mutations in the same region are separated with commas. Mutation entries are stored as follows,

```
$type$position:$reference>$query
```

where ``$type`` is ``S`` for substitution, ``D`` for deletion or ``I`` for indel. Position field ``$position`` marks either the substituted base, the first deleted base or the first base after insertion. Mutation positions are provided in Variable segment coordinates with the first Variable segment germline nucleotide having position of ``0`` (in contrast to BLAST output which is 1-based). Reference and query bases are provided for substitution, deleted and inserted bases are provided for deletions and insertions (omitting ``>``).

Amino-acid level mutations are provided as translations of codons adjacent to the mutated position. Thus, cumulative effect of mutations on the amino acid sequence is shown for mutation sets. The effect is however reported so that there is a one-to-one correspondence between nucleotide and amino acid mutation entries. For example, if ``S90:T>C`` and ``S91:C>A`` together lead to ``S30:S>H``, the ``S30:S>H,S30:S>H`` is reported, not ``S30:S>P,S30:S>Y``.

Output format for assembled clonotypes is the following:

Column              | Definition
--------------------|------------------------------------------------------------------------
freq                | clonotype frequency in (0,1]
count               | number of reads
v                   | Variable segment (top hit only)
d                   | Diversity segment (top hit only)
j                   | Joining segment (top hit only)
cdr3nt              | CDR3 nucleotide sequence
cdr3aa              | CDR3 amino acid sequence
mutations.nt.FR1    | Mutations in FR1 region, nucleotide level
mutations.nt.CDR1   | Mutations in CDR1 region, nucleotide level
mutations.nt.FR2    | Mutations in FR2 region, nucleotide level
mutations.nt.CDR2   | Mutations in CDR2 region, nucleotide level
mutations.nt.FR3    | Mutations in FR3 region, nucleotide level
mutations.nt.CDR3   | Mutations in V/D/J germline sequence in CDR3 region, nucleotide level
mutations.nt.FR4    | Mutations in FR4 region, nucleotide level
mutations.aa.FR1    | Mutations in FR1 region, amino-acid level
mutations.aa.CDR1   | Mutations in CDR1 region, amino-acid level
mutations.aa.FR2    | Mutations in FR2 region, amino-acid level
mutations.aa.CDR2   | Mutations in CDR2 region, amino-acid level
mutations.aa.FR3    | Mutations in FR3 region, amino-acid level
mutations.aa.CDR3   | Mutations in V/D/J germline sequence in CDR3 region, amino-acid level
mutations.aa.FR4    | Mutations in FR4 region, amino-acid level
cdr.insert.qual     | quality string N-nucleotides in CDR3 region
mutations.qual      | mutation quality string
v.end.in.cdr3       | V segment end (exclusive) in CDR3 region
d.start.in.cdr3     | D segment start in CDR3 region or *-1* if D segment is not defined
d.end.in.cdr3       | D segment end (exclusive) in CDR3 region or *-1* if D segment is not defined
j.start.in.cdr3     | J segment start in CDR3 region or *-1* if J segment is not defined
v.del               | Number of nucleotides deleted from V segment 3' end
d.del.5             | Number of nucleotides deleted from D segment 5' end or *-1* if D segment is not defined
d.del.3             | Number of nucleotides deleted from D segment 3' end or *-1* if D segment is not defined
j.del               | Number of nucleotides deleted from J segment 5' end or *-1* if J segment is not defined
pol.v               | Position of last nucleotide of V segment's P segment or *-1* if P segment was not found
pol.d.5             | Position of first nucleotide of D segment's 5' P segment or *-1* if P segment was not found
pol.d.3             | Position of last nucleotide of D segment's 3' P segment or *-1* if P segment was not found
pol.j               | Position of first nucleotide of J segment's P segment or *-1* if P segment was not found
has.cdr3            | *true* if CDR3 region is present (both V segment conserved residue is present)
in.frame            | *true* if receptor has no frameshifts
no.stop             | *true* if receptor contains no stop codons
complete            | *true* if CDR3 region is fully defined (both V and J conserved residues are present)
canonical           | *true* if CDR3 region starts with C residue and ends with F/W residue

Note that all coordinates are 0-based.

In case the ``--details ...`` option is specified, corresponding columns will be added to output. E.g. ``--details cdr1nt,contigaa`` will add CDR1 nucleotide sequence and translated complete receptor sequence to the table. By default, all nucleotide and amino acid sequences for FR/CDR regions and entire contigs are added.
