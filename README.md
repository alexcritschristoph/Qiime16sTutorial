QIIME 16s Tutorial - DiRuggiero Lab
===================
Methods in analysis of 16S sequencing data have evolved and changed rapidly over the past few years, making most available online tutorials for QIIME out of date with respect to either sequencing technology, QIIME syntax and scripts, or best practices in statistical analysis. On top of that, the sheer *number* of scripts and methods packaged with QIIME (or other microbial analysis pipelines!) might be overwhelming to someone new to the field. This page intends to be an accessible and all-inclusive introduction to how to analyze 16S sequencing data using statistically robust and current methods. 

----------

## Links and Tools
 - [QIIME Scripts](http://qiime.org/scripts/)
 - [Phyloseq](https://joey711.github.io/phyloseq/)
 - [Stability of operational taxonomic units: an important but neglected property for analyzing microbial diversity](http://www.microbiomejournal.com/content/3/1/20). He Y, Caporaso JG, Jiang XT, Sheng HF, Huse SM, Rideout JR, Edgar RC, Kopylova E et al. Microbiome 2015, 3:20 (20 May 2015)
 - McMurdie, P. J., & Holmes, S. (2014). [Waste not, want not: why rarefying microbiome data is inadmissible](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531). PLoS Comput Biol, 10(4), e1003531.
 - [Robust methods for differential abundance analysis in marker gene surveys](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010126/). JN Paulson, 2013.
 - [DADA2: High resolution sample inference from amplicon data](http://biorxiv.org/content/early/2015/08/06/024034)
 - [PiCRUST](http://picrust.github.io/picrust/)

## Generating Your Data

For the purpose of this tutorial, we'll be using a small dataset of 10 samples of 16S sequences from microbial communities inhabiting 3 different rock/soil environments: Luna (Calcite rock), Ignimbrite rock, and Soil (SAT) environments. All sequences from these samples have been combined into a single seqs.fna FASTA file for our analysis.

### Creating a mapping file

Before analyzing a set of samples, creating a mapping file is useful for the purpose of thinking about experimental design and hypothesis testing. The mapping file for QIIME includes information about your sequencing files and their associated metadata. It should be a tab-delimited text file - you can make it in Excel.
The columns SampleID, BarcodeSequence, LinkerSequence, and Description are required for each sample. SampleIDs should refer to the sequence headers used in your FASTA files. You can add other columns of metadata as needed - Description should always be the last column. The mapping file for this example is saved as `example_map.txt` in the repository.

### Picking Open Reference OTUs

The initial step we will perform in our analysis is OTU-picking. There are several OTU picking strategies in QIIME, and for almost all single-experiment analyses, it would be best to use [open reference OTU picking](http://qiime.org/scripts/pick_open_reference_otus.html). OTUs are Operating Taxonomic Units, clusters of sequences that are at least X% identical, where X is generally 97%. OTUs are not a perfect method of describing the data, but are a very widely used one. "Open reference" picking will use a database of known 16S genes to create OTU clusters while also allowing for the formation of OTUs which have sequences sufficiently different from the references.
 
> **Note:** Why do we use open-reference picking? See [Stability of operational taxonomic units: an important but neglected property for analyzing microbial diversity](http://www.microbiomejournal.com/content/3/1/20). Also check out the [QIIME page](http://qiime.org/tutorials/otu_picking.html) on OTU picking for cases in which you'd want to use closed reference OTU picking. I highly doubt that you'll be in a situation where you have to, which is mainly restricted to comparisons between difference sequencing regions.
> 

You'll need to download the [GreenGenes](http://greengenes.secondgenome.com/) database (ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz) of 16S sequences for this step- it is not the most recently updated 16S database, but it is considered to be high quality. Also, make sure [usearch](http://www.drive5.com/usearch/download.html) **version 6.1.554** is also installed in your path, by requesting a download from the link above and moving the executable to `/usr/bin`. Another version may cause an error.

The command to run is below. Remember to use \$PWD - the QIIME notes make a big deal out of using absolute paths here. 
> **Note:** A common mistake is to forget to run `sudo chmod 777 /usr/bin/usearch` and/or `sudo chmod 777 /usr/bin/usearch61` before trying to run usearch for the first time.
> 
```
pick_open_reference_otus.py -i $PWD/seqs.fna -r $PWD/97_otus.fasta -o $PWD/otu_output/ -s 0.3 -m usearch61
```

Explaining this command: -i seqs.fna is our sequences file, -r 97_otus.fasta is the reference OTU file from Greegenes, -o otu_output/ is our output directory, and -m usearch61 is our clustering algorithm (usearch).

This step has accomplished many things we used to do in QIIME with many steps in a single step: it has (1) *picked our OTUs*, (2) *generated a representative sequence for each OTU*, (3) *assigned known taxonomy to those OTUs*, (4) *created a phylogenetic tree*, and (5) *created an OTU table*. Let's now explore this output.

### The OTU Table

The important output of the last command are two files, `./otu_output/rep_set.tre`, `./otu_output/rep_set.fna`, and `./otu_output/otu_table_mc2_w_tax.biom`. The first is the phylogenetic tree describing the relationship of all of our sequences, the second is the list of representative sequences for each OTU (one sequence describing the most likely sequence for each OTU), and the third is our OTU table.

Now that we have clusters, abundances, and taxonomic labels for our sequencing data, so we can begin to interrogate it with statistical analysis. 

What does this OTU table look like? It's currently in a binary format called .biom, and needs to be converted to text using the command:

```
biom convert -i ./otu_output/otu_table_mc2_w_tax.biom -o ./otu_table.txt --to-tsv --header-key taxonomy
```
 You can then open otu_table.txt in your text editor of choice, and see the list of OTUs, their counts in each sample, and their assigned taxonomy, if any.
Let's move the OTU table to our working directory with the command `cp ./otu_output/otu_table_mc2_w_tax.biom ./otu_table.biom` to begin working with it.

##Data Analysis

### Summarize and plot taxonomy

We can begin our analysis by looking at the relative abundances of taxa per sample in the OTU table. It is useful to visualize relative abundance at specific taxonomic ranks for visual comparison and to obtain a basic overview of the members of the community. Briefly, let's see what the composition of our communities is at the phylum level. The level specified to passed to [summarize_taxa.py](http://qiime.org/scripts/summarize_taxa.html) (phylum is 2), and then some quick (hardly publication ready) plots are generated by `plot_taxa_summary.py`.

```
summarize_taxa.py -i otu_table.biom -L 2 -o ./phylum/
```

```
plot_taxa_summary.py -i phylum/otu_table_L2.txt -l phylum -c pie,bar,area -o phylum_charts/
```

Our output looks like:

![enter image description here](http://i.imgur.com/q6hs8Zf.png)

We can immediately see the four soil samples have a drastically different composition than the rock samples. 

### Calculating alpha diversity


> **Note:** I think the only place where you might want to rarefy your data to an even depth is when calculating alpha diversity statistics. 

### Identifying differentially abundant OTUs

> **Note:** DESeq2 always takes a non-normalized, raw OTU table as input.
>

A major question you will be interested in answering is: *What sorts of microbes are significantly different between two sample groupings?* These groupings could be different treatments, or different environments, and understanding how they differ in terms of microbial composition in a statistically significant way is critical. 
We'll want to use [DESeq2](http://qiime.org/scripts/differential_abundance.html), an R package which uses univariate tests to identify statistically significant features between two different groups. DESeq2 is a Bioconductor package in R, so to install it, you'll have to open up R (simply run the command "R" if it's installed), and run:
```
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

If that returns successfully, you can quit R and now be able to run the QIIME script which depends on DESeq2. For this example, let's see if there are any OTUs which are significantly more abundant in rock samples than in soil samples:

```
differential_abundance.py -i otu_table.biom -o diff_otus.txt -m example_map.txt -a DESeq2_nbinom -c Substrate -x Rock -y Soil -d
```


### Normalizing the OTU Table

An important confounding factor that must be considered when analyzing microbial data is that library sizes across an experiment are going to be of uneven depth. This can cause spurious correlations and cases where what appear to be significant differences in beta diversity analysis are actually caused by differences in library size.
The ways to deal with this are a fast-moving area of research, however, **rarefaction or cutting out all sequences at an event depth is probably a bad idea**. There are normalization techniques we will use instead. This normalization technique is called CSS, and is from one of the papers linked to above. 
```
normalize_table.py -i otu_table.biom -a CSS -o CSS_normalized_otu_table.biom
```

### Beta-diversity and PCoA

Beta diversity is a metric of diversity that describes how different the species composition of different habitats is. Beta diversity therefore tells us important information about how different every sample is from all of the rest, and whether any grouping of samples are more similar in composition than the average. 
How do we measure the difference between 2 samples? There are [mathematical and phylogenetic metrics](http://qiime.org/scripts/beta_diversity.html) which can be used. The full list is: 
```
Known metrics are: abund_jaccard, binary_chisq, binary_chord, binary_euclidean, binary_hamming, binary_jaccard, binary_lennon, binary_ochiai, binary_otu_gain, binary_pearson, binary_sorensen_dice, bray_curtis, bray_curtis_faith, bray_curtis_magurran, canberra, chisq, chord, euclidean, gower, hellinger, kulczynski, manhattan, morisita_horn, pearson, soergel, spearman_approx, specprof, unifrac, unifrac_g, unifrac_g_full_tree, unweighted_unifrac, unweighted_unifrac_full_tree, weighted_normalized_unifrac, weighted_unifrac
```

Two commonly used metrics are bray_curtis and weighted_unifrac, the latter of which is a ensemble phylogenetic distance. 
```
beta_diversity.py -i CSS_normalized_otu_table.biom -m bray_curtis,weighted_unifrac -t ./otu_output/rep_set.tre -o beta_div
```

The output of this command is a distance matrix with a value which essentially is the metric difference between every pair of samples. We can visualize this data in a Principle Coordinate plot (PCoA): 

```
principal_coordinates.py -i ./beta_div/bray_curtis_CSS_normalized_otu_table.txt -o ./beta_div_coords.txt
```
And produce a plot of the PCoA:
```
make_2d_plots.py -i beta_div_coords.txt -m example_map.txt
```
![enter image description here](http://i.imgur.com/lHsqtMd.png)

### Comparing categories

### Correlating with metadata






