# *Sum*mary Statistics for BCR and TCR *Rep*ertoires

## Installation Notes
### Installing `sumrep`
* Although `sumrep` is still in development, we can use the `devtools` package to load `sumrep` in R as a makeshift pacakge. Follow the instructions from the [`devtools` repo](https://github.com/hadley/devtools) for installation instructions.

* The following `R` command will install almost all of the remaining `R` package dependencies:
  ```
  install.packages(c("alakazam", "ape", "CollessLike", "data.table", "dplyr", "entropy", "HDMD", "jsonlite", "magrittr", "pegas", "Peptides", "RecordLinkage", "shazam", "seqinr", "stringdist", "stringr", "testthat", "textmineR", "yaml"))
  ```
  
* `Biostrings` needs to be installed in a special way:
    ```
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    ```
  (try `http://` if `https://` URLs are not supported)

* Once all of the above are installed, the following lines will load `sumrep` into your workspace:
  ```
  library(devtools)
  devtools::load_all("/path/to/sumrep/")
  ```
  where `/path/to/sumrep/` is the directory of sumrep from `git clone`.
  
* If you only want to use the `partis` functionality, only the following `R` packages need to be installed:
  ```
   install.packages(c("data.table", "dplyr", "magrittr", "stringr"))
  ```
  The section below discusses installing `partis` in more detail.


### Installing annotation/simulation tools (optional)
* Although `sumrep` can be run with any annotated dataset, several annotation and simulation tools can be invoked within `sumrep`. To use these features, follow the corresponding installation steps for each desired tool.

#### Installing partis
* [`partis`](https://github.com/psathyrella/partis) can be used for annotation and clonal family partitioning within `sumrep`; follow the link to its respository for specific installation directions. 
* We suggest running the following command after the regular `partis` install:
  ```
  conda update -y -c bioconda pysam
  ```
  Otherwise, you might encounter the following (or similar) error:
  ```
  AttributeError: 'csamtools.AlignedRead' object has no attribute 'cigarstring'
  ```
  This issue arises because `conda pysam` doesn't have a recent enough verion of `pysam`.

* Once `partis` is installed, you will need to either set an environmental variable `'PARTIS_PATH'` to the path of the partis   executable, or supply this path to any function which calls `partis` (via the `partis_path` argument).
  To set the variable, add the following line to your `~/.bash_profile`:
  ```
  export PARTIS_PATH='/path/to/partis/bin/partis'
  ```
  and make sure to either restart your terminal session or re-source the profile script, e.g. via the command
  ```
  . ~/.bash_profile
  ```
  
  NOTE: this needs to be the path to the actual partis executable, and *not* the parent folder, which is also named `partis`.
  
#### Installing IgBlast + Change-O
* [IgBlast](https://www.ncbi.nlm.nih.gov/igblast/) can also be used for annotation within `sumrep`; instructions for installing standalone `igblast` can be found [here](https://ncbi.github.io/igblast/cook/How-to-set-up.html).

* `sumrep` uses [`changeo`](https://changeo.readthedocs.io/en/version-0.4.4/index.html) to process the raw `igblast` output. Installation instructions can be found [here](https://changeo.readthedocs.io/en/version-0.4.4/install.html).

* Use the `getIgBlastAnnotations` function to obtain a `data.table` with annotations from a given input fasta file. You will need to include the binary directories for `igblast` and `changeo` via the arguments `igblast_dir` and `changeo_dir`, respectively.

## Usage
### Data structures
Most functions to retrieve and compare distributions between repertoires expect `data.table` objects as input.
For example, the `data` folder contains an annotations dataset (obtained from `getPartisAnnotations("data/test_data.fa") %$% annotations %>% fwrite("data/test_annotations.csv")`.
We can read this in as a `data.table` as follows:
```
dat <- data.table::fread("data/test_annotations.csv")
```
There are many other helper functions which take other types of data structures as input.
These are of course available to the user but are not as polished or standardized.
For example, the `getDistanceMatrix` function returns the pairwise distance matrix of a `list` or `vector` of input sequences.
This object may be of auxiliary interest to the user but is not directly useful for plotting or comparison functions within the package.

### Retrieving distributions
Functions for retrieving distributions are generally of the form `getXDistribution`.
For example, the pairwise distance distribution of `dat` can be obtained via
```
pairwise_distances <- getPairwiseDistanceDistribution(dat)
```
This returns a vector of pairwise distances rather than a matrix, which is more practical for plotting and comparison.
This function will by default compute this distribution on the `sequence` column.
To specify a different, column, say the `junction` column, if present, you would instead want:
```
junction_pairwise_distances <- getPairwiseDistanceDistribution(dat, column="junction")
```
Column defaults are chosen to agree with the [AIRR standard](http://docs.airr-community.org/en/latest/datarep/rearrangements.html#fields), but can be manually set via the `column` argument to many of these functions.

### Comparing distributions
Functions to compare distributions of two annotations datasets, say `dat_a` and `dat_b`, are in general of the form `compareXDistributions`, and expect two `data.tables` as input.
For example, to compare the pairwise distance distributions of `dat_a` and `dat_b`, we would have
```
divergence <- comparePairwiseDistanceDistributions(dat_a, dat_b)
```
In this case, the output is the scalar JS-divergence between the pairwise distance distributions of the two distributions.

### Running a full comparison
The function `compareRepertoires` performs the full suite of summary comparisons between two repertoires.
The main inputs of this function, `repertoire_1` and `repertoire_2` are `list`s, rather than `data.table`s.
These lists should include an `annotations` field corresponding to the annotations `data.table` objects discussed above, and optionally a `mutation_rates` field, which is described below.
For example, the following would do the trick:
```
repertoire_a <- list(annotations=dat_a)
repertoire_b <- list(annotations=dat_b)
compareRepertoires(repertoire_a, repertoire_b)
```
This would perform every comparison sans `comparePerGeneMutationRates` and `comparePerGenePerPositionMutationRates`.
Note that by default `getPartisAnnotations` and `getIgBlastAnnotations` return lists of this sort, although only `getPartisAnnotations` includes the `mutation_rates` object.

The `mutation_rates` object should be equivalent to a data structure returned by `getMutationInfo`, which is called within `getPartisAnnotations`.
This structure includes a field for each gene name (e.g. `` `IGHD2-21\*01` ``), which then includes subfields `overall_mut_rate` and `mut_rate_by_position`.

##
Please let us know if you run into any other installation issues not addressed above!
