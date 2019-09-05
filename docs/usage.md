
## Usage

### Data structures
Most functions to retrieve and compare distributions between repertoires expect `data.table` objects as input.
For example, the `data` folder contains an annotations dataset (which was obtained via `getPartisAnnotations("data/test_data.fa") %$% annotations %>% fwrite("data/test_annotations.csv")`.
We can read this in as a `data.table` as follows:
```
dat_a <- data.table::fread("data/test_annotations.csv")
```

While `sumrep` is able to handle rather general annotations datasets, things work best when the annotations dataset is a `data.table` object with specific defaults chosen by members of the [AIRR software working group](http://airr.irmacs.sfu.ca/node/35).
In particular, by default, `sumrep` adheres to the [AIRR standard](http://docs.airr-community.org/en/latest/datarep/rearrangements.html#fields), with specific default column names and definitions laid out in the [extended documention](extended_documentation.md).
Columns can be manually set via the `column` argument to any function which acts on a single column; functions acting on more than one column have similar arguments for each column involved.
An example of this is shown [below](#retrieving-distributions).

To encourage principled analyses, we recommend the following specification of `sequence_alignment` and `germline_alignment` (which have somewhat open-ended definitions in the general AIRR schema):

1. The recommended specification of `sequence_alignment` will be the aligned portion of the query sequence with IMGT gaps, constrained to the V(D)J region. Similarly, the recommended specification of `germline_alignment` will be the assembled, aligned, fully length inferred germline sequence spanning the same region as the sequence_alignment field, with IMGT gaps.
2. If IMGT gaps are not available, we will still recommend that the sequences in `sequence_alignment` be pairwise aligned and constrained to the V(D)J region, and that `germline_alignment` adhere to the same properties as above, just without the IMGT gaps.
3. The user may specify custom columns to any function which uses `sequence_alignment` and/or `germline_alignment` as defaults. In this case, we note that the results will be sensible insofar as these sequences are well-defined and compared over a span that is covered by the reads in both sets, so please take caution.

Moreover, any function which acts on a column of sequences in an annotations object by default performs the following filters:

1. Any gaps/spacers (defined here as `-` or `.` characteres) are removed for distance-based functions.
2. Any sequences with a stop codon (with respect to the start of the germline V sequence ), as determined by the `stop_codon` column, are excluded.
3. Any sequences that are out-of-frame with respect to the start of the germline V sequence, as determined by the `vj_in_frame` column, are excluded.

<!--- TODO: discuss dropping sequences shorter than some prespecified length -->

If the `stop_codon` and/or `vj_in_frame` columns are missing, the function will output a warning and return the values without applying the respective filter.
The user can adjust one or more of these flags for any summary that operates on sequences from a column in a `data.table`.

These defaults attempt to encourage sensible analyses that are not biased from sequencing artifacts or otherwise noisy data.
Nonetheless, we have tried to keep `sumrep` as flexible as possible to accomodate a wide range of use cases.

There are many helper functions which take other types of data structures as input.
These are available to the user but are not as polished or standardized.
For example, the `getDistanceMatrix` function returns the pairwise distance matrix of a `list` or `vector` of input sequences.
This object may be of auxiliary interest to the user but is not directly useful for plotting or comparison functions within the package.

### Retrieving distributions
Functions for retrieving distributions are generally of the form `getXDistribution`.
For example, the pairwise distance distribution of `dat` can be obtained via
```
pairwise_distances <- getPairwiseDistanceDistribution(dat_a)
```
This returns a vector of pairwise distances rather than a matrix, which is more practical for plotting and comparison.
This function will by default compute this distribution on the `sequence` column.
To specify a different, column, say the `junction` column, if present, you would instead want:
```
junction_pairwise_distances <- getPairwiseDistanceDistribution(dat_a, column="junction")
```
A complete table of available summary functions can be found in the [extended documentation](extended_documentation.md).

### Comparing distributions
Functions to compare distributions of two annotations datasets, say `dat_a` and `dat_b`, are in general of the form `compareXDistributions`, and expect two `data.tables` as input.
Let's read in another dataset, this time a post-processed annotations dataset in a previously-saved RDS file:

```
dat_b <- readRDS("data/test_dat_boot.rds")
```

Then, to compare the pairwise distance distributions of `dat_a` and `dat_b`, we simply run
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
compareRepertoires(repertoire_a, repertoire_b, locus="igh")
```
This would perform every comparison sans `comparePerGeneMutationRates` and `comparePerGenePerPositionMutationRates`.
Note that by default `getPartisAnnotations` and `getIgBlastAnnotations` return lists of this sort, although only `getPartisAnnotations` includes the `mutation_rates` object.

The `mutation_rates` object should be equivalent to a data structure returned by `getMutationInfo`, which is called within `getPartisAnnotations`.
This structure includes a field for each gene name (e.g. `` `IGHD2-21\*01` ``), which then includes subfields `overall_mut_rate` and `mut_rate_by_position`.

### Plotting distributions
The function `plotUnivariateDistributions` takes in a list of annotations datasets as well as a locus, and returns a plot of as many summaries are as relevant to the locus and data. For example, using our datasets above, we can run

```
plotUnivariateDistributions(list(dat_a, dat_b), locus="igh")
```

By default, this creates two plots: a frequency polygon of each distribution, as well as an empirical CDF of each distribution.

This function is easy to modify; for example, if we wish to only plot frequency polygons of the pairwise distance distribution and aromaticity distribution, we could run

```
plotUnivariateDistributions(list(dat_a, dat_b), locus="igh", plot_types="freqpoly", plot_function_strings=c("getPairwiseDistanceDistribution", "getAromaticityDistribution"))
```

### Examples
The `Examples` folder includes a few scripts that demonstrate basic `sumrep` usage:

* `ExampleComparisonUsingPartis.R` shows how to obtain `partis` annotations and parameters from `sumrep`; how to simulate from these parameters using `partis` from `sumrep`; and how to compare these observed and simulated annotations datasets with the `compareRepertoires` function.

* `ExampleComparisonWithoutPartis.R` loads pre-computed annotations, so that `partis` need not be installed, and shows how to compare these observed and simulated annotations datasets with the `compareRepertoires` function.

* `MDS.R` downloads many post-processed IgH datasets from [Gupta et al. (2017)](https://www.jimmunol.org/content/198/6/2489), computes divergences of each pairwise CDR3 length distribution, performs an multidimensional scaling of these divergences, and plots the first 2 coordinates.
