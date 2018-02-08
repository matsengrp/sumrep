# *Sum*mary Statistics for BCR and TCR *Rep*ertoires

## Installation Notes
* The following `R` command will install almost all of the `R` package dependencies:
  ```
  install.packages(c("alakazam", "ape", "data.table", "dplyr", "HDMD", "jsonlite", "magrittr", "pegas", "Peptides", "RecordLinkage", "shazam", "seqinr", "stringdist", "stringr", "textmineR", "yaml"))
  ```

* `Biostrings` needs to be installed in a special way:
    ```
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    ```
  (try `http://` if `https://` URLs are not supported)
  
* If you only want to use the partis functionality, only the following `R` packages need to be installed:
  ```
   install.packages(c("data.table", "dplyr", "magrittr", "stringr"))
  ```
  The section below discusses installing `partis` in more detail.

### Installing partis
* [`partis`](https://github.com/psathyrella/partis) is used for default annotation and clonal family partitioning;
  follow the link to its respository for specific installation directions. 
* We adivse to run the following command after the regular `partis` install:
  ```
  conda update -y -c bioconda pysam
  ```
  Otherwise, you might encounter the following (or similar) error:
  ```
  AttributeError: 'csamtools.AlignedRead' object has no attribute 'cigarstring'
  ```
  This issue arises because `conda pysam` doesn't have a recent enough verion of `pysam`.

* Once `partis` is installed, you will need to either set an environmental variable `'PARTIS_PATH'` to the path of the partis   executable,
  or supply this path to any function which calls partis (via the `partis_path` argument).
  To set the variable, add the following line to your `~/.bash_profile`:
  ```
  export PARTIS_PATH='/path/to/partis/bin/partis'
  ```
  and make sure to either restart your terminal session or re-source the profile script, e.g. via the command
  ```
  . ~/.bash_profile
  ```
  
  NOTE: this needs to be the path to the actual partis executable, and *not* the parent folder, which is also named `partis`.

## Running sumrep
* `sumrep` is still in development and is not yet a complete R package. However, you can use the `devtools` library to load it   as if it were a package.
* Follow the instructions from the [`devtools` repo](https://github.com/hadley/devtools) for installation instructions.
* Once installed, you can use the following lines will load `sumrep` into your workspace:
  ```
  library(devtools)
  devtools::load_all("/path/to/sumrep/")
  ```
  where `/path/to/sumrep/` is the directory of sumrep from `git clone`.

##
Please let us know if you run into any other installation issues not addressed above!
