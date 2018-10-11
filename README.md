# *Sum*mary Statistics for BCR and TCR *Rep*ertoires

## Installation Notes
### Installing `sumrep`
* The following `R` command will install almost all of the `R` package dependencies:
  ```
  install.packages(c("alakazam", "ape", "CollessLike", "data.table", "dplyr", "entropy", "HDMD", "jsonlite", "magrittr", "pegas", "Peptides", "phytools", "RecordLinkage", "shazam", "seqinr", "stringdist", "stringr", "textmineR", "yaml"))
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
  
* `sumrep` is still in development and is not yet a complete R package. However, you can use the `devtools` library to load it   as if it were a package. Follow the instructions from the [`devtools` repo](https://github.com/hadley/devtools) for installation instructions.
* Once `devtools` is installed, the following lines will load `sumrep` into your workspace:
  ```
  library(devtools)
  devtools::load_all("/path/to/sumrep/")
  ```
  where `/path/to/sumrep/` is the directory of sumrep from `git clone`.


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
  
#### Installing IgBlast
* [IgBlast](https://www.ncbi.nlm.nih.gov/igblast/) can also be used for annotation within `sumrep`; instructions for installing standalone `igblast` can be found [here](https://ncbi.github.io/igblast/cook/How-to-set-up.html).

* Use the `getIgBlastAnnotations` function to obtain a `data.table` with annotations from a given input fasta file.

##
Please let us know if you run into any other installation issues not addressed above!
