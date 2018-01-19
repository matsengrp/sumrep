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

* [`partis`](https://github.com/psathyrella/partis) is used for default annotation and clonal family partitioning. 
  If you only want to use the partis functionality, the following `R` command will suffice:
  ```
   install.packages(c("data.table", "dplyr", "magrittr", "stringr"))
  ```
  Once installed, you will need to either set an environmental variable `'PARTIS_PATH'` to the path of the partis executable,
  or supply this path into any function which calls partis.
  To set the variable, add the following line to your `~/.bash_profile`:
  ```
  export PARTIS_PATH='/path/to/partis/bin/partis'
  ```
  Note 1: this needs to be the path to the actual partis executable, and *not* the parent folder, which is also named `partis`.
  
  Note 2: In order for `PARTIS_PATH` to be seen in `.bash_profile`, you will need to either restart your terminal session or re-source the profile script, e.g. via the command `. ~/.bash_profile`.
  
* Please let us know if you run into any other installation issues not addressed above!
