## Installation Notes
### Installing `sumrep`
* Although `sumrep` is still in development, we can use the `devtools` package to load `sumrep` in R as a makeshift pacakge. Follow the instructions from the [`devtools` repo](https://github.com/hadley/devtools) for installation instructions.

* The following `R` command will install almost all of the remaining `R` package dependencies:
  ```
  install.packages(c("alakazam", "ape", "CollessLike", "data.table", "dplyr", "entropy", "jsonlite", "magrittr", "Peptides", "RecordLinkage", "shazam", "seqinr", "stringdist", "stringr", "testthat", "textmineR", "yaml"))
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
* Although `sumrep` can be run with any annotated dataset, several annotation and simulation tools can be invoked within `sumrep`.
The `sumrep` methods that run these tools are built to output `data.table`s which match the specification in the [usage](#usage) section.
To use these features, follow the corresponding installation steps for each desired tool.

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

#### Installing IGoR and pygor
* [IGoR](https://github.com/qmarcou/IGoR) is an annotation and simulation tool that can be used within `sumrep`.
Please follow the installation instructions if you wish to access these methods.
In particular, you must also install pygor, a Python module which is shipped with IGoR.
Once IGoR and [pip](https://pypi.org/project/pip/) are installed, execute the command
```
pip install ./pygor
```
in IGoR's root directory.
To get this working, you may or may not need to install [conda](https://conda.io/docs/) (if you haven't already), and create an environment as so:
```
conda create -n "pygor" python=3 pandas biopython matplotlib numpy scipy
```
Then, anytime you wish to run `getIgorAnnotations` you must execute the command
```
source activate pygor
```
beforehand (i.e., before running R).
When you're finished, run
```
source deactivate
```
to exit the custom-made pygor Python environment.
This is particuarly important if you are using other `conda` environments with `sumrep`, such as for `partis`.

Finally, you will need to an environmental variable `SUMREP_PATH` specifying the parent `sumrep` folder (e.g. `export SUMREP_PATH="/path/to/sumrep/"`), as well as an environmental variable `PYGOR_PATH` specifying the `pygor` directory (e.g., `export PYGOR_PATH="/path/to/IGoR/pygor"`).
You can add the above lines to your `.bash_profile` if you wish, and may need to resource this file or reset the terminal.

