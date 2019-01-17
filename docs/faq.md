## Frequently Asked Questions

### Does `sumrep` perform annotations, clonal family clustering, or phylogenetic inference on my data?

No -- at its core, `sumrep` does not perform these sorts of inferences on your dataset.
It does, however, include helper functions to run several external annotation tools through R, including [`partis`](https://github.com/psathyrella/partis), [`igblast`](https://www.ncbi.nlm.nih.gov/igblast/), and [`IGoR`](https://github.com/qmarcou/IGoR).
See the [Installation](docs/installation.md) page for more information.

### Can `sumrep` process paired-end sequencing reads?

Unfortunately, there is no current plan to support this sort of data directly.
Please see the [data specification](docs/extended_documentation.md#data-specification] for more details.

### Which operating systems is `sumrep` compatible with?

`sumrep` should be compatible with contemporary versions of OS X and Linux.
`sumrep` has not been tested on Windows.
