## Extended documentation

### Data specification

The following table details the expected columns in an annotations `data.table`:

| Name               | Type    | Description                                                                                                                                                       |
|--------------------|---------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------
| `sequence`           | string  | The unmodified query nucleotide sequence.                                                                                                                          |
| `vj_in_frame`        | boolean | True if the V and J segment alignments are in-frame.                                                                                                              |
| `v_call`             | string  | V gene with or without allele. For example, IGHV4-59\*01.                                                                                                          |
| `d_call`             | string  | D gene with or without allele. For example, IGHD3-10\*01.                                                                                                          |
| `j_call`             | string  | J gene with or without allele. For example, IGHJ4\*02.                                                                                                             |
| `sequence_alignment` | string  | Aligned portion of query sequence. By default constrained to variable region, but not required. Synonymous with "mature" sequence in sumrep.                      |
| `germline_alignment` | string  | Assembled, aligned, fully length inferred germline sequence spanning the same region as the sequence_alignment field. Synonymous with "naive" sequence in sumrep. |
| `junction`           | string  | Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons.                                            |
| `junction_aa`        | string  | Junction region amino acid sequence.                                                                                                                              |
| `v_3p_del`           | integer | Number of nucleotides in the V 3' deletion.                                                                                                                        |
| `d_5p_del`           | integer | Number of nucleotides in the D 5' deletion.                                                                                                                        |
| `d_3p_del`           | integer | Number of nucleotides in the D 3' deletion.                                                                                                                        |
| `j_5p_del`           | integer | Number of nucleotides in the J 5' deletion.                                                                                                                        |
| `vd_insertion`       | string  | Sequence of the insertion between the V and D segments (for heavy/beta chains).                                                                                    |
| `dj_insertion`       | string  | Sequence of the insertion between the D and J segments (for heavy/beta chains).                                                                                    |
| `vj_insertion`       | string  | Sequence of the insertion between the V and J segments (for light/alpha chains).                                                                                   |
| `np1_length`         | integer | Number of nucleotides between the V and D segments or V and J segments.                                                                                    |
| `np2_length`         | integer | Number of nucleotides between the D and J segments (for heavy/beta chains).                                                                                    |
| `clone_id`           | integer | Clonal familiy cluster assignment for the query sequence.                                                                                     |

Most of these names and definitions come directly from the [AIRR standard](http://docs.airr-community.org/en/latest/datarep/rearrangements.html#fields), with some exceptions and modifications.
Not every column is strictly required for sumrep to work (e.g., TCR datasets to not need a `clone_id`), but you will only be able to use functions for which the required columns are present.
See specific function man pages for more details.

### Summary functions

The following table details the available distribution retrieval functions in `sumrep`:

| sumrep function                                   | Summary statistic                                                                       | Assumptions              | Packages used | Comments                                  |
|---------------------------------------------------|-----------------------------------------------------------------------------------------|--------------------------|---------------|-------------------------------------------|
| `getPairwiseDistanceDistribution`                   | Vector of Levenshtein distances of each sequence to each other sequence                                              | None                     | `stringdist`    |                                           |
| `getNearestNeighborDistribution`                    | Vector of nearest neighbor (NN) distances, where the NN distance of a sequence is the minimum Levenshtein distance to each other sequence                                      | None                     | `stringdist`    |                                           |
| `getGCContentDistribution`                          | Vector of sequence-wise GC contents                                                               | None                     | `ape`           |                                           |
| `getHotspotCountDistribution`                       | Vector of sequence-wise hotspot counts                                                            | None                     | `Biostrings`    | "WRC" and "WA" are default hotspot motifs |
| `getColdspotCountDistribution`                      | Vector of sequence-wise coldspot counts                                                           | None                     | `Biostrings`    | "SYC" are default coldspot motifs         |
| `getDistanceFromGermlineToSequenceDistribution`     | Vector of Levenshtein distances from `germline_alignment` to `sequence_alignment` | Annotations              | `stringdist`    |                                           |
| `getCDR3LengthDistribution`                         | Vector of CDR3 lengths, including conserved CDR3 anchors                          | Annotations              |               |                                           |
| `getGRAVYDistribution`                              | Vector of [GRAVY](https://www.ncbi.nlm.nih.gov/pubmed/7108955) indices                                                           | Annotations              | `alakazam`      |                                           |
| `getKideraFactorDistributions`                      | Vector of each of the ten [Kidera factors](https://link.springer.com/article/10.1007/BF01025492)                                         | Annotations              | `Peptides`      |                                           |
| `getAtchleyFactorDistributions`                     | Vector of each of the five [Atchley factors](https://www.pnas.org/content/102/18/6395)                                       | Annotations              | `HDMD`          |               
| `getAliphaticIndexDistribution`                     | Vector of sequence-wise [aliphatic indices](https://www.cymobase.org/cymobase/help_popup?for_window=true&topic=aliphatic_index/)                                         | Annotations              | `Peptides`          |               |
| `getPolarityDistribution`                           | Vector of sequence-wise [polarity](https://www.ncbi.nlm.nih.gov/pubmed/4843792) values                                           | Annotations              | `alakazam`      |                                           |
| `getChargeDistribution`                             | Vector of sequence-wise charge values                                             | Annotations              | `alakazam`      |                                           |
| `getBasicityDistribution`                           | Vector of sequence-wise basicity values                                           | Annotations              | `alakazam`      |                                           |
| `getAcidityDistribution`                            | Vector of sequence-wise acidity values                                            | Annotations              | `alakazam`      |                                           |
| `getAromaticityDistribution`                        | Vector of sequence-wise aromaticity values                                        | Annotations              | `alakazam`      |                                           |
| `getBulkinessDistribution`                          | Vector of sequence-wise [bulkiness](https://www.sciencedirect.com/science/article/pii/0022519368900696) values                                          | Annotations              | `alakazam`      |                                           |
| `getCDR3PairwiseDistanceDistribution`               | Vector of pairwise Levenshtein distances of CDR3 sequences                                         | Annotations              |               |                                           |
| `getPerGeneMutationRates`                           | List of mutation rates of each observed germline gene                                           | Annotations              |               |                                           |
| `getPerGenePerPositionMutationRates`                | List of mutation rate vectors over each position, over each observed germline gene                          | Annotations              |               |                                           |
| `getSubstitutionModel`                              | Inferred [substitution matrix](https://shazam.readthedocs.io/en/version-0.1.7---parallel-brown-bag/topics/createSubstitutionMatrix/) for somatically hypermutated sequences                     | Annotations              | `shazam`        |                                           |
| `getMutabilityModel`                                | Inferred [mutability matrix](https://shazam.readthedocs.io/en/version-0.1.7---parallel-brown-bag/topics/createMutabilityMatrix/) for somatically hypermutated sequences                       | Annotations              | `shazam`        |                                           |
| `getPositionalDistanceBetweenMutationsDistribution` | Vector of positional distances between mutations over all sequences               | Annotations              |               | Defined only for sequence reads with two or more mutations from the inferred germline ancestor                                          |
| `getVGene3PrimeDeletionLengthDistribution`          | Vector of V 3' intron lengths                                                     | Annotations              |               |                                           |
| `getVGene5PrimeDeletionLengthDistribution`          | Vector of V 5' intron lengths                                                     | Annotations              |               |                                           |
| `getDGene3PrimeDeletionLengthDistribution`          | Vector of D 3' intron lengths                                                     | Annotations              |               |                                           |
| `getDGene5PrimeDeletionLengthDistribution`          | Vector of D 5' intron lengths                                                     | Annotations              |               |                                           |
| `getJGene3PrimeDeletionLengthDistribution`          | Vector of J 3' intron lengths                                                     | Annotations              |               |                                           |
| `getJGene5PrimeDeletionLengthDistribution`          | Vector of J 5' intron lengths                                                     | Annotations              |               |                                           |
| `getVDInsertionLengthDistribution`                  | Vector of VD exon lengths                                                         | Annotations              |               |                                           |
| `getDJInsertionLengthDistribution`                  | Vector of DJ exon lengths                                                         | Annotations              |               |                                           |
| `getVJInsertionLengthDistribution`                  | Vector of VJ exon lengths                                                         | Annotations              |               |                                           |
| `getVDInsertionMatrix`                              | Empirical transition matrix for VD exons                                                | Annotations              |               |                                           |
| `getDJInsertionMatrix`                              | Empirical transition matrix for DJ exons                                                | Annotations              |               |                                           |
| `getVJInsertionMatrix`                              | Empirical transition matrix for VJ exons                                                | Annotations              |               |                                           |
| `getInFramePercentage`                              | Percentage of sequences whose V and J regions are in-frame                              | Annotations              |               |                                           |
| `getClusterSizeDistribution`                        | Vector of clonal family cluster sizes                                             | Clonal family clustering |               |                                           |
| `getHillNumbers`                                    | Vector of Hill numbers of the supplied diversity orders of clonal family clusters                                                  | Clonal family clustering | `alakazam`      |                                           |
| `getSelectionEstimate`                              | Vector of estimated selection strengths of clonal family clusters                 | Clonal family clustering | `shazam`        |  This method uses [`shazam::calcBaseline`](https://shazam.readthedocs.io/en/version-0.1.7---parallel-brown-bag/topics/calcBaseline/) to compute the BASELINe posterior density for estimating selection                                      |

