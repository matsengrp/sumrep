## Extended documentation

### Data specification

The following table details the expected columns in an annotations `data.table`.
Note that 

| Name               | Type    | Description                                                                                                                                                       |
|--------------------|---------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------
| `sequence_alignment` | string  | Aligned portion of query sequence. By default constrained to variable region, but not required. Synonymous with "mature" sequence in sumrep.                      |
| `germline_alignment` | string  | Assembled, aligned, fully length inferred germline sequence spanning the same region as the sequence_alignment field. Synonymous with "naive" sequence in sumrep. |
| `v_call`             | string  | V gene with or without allele. For example, IGHV4-59\*01.                                                                                                          |
| `d_call`             | string  | D gene with or without allele. For example, IGHD3-10\*01.                                                                                                          |
| `j_call`             | string  | J gene with or without allele. For example, IGHJ4\*02.                                                                                                             |
| `junction`           | string  | Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons.                                            |
| `junction_aa`        | string  | Junction region amino acid sequence.                                                                                                                              |
| `vj_in_frame`        | boolean | True if the V and J segment alignments are in-frame.                                                                                                              |
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

| sumrep function                                   | Summary statistic                                                                       | Default column(s)              | Packages used | Comments                                  |
|---------------------------------------------------|-----------------------------------------------------------------------------------------|--------------------------|---------------|-------------------------------------------|
| `getPairwiseDistanceDistribution`                   | Vector of Levenshtein distances of each sequence to each other sequence                                              | `sequence_alignment`                     | `stringdist`    |                                           |
| `getNearestNeighborDistribution`                    | Vector of nearest neighbor (NN) distances, where the NN distance of a sequence is the minimum Levenshtein distance to each other sequence                                      | `sequence_alignment`                     | `stringdist`    |                                           |
| `getGCContentDistribution`                          | Vector of sequence-wise GC contents                                                               | `sequence_alignment`                     | `ape`           |                                           |
| `getHotspotCountDistribution`                       | Vector of sequence-wise hotspot counts                                                            | `sequence_alignment`                     | `Biostrings`    | "WRC" and "WA" are default hotspot motifs |
| `getColdspotCountDistribution`                      | Vector of sequence-wise coldspot counts                                                           | `sequence_alignment`                     | `Biostrings`    | "SYC" are default coldspot motifs         |
| `getDistanceFromGermlineToSequenceDistribution`     | Vector of Levenshtein distances from `germline_alignment` to `sequence_alignment` | `sequence_alignment`, `germline_alignment`              | `stringdist`    |                                           |
| `getCDR3LengthDistribution`                         | Vector of CDR3 lengths, including conserved CDR3 anchors                          | `junction_aa`, `junction`, or `junction_length`              |               |                                           |
| `getGRAVYDistribution`                              | Vector of [GRAVY](https://www.ncbi.nlm.nih.gov/pubmed/7108955) indices                                                           | `junction_aa`              | `alakazam`      |                                           |
| `getKideraFactorDistributions`                      | Vector of each of the ten [Kidera factors](https://link.springer.com/article/10.1007/BF01025492)                                         | `junction_aa`              | `Peptides`      |                                           |
| `getAtchleyFactorDistributions`                     | Vector of each of the five [Atchley factors](https://www.pnas.org/content/102/18/6395)                                       | `junction_aa`              | `HDMD`          |               
| `getAliphaticIndexDistribution`                     | Vector of sequence-wise [aliphatic indices](https://www.cymobase.org/cymobase/help_popup?for_window=true&topic=aliphatic_index/)                                         | `junction_aa`              | `Peptides`          |               |
| `getPolarityDistribution`                           | Vector of sequence-wise [polarity](https://www.ncbi.nlm.nih.gov/pubmed/4843792) values                                           | `junction_aa`              | `alakazam`      |                                           |
| `getChargeDistribution`                             | Vector of sequence-wise charge values                                             | `junction_aa`              | `alakazam`      |                                           |
| `getBasicityDistribution`                           | Vector of sequence-wise basicity values                                           | `junction_aa`              | `alakazam`      |                                           |
| `getAcidityDistribution`                            | Vector of sequence-wise acidity values                                            | `junction_aa`              | `alakazam`      |                                           |
| `getAromaticityDistribution`                        | Vector of sequence-wise aromaticity values                                        | `junction_aa`              | `alakazam`      |                                           |
| `getBulkinessDistribution`                          | Vector of sequence-wise [bulkiness](https://www.sciencedirect.com/science/article/pii/0022519368900696) values                                          | `junction_aa`              | `alakazam`      |                                           |
| `getCDR3PairwiseDistanceDistribution`               | Vector of pairwise Levenshtein distances of CDR3 sequences                                         | `junction_aa`              |               |                                           |
| `getPerGeneMutationRates`                           | List of mutation rates of each observed germline gene                                           | N/A              |               |                                           |
| `getPerGenePerPositionMutationRates`                | List of mutation rate vectors over each position, over each observed germline gene                          | N/A              |               |                                           |
| `getSubstitutionModel`                              | Inferred [substitution matrix](https://shazam.readthedocs.io/en/version-0.1.7---parallel-brown-bag/topics/createSubstitutionMatrix/) for somatically hypermutated sequences                     | `sequence_alignment`, `germline_alignment`, `v_call`              | `shazam`        |                                           |
| `getMutabilityModel`                                | Inferred [mutability matrix](https://shazam.readthedocs.io/en/version-0.1.7---parallel-brown-bag/topics/createMutabilityMatrix/) for somatically hypermutated sequences                       | `sequence_alignment`, `germline_alignment`, `v_call`              | `shazam`        |                                           |
| `getPositionalDistanceBetweenMutationsDistribution` | Vector of positional distances between mutations over all sequences               | `sequence_alignment`, `germline_alignment`              |               | Defined only for sequence reads with two or more mutations from the inferred germline ancestor                                          |
| `getVGene3PrimeDeletionLengthDistribution`          | Vector of V 3' intron lengths                                                     | `v_3p_del`              |               |                                           |
| `getVGene5PrimeDeletionLengthDistribution`          | Vector of V 5' intron lengths                                                     | `v_5p_del`              |               |                                           |
| `getDGene3PrimeDeletionLengthDistribution`          | Vector of D 3' intron lengths                                                     | `d_3p_del`              |               |                                           |
| `getDGene5PrimeDeletionLengthDistribution`          | Vector of D 5' intron lengths                                                     | `d_5p_del`              |               |                                           |
| `getJGene3PrimeDeletionLengthDistribution`          | Vector of J 3' intron lengths                                                     | `j_3p_del`              |               |                                           |
| `getJGene5PrimeDeletionLengthDistribution`          | Vector of J 5' intron lengths                                                     | `j_5p_del`              |               |                                           |
| `getVDInsertionLengthDistribution`                  | Vector of VD exon lengths                                                         | `np1_length`              |               |                                           |
| `getDJInsertionLengthDistribution`                  | Vector of DJ exon lengths                                                         | `np2_length`              |               |                                           |
| `getVJInsertionLengthDistribution`                  | Vector of VJ exon lengths                                                         | `np1_length`              |               |                                           |
| `getVDInsertionMatrix`                              | Empirical transition matrix for VD exons                                                | `vd_insertion`              |               |                                           |
| `getDJInsertionMatrix`                              | Empirical transition matrix for DJ exons                                                | `dj_insertion`              |               |                                           |
| `getVJInsertionMatrix`                              | Empirical transition matrix for VJ exons                                                | `vj_insertion`              |               |                                           |
| `getInFramePercentage`                              | Percentage of sequences whose V and J regions are in-frame                              | `vj_in_frame`              |               |                                           |
| `getClusterSizeDistribution`                        | Vector of clonal family cluster sizes                                             | `clone_id` |               |                                           |
| `getHillNumbers`                                    | Vector of Hill numbers of the supplied diversity orders of clonal family clusters                                                  | `clone_id` | `alakazam`      |                                           |
| `getSelectionEstimate`                              | Vector of estimated selection strengths of clonal family clusters                 | `sequence_alignment`, `germline_alignment` | `shazam`        |  This method uses [`shazam::calcBaseline`](https://shazam.readthedocs.io/en/version-0.1.7---parallel-brown-bag/topics/calcBaseline/) to compute the BASELINe posterior density for estimating selection                                      |


In general, these summaries are grouped into hierarchical levels, whose assumptions are described in the following table.

| Level | Assumptions | Main expected column(s) |
| ----- | ----------- | ------- |
| 0     | None        | `sequence` |
| 1     | Pairwise alignment | `sequence_alignment` |
| 2     | Annotations | `germline_alignment`, `junction_aa` |
| 3     | Clonal family clustering | `clone_id` |
| 4     | Phylogenies (BCR only) | (Currently unavailable) |

The various indel statistics tabulated above depend on the locus (e.g., whether or not D gene statistics are relevant) and fall into the Level 2 category.
These levels are hierarchical in that each level depends on the assumptions of the previous level.
For example, to obtain annotations, you would have first needed to pairwise align the sequences.
Level 0 includes raw query sequences and, while supported, does not comprise the default level for any sumrep function.
Level 4 includes tree statistics for clonal family trees of BCR sequences; while sumrep currently contains some tree functions, these functions have not yet been tested on experimental data or incorporated in comparison routines.
