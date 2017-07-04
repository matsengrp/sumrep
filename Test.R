source("SummaryStats.R")

db <- read.csv("Data/test_data.csv")

seqs <- db$SEQUENCE_VDJ %>% sapply(toString)

sequence.1 <- c(A="ATGGC", B="ATGGC", C="ATGGG", D="AT--C")
sequence.2 <- c(A="TTGGC", B="ATGGC", C="ATAGG", D="ATGAC")
sequence.3 <- c(A="CGAAC", B="GGAAT", C="TACCC", D="ATAGCAC")

c1 <- compare.pairwise.distance.distribution(sequence.1, sequence.2)
c2 <- compare.pairwise.distance.distribution(sequence.1, sequence.3)
c3 <- compare.pairwise.distance.distribution(sequence.2, sequence.3)

cat(c1, c2, c3, '\n')

d1 <- compare.NN.distance.distribution(sequence.1, sequence.2)
d2 <- compare.NN.distance.distribution(sequence.1, sequence.3)
d3 <- compare.NN.distance.distribution(sequence.2, sequence.3)
cat(d1, d2, d3, '\n')
