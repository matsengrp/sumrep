source("SummaryStats.R")

db <- read.csv("Data/test_data.csv")

seqs <- db$SEQUENCE_VDJ %>% sapply(toString)

sequence.1 <- c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C")
sequence.2 <- c(A="TTGGC", B="ATGGG", C="ATGGG", D="ATGAC")
sequence.3 <- c(A="CGAAC", B="GGAAT", C="TACCC", D="ATAGCAC")

c1 <- compare.pairwise.distance.distribution(sequence.1, sequence.2)
c2 <- compare.pairwise.distance.distribution(sequence.1, sequence.3)
c3 <- compare.pairwise.distance.distribution(sequence.2, sequence.3)

cat(c1, c2, c3, '\n')
