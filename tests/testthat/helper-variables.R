seq_1 <- c(A="AT", B="AT")
seq_2 <- list(A="AA", B="AC")
seq_3 <- c(A="", B="C-TA")
seq_4 <- c(A="", B="")
seq_5 <- c("AAA", "AAT", "TTT")
seq_6 <- c("A-T-", "GCAC", "CGGC", "ACGT")
seq_7 <- list(c("G", "T"), c("T", "A"))
seq_8 <- c("AAAAAA", "TTTTTT", "AAAA", "AAAAACC", "AAAAAAATTTT", "GGGCCCCCGG")
seq_9 <- list("AAAAAC", c("T", "T", "T", "C", "C", "T"), "AAAAGGGG", "AAAAACCCC", "CCCCTTTC", 
              list("G", "G", "G", "C", "C", "C", "C", "C", "G", "G"))

dat_1 <- data.frame(sequence_alignment=seq_5)
dat_2 <- data.frame(sequence_alignment=seq_6)
