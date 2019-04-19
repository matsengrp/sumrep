source("Examples/RunAnnotations.R")

tcr_dir <- "/fh/fast/matsen_e/bolson2/mike_aging/converted"

writeAnnotations(file.path(tcr_dir, "A4-i194_R2.fa"),
                 "A4_i194",
                 "igor",
                 locus="trb"
                )
writeAnnotations(file.path(tcr_dir, "A5-S22_R2.fa"),
                 "A5_S22",
                 "igor",
                 locus="trb"
                )
writeAnnotations(file.path(tcr_dir, "A5-S9_R2.fa"),
                 "A5_S9",
                 "igor",
                 locus="trb"
                )
writeAnnotations(file.path(tcr_dir, "A5-S10_R2.fa"),
                 "A5_S10",
                 "igor",
                 locus="trb"
                )
writeAnnotations(file.path(tcr_dir, "A5-S15_R2.fa"),
                 "A5_S15",
                 "igor",
                 locus="trb"
                )
writeAnnotations(file.path(tcr_dir, "A4-i107_R2.fa"),
                 "A4_i107",
                 "igor",
                 locus="trb"
                )
