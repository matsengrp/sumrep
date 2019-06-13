source("Examples/RunAnnotations.R")

igb_germline_dir <- "~/Software/igblast/partis_friendly_bin"

writeAnnotations("~/Data/FV-igh-m1h.fa", 
                 "pi_f1",
                 "partis",
                 germline_dir=igb_germline_dir)
writeAnnotations("~/Data/FV-igh-m8d.fa", 
                 "pi_f2", 
                 "partis",
                 germline_dir=igb_germline_dir)
writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                 "pi_g1", 
                 "partis",
                 germline_dir=igb_germline_dir)
writeAnnotations("~/Data/GMC-igh-m8d.fa", 
                 "pi_g2", 
                 "partis",
                 germline_dir=igb_germline_dir)
writeAnnotations("~/Data/IB-igh-m1h.fa", 
                 "pi_i1", 
                 "partis",
                 germline_dir=igb_germline_dir)
writeAnnotations("~/Data/IB-igh-m8d.fa", 
                 "pi_i2", 
                 "partis",
                 germline_dir=igb_germline_dir)
