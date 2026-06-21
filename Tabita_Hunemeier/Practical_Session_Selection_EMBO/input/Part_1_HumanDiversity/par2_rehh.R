## Load libraries
library("vroom")
library("rehh")

## Transform the data into haplohh format
hh <- data2haplohh(hap_file = "Chr2_EDAR_CHS_500K.recode.vcf",
                   polarize_vcf = FALSE)

## Calculate ehh
res <- calc_ehh(hh, mrk = "rs3827760")

## plot rees
plot(res)

##
hh2 <- data2haplohh(hap_file = "Chr2_EDAR_LWK_500K.recode.vcf",
                    polarize_vcf = FALSE)

## Calculate ehh
res2 <- calc_ehh(hh2, mrk = "rs3827760")

## plot rees
plot(res2)
