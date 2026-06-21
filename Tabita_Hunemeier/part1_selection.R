### Load libraries
library("dplyr")
library("vroom")
library("ggplot2")

## Load files
eas_eur <- vroom("EAS_EUR.weir.fst")
afr_eas <- vroom("AFR_EAS.weir.fst")
afr_eur <- vroom("AFR_EUR.weir.fst")

## Actually, no duplicates
## remove NAS
eas_eur <- na.omit(eas_eur)
afr_eas <- na.omit(afr_eas)
afr_eur <- na.omit(afr_eur)

## Set negative values to zero
eas_eur$WEIR_AND_COCKERHAM_FST[eas_eur$WEIR_AND_COCKERHAM_FST < 0] <- 0
afr_eas$WEIR_AND_COCKERHAM_FST[afr_eas$WEIR_AND_COCKERHAM_FST < 0] <- 0
afr_eur$WEIR_AND_COCKERHAM_FST[afr_eur$WEIR_AND_COCKERHAM_FST < 0] <- 0

## Combine datasets
final_merged <- eas_eur %>%
  inner_join(afr_eas, by = c("CHROM", "POS")) %>%
  inner_join(afr_eur, by = c("CHROM", "POS"))

## Filter SNP
SNP_df <- final_merged %>% 
  filter(POS == 109513601)

## Plot pairwise Fst around 109513601 in a 10kb window, highlighting the candidate SNP
## filter to keep 10kb first
SNP_df <- final_merged %>% 
  filter(POS >= 109503601 & POS <= 109523601)

## Plotting per pair-wise plot
p1 <- ggplot(SNP_df, aes(x = POS, y = WEIR_AND_COCKERHAM_FST.x)) +
  geom_point() +
  geom_vline(xintercept = 109513601, color = "pink") +
  labs(title = "EAS vs EUR Fst", x = "Position", y = "Fst")
p1

## Second plot
p2 <- ggplot(SNP_df, aes(x = POS, y = WEIR_AND_COCKERHAM_FST.y)) +
  geom_point() +
  geom_vline(xintercept = 109513601, color = "purple") +
  labs(title = "AFR vs EAS Fst", x = "Position", y = "Fst")
p2

## Third plot
p3 <- ggplot(SNP_df, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point() +
  geom_vline(xintercept = 109513601, color = "blue") +
  labs(title = "AFR vs EUR Fst", x = "Position", y = "Fst")
p3

## merge the plots
library("cowplot")
all_plots <- plot_grid(p1, p2, p3, ncol = 1)
all_plots

## Code for PBS for this region
pbs_snp_region <- SNP_df %>% 
  mutate(
    PBS = (WEIR_AND_COCKERHAM_FST.y +
             WEIR_AND_COCKERHAM_FST.x -
             WEIR_AND_COCKERHAM_FST / 2
  ))

## Convert branches to zero
pbs_snp_region <- pbs_snp_region %>% 
  mutate(PBS = ifelse(PBS < 0, 0, PBS))

## plot
p3 <- ggplot(pbs_snp_region, aes(x = POS, y = PBS)) +
  geom_point() +
  geom_vline(xintercept = 109513601, color = "pink") +
  labs(title = "PBS", x = "Position", y = "PBS") +
  theme_bw()
p3

#####################
## Second part 
eigenval <- read.table("plink_pca.eigenval")
eigenvec <- read.table("plink_pca.eigenvec")

## rename
colnames(eigenvec) <- c("sampleName","sampleName2", paste0("PC", 1:(ncol(eigenvec)-2)))

## merge with population information
info <- vroom("sample_info.txt")

## merge with breed
pca <- eigenvec %>% 
  inner_join(info, by = c("sampleName"))

## plot PCA
pca_points <-
  ggplot(pca, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = paste0("PC1: " ,eigenval[1,]), y = paste0("PC2: " ,eigenval[2,]))
pca_points


## run pcaadapt
library("pcadapt")

## how to run
bed_file <- read.pcadapt("subset_chr15.bed", type = "bed")
pcadapt_results <- pcadapt(input = bed_file, K = 2, LD.clumping = list(size = 500, thr = 0.1))

## plotting
plot(pcadapt_results , option = "manhattan")
