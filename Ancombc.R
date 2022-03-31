library("devtools")
install_github("opisthokonta/tsnemicrobiota")
library(phyloseq)
library(tsnemicrobiota)
library(ggplot2)
library(microbiome)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("BDMMAcorrect"))
phy <- readRDS("/media/edwin/Transcend/datasets/dada2/meta_analysis/meta_analysis_phyloseq.RDS")
TREE = read_tree("/media/edwin/Transcend/datasets/dada2/meta_analysis/picrust/meta.tre")
phy_tree(phy) <- TREE

# Check the number of reads per sample
reads <- sample_sums(phy) 
head(reads)
table(tax_table(phy)[, "Kingdom"], exclude = NULL)
ps <- subset_taxa(phy, !is.na(Species) & !Phylum %in% c("", "uncharacterized"))

filterKingdom = c("Archaea", "NA")
phy = subset_taxa(ps, !Kingdom %in% filterKingdom)


# Check the number of reads per sample
reads <- sample_sums(phy) 
reads

# Standardize sample read count to compare between samples
## Set (arbitrary) cutoff for number of acceptable reads/sample
length(which(reads<5000))
## Find median sample read count
total = median(sample_sums(phy))
## Function to standardize to median sample read count
stand_f = function(x, t=total) round(t * (x / sum(x)))

# Apply to phyloseq object
phy_std = transform_sample_counts(phy, stand_f)
ntaxa(phy_std)
ps.taxa <- tax_glom(phy_std, taxrank = 'Species', NArm = TRUE)
M.f = filter_taxa(ps.taxa,function(x) sum(x > 100) > (0.2*length(x)) | sum(x) > 0.01*total, TRUE)
ntaxa(M.f)

tax_mat = as(tax_table(phylum_data), "matrix")

M.f = filter_taxa(ps.taxa,function(x) sum(x > 100) > (0.2*length(x)) | sum(x) > 0.01*total, TRUE)
ntaxa(M.f)
ps.taxa.pse.sub <- subset_samples(M.f, Conditions %in% c("Healthy_control", "Gastritis"))
# Run ancombc function
out = ancombc(phyloseq = ps.taxa.pse.sub, formula = "Conditions",
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
              group = "Conditions", struc_zero = TRUE, neg_lb = TRUE,
              tol = 1e-5, max_iter = 1000, conserve = TRUE,
              alpha = 0.05, global = TRUE)
out$delta_em
out$delta_wls
head(out$res$diff_abn)
diffabun <- out$res$diff_abn
diffabun <- rename(diffabun, c("diffabun" = "ConditionsHealthy_control"))
diffabun <- subset(diffabun, diffabun == 'TRUE')
log2foldchan = cbind(as(diffabun, "data.frame"), 
                     out$res$beta[rownames(diffabun), ])
log2foldchan <- rename(log2foldchan, c("ConditionsHealthy_control" = "out$res$beta[rownames(diffabun), ]"))
head(log2foldchan)
length(log2foldchan)
log2foldchan1 <- subset(log2foldchan, ConditionsHealthy_control >= 1 | ConditionsHealthy_control <= -1 ,)
head(log2foldchan1)
sigtab = cbind(as(log2foldchan1, "data.frame"), as(tax_table(ps.taxa)[rownames(log2foldchan1), ], "matrix"))

out$res$p_val
sigtab1 = cbind(as(log2foldchan1, "data.frame"), 
                out$res$q_val[rownames(log2foldchan1), ])
sigtab1<- rename(sigtab1, c("adjust_P" = "out$res$q_val[rownames(log2foldchan1), ]"))
sigtab2 = cbind(as(sigtab1, "data.frame"), as(tax_table(ps.taxa)[rownames(sigtab1), ], "matrix"))
head(sigtab2)

samp_frac = out$samp_frac
length(samp_frac)
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
rownames(samp_frac)
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(abundances(ps.taxa.pse.sub) + 1) 
# Adjust the log observed abundances
log_obs_abn_adj = as.data.frame(t(t(log_obs_abn) - samp_frac))
max(log_obs_abn_adj$ERR2014724)
abun = cbind(as(sigtab2, "data.frame"), 
             log_obs_abn_adj[rownames(sigtab2), ])
abun[1:5,1:15]
setwd("/media/edwin/Transcend/datasets/dada2/meta_analysis/DA")
write.table(abun, "Gastritis Vs Cancer", sep = "\t")

# Show the first 6 samples
log_obs_abn_adj[1:5, 1:5]


# Look at the ASVs that were significantly different between the two BV groups
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x = tapply(sigtab2$ConditionsHealthy_control, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Phylum = factor(as.character(sigtab2$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab2$ConditionsHealthy_control, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)

sigtab$Species = factor(as.character(sigtab2$Species), levels=names(x))
ggplot(log2foldchan, aes(x=ConditionsHealthy_control, y= Species, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ggtitle("Differential Abundance of Metaplasia and Cancer Individuals")

