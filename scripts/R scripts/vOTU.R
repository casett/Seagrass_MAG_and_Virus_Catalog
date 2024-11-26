library(phyloseq)
library(vroom)

tpmean.75 <- vroom("data/coverm_relative_abundance/votu.tmean.tsv")
meta <- read.csv("data/vOTU_metadata.csv")
tax <- readRDS("results/tax_vc_phyloseq_ictv.RDS") 

rownames(meta) <- meta$SampleID

tpmean.75 <- tpmean.75 %>% arrange(Contig)
tpmean.75 <- as.data.frame(tpmean.75)
rownames(tpmean.75) <- tpmean.75$Contig

tax <- tax %>% arrange(Genome)
rownames(tax) <- tax$Genome

tax <- as.matrix(mixed.tax[-c(1)])

#import into phyloseq
otu_tab = otu_table(tpmean.75[-c(1)], taxa_are_rows = TRUE)
mapping_file = sample_data(meta)
tax_file = tax_table(tax)

vOTU.names <- as.data.frame(rownames(otu_tab))
vOTU.names$votu.id <- paste0("vOTU", 1:nrow(vOTU.names))

#rename votus
rownames(otu_tab) <- paste0("vOTU", 1:nrow(otu_tab))
rownames(tax_file) <- paste0("vOTU", 1:nrow(tax_file))




ps <- phyloseq(otu_tab, tax_file, mapping_file)

ps <-prune_samples(sample_sums(ps) > 0, ps) 
#remove samples with no viruses after removing low-quality


#beta div
library(microbiome)

#votus.high <- unique(na.omit(checkv_dRep_tax$votu.id[checkv_dRep_tax$checkv_quality %in% c("High-quality", "Complete", "Medium-quality") ]))

votus.high <- unique(na.omit(checkv_dRep_tax_high$votu.id))


ps.HQ <- prune_taxa(votus.high, ps)
ps.HQ <-prune_samples(sample_sums(ps.HQ) > 0, ps.HQ) 

ps.host <- subset_samples(ps.HQ, Tissue != 'sediment')

ps_hell <- transform(ps.HQ, "hellinger")

ps_hell_ord <- ordinate(physeq = ps_hell, method = "PCoA", distance = "euclidean")

plot_ordination(physeq = ps_hell, 
                ordination = ps_hell_ord, color = "SeagrassSpecies", shape = "TissueLocation")

plot_ordination(physeq = ps_hell, 
                ordination = ps_hell_ord, color = "SeagrassSpecies", shape = "MetaType")

plot_ordination(physeq = ps_hell, 
                ordination = ps_hell_ord, color = "SeagrassSpecies", shape = "Tissue")

plot_ordination(physeq = ps_hell, 
                ordination = ps_hell_ord, color = "SeagrassSpecies", shape = "StudyLastAuthor")

#alpha div 
plot_richness(ps, measures = c("Shannon"), color="SeagrassSpecies", shape="Tissue")



#relative abundance
se <- function(x) sqrt(var(x)/length(x))


ps.RA <- transform_sample_counts(ps.HQ,  function(x) x/sum(x))

# ps.RA.filt = filter_taxa(ps.RA, function(x) mean(x) > 
#                                    0.001, TRUE)

df_ps.RA.filt <- psmelt(ps.RA)

df_ps.RA.filt <- df_ps.RA.filt %>% 
  mutate(Class = ifelse(is.na(Class), "Unclassified", as.character(Class))) %>%
  mutate(Class2 = ifelse(Class == "Unassigned", ifelse(ClusterStatus == "Clustered", "Unclassified", "Unclustered"), as.character(Class)))

grouped <- group_by(df_ps.RA.filt, SampleID, Class2)


avgs_grouped <- summarise(grouped, mean = 100 * mean(Abundance), 
                          sd = 100 * sd(Abundance), se = 100 * se(Abundance))

plot = ggplot(grouped, aes(x = SampleID, y = Abundance, fill = Class2)) + 
  geom_bar(stat = "identity", position = "stack") 

plot + facet_wrap(TissueLocation~MetaType, scales = "free") + 
  theme(axis.text.x = element_text(angle = -70, hjust = 0, vjust = 0.5), text = element_text(size = 18)) + 
  ylab("Mean Relative Abundance") + xlab("vOTU") 





