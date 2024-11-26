#vcontact2 processing

#load libraries we are using
library(GGally)
library(tidyverse)
library(network)
library(vroom)
#this script will take in the vcontact2 results and connect them back to our samples/summarize them

#load in data
genome.ov <- read.table("data/vcontact2/genome_by_genome_overview.csv", header = T, sep = ",") #results from vcontact2
ntwk <- read.table("data/vcontact2/c1.ntw", header = F, sep = " ", col.names = c("OTU1", "OTU2", "Score")) #network from vcontact2
ictv_names <- vroom("data/vcontact2/1Aug2022_data.tsv")

samples <- c("Ettinger","Crump", "Cucio", "Fraser", "Fraser2023","Miranda", "Mohr","Rubio_Protillo","Schorn", "Sogin", "SoginMETA")

#define the source of each vOTU, either SG or from the reference database refseq
genome.ov <- genome.ov %>% 
  mutate(Source = ifelse(str_detect(Genome, paste(samples, collapse = "|")), "SG", "refseq"))

genome.ov <- left_join(genome.ov, ictv_names, by=c("Genome"="Accession"))

genome.ov <- genome.ov %>%
  mutate(Order = ifelse(is.na(Order.y), Order.x, Order.y)) %>%
  mutate(Family = ifelse(is.na(Family.y), Family.x, Family.y)) %>%
  mutate(Genus = ifelse(is.na(Genus.y), Genus.x, Genus.y))



#check that # contigs matches number expected vOTUs (3633)
genome.ov %>% 
  filter(Source == "SG") %>% 
  group_by(Genome) %>% 
  count()

#Generate a data frame that specifies the composition 
#of each viral cluster in terms of source of network 
#nodes (ie. does the cluster have sequences from seagrass?)
#ClstrComp = cluster composition, either from seagrass (SG), reference data(refseq) or both 
clstr.source <- genome.ov %>% 
  filter(VC.Status == "Clustered") %>% 
  filter(VC != "nan") %>% 
  mutate(SG = str_detect(Genome,  paste(samples, collapse = "|"))) %>% 
  group_by(VC, Size) %>% 
  summarise(pSG = sum(SG)/n()) %>% 
  mutate(ClstrComp = case_when(pSG == 0 ~ "refseq",
                               pSG == 1 ~ "SG",
                               TRUE ~ "both" )) 

#Let's create a data frame with the consensus viral taxonomy for each cluster. 
#the strategy is to check each taxonomic rank and see if there's only one assignment or a mix. 
#Some refseq entrees have unassigned ranks so we are ignoring those 
#and deciding the proper classification based on the rest of the VC members 

clstr.order <- genome.ov %>% 
  filter(VC.Status == "Clustered") %>% 
  filter(VC != "nan") %>% 
  filter(ICTV_Order != "Unclassified") %>%
  filter(!is.na(ICTV_Order)) %>% 
  group_by(VC, ICTV_Order) %>% 
  count() %>% 
  group_by(VC) %>% 
  mutate(Duplicates = n()) %>% 
  mutate(Order = ifelse(Duplicates > 1, "Unclassified", as.character(ICTV_Order))) %>% 
  group_by(VC, Order) %>% 
  count() %>% 
  select(-n)

clstr.family <- genome.ov %>% 
  filter(VC.Status == "Clustered") %>% 
  filter(VC != "nan") %>% 
  filter(ICTV_Family != "Unclassified") %>% 
  filter(!is.na(ICTV_Family)) %>% 
  group_by(VC, ICTV_Family) %>% 
  count() %>% 
  group_by(VC) %>% 
  mutate(Duplicates = n()) %>% 
  mutate(Family = ifelse(Duplicates > 1, "Unclassified", as.character(ICTV_Family))) %>% 
  group_by(VC, Family) %>% 
  count() %>% 
  select(-n)

clstr.genus <-  genome.ov %>% 
  filter(VC.Status == "Clustered") %>% 
  filter(VC != "nan") %>% 
  filter(ICTV_Genus != "Unclassified") %>% 
  filter(!is.na(ICTV_Genus)) %>% 
  group_by(VC, ICTV_Genus) %>% 
  count() %>% 
  group_by(VC) %>% 
  mutate(Duplicates = n()) %>% 
  mutate(Genus = ifelse(Duplicates > 1, "Unclassified", as.character(ICTV_Genus))) %>% 
  group_by(VC, Genus) %>% 
  count() %>% 
  select(-n)

clstr.phyla <-  genome.ov %>% 
  filter(VC.Status == "Clustered") %>% 
  filter(VC != "nan") %>% 
  filter(ICTV_Phylum != "Unclassified") %>% 
  filter(!is.na(ICTV_Phylum)) %>% 
  group_by(VC, ICTV_Phylum) %>% 
  count() %>% 
  group_by(VC) %>% 
  mutate(Duplicates = n()) %>% 
  mutate(Phylum = ifelse(Duplicates > 1, "Unclassified", as.character(ICTV_Phylum))) %>% 
  group_by(VC, Phylum) %>% 
  count() %>% 
  select(-n)

clstr.class <-  genome.ov %>% 
  filter(VC.Status == "Clustered") %>% 
  filter(VC != "nan") %>% 
  filter(ICTV_Class != "Unclassified") %>% 
  filter(!is.na(ICTV_Class)) %>% 
  group_by(VC, ICTV_Class) %>% 
  count() %>% 
  group_by(VC) %>% 
  mutate(Duplicates = n()) %>% 
  mutate(Class = ifelse(Duplicates > 1, "Unclassified", as.character(ICTV_Class))) %>% 
  group_by(VC, Class) %>% 
  count() %>% 
  select(-n)

#combine
clstr.master <- clstr.source %>% 
  left_join(clstr.order, by = "VC") %>% 
  left_join(clstr.family, by = "VC") %>% 
  left_join(clstr.genus, by = "VC") %>% 
  left_join(clstr.phyla, by = "VC") %>% 
  left_join(clstr.class, by = "VC") %>% 
  mutate(Order = ifelse(is.na(Order), "Unassigned", Order)) %>% 
  mutate(Family = ifelse(is.na(Family), "Unassigned", Family)) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Unassigned", Genus)) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unassigned", Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), "Unassigned", Class)) 

#make a taxonomy table for phyloseq
genome.SG <- genome.ov[which(genome.ov$Source=='SG'),]

genome.SG<- genome.SG[-c(1,3:5,19:59)]

#join the taxomy back to the main cluster dataset
genome.SG.tax <- left_join(genome.SG, clstr.master, by="VC")

#edit taxonomy a bit more for use in plots
genome.SG.tax <- genome.SG.tax %>%
  mutate(ClusterStatus = ifelse(is.na(Size.x), "Unclustered", "Clustered"))%>%
  mutate(VCStatus = ifelse(VC == "", "Unclustered", VC)) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unassigned", Phylum)) %>% 
  mutate(Class = ifelse(is.na(Class), "Unassigned", Class)) %>% 
  mutate(Order = ifelse(is.na(Order), "Unassigned", Order)) %>% 
  mutate(Family = ifelse(is.na(Family), "Unassigned", Family)) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Unassigned", Genus))

genome.SG.justtax <- genome.SG.tax[-c(2:17)]


genome.SG.justtax <- genome.SG.justtax %>%
  select(Genome, ClusterStatus, Phylum, Class, Order, Family, Genus, VCStatus)

#vcontact2 also makes a network of protein similarity
#it uses this to make VCs
#we can use this to plot viruses and color by taxonomy
#or host
#or genome quality
#first we need to import network info into R and save as R-readable files
#get network information - too big these steps done on cluster
# nodes <- ggnet2(ntwk[,-3], mode = "fruchtermanreingold", layout.par = list(list=(niter=2000)))$data %>%
#   rename("Genome" = "label")
# 
# edges <- ntwk %>%
#   mutate(Pair = paste(OTU1, OTU2, sep = ".")) %>%
#   gather(key = "Member", value = "Genome", -Pair, -Score) %>%
#   inner_join(nodes, by = "Genome")
# 
# #save files for later
# saveRDS(nodes, "results/ntwk_nodes.RDS")
# saveRDS(edges, "results/ntwk_edges.RDS")
# 

#save all the files we made above
saveRDS(genome.ov, "results/genome_vc_master_ictv.RDS")
write_csv(genome.ov, "results/genome_vc_master_ictv.csv")

saveRDS(clstr.master, "results/cluster_vc_master_ictv.RDS")
write_csv(clstr.master, "results/cluster_vc_master_ictv.csv")

saveRDS(genome.SG.justtax, "results/tax_vc_phyloseq_ictv.RDS")
write_csv(genome.SG.justtax, "results/tax_vc_phyloseq_ictv.csv")

