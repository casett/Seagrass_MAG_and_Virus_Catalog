#network plotting

#load in network files, made in vcontact2_processing.R
nodes <- readRDS("results/ntwk_nodes.RDS")
edges <- readRDS("results/ntwk_edges.RDS")
genome.ov <- readRDS("results/genome_vc_master_ictv.RDS")
genome.SG.justtax <- readRDS("results/tax_vc_phyloseq_ictv.RDS")

clstrd.nodes <- nodes %>% 
  left_join(genome.ov, by = "Genome") %>% 
  filter(VC.Status == "Clustered") 

#can plot networks based on vcontact2

v.ids <- genome.SG.justtax$Genome

ntwk.p <- clstrd.nodes %>% 
  ggplot(aes(x, y)) +
  geom_line(data = filter(edges, Genome %in% clstrd.nodes$Genome), aes(group = Pair), color = "black", size = 0.5, alpha = 0.1) +
  geom_point(data = filter(clstrd.nodes, !Genome %in% c(v.ids)), alpha = 0.8, shape = 16, size = 2, color = "gray") +
  geom_point(data = filter(clstrd.nodes, Genome %in% v.ids), alpha = 0.8, shape = 16, size = 2, color = RColorBrewer::brewer.pal(9, "Blues")[6]) +
  theme_minimal() +
  theme(text = element_text(size = 15), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

#legend
legend.p <- data.frame(Type = c("RefSeq", "SG"),
                       Value = 1:2) %>% 
  ggplot(aes(1, Value, color = Type)) +
  geom_point(size = 5) +
  scale_color_manual(name = "Source", values = c("gray", RColorBrewer::brewer.pal(9, "Blues")[6])) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = "right")

ntwk.p

#plot by fam
genome.ov.SG <- left_join(genome.ov, genome.SG.justtax)

vir.fam.clust <- genome.ov.SG %>% 
  filter(Source == "SG") %>%  
  select(Genome, VC, VC.Status) %>%
  inner_join(clstr.master, by = "VC") %>% 
  filter(ClstrComp == "both") %>% 
  filter(Class != "Unassigned") %>% 
  group_by(Class, Family) %>% 
  count() %>% 
  group_by(Class) %>% 
  mutate(nOrd = sum(n)) %>% 
  arrange(nOrd, n) %>% 
  ungroup() %>% 
  mutate(Rank = 1:n()) 

vir.fam.nodes <- clstrd.nodes %>% 
  left_join(select(clstr.master, VC, ClstrComp), by = "VC") %>% 
  mutate(Family2 = case_when(Source == "SG" ~ "This Study",
                             Family %in% vir.fam.clust$Family ~ as.character(Family),
                             TRUE ~ "Other")) %>% 
  left_join(vir.fam.clust, by = c("Family2" = "Family")) %>% 
  filter(Source == "refseq" | ClstrComp == "both") 

# Plot
ntwk.vir <-  vir.fam.nodes %>% 
  ggplot(aes(x, y)) +
  geom_line(data = filter(edges, Genome %in% vir.fam.nodes$Genome), aes(group = Pair), alpha = 0.1, color = "gray25", size = 0.5) +
  geom_point(alpha = 0.8, size = 2, shape = 16, aes(color = Family2)) +
  scale_color_manual(name = "Virus Family",
                     values = c(RColorBrewer::brewer.pal(8, "Set2")[c(1:5)], "slateblue4","gray10" , "gray75", "orange")) +
  theme_minimal() +
  theme(text = element_text(size = 15), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")

ntwk.vir
