#pre-processing checkv / dram results to make summary files

#Connect CheckV results for each viral genome back to vOTUs and to VCs 
#(clustered levels based on similarity)

#combine individual CheckV file and connect to 'Genome' name used in VContact2 results
CV_data_path <- "data/checkV_results/"   # path to the data
CV_files <- dir(CV_data_path, pattern = "*.tsv") # get file names


CV_data <- data_frame(filename = CV_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(CV_data_path, .))) # a new data column
  )  
CV_data

CV_data_un <- unnest(CV_data, cols = c(file_contents))

CV_data_un <- CV_data_un %>% 
  mutate_at("filename", str_replace, "_megahit_coassembly.checkv.quality_summary.tsv", "") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly.checkv.quality_summary.tsv", "Ettinger")

CheckV_results <- CV_data_un %>%
  mutate(Genome = paste0(filename, "_", contig_id)) %>%
  separate_wider_delim(Genome, delim = "|", names = c("GenomeID"), too_many = "drop")


saveRDS(CheckV_results, "results/checkv_results.RDS")
write_csv(CheckV_results, "results/checkv_results.csv")


#Get cluster membership information from dRep, to connect to vOTUs
drep <- vroom("data/dRep_votus/Cdb.csv")
drep_reps <- vroom("data/dRep_votus/Widb.csv")


drep <- drep %>%
  mutate_at("genome", str_replace, ".fa", "") %>%
  select(genome, secondary_cluster)

drep_reps <- drep_reps %>%
  mutate_at("genome", str_replace, ".fa", "") %>%
  mutate(Representative = TRUE) %>%
  select(genome, Representative)

checkv_dRep <- left_join(CheckV_results, drep, by = c("GenomeID" ="genome")) %>%
  left_join(drep_reps, by = c("GenomeID" ="genome"))

saveRDS(checkv_cdhit, "results/checkv_drep_results.RDS")
write_csv(checkv_cdhit, "results/checkv_drep_results.csv")

#subset to get rep seqs for each cluster, 
#combine with taxonomy VC 
#and then remove genome ID and extra info
#then combine back based on clstr # to get taxonomy for all 

checkv_dRep_rep <- checkv_dRep %>% 
  filter(Representative == TRUE)

# genome.SG.justtax<- genome.SG.justtax %>%
#   mutate_at("Genome", str_replace, "RubioProtillo", "Rubio_Protillo")
# 

checkv_drep_rep_tax <- left_join(checkv_dRep_rep, genome.SG.justtax, by=c("GenomeID" = "Genome"))



## ADD votu ids

#load from phyloseq file

vOTU.names


vOTU.names$Genome = vOTU.names$`rownames(otu_tab)`

vOTUs.ids.for.join <- vOTU.names[-c(1)]

checkv_drep_rep_tax_votu <- left_join(checkv_drep_rep_tax, vOTUs.ids.for.join, by=c("GenomeID" = "Genome"))



checkv_drep_rep_tax_cltr <- checkv_drep_rep_tax_votu %>% 
  select(secondary_cluster, votu.id)

checkv_dRep_tax <- left_join(checkv_dRep, checkv_drep_rep_tax_cltr)

saveRDS(checkv_dRep_tax, "results/checkv_drep_results_with_tax_ictv.RDS")
write_csv(checkv_dRep_tax, "results/checkv_drep_results_with_tax_ictv.csv")










##note DRAM files have different viral contig names and need a little love to connect back to checkv/vcontact2 reuslts

#Combine dram distill quality info 
DD_data_path <- "data/dramv_results/"   # path to the data
DMAG_files <- dir(DD_data_path, pattern = "*.vMAG_stats.tsv") # get file names


DMAG_data <- data_frame(filename = DMAG_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(DD_data_path, .))) # a new data column
  )  

DMAG_data_un <- unnest(DMAG_data, cols = c(file_contents))

DMAG_data_un <- DMAG_data_un %>% 
  mutate_at("filename", str_replace, "_megahit_coassembly_dramv-distill_vMAG_stats.tsv", "") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_1000_vMAG_stats.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_2000_vMAG_stats.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_3000_vMAG_stats.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_4000_vMAG_stats.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_4272_vMAG_stats.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Schorn_megahit_coassembly_dramv-distill_0_vMAG_stats.tsv", "Schorn") %>%
  mutate_at("filename", str_replace, "Schorn_megahit_coassembly_dramv-distill_1000_vMAG_stats.tsv", "Schorn") %>%
  mutate_at("filename", str_replace, "Schorn_megahit_coassembly_dramv-distill_2000_vMAG_stats.tsv", "Schorn") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_0_vMAG_stats.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_1000_vMAG_stats.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_2000_vMAG_stats.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_3000_vMAG_stats.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_4000_vMAG_stats.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_5000_vMAG_stats.tsv", "SoginMETA") 
  


DMAG_data_un <- DMAG_data_un %>% 
  # mutate_at("...1", str_replace, "-cat_1", "") %>%
  # mutate_at("...1", str_replace, "-cat_2", "") %>%
  # mutate_at("...1", str_replace, "-cat_3", "") %>%
  # mutate_at("...1", str_replace, "__", "||") %>%
  # mutate(Genome_old = paste0(filename, "_", ...1)) %>%
  separate_wider_delim("...1", delim = "__", names = c("Genome"), too_many = "drop") %>%
  mutate(Genome = paste0(filename, '_', Genome))




saveRDS(DMAG_data_un[-c(1)], "results/dram_vmag_stats.RDS")
write_csv(DMAG_data_un[-c(1)], "results/dram_vmag_stats.csv")





dram_vmag_stats <- left_join(checkv_drep_tax, DMAG_data_un[-c(1)])

saveRDS(dram_vmag_stats, "results/checkv_dram_vmag_stats_with_tax.RDS")
write_csv(dram_vmag_stats, "results/checkv_dram_vmag_stats_with_tax.csv")


#combine dram AMG gene info

DD_data_path <- "data/dramv_results/"   # path to the data
DAMG_files <- dir(DD_data_path, pattern = "*amg_summary.tsv") # get file names


DAMG_data <- data_frame(filename = DAMG_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(DD_data_path, .)) %>% mutate(auxiliary_score = as.character(auxiliary_score)))) 



DAMG_data_un <- unnest(DAMG_data, cols = c(file_contents))

DAMG_data_un <- DAMG_data_un %>% 
  mutate_at("filename", str_replace, "_megahit_coassembly_dramv-distill_amg_summary.tsv", "") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_1000_amg_summary.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_2000_amg_summary.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_3000_amg_summary.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_4000_amg_summary.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_dramv-distill_4272_amg_summary.tsv", "Ettinger") %>%
  mutate_at("filename", str_replace, "Schorn_megahit_coassembly_dramv-distill_0_amg_summary.tsv", "Schorn") %>%
  mutate_at("filename", str_replace, "Schorn_megahit_coassembly_dramv-distill_1000_amg_summary.tsv", "Schorn") %>%
  mutate_at("filename", str_replace, "Schorn_megahit_coassembly_dramv-distill_2000_amg_summary.tsv", "Schorn") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_0_amg_summary.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_1000_amg_summary.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_2000_amg_summary.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_3000_amg_summary.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_4000_amg_summary.tsv", "SoginMETA") %>%
  mutate_at("filename", str_replace, "SoginMETA_megahit_coassembly_dramv-distill_5000_amg_summary.tsv", "SoginMETA") 


DAMG_data_un <- DAMG_data_un %>% 
  separate_wider_delim("scaffold", delim = "__", names = c("Genome"), too_many = "drop") %>%
  mutate(Genome = paste0(filename, '_', Genome))


  
saveRDS(DAMG_data_un, "results/dram_amg_results.RDS")
write_csv(DAMG_data_un, "results/dram_amg_results.csv")





### virsorter boundaries

VS_data_path <- "data/virsorter2/"   # path to the data
VS_files <- dir(VS_data_path, pattern = "*.final-viral-score.tsv") # get file names


VS_data <- data_frame(filename = VS_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(VS_data_path, .))))


VS_data_un <- unnest(VS_data, cols = c(file_contents))

VS_data_un <- VS_data_un %>%
  mutate_at("filename", str_replace, ".final-viral-score.tsv", "") %>%
  mutate_at("filename", str_replace, "_megahit_coassembly_vs2", "") %>%
  mutate_at("filename", str_replace, "Cassie_seagrass_metagenome_chytrids_coassembly_vs2", "Ettinger")


VS_data_un <- VS_data_un %>%
  separate_wider_delim("seqname", delim = "||", names = c("Genome"), too_many = "drop") %>%
  mutate(Genome = paste0(filename, "_", Genome))

saveRDS(VS_data_un, "results/virsorter_score_results.RDS")
write_csv(VS_data_un, "results/virsorter_score_results.csv")





checkv_dRep_tax_high <- checkv_dRep_tax %>% 
  left_join(VS_data_un, by=c("GenomeID" = "Genome")) %>%
  #filter(miuvig_quality == "High-quality") %>%
  filter(checkv_quality %in% c("Medium-quality", "High-quality", "Complete")) %>%
  filter(contig_length >= 10000) %>%
  filter(max_score_group != "NCLDV")


saveRDS(checkv_cdhit_tax_high, "results/checkv_cdhit_results_with_tax_highquality.RDS")
write_csv(checkv_cdhit_tax_high, "results/checkv_cdhit_results_with_tax_highquality.csv")




host_predictions <- vroom("data/iHop_predictions/Host_prediction_to_genome_m90.csv")
host_predictions_genus <- vroom("data/iHop_predictions/Host_prediction_to_genus_m90.csv")

host_predictions_genus <- host_predictions_genus %>%
  separate_wider_delim(cols= `Host genus`, delim = ";", names = c("Host_Domain", "Host_Phylum","Host_Class","Host_Order","Host_Family","Host_Genus"), too_few = "align_start")

checkv_dRep_tax_high %>% 
  left_join(host_predictions_genus, by=c("GenomeID"= "Virus")) %>%
  left_join(host_predictions, by=c("GenomeID" = "Virus")) %>%
  filter(!is.na(Host_Class)) %>%
  group_by(Host_Class) %>%
  tally() %>%
  arrange(desc(n)) %>%
  group_by(Host_Class) %>%
  mutate(n2 = sum(n)) %>%
  mutate(lab_ypos = n2 + 3) %>%
  ggplot(aes(x = Host_Class, fill = Host_Class, y = n)) + theme_bw() +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_viridis_d(option = "B") + ylab("Number of viral sequences") +
  xlab("") + theme(legend.position = "none") + coord_flip() +
  scale_x_discrete(limits = rev) 



