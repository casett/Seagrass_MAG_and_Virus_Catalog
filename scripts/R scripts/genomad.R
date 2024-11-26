

genomad <- vroom("data/genomad_tax/dRep_clustered_taxonomy.tsv")

genomad <- genomad %>%
  separate_wider_delim(cols= lineage, delim = ";", names = c("Domain","Realm", "Kingdom", "Phylum","Class","Order","Family"), too_few = "align_start")

genomad.tax <- genomad[-c(2,3,4)]

genomad.tax <- as.data.frame(genomad.tax)

genomad.tax <- genomad.tax  %>% arrange(seq_name)

rownames(genomad.tax) <- genomad.tax$seq_name

genomad.tax.phylo <- as.matrix(genomad.tax[-c(1)])


mixed.tax <- left_join(tax[-c(3:7)], genomad.tax, by=c("Genome" = "seq_name"))

rownames(mixed.tax) <- mixed.tax$Genome


mixed.tax <- mixed.tax[-c(2,4)] %>%
  select(Genome, Realm, Kingdom, Phylum, Class, Order, Family, VCStatus) 


