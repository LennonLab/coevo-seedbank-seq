library(here)
library(tidyverse)
library(cowplot)


# Teichoic acid synthesis genes in delta6 ---------------------------------

delta6_168 <- read_csv(here("data/teichoic_acid","delta6_168_cds_matched.csv"),trim_ws = T, name_repair = "universal")
categories_168 <- read_csv(here("data/teichoic_acid","geneCategories-2022-06-27.csv"),trim_ws = T, name_repair = "universal")
SW.export_168 <- read_csv(here("data/teichoic_acid","subtiwiki.gene.export.2022-06-27.csv"),trim_ws = T)

# bsu tags of Teichoic acid biosynthesis genes (TA)
TA_bsu <- categories_168 %>% 
  filter(str_detect(category, "Biosynthesis of teichoic acid")) %>% 
  filter(!duplicated(gene.locus)) %>% 
  pull(gene.locus)


# delta6 tags of Teichoic acid biosynthesis genes (TA)
TA_delta6 <- delta6_168 %>% 
  filter(locus_tag.168 %in% TA_bsu) %>% 
  # add gene tags
  left_join(., SW.export_168 %>%
              filter(locus %in% TA_bsu), by = c("locus_tag.168"="locus"))


# host gene multiplicity --------------------------------------------------

gene_mult <- read_csv(here("data/mult_host.csv")) %>% 
  rename(trt = ...1) %>% 
  # long format
  pivot_longer(-1, names_to = "locus_tag.d6", values_to = "mult")

gene_mult_TA <- gene_mult %>% 
  filter(locus_tag.d6 %in% TA_delta6$locus_tag.d6)

# add gene tags
gene_mult_TA <- TA_delta6 %>%
  select(locus_tag.d6, gene = title) %>% 
  left_join(gene_mult_TA,.)

# parse treatments
gene_mult_TA <- gene_mult_TA %>% 
  mutate(phage = if_else(str_detect(trt, "SPO1"), "SPO1", "no Phage") %>% fct_rev(),
         seed.bank =if_else(str_detect(trt, "long"), "with-seed-bank", "no-seed-bank")%>% fct_rev(),
         replicate_pop = str_remove(trt, ".*_"))

# arrange genes
gene_mult_TA %>% 
  mutate(gene = as_factor(gene) %>% fct_relevel("tagF", after = 3))
  
  filter(phage == "SPO1") %>% 
  group_by(gene) %>% 
  summarise(m = median(mult)) %>% 
  arrange(m) %>% 
  pull(gene)


p <- gene_mult_TA %>% 
  mutate(`log10(m)` = if_else(mult>0,mult,  NA_real_) %>% log10()) %>% 
  mutate(gene = as_factor(gene) %>% fct_relevel("tagF", after = 3) %>% fct_rev()) %>% 
  ggplot(aes(replicate_pop, gene)) +
  geom_tile(aes(fill = `log10(m)`), color = "white", size=1)+
  scale_fill_gradient(na.value = "white", low = "grey80", high = "black")+
  facet_grid(phage ~ seed.bank)+
  theme_classic(base_size = 12)+
  panel_border(color = "black")+
  # guides(fill = guide_legend(title = "log10(multiplicity)"))+
  theme(legend.position = "bottom", 
         strip.background = element_blank())

ggsave(here("analysis/teichoic_acid.png"), p , width = 3, height = 4)

