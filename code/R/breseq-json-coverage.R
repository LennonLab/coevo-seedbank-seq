# parse coverage from breseq Json

setwd("/N/slate/danschw/coevo-seedbank-seq/")
library(here)
library(tidyverse)
library(cowplot)
library(jsonlite)

pops <- list.dirs(here("data/map-EVOL/host/breseq2/"), recursive = F)

# prep tibble
d <- tibble()

for (p in pops){
  
  js <- fromJSON(here(p, "/output/summary.json")) 
  
  if(str_detect(p, "SN")){
    d <- tibble(pop = str_remove(p, ".*/"),
           coverage = js$references$reference$Exported$coverage_average,
           dispersion = js$references$reference$Exported$coverage_dispersion,
           total_fraction_aligned_reads = js$reads$total_fraction_aligned_reads ) %>% 
      bind_rows(d, .)
    
  } else{
    d <- tibble(pop = str_remove(p, ".*/"),
           coverage = js$references$reference$NZ_CP015975$coverage_average,
           dispersion = js$references$reference$NZ_CP015975$coverage_dispersion,
           total_fraction_aligned_reads = js$reads$total_fraction_aligned_reads) %>% 
      bind_rows(d, .)
  }

}


d <- d %>% 
  mutate(subpop = case_when(str_detect(pop,"pl") ~ "pellet",
                            str_detect(pop,"rS") ~ "revived_spore",
                            str_detect(pop, "rV") ~ "revived_total")) %>% 
  mutate(subpop=as.factor(subpop))%>%
  mutate(subpop = fct_relevel(subpop, "pellet","revived_total","revived_spore")) %>% 
  
  #ajust factor order for seed bank
  mutate (seed.bank = case_when(str_starts(pop, "WL") ~ "long seed bank",
                                str_starts(pop, "WS") ~ "short seed bank",
                                str_starts(pop, "SN") ~ "no seed bank")) %>%
  mutate(seed.bank=as.factor(seed.bank))%>%
  mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank")) %>% 
  #ajust factor order for seed bank
  mutate (phage = case_when(str_detect(pop, "O") ~ "SPO1",
                            str_detect(pop, "Ct") ~ "no phage")) %>% 
  mutate(line = str_extract(pop, "-L.-") %>%  str_remove_all("-"))

d %>% 
  ggplot(aes(seed.bank, coverage))+
  geom_pointrange(aes(ymin = coverage-dispersion, ymax = coverage+dispersion), shape=21,
                  position = position_jitter(width = .1, height = 0))+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  facet_grid(phage ~ subpop)+
  ylim(0, NA)+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 35, hjust = 1))+
  ggsave(here("plots","host_coverage_total.png"))

d %>% 
  ggplot(aes(seed.bank, total_fraction_aligned_reads))+
  geom_point( size=2, shape=21,position = position_jitter(width = .1, height = 0))+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  facet_grid(phage ~ subpop)+
  ylim(0, 1)+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 35, hjust = 1))+
  ggsave(here("plots","host_fract_aligned.png"))
