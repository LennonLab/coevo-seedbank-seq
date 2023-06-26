library(here)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggnewscale)


# host annotation --------------------------------------------------------

delta6_168 <- read_csv(here("data/teichoic_acid","delta6_168_cds_matched.csv"),trim_ws = T, name_repair = "universal")
categories_168 <- read_csv(here("data/teichoic_acid","geneCategories-2022-06-27.csv"),trim_ws = T, name_repair = "universal")
SW.export_168 <- read_csv(here("data/teichoic_acid","subtiwiki.gene.export.2022-06-27.csv"),trim_ws = T)


# Find hi-freq loci --------------------------------------------------------

# bottom threshold for mutation to be high frequency
freq_cutoff <- 0.3

f_freq <- list.files(here("data/timecourse_final_breseq/"))

f_freq <- f_freq[str_detect(f_freq, "host")]
f_freq <- f_freq[str_detect(f_freq, "total")]
f_freq <- f_freq[str_detect(f_freq, "WL|SN")]

hi_freq <- tibble()
d.plot <- tibble()

for(f in f_freq){
  d <- read_csv((here("data/timecourse_final_breseq/", f))) %>% 
    filter(Gene != "intergenic") %>%
    filter(Annotation != "synonymous") %>% 
    filter(Annotation != "noncoding") %>% 
    # unified label for coding and promoter regions
    mutate(full_gene = Gene, 
           Gene = str_remove(full_gene, "possible promotor region for"))
  
  # Mutations detected > twice
  mut_3 <- d %>% 
    #make uniqe id
    mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
    # pull(mut_id) %>%  anyDuplicated(.) # it is sufficient
    mutate(`Freq:0` = 0) %>% 
    select(mut_id, contains("Freq")) %>% 
    pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>% 
    mutate(t_sample = parse_number(t_sample))  %>%
    mutate(detected = freq > 0) %>%
    group_by(mut_id) %>% 
    summarise(n = sum(detected)) %>% 
    filter (n > 0) %>% 
    pull(mut_id)
  
  # Mutations fixed
  fixed <- d %>% 
    #make uniqe id
    mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
    # pull(mut_id) %>%  anyDuplicated(.) # it is sufficient
    mutate(`Freq:0` = 0) %>% 
    select(mut_id, contains("Freq")) %>% 
    # slice_head(n = 100) %>% 
    pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>% 
    mutate(t_sample = parse_number(t_sample))  %>%
    filter(freq > freq_cutoff) %>% 
    pull(mut_id)
  
  to_pull = intersect(mut_3, fixed)
  
  # if(length(to_pull) < 1) next
  
  hi_freq <- d %>% 
    #make uniqe id
    mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
    filter(mut_id %in% to_pull) %>% 
    mutate(f = f) %>% 
    bind_rows(hi_freq, .)
  
  #keep all mutations to plot
  d.plot <- d %>% 
    #make uniqe id
    mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
    mutate(observed_3times_or_more =(mut_id %in% mut_3)) %>%
    mutate(reached_hi_freq = (mut_id %in% fixed)) %>% 
    mutate(f = f) %>% 
    bind_rows(d.plot, .)
}


# gene name per locus ---------------------------------------------

  
gene_tags <-
  delta6_168 %>% 
  filter(locus_tag.d6 %in% unique(d.plot$Gene)) %>% 
  select(locus_tag.d6,locus_tag.168)


all_freq <-
  SW.export_168 %>% 
  filter(locus %in% gene_tags$locus_tag.168) %>% 
  left_join(. ,gene_tags, by = c("locus" = "locus_tag.168")) %>% 
  rename(Gene = locus_tag.d6, Name = title, Function = `function`) %>% 
  left_join(d.plot,., by = "Gene") 


# prepare table  --------------------------------------------------------------------


d <- all_freq %>%
  # parse treatments
  mutate(trt = str_extract(f, "^.*_L."),
         seed.bank = if_else(str_detect(trt,"WL."), 
                             "seed bank +","seed bank -"),
         phage=if_else(str_detect(trt,"O"), 
                       "phage +","phage -"),
         pop = parse_number(trt)) %>% 
  # make codon position 1-3 (rather than 0-2) 
  mutate(`Position in codon` = parse_number(`Position in codon`)+1) %>%
  select(seed.bank,phage,pop,delta6_locus = full_gene, Bs168_locus=locus, gene_tag=Name, description, Function,
         genome_position=Position, `Mutation type`,  `Position in gene`,`Codon position in gene`,
         Codon, `Position in codon`, mutation_effect=Annotation, new_allele=Allele, reached_hi_freq,observed_3times_or_more,
         contains("Freq")) %>%
    # slice_head(n = 100) %>%
    pivot_longer(contains("Freq:"), names_to = "t_sample", values_to = "freq") %>%
    mutate(t_sample = parse_number(t_sample))

write_csv(d, here("data/host_observed_mutations.csv"))
