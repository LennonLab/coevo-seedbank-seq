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
    filter(Annotation != "synonymous")
  
  # Mutations detected > twice
  mut_3 <- d %>% 
    #make uniqe id
    mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
    # pull(mut_id) %>%  anyDuplicated(.) # it is sufficient
    mutate(`Freq:0` = 0) %>% 
    select(mut_id, contains("Freq")) %>% 
    # slice_head(n = 100) %>% 
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
    filter(mut_id %in% mut_3) %>% 
    mutate(f = f) %>% 
    bind_rows(d.plot, .)
}


# gene name per hi-freq locus ---------------------------------------------

  
hi_freq_tags <-
  delta6_168 %>% 
  filter(locus_tag.d6 %in% unique(hi_freq$Gene)) %>% 
  select(locus_tag.d6,locus_tag.168)


hi_freq <-
  SW.export_168 %>% 
  filter(locus %in% hi_freq_tags$locus_tag.168) %>% 
  left_join(. ,hi_freq_tags, by = c("locus" = "locus_tag.168")) %>% 
  rename(Gene = locus_tag.d6, Name = title, Function = `function`) %>% 
  left_join(hi_freq,., by = "Gene") %>% 
  filter(if_any(.cols = contains("Freq"), ~ . > 0.4))
# categories_168 %>% 
#   filter(gene.locus %in% hi_freq_tags) %>% view


# Plot --------------------------------------------------------------------


d <- d.plot %>%
  # parse treatments
  mutate(trt = str_extract(f, "^.*_L."),
         seed.bank = if_else(str_detect(trt,"WL."), 
                             "with seed bank","w/o seed bank"),
         phage=if_else(str_detect(trt,"O"), 
                       "with phage","w/o phage"),
         pop = parse_number(trt)) %>% 
  #add T = 0  
  mutate(`Freq:0` = 0) %>%
  select(mut_id,trt, seed.bank,phage,pop, contains("Freq")) %>%
    # slice_head(n = 100) %>%
    pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>%
    mutate(t_sample = parse_number(t_sample))

# hi freq trajectories
d.hi <- 
  d %>% 
  filter(mut_id %in% hi_freq$mut_id) %>% 
  # check that hi freq
  mutate(hi = freq>freq_cutoff) %>% 
  group_by(mut_id, trt) %>% 
  filter(any(hi)) %>% 
  
  left_join(., select(hi_freq, mut_id, Name), by = "mut_id") %>% 
  mutate(Name = if_else(str_detect(mut_id,"intergenic"), "intergenic", Name)) %>% 
  #adjust legend order
  arrange(phage, seed.bank) %>% 
  mutate(Name = fct_inorder(Name)) %>% 
  mutate(phage = as_factor(phage) %>% fct_rev(),
         seed.bank = as_factor(seed.bank) %>% fct_rev()) 

# hi freq names
lab.hi <- d.hi %>% 
  group_by(seed.bank, phage, pop) %>% 
  summarise(Name) %>% 
  distinct() %>% 
  mutate(y.idx = row_number(),
         y.pos = 0.98- 0.12*(y.idx-1)) %>% 
  mutate(phage = as_factor(phage) %>% fct_rev(),
         seed.bank = as_factor(seed.bank) %>% fct_rev())
# plot

  p <- d %>%
    mutate(phage = as_factor(phage) %>% fct_rev(),
           seed.bank = as_factor(seed.bank) %>% fct_rev()) %>% 
    ggplot(aes(t_sample, freq))+
    geom_line(aes(color = mut_id), size=.8, show.legend = F)+
    scale_color_grey(guide = "none")+
    
    new_scale_color() +
    geom_line(data = d.hi,
              aes(group = mut_id,color = Name), size=.8,alpha = 0.7 ,show.legend = T)+
    geom_point(data = d.hi,
              aes(group = mut_id, color = Name), shape=21, fill="white",
              alpha = 0.7 ,size=.8, show.legend = T)+
    scale_color_viridis_d()+
    
    geom_text(data = lab.hi, 
              aes(label = Name, color = Name, y=  y.pos),
              x=0.1, hjust = 0, size=2)+
    
    theme_classic()+
    theme(legend.position = "none")+
    facet_grid(pop ~ phage + seed.bank)+
    panel_border(color = "black")+
    ylim(0,1)+
    ylab(expression("Allele frequency,"~italic("f(t)")))+
    xlab(expression("Transfer,"~italic("t")))
    # ggtitle("B. subtilis mutation frequencies")
s = 0.6
ggsave(here("analysis/host_mutation_trajectories2.png"),p, height = s* 6, width = s* 8, units = "in")
  

# # distribution ------------------------------------------------------------
# 
# d.plot %>% 
#   filter(!Annotation %in% c("noncoding", "unknown")) %>% 
#   mutate(trt = str_extract(f, "^.*_L."),
#          seed.bank = if_else(str_detect(trt,"WL."), 
#                              "with seed bank","without seed bank"),
#          phage=if_else(str_detect(trt,"O"), 
#                        "with phage","without phage"),
#          pop = parse_number(trt) %>% as.character()) %>% 
#   select(mut_id,trt, seed.bank,phage,pop, contains("Freq")) %>%
#   pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>%
#   
#   filter(freq > 0) %>% 
#   
#   mutate(t_sample = parse_number(t_sample)) %>% 
#   arrange(desc(freq)) %>% 
#   group_by(t_sample, trt, group=pop) %>% 
#   mutate(rank = row_number()) %>% 
#   
#   ggplot(aes(rank, freq, group = interaction(seed.bank,pop)))+
#   geom_line(aes( color = seed.bank), size = 0.8)+
#   facet_grid(phage ~ t_sample)+
#   theme_classic(base_size = 18)+
#   panel_border(color = "black", size = 1.5)+
#   scale_y_log10(labels = trans_format("log10", math_format(10^.x)), 
#                 # limits = c(1e-5, NA)
#   )+
#   scale_color_grey(start = 0, end = 0.7)+
#   annotation_logticks(sides = "l")+
#   ylab("frequency")+
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold"))
# 
