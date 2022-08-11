library(here)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggnewscale)


# Phage annotation --------------------------------------------------------
d.spo1 <- 
  read_delim(here("data/SPO1-ANC_no_fasta.gff3"),
             col_names = c(letters[1:9]),
             trim_ws = T,
             comment = "#") %>% 
  select(seqname =a, feature =c, start = d, end = e,
         strand = g, frame = h, attribute =i)

d.spo1_genes <- 
  d.spo1 %>% 
  filter(feature == "gene") %>% 
  separate(attribute, into = c("Alias", "ID", "Name", "Pseudo"), sep = ";") %>% 
  mutate(Alias = str_remove(Alias, ".*="),
         ID = str_remove(ID, ".*="),
         Name = str_remove(Name, ".*="),
         Pseudo = str_remove(Pseudo, ".*="))
  


# Find hi-freq loci --------------------------------------------------------

# bottom threshold for mutation to be high frequency
freq_cutoff <- 0.4

f_freq <- list.files(here("data/timecourse_final_breseq/"))

f_freq <- f_freq[str_detect(f_freq, "phage")]
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
    filter (n > 2) %>% 
    pull(mut_id)
  
  # Mutations (nearly) fixed
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
  
  if(length(to_pull) < 1) next
  
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
hi_freq <- 
  d.spo1_genes %>% 
  select(Gene = ID, Name) %>% 
  left_join(hi_freq,., by = "Gene") %>% 
  mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) 
  


# Plot --------------------------------------------------------------------

# d <- read_csv((here("data/timecourse_final_breseq/SNO_L1_host_no_seed_bank_SPO1_revived_total_annotated_timecourse.txt")))

d <- d.plot %>%
  # parse treatments
  mutate(trt = str_extract(f, "^..._L."),
         seed.bank = if_else(str_detect(trt,"WL."), 
                             "with seed bank","without seed bank"),
         pop = parse_number(trt)) %>% 
  #add T = 0  
  mutate(`Freq:0` = 0) %>%
  select(mut_id,trt, seed.bank,pop, contains("Freq")) %>%
    # slice_head(n = 100) %>%
    pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>%
    mutate(t_sample = parse_number(t_sample))


# hi freq trajectories
d.hi <- 
  d %>% 
  filter(mut_id %in% hi_freq$mut_id) %>% 
  left_join(., select(hi_freq, mut_id, Name), by = "mut_id") %>% 
  mutate(Name = if_else(str_detect(mut_id,"intergenic"), "intergenic", Name)) %>%
  # filter(!str_detect(mut_id,"intergenic")) %>% 
  #adjust legend order
  arrange(Name) %>% 
  mutate(Name = fct_inorder(Name))

# hi freq names
lab.hi <- d.hi %>% 
  group_by(seed.bank, pop) %>% 
  summarise(Name) %>% 
  distinct() %>% 
  mutate(y.idx = row_number(),
         y.pos = 0.98- 0.12*(y.idx-1)) 

# plot

  p <- d %>%
    ggplot(aes(t_sample, freq))+
    geom_line(aes(color = mut_id), size=.8, show.legend = F)+
    scale_color_grey(guide = "none")+
    
    new_scale_color() +
    geom_line(data = d.hi,
              aes(group = mut_id, color = Name), size=.8, show.legend = T)+
    geom_point(data = d.hi,
              aes(group = mut_id, color = Name), shape=21, fill="white",
              size=.8, show.legend = T)+
    scale_color_viridis_d(guide = guide_legend(title = "SPO1 gene"))+
    
    geom_text(data = lab.hi, 
              aes(label = Name, color = Name, y=  y.pos),
              x=-0.5, hjust = 0)+
    
    theme_classic()+
    theme(legend.position = "none")+
    facet_grid(pop ~ seed.bank)+
    panel_border(color = "black")+
    ylim(0,1)+
    ylab(expression("Allele frequency,"~italic("f(t)")))+
    xlab(expression("transfer,"~italic("t")))
    # ggtitle("Phage SPO1 mutation frequencies")

ggsave(here("analysis/phage_mutation_trajectories.png"),p, height = 5.52, width = 4.25, units = "in")
  

