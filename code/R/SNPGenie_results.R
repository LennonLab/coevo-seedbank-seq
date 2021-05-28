
setwd("/N/slate/danschw/coevo-seedbank-seq/")
library(here)
library(tidyverse)
library(cowplot)


host_pops <- list.dirs(here("data/SnpGenie/host"), full.names = T)[-1]
phage_pops <- list.dirs(here("data/SnpGenie/phage"), full.names = T)[-1]


d <- tibble() 

for (f in host_pops){
  cur <-  str_remove(f, ".*SnpGenie/host/")
    d <- 
      read_tsv(paste0(f,"/population_summary.txt"),
               col_types = cols(
                 .default = col_double(),
                 file = col_character(),
                 mean_S_gdiv = col_character()
               ) ) %>%
      mutate(pop = cur) %>%
      bind_rows(d, .)
}

d <- d%>%
  separate(pop, into = c("trt", "Line", "t.exp", "subpop")) %>%
  filter(subpop != "pl")
  

# d %>% 
#   select(where(is.numeric)) %>% 
#   as.matrix() %>% pairs()

#overview
d %>%
  select(-mean_S_gdiv)%>%
  pivot_longer(cols = is_numeric,values_to = "value", names_to = "param")%>%
  mutate(transfer = parse_number(t.exp)) %>% 
  ggplot(aes(x=subpop, y = value))+
  geom_point(aes(color = trt, shape = Line), alpha = 0.5)+
  facet_wrap(~param, scales = "free_y")+
  theme_classic()
  # ggsave(here("snpGenie_pop_params.png"), width = 8, height = 5)

d %>%
  select(trt, Line, transfer=t.exp ,subpop, pi, piN, piS) %>% 
  pivot_longer(cols=c("pi", "piN", "piS")) %>% 
  group_by(trt, transfer, subpop, name) %>% 
  summarise(m=mean(value), v=sd(value))  %>% 
  ggplot(aes(x=trt, y = m))+
  geom_line()+
  geom_pointrange( aes(ymin = m-v, ymax = m+v), fill="white")+
  theme_classic(base_size = 16)+
  # scale_y_log10()+
  facet_grid(subpop~name)+
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", size = 1.5, fill = NA))+
  labs(caption = " allele freq. > 1%")+
  coord_flip()


  d %>%
    mutate(transfer = parse_number(t.exp)) %>% 
    select(trt, Line, transfer ,trt, pi, piN, piS) %>% 
    pivot_longer(cols=c("pi", "piN", "piS")) %>% 
  group_by(trt,transfer, name) %>% 
  summarise(m = mean(value), v = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x=transfer, y = m, color = trt))+
  geom_line()+
  geom_pointrange(aes(ymin = m-v, ymax = m+v), shape = 21, fill="white")+
  theme_classic(base_size = 16)+
  # scale_y_log10()+
    facet_wrap(~name)+
    scale_x_continuous(breaks = c(1,4,7,10,14))+
    theme(legend.position = "bottom",
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.border = element_rect(colour = "black", size = 1.5, fill = NA))+
    labs(caption = " allele freq. > 1%\t mean±sem (n=3)")+
  ggsave(here("pi.png"), width = 6, height = 4)
  
  for (i in unique(d$t.exp)){
    print(i)
    pairwise.t.test(log10(d$pi[d$t.exp==i]), d$trt[d$t.exp==i], p.adjust.method="none") %>% print()
  }
       
  #####################################
  #time as muneric
 d <- d %>% 
    mutate(transfer  = parse_number(t.exp)) %>% 
    #ajust factor order for seed bank
    mutate (seed.bank = case_when(trt == "WLO" ~ "long seed bank",
                                  trt ==  "WSO" ~ "short seed bank",
                                  trt ==  "SNO" ~ "no seed bank")) %>%
    mutate(seed.bank=as.factor(seed.bank))%>%
    mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank"))
  
  d.sum <- d %>%
    mutate(transfer = parse_number(t.exp)) %>% 
    select(seed.bank, Line, transfer ,trt, pi, piN, piS) %>% 
    pivot_longer(cols=c("pi", "piN", "piS")) %>% 
    group_by(seed.bank,transfer, name) %>% 
    summarise(m = mean(value), v = sd(value)/sqrt(n())) 
  
  p1 <- d %>% 
    ggplot(aes(seed.bank, pi, color = seed.bank))+
    geom_crossbar(data = filter(d.sum,name == "pi"), aes(y = m, ymin = m-v, ymax = m+v))+
    geom_point(aes(fill = seed.bank), color="black", shape=21, size=3)+
    facet_wrap(~transfer, nrow = 1)+
    theme_cowplot()+
    panel_border(size = 1.5, color = "black")+
    theme(strip.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "white"),
          axis.text.x = element_blank())+
    ylab("Nucleotide Diversity (pi)")+
    labs(caption = "allele freq. > 1%\tbox depicts mean±sem")+
    ggsave(here("pi.png"), width = 7, height = 4)
  
  p2 <- d %>% 
    ggplot(aes(seed.bank, piN, color = seed.bank))+
    geom_crossbar(data = filter(d.sum,name == "piN"), aes(y = m, ymin = m-v, ymax = m+v))+
    geom_point(aes(fill = seed.bank), color="black", shape=21, size=3)+
    facet_wrap(~transfer, nrow = 1)+
    theme_cowplot()+
    panel_border(size = 1.5, color = "black")+
    theme(strip.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "white"),
          axis.text.x = element_blank())+
    ylab("Nucleotide Diversity (piN)")+
    labs(caption = "allele freq. > 1%\tbox depicts mean±sem")+
  ggsave(here("piN.png"), width = 7, height = 4)
  
 p3 <-  d %>% 
    ggplot(aes(seed.bank, piS, color = seed.bank))+
    geom_crossbar(data = filter(d.sum,name == "piS"), aes(y = m, ymin = m-v, ymax = m+v))+
    geom_point(aes(fill = seed.bank), color="black", shape=21, size=3)+
    facet_wrap(~transfer, nrow = 1)+
    theme_cowplot()+
    panel_border(size = 1.5, color = "black")+
    theme(strip.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "white"),
          axis.text.x = element_blank())+
    ylab("Nucleotide Diversity (piS)")+
    labs(caption = "allele freq. > 1%\tbox depicts mean±sem")+
  ggsave(here("piS.png"), width = 6, height = 4)

 plot_grid(p1,p2,p3, ncol=1)+
   ggsave(here("pi-all.png"), width = 7, height = 9)
 
 #####################################

d.product <- tibble() 

for (f in pops){
  cur <-  str_remove(f, ".*g_SnpGenie//")
  d.product <- 
    read_tsv(paste0(f,"/product_results.txt")) %>%
    mutate(pop = cur) %>%
    bind_rows(d.product, .)
}

d.product <- d.product%>%
  separate(pop, into = c("trt", "Line", "t.exp"))%>%
  mutate(transfer = parse_number(t.exp)) %>% 
  mutate(gene =str_remove(product, "SPO1_") %>% as.integer())%>%
  mutate(ns = piN/piS)

d.product%>%
  filter(!is.infinite(ns)) %>% 
  filter(!is.na(ns)) %>% 
  ggplot(aes(x=gene, y=ns))+
  geom_point(aes(color=trt))+
  facet_grid(transfer~.)+
  theme_classic(base_size = 16)
