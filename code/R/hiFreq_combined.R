library(here)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggnewscale)
library(ggh4x)


# load data ---------------------------------------------------------------

# host, from "hiFreq_host.R"
load(here("data","host_trajectories.Rdata"))
d_host <- d %>% mutate(who = "host"); rm(d)
d.hi_host <- d.hi %>% mutate(who = "host"); rm(d.hi)
lab.hi_host <- lab.hi %>% mutate(who = "host"); rm(lab.hi)

# phage, from "hiFreq_phage.R"
load(here("data","phage_trajectories.Rdata"))
d_phage <- d %>% mutate(phage = "with phage", who = "phage"); rm(d)
d.hi_phage <- d.hi %>% mutate(phage = "with phage", who = "phage"); rm(d.hi)
lab.hi_phage <- lab.hi %>% mutate(phage = "with phage", who = "phage"); rm(lab.hi)



# plot infected populations without seed bank ------------------------------

# all mutations
d <- bind_rows(d_host,d_phage) %>% 
  filter(phage == "with phage") %>% 
  filter(seed.bank == "without seed bank")

# high frequency mutations
d.hi <- bind_rows(d.hi_host,d.hi_phage)%>%
  ungroup() %>% 
  filter(phage == "with phage") %>% 
  # needed to make gene colors conistent across plots
  mutate(Name = droplevels(Name) %>% fct_inseq()) %>% 
  filter(seed.bank == "without seed bank")

# high frequency gene labels
lab.hi <- bind_rows(lab.hi_host,lab.hi_phage)%>% 
  ungroup() %>% 
  filter(phage == "with phage") %>% 
  # needed to make gene colors conistent across plots
  mutate(Name = droplevels(Name) %>% fct_inseq()) %>% 
  filter(seed.bank == "without seed bank")

# plot
p_noSB <- d %>%
  ggplot(aes(t_sample, freq))+
  
  #all mutations
  geom_line(aes(color = mut_id), size=.8, show.legend = F)+
  scale_color_grey(guide = "none")+
  
  # high frequency mutations
  new_scale_color() +
  geom_line(data = d.hi,
            aes(group = mut_id,color = Name), size=.8,alpha = 0.7 ,show.legend = T)+
  geom_point(data = d.hi,
             aes(group = mut_id, color = Name), shape=21, fill="white",
             alpha = 0.7 ,size=.8, show.legend = T)+
  scale_color_viridis_d(drop=F)+
  
  # high frequency gene labels
  geom_text(data = lab.hi, 
            aes(label = Name, color = Name, y=  y.pos),
            x=0.1, hjust = 0, size = 3)+
  
  # general design
  theme_classic()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        strip.background.x = element_blank())+
  facet_nested(pop ~ seed.bank + who, nest_line = element_line())+
  panel_border(color = "black")+
  ylim(0,1)+
  ylab(expression("Allele frequency,"~italic("f(t)")))+
  xlab(expression("transfer,"~italic("t")))

# plot infected populations with a seed bank ------------------------------

# all mutations
d <- bind_rows(d_host,d_phage) %>% 
  filter(phage == "with phage") %>% 
  filter(seed.bank == "with seed bank")

# high frequency mutations
d.hi <- bind_rows(d.hi_host,d.hi_phage)%>% 
  ungroup() %>% 
  filter(phage == "with phage") %>% 
  # needed to make gene colors conistent across plots
  mutate(Name = droplevels(Name) %>% fct_inseq()) %>% 
  filter(seed.bank == "with seed bank")

# high frequency gene labels
lab.hi <- bind_rows(lab.hi_host,lab.hi_phage)%>% 
  ungroup() %>% 
  filter(phage == "with phage") %>% 
  # needed to make gene colors conistent across plots
  mutate(Name = droplevels(Name) %>% fct_inseq()) %>% 
  filter(seed.bank == "with seed bank")


# plot
p_wSB <- d %>%
  ggplot(aes(t_sample, freq))+
  #all mutations
  geom_line(aes(color = mut_id), size=.8, show.legend = F)+
  scale_color_grey(guide = "none")+
  
  # high frequency mutations
  new_scale_color() +
  geom_line(data = d.hi,
            aes(group = mut_id,color = Name), size=.8,alpha = 0.7 ,show.legend = T)+
  geom_point(data = d.hi,
             aes(group = mut_id, color = Name), shape=21, fill="white",
             alpha = 0.7 ,size=.8, show.legend = T)+
  scale_color_viridis_d(drop=F)+
  
  # high frequency gene labels
  geom_text(data = lab.hi, 
            aes(label = Name, color = Name, y=  y.pos),
            x=0.1, hjust = 0, size = 3)+
  
  # general design  
  theme_classic()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        strip.background.x = element_blank())+
  facet_nested(pop ~ seed.bank + who, nest_line = element_line())+
  panel_border(color = "black")+
  ylim(0,1)+
  ylab(expression("Allele frequency,"~italic("f(t)")))+
  xlab(expression("transfer,"~italic("t")))


# combine plots and save --------------------------------------------------

ggsave(here("analysis/main_mutation_trajectories.png"),
       plot_grid(p_noSB, p_wSB, labels = c("(a)","(b)")),
       height = 5, width =8, units = "in")




# no phage host plot - SUPPL ------------------------------------------------

# all mutations
d <- bind_rows(d_host,d_phage) %>%
  filter(phage == "without phage")

# high frequency mutations
d.hi <- bind_rows(d.hi_host,d.hi_phage)%>%
  ungroup() %>%
  filter(phage == "without phage") %>%
  mutate(Name = droplevels(Name))

# high frequency gene labels
lab.hi <- bind_rows(lab.hi_host,lab.hi_phage)%>%
  ungroup() %>%
  filter(phage == "without phage")%>%
  mutate(Name = droplevels(Name))



p_sup <- d %>%
  ggplot(aes(t_sample, freq))+
  #all mutations
  geom_line(aes(color = mut_id), size=.8, show.legend = F)+
  scale_color_grey(guide = "none")+
  
  # high frequency mutations
  new_scale_color() +
  geom_line(data = d.hi,
            aes(group = mut_id,color = Name), size=.8,alpha = 0.7 ,show.legend = T)+
  geom_point(data = d.hi,
             aes(group = mut_id, color = Name), shape=21, fill="white",
             alpha = 0.7 ,size=.8, show.legend = T)+
  scale_color_viridis_d(drop=F)+
  
  # high frequency gene labels
  geom_text(data = lab.hi, 
            aes(label = Name, color = Name, y=  y.pos),
            x=0.1, hjust = 0, size = 3)+
  
  # general design 
  theme_classic()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        strip.background.x = element_blank())+
  facet_nested(pop ~ phage + fct_rev(seed.bank), nest_line = element_line())+
  panel_border(color = "black")+
  ylim(0,1)+
  ylab(expression("Allele frequency,"~italic("f(t)")))+
  xlab(expression("transfer,"~italic("t")))


ggsave(here("analysis/noPhage_mutation_trajectories.png"),p_sup,
       height = 5, width =5, units = "in")


# # Single plot -------------------------------------------------------------
# 
# 
# 
# d <- bind_rows(d_host,d_phage) %>% 
#   filter(phage == "with phage")
# 
# d.hi <- bind_rows(d.hi_host,d.hi_phage)%>%
#   ungroup() %>% 
#   filter(phage == "with phage") %>%   
#   mutate(Name = droplevels(Name) %>% fct_inseq())
# 
# lab.hi <- bind_rows(lab.hi_host,lab.hi_phage)%>% 
#   ungroup() %>% 
#   filter(phage == "with phage")%>%  
#   mutate(Name = droplevels(Name) %>% fct_inseq())
# 
# 
# 
# p <- d %>%
#   ggplot(aes(t_sample, freq))+
#   geom_line(aes(color = mut_id), size=.8, show.legend = F)+
#   scale_color_grey(guide = "none")+
#   
#   new_scale_color() +
#   geom_line(data = d.hi,
#             aes(group = mut_id,color = Name), size=.8,alpha = 0.7 ,show.legend = T)+
#   geom_point(data = d.hi,
#              aes(group = mut_id, color = Name), shape=21, fill="white",
#              alpha = 0.7 ,size=.8, show.legend = T)+
#   scale_color_viridis_d(drop=F)+
#   
#   geom_text(data = lab.hi, 
#             aes(label = Name, color = Name, y=  y.pos),
#             x=0.1, hjust = 0)+
#   
#   theme_classic()+
#   theme(legend.position = "none",
#         strip.text = element_text(face = "bold"),
#         strip.background.x = element_blank())+
#   facet_nested(pop ~ fct_rev(seed.bank) + who, nest_line = element_line())+
#   panel_border(color = "black")+
#   ylim(0,1)+
#   ylab(expression("Allele frequency,"~italic("f(t)")))+
#   xlab(expression("transfer,"~italic("t")))
# 
# ggsave(here("analysis/main2_mutation_trajectories.png"),p,
#        height = 5, width =8, units = "in")
