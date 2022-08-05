library(here)
library(tidyverse)
library(cowplot)


# d <- read_csv((here("data/timecourse_final/SNO_L1_host_no_seed_bank_SPO1_revived_total_annotated_timecourse.txt")))

# d <- d %>%
#   rowwise() %>%
#   mutate(f_max=max(c_across(contains("Freq"))))
# 
# filter(d,f_max> 0.5) %>% view
#   pull(f_max) %>% 
#   hist(.)
# (d$f_max>0.4)

  

# plot freq ---------------------------------------------------------------
d <- read_csv((here("data/timecourse_final/SNO_L1_host_no_seed_bank_SPO1_revived_total_annotated_timecourse.txt")))

p1 <- d %>% 
    #make uniqe id
    mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
    # pull(mut_id) %>%  anyDuplicated(.) # it is sufficient
    mutate(`Freq:0` = 0) %>% 
    select(mut_id, contains("Freq")) %>% 
    # slice_head(n = 100) %>% 
    pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>% 
    mutate(t_sample = parse_number(t_sample)) %>% 
    
    
    ggplot(aes(t_sample, freq))+
    geom_line(aes(color = mut_id), size=.8, show.legend = F)+
    scale_color_viridis_d()+
    theme_classic()+
    panel_border(color = "black")+
  ylim(0,1)+
  labs(title = "Replicate 1",
       caption = "File: SNO_L1_host_no_seed_bank_SPO1\nrevived_total_annotated_timecourse.txt")

###
d <- read_csv((here("data/timecourse_final/SNO_L2_host_no_seed_bank_SPO1_revived_total_annotated_timecourse.txt")))

p2 <- d %>% 
  #make uniqe id
  mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
  # pull(mut_id) %>%  anyDuplicated(.) # it is sufficient
  mutate(`Freq:0` = 0) %>% 
  select(mut_id, contains("Freq")) %>% 
  # slice_head(n = 100) %>% 
  pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>% 
  mutate(t_sample = parse_number(t_sample)) %>% 
  
  
  ggplot(aes(t_sample, freq))+
  geom_line(aes(color = mut_id), size=.8, show.legend = F)+
  scale_color_viridis_d()+
  theme_classic()+
  panel_border(color = "black")+
  ylim(0,1)+
  labs(title = "Replicate 2",
       caption = "File: SNO_L2_host_no_seed_bank_SPO1\nrevived_total_annotated_timecourse.txt")
###
d <- read_csv((here("data/timecourse_final/SNO_L3_host_no_seed_bank_SPO1_revived_total_annotated_timecourse.txt")))

p3 <- d %>% 
  #make uniqe id
  mutate(mut_id = str_c(Position, Gene, Allele, sep = "_")) %>% 
  # pull(mut_id) %>%  anyDuplicated(.) # it is sufficient
  mutate(`Freq:0` = 0) %>% 
  select(mut_id, contains("Freq")) %>% 
  # slice_head(n = 100) %>% 
  pivot_longer(contains("Freq"), names_to = "t_sample", values_to = "freq") %>% 
  mutate(t_sample = parse_number(t_sample)) %>% 
  
  
  ggplot(aes(t_sample, freq))+
  geom_line(aes(color = mut_id), size=.8, show.legend = F)+
  scale_color_viridis_d()+
  theme_classic()+
  panel_border(color = "black")+
  ylim(0,1)+
  labs(title = "Replicate 3",
       caption = "File: SNO_L3_host_no_seed_bank_SPO1\nrevived_total_annotated_timecourse.txt")
###

ggsave(here("analysis/SNO.png"),plot_grid(p1,p2,p3, nrow = 1))  