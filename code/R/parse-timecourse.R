# parse mutation time course

setwd("/N/slate/danschw/coevo-seedbank-seq/")
library(here)
library(tidyverse)

pop <- "SNO-L1"

d.raw <-  read_csv(here("data/timecourses/phage/timecourse_merged/SNO-L2_timecourse.txt"),
                    col_names = c("chromosome", "position", "alt_allele", "times", "n.allele", "n.depth"))

time.key <- strsplit(d.raw$times[1], " ", fixed = T) %>%
  unlist() 

d <- d.raw %>% 
  mutate(mut = paste0(as.character(position), alt_allele)) %>% 
  separate(col=n.allele, into = paste0("At-",time.key), sep = " ") %>% 
  separate(col=n.depth, into = paste0("Dt-",time.key), sep = " ") %>% 
  # select()
  pivot_longer(cols = contains("At")|contains("Dt"),
               names_to = c(".value", "transfer"),
               names_sep = "-") %>% 
  mutate(At = as.integer(At), Dt = as.integer(Dt),
         transfer = as.integer(transfer)) %>% 
  mutate(ft =  At/Dt)


d %>% 
  ggplot(aes(transfer, ft, color = mut))+
  geom_line(show.legend = F)+
  theme_classic()+
  scale_colour_viridis_d()
  

