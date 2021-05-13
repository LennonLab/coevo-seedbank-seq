# parse mutation time course

setwd("/N/slate/danschw/coevo-seedbank-seq/")
library(here)
library(tidyverse)
library(cowplot)

# extract bz file
xtr.bz <- function(file.name){
  new.file <- gsub(".bz", ".txt", file.name)
  if (file.exists(new.file))  file.remove(new.file)
  system(paste("bzip2 -dkc", file.name, ">", new.file))
  return(new.file)
}


#-------------------#
# files to read
time_course_files <- list.files(here("data/CompSubPops/host/CompSubPop_merged2"),
                                pattern = ".bz",full.names = TRUE, recursive = TRUE)
# time_course_files <- list.files(here("data/CompSubPops/host/"),
                                # pattern = ".txt",full.names = TRUE, recursive = TRUE)

d.raw <- tibble()

for (f in time_course_files){
  
  #get pop name
  pop_name <- gsub(".*_merged2/","", f) %>% gsub("_merged_Comp.*","", .)
  # pop_name <- "WLCt-L1"
  
  #extract and read file
  txt.file <- xtr.bz(f)
  
  d.raw <-  read_csv(txt.file, 
                     col_names = c("chromosome", "position", "alt_allele", "times", "n.allele", "n.depth")) %>% 
    bind_cols(pop = pop_name, .) %>% 
    bind_rows(d.raw,.)
}



########################
#Do separetly for WT and mut
########################
time.key <- strsplit(d.raw$times[nchar(d.raw$times) %>% which.max()],
                     " ", fixed = T) %>%
  unlist() 

d.wt <- d.raw %>% 
  filter(str_starts(pop, "W")) %>% 
  mutate(mut = paste0(as.character(position), alt_allele)) %>% 
  separate(pop, into = c("trt", "rep")) %>% 
  separate(col=n.allele, into = paste0("At-",time.key), sep = " ") %>% 
  separate(col=n.depth, into = paste0("Dt-",time.key), sep = " ") %>% 
  # select()
  pivot_longer(cols = contains("At")|contains("Dt"),
               names_to = c(".value", "transfer"),
               names_sep = "-") %>% 
  mutate(At = as.integer(At), Dt = as.integer(Dt),
         transfer = as.integer(transfer)) %>% 
  mutate(ft =  At/Dt)

time.key <- strsplit(d.raw$times[nchar(d.raw$times) %>% which.min()],
                     " ", fixed = T) %>%
  unlist() 

d.mut <- d.raw %>% 
  filter(str_starts(pop, "S")) %>% 
  mutate(mut = paste0(as.character(position), alt_allele)) %>% 
  separate(pop, into = c("trt", "rep")) %>% 
  separate(col=n.allele, into = paste0("At-",time.key), sep = " ") %>% 
  separate(col=n.depth, into = paste0("Dt-",time.key), sep = " ") %>% 
  # select()
  pivot_longer(cols = contains("At")|contains("Dt"),
               names_to = c(".value", "transfer"),
               names_sep = "-") %>% 
  mutate(At = as.integer(At), Dt = as.integer(Dt),
         transfer = as.integer(transfer)) %>% 
  mutate(ft =  At/Dt)  



d <- bind_rows(d.wt, d.mut) %>% 
  mutate(subpop = case_when(transfer == 100 ~ "pellet",
                            transfer == 200 ~ "revived_spore",
                            transfer == 300 ~ "revived_total")) %>% 
  #ajust factor order for seed bank
  mutate (seed.bank = case_when(str_starts(trt, "WL") ~ "long seed bank",
                                str_starts(trt, "WS") ~ "short seed bank",
                                str_starts(trt, "SN") ~ "no seed bank")) %>%
  mutate(seed.bank=as.factor(seed.bank))%>%
  mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank"))

d1 <- d %>% 
  filter(Dt > 50) %>%
  select(trt, rep, mut, subpop, ft) %>% 
  pivot_wider(names_from = "subpop", values_from = ft)

d1 %>% 
  ggplot(aes(revived_spore, revived_total))+
  geom_abline(slope = 1, intercept = 0, color = "grey")+
  geom_point(shape = 21)+
  theme_classic()+
  ylim(0,1)+
  xlim(0,1)+
  facet_grid(trt~rep)+
  panel_border()

d1 %>% 
  ggplot(aes(pellet, revived_total))+
  geom_abline(slope = 1, intercept = 0, color = "grey")+
  geom_point(shape = 21)+
  theme_classic()+
  ylim(0,1)+
  xlim(0,1)+
  facet_grid(trt~rep)+
  panel_border()

d1 %>% 
  ggplot(aes(pellet, revived_spore))+
  geom_abline(slope = 1, intercept = 0, color = "grey")+
  geom_point(shape = 21)+
  theme_classic()+
  ylim(0,1)+
  xlim(0,1)+
  facet_grid(trt~rep)+
  panel_border()
  