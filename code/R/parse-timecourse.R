# parse mutation time course

setwd("/N/slate/danschw/coevo-seedbank-seq/")
library(here)
library(tidyverse)

xtr.bz <- function(file.name){
  new.file <- gsub(".bz", ".txt", file.name)
  if (file.exists(new.file))  file.remove(new.file)
  system(paste("bzip2 -dkc", file.name, ">", new.file))
  return(new.file)
}

pop <- "WSO-L3"

# first time course
# d.raw <-  read_csv(paste0(here("data/timecourses/phage/timecourse_merged"),"/",
#                           pop, "_timecourse.txt"),
#                     col_names = c("chromosome", "position", "alt_allele", "times", "n.allele", "n.depth"))
# second time course
txt.file <- 
  paste0(here("data/timecourses/phage/timecourse_merged2"),"/",
                   pop, "_merged_timecourse.bz") %>%
  xtr.bz() 
d.raw <-  read_csv(txt.file, col_names = c("chromosome", "position", "alt_allele", "times", "n.allele", "n.depth"))
if (file.exists(txt.file))  file.remove(txt.file)

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

# add freq = 0 at T = 0
d0 <- 
  d.raw %>% 
  select(position, alt_allele ) %>% 
  mutate(mut = paste0(as.character(position), alt_allele)) %>% 
  mutate(transfer = 0, ft = 0)

d %>% 
  bind_rows(d, d0) %>% 
  ggplot(aes(transfer, ft, color = mut))+
  geom_line(show.legend = F)+
  theme_classic()+
  scale_colour_viridis_d()+
  ggtitle(pop)
  

