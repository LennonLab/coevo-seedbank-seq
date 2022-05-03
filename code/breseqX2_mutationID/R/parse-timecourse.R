# parse mutation time course

setwd("/N/slate/danschw/coevo-seedbank-seq/")
library(here)
library(tidyverse)

# extract bz file
xtr.bz <- function(file.name){
  new.file <- gsub(".bz", ".txt", file.name)
  if (file.exists(new.file))  file.remove(new.file)
  system(paste("bzip2 -dkc", file.name, ">", new.file))
  return(new.file)
}


#-------------------#
# files to read
time_course_files <- list.files(here("data/timecourses/phage/timecourse_merged2"),
                                pattern = ".bz",full.names = TRUE, recursive = TRUE)

d.raw <- tibble()

for (f in time_course_files){
  
  #get pop name
  pop_name <- gsub(".*merged2/","", f) %>% gsub("_merged.*","", .)
  
  #extract and read file
  txt.file <- xtr.bz(f)
  
  d.raw <-  read_csv(txt.file, 
                     col_names = c("chromosome", "position", "alt_allele", "times", "n.allele", "n.depth")) %>% 
    bind_cols(pop = pop_name, .) %>% 
    bind_rows(d.raw,.)
  if (file.exists(txt.file))  file.remove(txt.file)
  
}



# first time course
# d.raw <-  read_csv(paste0(here("data/timecourses/phage/timecourse_merged"),"/",
#                           pop, "_timecourse.txt"),
#                     col_names = c("chromosome", "position", "alt_allele", "times", "n.allele", "n.depth"))
# second time course

time.key <- strsplit(d.raw$times[1], " ", fixed = T) %>%
  unlist() 

d <- d.raw %>% 
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
  mutate(ft =  At/Dt) %>% 
  #ajust factor order for seed bank
  mutate (seed.bank = case_when(trt == "WLO" ~ "long seed bank",
                                trt ==  "WSO" ~ "short seed bank",
                                trt ==  "SNO" ~ "no seed bank")) %>%
  mutate(seed.bank=as.factor(seed.bank))%>%
  mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank"))

# add freq = 0 at T = 0
d0 <- d %>% 
  group_by(seed.bank, rep, mut) %>% 
  summarise()%>% 
  mutate(transfer = 0, ft = 0)



d %>% 
  bind_rows(d, d0) %>% 
  ggplot(aes(transfer, ft, color = mut))+
  geom_line(show.legend = F)+
  theme_classic()+
  scale_colour_viridis_d()+
  facet_wrap(seed.bank ~ rep)+
  scale_x_continuous(breaks = c(1,4,7,10,14))+
  ylab("mutation freq (A/D)") +
  ggsave(filename = here("plots","merged_mut_freq.png"), 
        width = 10, height = 7, units = "in")
  

d %>% 
  bind_rows(d, d0) %>% 
  filter(position > 47500 & position < 74316) %>% 
  ggplot(aes(transfer, ft, color = mut))+
  geom_line(show.legend = F)+
  theme_classic()+
  scale_colour_viridis_d()+
  facet_wrap(seed.bank ~ rep)+
  scale_x_continuous(breaks = c(1,4,7,10,14))+
  ylab("mutation freq (A/D)") +
  ggsave(filename = here("plots","merged_mut_freq_TAIL_REGION.png"), 
         width = 10, height = 7, units = "in")


d %>% 
  mutate(transfer = paste0("T", transfer)) %>% 
  ggplot(aes(Dt, At))+
  geom_point(aes(color=transfer), shape=21, alpha = 0.5)+
  theme_classic()+
  scale_colour_viridis_d()+
  facet_wrap(seed.bank ~ rep)+
  ylab("Allele (#reads)") +
  xlab("Depth (#reads)") + 
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks()+
  ggsave(filename = here("plots","phage_A-D.png"), 
         width = 10, height = 7, units = "in")
