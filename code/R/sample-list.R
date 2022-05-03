
# purpose -----------------------------------------------------------------


# to make a master list of data for facilitating the downstream analysis.
# I will make a master list of all the samples, trier attributes, 
# and the corresponding data


# dependencies ----------------------------------------------------------------
library(here, quietly = T, verbose = F)
library(tidyverse, quietly = T, verbose = F)



# List data files -------------------------------------------------------------

# > raw sequencing output -----------------------------------------------------
# (= analysis input) 
data_dir <- "data/input/fastq"

# phage
phage.raw <-
  tibble(phageORhost = "phage",
         path = list.files(here(data_dir, "phage"), 
                           pattern = "fastq", full.names = TRUE) ,
         fq = list.files(here(data_dir, "phage"), 
                        pattern = "fastq", full.names = FALSE) )%>%
  separate(fq, into=c("seq_run","sample_num","trt", "line", "transfer", "NexteraBC", "read", "sufx"),
           sep=("-|_"), remove = F) %>% 
  # subpop for compatability with host
  mutate(subpop = "filtered_phage") %>% 
  select(-sample_num) %>% 
  # sample name
  mutate(unq.sample = str_c(trt, line, transfer, sep = '-'))

# host
host.raw <-
  tibble(phageORhost = "host",
         path = list.files(here(data_dir, "host"), 
                           pattern = "fastq", full.names = TRUE) ,
         fq = list.files(here(data_dir, "host"), 
                        pattern = "fastq", full.names = FALSE) )%>%
  separate(fq,into=c("seq_run","trt", "line", "transfer", "subpop", "NexteraBC", "read", "sufx"),
           sep=("-|_"), remove = F) %>%   
  #fix extraction names to unify across sequencing runs
  mutate(subpop = case_when(
    subpop %in% c("SB","rS") ~ "revived_spore",
    subpop %in% c("veg","rV") ~ "revived_total",
    subpop %in% c("pl") ~ "frozen_pellet")) %>% 
  # sample name
  mutate(unq.sample = str_c(trt, line, transfer, sep = '-'))

d.raw <- bind_rows(host.raw, phage.raw)

d.raw <-  bind_rows(host.raw, phage.raw) %>% 
  select("unq.sample", "phageORhost", "trt", "line", "transfer", "subpop",
         "seq_run" , "read", "fq") %>% 
  mutate(fq = str_c(data_dir ,phageORhost ,fq, sep = "/"),
         read = str_c("raw_", read)) %>% 
  pivot_wider(names_from = read, values_from = fq)


# de-duplicated  reads ----------------------------------------------------

data_dir <- "data/ddup-fastq"
# phage
phage.ddup <-
  tibble(phageORhost = "phage",
         path = list.files(here(data_dir, "phage"), 
                           pattern = "fastq", full.names = TRUE) ,
         fq = list.files(here(data_dir, "phage"), 
                         pattern = "fastq", full.names = FALSE) )%>%
  separate(fq, into=c("stage","seq_run","sample_num","trt", "line", "transfer", "NexteraBC", "read", "sufx"),
           sep=("-|_"), remove = F) %>% 
  # subpop for compatability with host
  mutate(subpop = "filtered_phage") %>% 
  select(-sample_num) %>% 
  # sample name
  mutate(unq.sample = str_c(trt, line, transfer, sep = '-'))

# host
host.ddup <-
  tibble(phageORhost = "host",
         path = list.files(here(data_dir, "host"), 
                           pattern = "fastq", full.names = TRUE) ,
         fq = list.files(here(data_dir, "host"), 
                         pattern = "fastq", full.names = FALSE) )%>%
  separate(fq,into=c("stage","seq_run","trt", "line", "transfer", "subpop", "NexteraBC", "read", "sufx"),
           sep=("-|_"), remove = F) %>%   
  #fix extraction names to unify across sequencing runs
  mutate(subpop = case_when(
    subpop %in% c("SB","rS") ~ "revived_spore",
    subpop %in% c("veg","rV") ~ "revived_total",
    subpop %in% c("pl") ~ "frozen_pellet")) %>% 
  # sample name
  mutate(unq.sample = str_c(trt, line, transfer, sep = '-'))

d.ddup <- bind_rows(host.ddup, phage.ddup) %>% 
  select("unq.sample", "phageORhost", "trt", "line", "transfer", "subpop",
         "seq_run" , "read", "fq") %>% 
  mutate(fq = str_c(data_dir ,phageORhost ,fq, sep = "/"),
         read = str_c("ddup_", read)) %>% 
  pivot_wider(names_from = read, values_from = fq)


# merge raw and deduplicated ----------------------------------------------

# separate out ancestral samples
d.ancestral <-  
  full_join(
    d.raw %>% filter(!transfer %in% (c("T1","T4", "T7", "T10", "T14"))),
    d.ddup %>% filter(!transfer %in% (c("T1","T4", "T7", "T10", "T14"))),
    by = c("unq.sample", "phageORhost", "trt", "line", 
           "transfer", "subpop","seq_run" )
  )

# separate out evolved samples
d.evol <-  
  full_join(
    d.raw %>% filter(transfer %in% (c("T1","T4", "T7", "T10", "T14"))),
    d.ddup %>% filter(transfer %in% (c("T1","T4", "T7", "T10", "T14"))),
    by = c("unq.sample", "phageORhost", "trt", "line", 
           "transfer", "subpop","seq_run" )
  )

# parse trt

d.evol <- d.evol %>% 
  
  mutate (ancestor = case_when(phageORhost == "phage" ~ "SPO1-ANC",
                               str_starts(trt, "W") ~ "delta6-ANC",
                               str_starts(trt, "S") ~ "dspoIIE-ANC")) %>% 
  mutate (ancestor_fasta = str_c("data/",ancestor, ".fa", sep = "")) %>% 
  mutate (ancestor_gff = str_c("data/",ancestor, ".gff3", sep = "")) %>% 

  mutate (seed.bank = case_when(str_starts(trt, "WL") ~ "long_seed_bank",
                                str_starts(trt, "WS") ~ "short_seed_bank",
                                str_starts(trt, "SN") ~ "no_seed_bank")) %>% 
  
  mutate(phage_trt = case_when(str_ends(trt, "O") ~ "SPO1",
                               str_ends(trt, "Ct") ~ "noPhage")) %>% 
  
  relocate(ancestor, seed.bank, phage_trt, .before = subpop )


# Export data -------------------------------------------------------------

write_csv(d.evol, here("data", "samples_evolved.csv"))  
write_csv(d.ancestral, here("data", "samples_ancestors.csv")) 


