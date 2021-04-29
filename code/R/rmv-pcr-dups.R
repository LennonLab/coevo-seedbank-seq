# Removing PCR duplicates using FastUniq v1.1
# https://sourceforge.net/projects/fastuniq/files/
# installed in my home directory at 
# `/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source`. 
fastuniq.path <- "/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq"

# phage or host?
  input.arg <- NULL
  
  input.arg <- commandArgs(trailingOnly = TRUE)
  # test if there is an argument: if not, return an error
  if ( length(input.arg)!=1 ) 
    stop("one argument must be supplied: phage OR host", call.=FALSE)
  
  if (! input.arg %in% c("phage", "host"))
    stop("one argument must be supplied: phage OR host", call.=FALSE)


# load libraries
  library(here, quietly = T, verbose = F)
  library(tidyverse, quietly = T, verbose = F)


############################
# copy fastq files and unzip
############################
  
  #make directories
  if (!dir.exists(here("data", "ddup-fastq")))
    dir.create(here("data", "ddup-fastq"))
  if (!dir.exists(here("data", "ddup-fastq", input.arg)))
    dir.create(here("data", "ddup-fastq", input.arg))

  #copy fastq files
  fastq <- list.files(here("data/input/fastq/", input.arg), 
                      pattern = "fastq", full.names = T)           
  file.copy(fastq, here("data", "ddup-fastq", input.arg))
  
  # unzip
  fastq <- list.files(here("data", "ddup-fastq", input.arg), pattern = "fastq.gz", full.names = T)
  for (i in fastq) 
    system(paste("gunzip", i))

############################
# make input lists for FastUniq
############################
# Each list file should contain the file names of paired sequencing files
  
  # parse file name
  reads <-
    tibble(f = list.files(here("data", "ddup-fastq", input.arg), pattern = "fastq"))%>%
    separate(1,into=c("run","num","trt", "line", "sufx"), sep="-", remove = F) %>% 
    mutate(unq.sample = interaction(trt, line, sep = '-'))
  
  # make list files
  if (!dir.exists(here("data", "ddup-lists")))
    dir.create(here("data", "ddup-lists"))
  if (!dir.exists(here("data", "ddup-lists", input.arg)))
    dir.create(here("data", "ddup-lists", input.arg))
  
  for (i in unique(reads$unq.sample)){
    cur <- reads %>%
      filter(unq.sample == i)

    write_lines(cur$f, file = here("data", "ddup-lists", input.arg, paste0("input_list_",i,".txt")))
  }


############################
# Run FastUniq 
############################
  # save list of files to be removed (before duplication)
  fastq <- list.files(here("data", "ddup-fastq", input.arg), pattern = "fastq", full.names = T)
  
  #save the fastuniq commands as they are passed to system
  path.sh <- here("data/ddup-lists", paste0("deduplicate-", input.arg, ".sh"))
  

  setwd(here("data", "ddup-fastq", input.arg))

  inputs <-  list.files( here("data", "ddup-lists", input.arg), pattern = ".txt", full.names = T)

  for (i in 1:length(inputs)){
    out <- read_lines(inputs[i]) %>%
            paste0(here("data/ddup-fastq/"), input.arg, "/ddup-",.)

  sys.cmd <- paste(
     fastuniq.path,
    "-i", inputs[i], # input list
    "-t q", # Output sequence format: fastq
    "-o", out[1], # first output file
    "-p", out[2], # second output file
    "-c 1") # New serial numbers assigned by FastUniq

  write_lines(paste(sys.cmd,"\n"), path.sh, append = T)
  
  #run fastuniq
  system(sys.cmd)
  }

############################ 
# clean up
############################
  # remove the copies of the original fastq files
  file.remove(fastq)
  
  #zip dduplicated files
  system("ls | grep \"ddup\" | xargs gzip")

quit(save = "no")