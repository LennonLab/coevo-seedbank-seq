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


#------------------#
# Collect the data #
#------------------#

# files to read
time_course_files <- list.files(here("data/CompSubPops/host/CompSubPop_merged2"),
                                pattern = ".bz",full.names = TRUE, recursive = TRUE)


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

#------------------#
# Arrange the data #
#------------------#

# parameters needed to correct for spoIIE deletion
    del.start <- 70538
    del.wt.end <- 73021
    del.mut.end <- 70692
    del.adj <- del.wt.end-del.mut.end

# Separate out the different "time" points
# Done separetly for WT and mut

#### process WT ####
time.key <- strsplit(d.raw$times[nchar(d.raw$times) %>% which.max()],
                     " ", fixed = T) %>%
  unlist() 

d.wt <- d.raw %>% 
  filter(str_starts(pop, "W")) %>% 

  #correct for spoIIE deletion
  # remove mutations in deleted genes
  filter(position < del.start | position > del.wt.end) %>% 

  mutate(mut = paste0(as.character(position), alt_allele)) %>% 
  separate(pop, into = c("trt", "rep")) %>% 
  separate(col=n.allele, into = paste0("At-",time.key), sep = " ") %>% 
  separate(col=n.depth, into = paste0("Dt-",time.key), sep = " ") %>% 
  pivot_longer(cols = contains("At")|contains("Dt"),
               names_to = c(".value", "transfer"),
               names_sep = "-") %>% 
  mutate(At = as.integer(At), Dt = as.integer(Dt),
         transfer = as.integer(transfer)) %>% 
  mutate(ft =  At/Dt)

#### process mutant #### 
time.key <- strsplit(d.raw$times[nchar(d.raw$times) %>% which.min()],
                     " ", fixed = T) %>%
  unlist() 

d.mut <- d.raw %>% 
  filter(str_starts(pop, "S")) %>% 
  
  #correct for spoIIE deletion
  # 1. remove mutations in deletion scar
  filter(position < del.start | position > del.mut.end) %>% 
  # 2. adjust the position number of mutations left of deletion
  mutate(position = if_else(position > del.mut.end, position + del.adj, position)) %>% 
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

#### Put everything together #### 

d <- bind_rows(d.wt, d.mut) %>% 
  mutate(subpop = case_when(transfer == 100 ~ "pellet",
                            transfer == 200 ~ "revived_spore",
                            transfer == 300 ~ "revived_total")) %>% 
  mutate(subpop=as.factor(subpop))%>%
  mutate(subpop = fct_relevel(subpop, "pellet","revived_total","revived_spore")) %>% 
  
  #ajust factor order for seed bank
  mutate (seed.bank = case_when(str_starts(trt, "WL") ~ "long seed bank",
                                str_starts(trt, "WS") ~ "short seed bank",
                                str_starts(trt, "SN") ~ "no seed bank")) %>%
  mutate(seed.bank=as.factor(seed.bank))%>%
  mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank")) %>% 
  #ajust factor order for seed bank
  mutate (phage = case_when(str_detect(trt, "O") ~ "SPO1",
                            str_detect(trt, "Ct") ~ "no phage")) 

# distribution of depths
d %>% 
  ggplot(aes(position, Dt))+
  geom_point()+
  facet_wrap(~ interaction(trt, rep), scales = "free")
# there are a few loci tha consistenetly have much bigher depths.
# X500 seens to be a usefule cutoff

d %>% 
  filter(Dt < 500) %>% 
  ggplot(aes(position, Dt))+
  geom_point()+
  facet_wrap(~ interaction(trt, rep), scales = "free")


d %>% 
  filter(Dt < 500) %>% 
  ggplot(aes(position, ft))+
  geom_point(aes(fill = phage), alpha=0.5, shape=21)+
  facet_grid(seed.bank~ subpop)+
  theme_cowplot()+
  panel_border()+
  scale_fill_manual(values = c("white", "black"))+
  ggsave(here("plots/host.mut.freq.png"), width = 11, height = 8)

d %>% 
  filter(Dt < 500) %>% 
  ggplot(aes(Dt))+
  geom_histogram(aes(fill = phage), bins=100)+
  facet_grid(seed.bank ~ subpop, scales = "free")+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"))+
  ggtitle("Sequencing depth at mutated positions")+
  xlab("Depth")+
  ggsave(here("plots/host-depth-histogram.png"), width = 8, height = 5)

d %>% 
  filter(Dt < 500) %>% 
  filter(Dt > 5) %>% 
  filter(At > 5) %>% 
  filter(ft > 0.01) %>% 
  ggplot(aes(At))+
  geom_histogram(aes(fill = phage), bins=100)+
  facet_grid(seed.bank ~ subpop, scales = "free")+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"))+
  # ggtitle("Sequencing depth at mutated positions")+
  xlab("Allele reads")


#---------------------------------#
# compare pellet to revived total #
#---------------------------------# 

d %>% 
  filter (phage == "no phage") %>% 
  filter (subpop != "revived_spore") %>% 
  select(trt, rep, mut, seed.bank, subpop, ft) %>% 
  pivot_wider(names_from = "subpop", values_from = ft) %>% 
  
  ggplot(aes(pellet, revived_total))+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype=2)+
  geom_point(shape=21, size = 1)+
  facet_grid(rep~seed.bank)+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"))+
  # scale_x_log10()+scale_y_log10()+
  ggtitle("Control populations (no phage)")+
  ggsave(here("plots/pellet-vs-revived.png"), width = 8, height = 5)

d %>% 
  filter (phage == "no phage") %>% 
  filter (subpop != "pellet") %>%
  filter (seed.bank != "no seed bank") %>% 
  select(trt, rep, mut, seed.bank, subpop, ft) %>% 
  pivot_wider(names_from = "subpop", values_from = ft) %>% 
  
  ggplot(aes(revived_spore, revived_total))+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype=2)+
  geom_point(shape=21, size = 1)+
  facet_grid(rep~seed.bank)+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"))+
  # scale_x_log10()+scale_y_log10()+
  ggtitle("Control populations (no phage)")+
  ggsave(here("plots/total-vs-spore.png"), width = 5, height = 5)

d %>% 
  filter (phage != "no phage") %>% 
  filter (subpop != "pellet") %>%
  filter (seed.bank != "no seed bank") %>% 
  select(trt, rep, mut, seed.bank, subpop, ft) %>% 
  pivot_wider(names_from = "subpop", values_from = ft) %>% 
  
  ggplot(aes(revived_spore, revived_total))+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype=2)+
  geom_point(shape=21, size = 1)+
  facet_grid(rep~seed.bank)+
  theme_classic()+
  panel_border(size = 1.5, color = "black")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"))+
  # scale_x_log10()+scale_y_log10()+
  ggtitle("Infected populations (SPO1 phage)")+
  ggsave(here("plots/total-vs-spore-infected.png"), width = 5, height = 5)



p <- d %>%
  filter(subpop == "revived_total") %>%
  ggplot(aes(x=position, y=ft))+
  geom_point(alpha=0.5)+
  geom_segment(aes(xend=position, yend=0),alpha=0.5)+
  facet_grid(phage + rep ~ seed.bank )+
  theme_bw()
save_plot(here("plots/host_all_muts_GoodPipe.png"),p, base_width = 15, base_height = 8)



#----------------------#
# distance to ancestor #
#----------------------# 
# Focusing on data from revived populations

# make unique identifers for mutations and populations
  d.mat <- d %>%
    filter(Dt < 500) %>% 
    filter(Dt > 5) %>% 
    filter(At > 5) %>% 
    filter(ft > 0.01) %>% 
      filter(subpop != "pellet") %>% 
  mutate(pop=interaction(trt, rep, subpop)) %>% 
  select(pop, mut, frequency=ft)

# verify that each mutation is unique
check <- d.mat %>% 
  group_by(pop, mut)%>%
  summarise(n=n(), .groups = "drop")%>%
  arrange(desc(n)) %>% 
  head()


#arrange as distance matrix
d.mat <- d.mat %>%
  select(pop, mut, frequency) %>%
  pivot_wider(values_from = frequency, names_from = mut, values_fill = 0)

# add ancestor - by definition has none of the mutations
d.mat <- d.mat %>%
  add_row(pop = "ANC-ANC-T0") 


d.mat[d.mat$pop=="ANC-ANC-T0",-1] <- 0.0


m <- as.matrix(d.mat[,-1])
dimnames(m) = list(d.mat$pop, colnames(d.mat)[-1])

dst.m <- dist(m, method = "euclidean")
# heatmap(m, Colv = NA, Rowv = NA)

require(ape)

tr <- bionj(dst.m)
tr <-root(tr, outgroup = "ANC-ANC-T0", resolve.root = T)

#make tip labels
tips <- d.mat %>%
  select(pop) %>%
  separate(pop, into = c("trt", "line", "subpop"), sep = "\\.", remove = F) %>%
  mutate(tip.color = case_when(str_starts(trt, "WLO") ~ rgb(1,0,0),
                               str_starts(trt, "WSO") ~ rgb(0,1,0),
                               str_starts(trt, "SNO") ~ rgb(0,0,1),
                               str_starts(trt, "WLCt") ~ rgb(1,0.5,0),
                               str_starts(trt, "WSCt") ~ rgb(0.5,1,0),
                               str_starts(trt, "SNCt") ~ rgb(0.5,0,1))) 


# mutate (seed.bank = case_when(str_starts(trt, "WL") ~ "long seed bank",
#                               str_starts(trt, "WS") ~ "short seed bank",
#                               str_starts(trt, "SN") ~ "no seed bank")) %>%
#   mutate(seed.bank=as.factor(seed.bank))%>%
#   mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank")) %>% 
#   #ajust factor order for seed bank
 

png(filename = here("plots/host_T14_bionj.png"),width = 12, height = 12, units = "in", res=150)
plot(tr, show.tip.label = T, type = "fan", root.edge = F, tip.color = tips$tip.color, font = 2)

title("BioNJ from Euclidean distance matrix ")
# tiplabels(tips$pop, cex = .8, frame = "none",  bg = NA, adj = 0, col = tips$tip.color)
dev.off()

###########
br.ln <- cophenetic.phylo(tr) %>% 
  as_tibble() %>% 
  mutate(pop = colnames(.))

d.anc <- select(br.ln, pop, dist2anc = `ANC-ANC-T0`) %>% 
  filter(!str_detect(pop, "ANC")) %>% 
  separate(pop, into = c("trt", "line", "subpop"), sep = "\\.") %>% 
  #ajust factor order for subpop
  mutate(subpop=as.factor(subpop))%>%
  mutate(subpop = fct_relevel(subpop, "pellet","revived_total","revived_spore")) %>% 
  
  #ajust factor order for seed bank
  mutate (seed.bank = case_when(str_starts(trt, "WL") ~ "long seed bank",
                                str_starts(trt, "WS") ~ "short seed bank",
                                str_starts(trt, "SN") ~ "no seed bank")) %>%
  mutate(seed.bank=as.factor(seed.bank))%>%
  mutate(seed.bank = fct_relevel(seed.bank, "long seed bank","short seed bank","no seed bank")) %>% 
  #ajust factor order for seed bank
  mutate (phage = case_when(str_detect(trt, "O") ~ "SPO1",
                            str_detect(trt, "Ct") ~ "no phage")) 
sum.d.anc <- d.anc %>% 
  group_by(seed.bank, subpop, phage) %>% 
  summarise(m=mean(dist2anc), v=sd(dist2anc)/sqrt(n()), n=n())


d.anc %>% 
  ggplot(aes(seed.bank, dist2anc, color = seed.bank))+
  geom_crossbar(data = sum.d.anc, aes(y = m, ymin = m-v, ymax = m+v))+
  geom_point(aes(fill = seed.bank), color="black", shape=21, size=3)+
  facet_grid(subpop ~ phage)+
  theme_cowplot()+
  panel_border(size = 1.5, color = "black")+
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white"),
        axis.text.x = element_blank())+
  ylab("Distance from Anc. host")+
  labs(caption = "box depicts mean±sem")+
  ggsave(here("plots/dist2anc.png"), width = 6, height = 4)

####################################
library(vegan)

nmds <- 
  br.ln %>% 
  select(-pop, -`ANC-ANC-T0`) %>% 
  metaMDS(., distance = "euclidean")
 
stressplot(nmds)

 sp <- nmds$species %>%
   tibble() %>% 
   rownames_to_column()

 ordiplot(nmds,type = "text")
################
d %>% 
  ggplot(aes(Dt))+
  geom_histogram(bins=300)+
  facet_wrap(~ interaction(trt, rep), scales = "free")


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
  