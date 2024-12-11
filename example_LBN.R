library("devtools")
devtools::install_github("r-lib/svglite")
library(tidyverse)
library(fs)

##import files
setwd=("/Volumes/Data/NGSraw/nkeith/DeepSeq6/pamStats/")
allfiles = setwd %>% 
  dir_ls(regexp = "\\.tsv$") %>% 
  map_dfr(read_tsv, .id = "source", na="noPhage") 

##getting set up
collapsed = allfiles %>%
  select(-score) %>%
  select(-name) %>%
  group_by(source, chrom, start, stop, NGG, strand, spacer) %>%
  mutate(unique = n()) %>%
  ungroup() %>%
  distinct_all()

##lastbasenorm 
fullDF = collapsed %>%
  ungroup() %>%
  filter(chrom == "Spyogenes_M1GAS_NC_002737.2") %>%
  mutate(normed = case_when(str_detect(spacer, "a$") ~ unique*1.812265748, 
                            str_detect(spacer, "c$") ~ unique*0.267975994, 
                            str_detect(spacer, "g$") ~ unique*5.640697298, 
                            str_detect(spacer, "t$") ~ unique*2.467159171))

P0_df = fullDF %>%
  filter(source == "/Volumes/Data/NGSraw/nkeith/DeepSeq6/pamStats/JW5507_P0.tsv") %>% 
  mutate(newTotal = sum(normed)) %>%
  mutate(RPM = case_when(str_detect(source, "JW5507_P0.tsv") ~ normed/(newTotal/1000000))) 


P8_df = fullDF %>%
  filter(source == "/Volumes/Data/NGSraw/nkeith/DeepSeq6/pamStats/JW5507_P8.tsv") %>% 
  mutate(newTotal = sum(normed)) %>%
  mutate(RPM = case_when(str_detect(source, "JW5507_P8.tsv") ~ normed/(newTotal/1000000))) 



P0ngg_df = fullDF %>%
  filter(source == "/Volumes/Data/NGSraw/nkeith/DeepSeq6/pamStats/JW5507_P0.tsv") %>% 
  mutate(newTotal = sum(normed)) %>%
  mutate(RPM = case_when(str_detect(source, "JW5507_P0.tsv") ~ normed/(newTotal/1000000))) %>%
  filter(NGG =="y") 
#filter(Phage == "noPhage") %>%

P8ngg_df = fullDF %>%
  filter(source == "/Volumes/Data/NGSraw/nkeith/DeepSeq6/pamStats/JW5507_P8.tsv") %>% 
  mutate(newTotal = sum(normed)) %>%
  mutate(RPM = case_when(str_detect(source, "JW5507_P8.tsv") ~ normed/(newTotal/1000000))) %>%
  filter(NGG =="y") 
#filter(Phage == "noPhage") %>%

forPlot = rbind(P0_df, P8_df)
forPlotngg = rbind(P0ngg_df, P8ngg_df)


ggplot(forPlot, aes(x=start, y=RPM, color = source)) + 
  geom_line(aes(y=RPM), stat="summary_bin", binwidth=10000, fun = "sum")  + 
  annotate("rect", xmin = 529587, xmax = 570505, ymin = 0, ymax = 50000, alpha = .3, fill = "yellow" ) +
  annotate("rect", xmin = 778520, xmax = 821005, ymin = 0, ymax = 50000, alpha = .2 , fill = "red") +
  annotate("rect", xmin = 1189120, xmax = 1222649, ymin = 0, ymax = 50000, alpha = .2 , fill = "orange") +
  annotate("rect", xmin = 1773339, xmax = 1786888, ymin = 0, ymax = 50000, alpha = .2 , fill = "green" ) +
  theme_classic() +
  theme(legend.position="none")

ggsave("/Volumes/Data/NGSraw/nkeith/DeepSeq6/Graphs/JW5507_P0_P8.pdf")


ggplot(forPlotngg, aes(x=start, y=RPM, color = source)) + 
  geom_line(aes(y=RPM), stat="summary_bin", binwidth=10000, fun = "sum")  + 
  annotate("rect", xmin = 529587, xmax = 570505, ymin = 0, ymax = 50000, alpha = .3, fill = "yellow" ) +
  annotate("rect", xmin = 778520, xmax = 821005, ymin = 0, ymax = 50000, alpha = .2 , fill = "red") +
  annotate("rect", xmin = 1189120, xmax = 1222649, ymin = 0, ymax = 50000, alpha = .2 , fill = "orange") +
  annotate("rect", xmin = 1773339, xmax = 1786888, ymin = 0, ymax = 50000, alpha = .2 , fill = "green" ) +
  theme_classic() +
  theme(legend.position="none")

ggsave("/Volumes/Data/NGSraw/nkeith/DeepSeq6/Graphs/JW5507_P0_P8_NGG.pdf")


write_delim(forPlot, "/Volumes/Data/NGSraw/nkeith/DeepSeq6/LBN/finalNumbers/JW5507_P0_P8_LBN.tsv", delim = "\t")







