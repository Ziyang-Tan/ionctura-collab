library(stats)
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(ggpubr)

info_path <- file.path('/Users/tan/Ionctura-collab/data/Sample_info_reformat.xlsx')
fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  #rbind(readxl::read_excel(info_path, sheet=2)) %>% # only use 37 and 56
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  mutate(timepoint = gsub('*', '', timepoint, fixed=TRUE)) %>%
  mutate(timepoint = factor(timepoint, levels = c('AC1D1', 'AC1D2', 'AC1D15', 'AC2D1', 'AC3D1', 'AC4D1',
                                                  'AC5D1','AC6D1','AC7D1','AC8D1','AC9D1','AC10D1','AC11D1',
                                                  'AC12D1','AC13D1','AC14D1','AC15D1','AC17D1','AC18D1')))

data_path <- '/Users/tan/ionctura-collab/data/Frequency tables'

#
data <- lapply(dir(data_path), function(x){
  readxl::read_excel(file.path(data_path, x)) %>%
    mutate(Sample_ID = gsub('(/.*)|(\\[ )|( \\])|(\\.fcs)', '', Name)) %>%
    filter(grepl('Tregs|Monocytes', Name)) %>%
    mutate(subpop = gsub('.*/', '', Name))
}) %>% do.call(what=rbind) %>%
  left_join(fileInfo, by = 'Sample_ID')

# plot 
sub_name = 'Monocytes'
g1 <- ggplot(data %>% filter(subpop == sub_name),
             aes(x = timepoint, y = Statistic, color = Cohort_and_dose)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = as.numeric(timepoint), y = Statistic), method='lm') +
  stat_cor(aes(x = as.numeric(timepoint), y = Statistic), 
           label.x = 3,
           label.x.npc = "right",
           method = 'spearman') +
  ylab(paste0(sub_name, ' frequency of parent subtype')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g2 <- ggplot(data %>% filter(subpop == sub_name),
       aes(x = timepoint, y = Statistic, color = Tumor_Groups)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = as.numeric(timepoint), y = Statistic), method='lm') +
  stat_cor(aes(x = as.numeric(timepoint), y = Statistic), 
           label.x = 3,
           label.x.npc = "right",
           method = 'spearman') +
  ylab(paste0(sub_name, ' frequency of parent subtype')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(g1, g2, ncol=1) %>%
  ggexport(filename = paste0('/Users/tan/Ionctura-collab/manual_results/', sub_name, ' changes.pdf'),
           width = 10, height = 10)
