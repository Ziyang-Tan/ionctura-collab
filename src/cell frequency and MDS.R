library(stats)
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(robCompositions)
library(ggpubr)

info_path <- file.path('/Users/tan/Ionctura-collab/data/Sample_info_reformat.xlsx')
fileInfo = readxl::read_excel(info_path, sheet=1) %>%
  rbind(readxl::read_excel(info_path, sheet=2)) %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  mutate(timepoint = gsub('*', '', timepoint, fixed=TRUE)) %>%
  mutate(timepoint = factor(timepoint, levels = c('AC1D1', 'AC1D2', 'AC1D15', 'AC2D1', 'AC3D1', 'AC4D1',
                                                  'AC5D1','AC6D1','AC7D1','AC8D1','AC9D1','AC10D1','AC11D1',
                                                  'AC12D1','AC13D1','AC14D1','AC15D1','AC17D1','AC18D1')))
projs<- unique(fileInfo$Expt_ID)
data_paths <- file.path('/Users/tan/cytof_data/', projs, 'classifiedV3', 'abundance.csv')
raw <- lapply(data_paths, function(x){
  read_csv(x, show_col_types = FALSE) %>%
    mutate(subpop = paste(...1,...2,sep='/')) %>%
    select(-c('...1', '...2'))
}) %>% Reduce(f = function(...){left_join(..., by = 'subpop')})

# quality control by cell counts?
#

tmp <- raw %>% select(-subpop) %>% t()
colnames(tmp) = raw$subpop
dat <- tmp %>% as_tibble() %>% add_column(Sample_ID = rownames(tmp)) %>%
  replace(is.na(.), 0) %>%
  mutate(`abT-cells/CD4 T/` = `abT-cells/CD4 T/` + `abT-cells/CD4 T`) %>% # fix the abundance matrix
  select(-c('abT-cells/CD4 T'))

subpopDat <- dat %>%
  select(contains('abT')) %>%
  select(contains('CD4')) %>%
  mutate(`T regs` = rowSums(select(.,contains('Tregs')))) %>%
  select(-contains('Tregs')) %>%
  mutate(`T regs of CD4` = `T regs`/rowSums(.)) %>%
  add_column(Sample_ID = dat$Sample_ID) %>%
  tidyr::pivot_longer(!Sample_ID, names_to = 'subpop', values_to = 'frequency') %>%
  left_join(fileInfo, by='Sample_ID')

# plot Treg freq changes
g1 <- ggplot(subpopDat %>% filter(subpop == 'T regs'),
             aes(x = timepoint, y = frequency, color = Cohort_and_dose)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = as.numeric(timepoint), y = frequency), method='lm') +
  stat_cor(aes(x = as.numeric(timepoint), y = frequency), 
           label.x = 3,
           label.x.npc = "right",
           method = 'spearman') +
  ylab('T reg frequency of total') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

g2 <- ggplot(subpopDat %>% filter(subpop == 'T regs of CD4'),
             aes(x = timepoint, y = frequency, color = `Tumor_Groups`)) +
  geom_point() +
  geom_line(aes(group = Subject_ID)) +
  geom_smooth(aes(x = as.numeric(timepoint), y = frequency), method='lm') +
  stat_cor(aes(x = as.numeric(timepoint), y = frequency), 
           label.x = 3, 
           label.x.npc = "left",
           method = 'spearman') +
  ylab('T reg frequency of CD4 T cells') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(g1, g2, ncol=1) %>%
  ggexport(filename = '/Users/tan/Ionctura-collab/results/Treg changes_dose.pdf',
           width = 10, height = 10)

#ggplot(d, aes(fill=subpop, y=frequency, x=timepoint)) + 
#  geom_bar(position="fill", stat="identity")

# MDS
m <- (dat %>% select(-Sample_ID) %>% replace(.==0, 1e-9) %>% as.matrix())
datDist <- aDist(m)
datMDS = cmdscale(datDist)
colnames(datMDS) = c('MDS1', 'MDS2')
datMDS <- datMDS %>% 
  as_tibble() %>% 
  add_column(Sample_ID = dat$Sample_ID) %>%
  left_join(fileInfo, by='Sample_ID') %>%
  arrange(`Subject_ID`, `timepoint`)

nameList <- c('Tumor_Groups', 'Cohort_and_dose', 'timepoint', 'Expt_ID')
gList <- lapply(nameList, function(x){
  g <- ggplot(datMDS, aes_string(x='MDS1',y='MDS2',colour = paste0('`', x, '`'))) + 
    geom_point() +
    geom_path(aes_string(group='`Subject_ID`'))
  return(g)
})
ggarrange(plotlist = gList, ncol = 2, nrow = 2, align = 'hv') %>%
  ggexport(filename = '/Users/tan/Ionctura-collab/results/cell_frequency_mds.pdf',
           width = 12, height = 6)


