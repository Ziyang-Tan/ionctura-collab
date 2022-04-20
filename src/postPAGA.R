library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tibble)
library(wesanderson)


data_dir = './PAGA_result_data_manual_gating_scale/'

for (tumor_group in c('Tumor Group 1', 'Tumor Group 2', 'Tumor Group 3')){
  sub_names <- sapply(dir(path=paste0(data_dir, tumor_group)), 
                      function(x){strsplit(x,split = '_')[[1]][3]})
  for (sub_name in sub_names){
    result_dir = file.path('./figures_manual_gating_scale', tumor_group, sub_name)
    file_path = Sys.glob(paste0(data_dir, tumor_group, '/*', sub_name, '*raw.h5ad'))
    datah5 = readH5AD(file = file_path)
    data = t(assay(datah5)) %>% as_tibble()
    meta = colData(datah5) %>% 
      as_tibble() %>%
      select(timepoint, Subject.ID, dose, tumor_group, leiden, Subject.ID) %>%
      mutate(timepoint = factor(timepoint, levels = c('AC1D1', 'AC1D2', 'AC1D15', 'AC2D1', 'AC3D1', 'AC4D1',
                                                      'AC5D1','AC6D1','AC7D1','AC8D1','AC9D1','AC10D1','AC11D1',
                                                      'AC12D1','AC13D1','AC14D1','AC15D1','AC17D1','AC18D1')))
    
    # subpop freq change
    data_freq = meta %>% 
      group_by(timepoint, tumor_group, leiden) %>%
      summarise(n = n()) %>%
      # tidyr::complete(timepoint, tumor_group, leiden, fill = list(n = 0)) %>%
      ungroup() %>%
      group_by(timepoint, tumor_group) %>% # for each timepoint x tumor_group, frequency add up to 1
      mutate(freq = n / sum(n))
    
    g_list <- lapply(sort(unique(meta$leiden)), function(x){
      ggplot(data = data_freq %>% filter(leiden == x), 
             aes(x = timepoint, y = freq, color = tumor_group, group = tumor_group)) +
        geom_line() +
        geom_point() +
        labs(title = paste0(sub_name, '_cluster_', x)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })
    ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
      ggexport(filename = file.path(result_dir, paste0(sub_name, '_cluster_changes_timepoint.pdf')),
               width = 20,
               height = 10)
    # subpop freq change by subject
    data_freq = meta %>% 
      group_by(timepoint, leiden, Subject.ID) %>%
      summarise(n = n()) %>%
      # tidyr::complete(timepoint, leiden, Subject.ID, fill = list(n = 0)) %>%
      ungroup() %>%
      group_by(timepoint, Subject.ID) %>% # for each timepoint x Subject.ID, frequency add up to 1
      mutate(freq = n / sum(n))
    
    g_list <- lapply(sort(unique(meta$leiden)), function(x){
      ggplot(data = data_freq %>% filter(leiden == x), 
             aes(x = timepoint, y = freq, color = Subject.ID, group = Subject.ID)) +
        geom_line() +
        geom_point() +
        labs(title = paste0(sub_name, '_cluster_', x)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })
    ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
      ggexport(filename = file.path(result_dir, paste0(sub_name, '_cluster_changes_timepoint_by_subject.pdf')),
               width = 20,
               height = 10)
    ################################
    
    # subpop marker expression
    data_anno <- data %>% add_column(leiden = meta$leiden)
    
    g_list <- lapply(colnames(data), function(x){
      d <- data_anno %>% select(x, leiden) %>% rename(marker_expression = x)
      ggplot(d, 
             aes(x = marker_expression, y = leiden, fill = leiden)) +
        geom_density_ridges(scale = 4, rel_min_height=.01) +
        scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$leiden), type = "continuous")) +
        theme_ridges() + 
        theme(legend.position = "none") +
        xlim(quantile(d$marker_expression, 0.01),
             quantile(d$marker_expression, 0.99)) +
        labs(title = x)
    })
    ggarrange(plotlist = g_list, ncol = 5, nrow = 9) %>% 
      ggexport(filename = file.path(result_dir, paste0(sub_name, '_cluster_marker_expressions.pdf')),
               width = 20,
               height = 40)
    
  }
}
