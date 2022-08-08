#
# set environment
#
# linux batch parameter
basename='batch5_EC'
data_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/data'
CPU=40

# windows batch parameter
base_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/whole_gut_all'
command_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/microbiome_batch'
manifest='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/whole_gut_all/manifest.tsv'
base_on='flt_unoise3_200_4'


#
# init env
#
working_path=paste(base_path, base_on, sep = '/')
prefix = paste0(basename, '_')

#
# init alpha data
#
data = read.csv(paste0(working_path, '/alpha_diversity/all_metrics_tab.tsv'), sep = '\t')
data = subset(data, select = -id)
manifest_tab = read.csv(manifest, sep = '\t')
data = cbind(manifest_tab, data)


#
# plotting box plot
#
library(ggplot2)
library(ggpubr)
library(ggsignif)

alpha_plot_path=paste(working_path, 'alpha_diversity', 'plot', sep = '/')
if(dir.exists(alpha_plot_path)){
  file.remove(dir(alpha_plot_path , pattern="*"))
}else{
  dir.create(alpha_plot_path)
}

metrics = c('ace', 'chao1', 'faith_pd', 'observed_features')

# plot(as.factor(data$strain), data$shannon_entropy)
# 
# ggplot(data,aes(strain, shannon_entropy,fill=strain)) +
#   geom_boxplot() +
#   geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
#   geom_signif(comparisons = list(c("Mutant", "Wildtype")), map_signif_level=T, textsize=6,test=t.test ,step_increase=0.2)

ace_plot = ggplot(data,aes(strain, ace ,fill=strain)) + 
  geom_boxplot()+
  # scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  stat_compare_means()+theme_bw()
ggsave(paste0(alpha_plot_path, '/ace_plot.png'), ace_plot)

chao1_plot = ggplot(data,aes(strain, chao1 ,fill=strain)) + 
  geom_boxplot()+
  # scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  stat_compare_means()+theme_bw()
ggsave(paste0(alpha_plot_path, '/chao1_plot.png'), chao1_plot)

faith_pd_plot = ggplot(data,aes(strain, faith_pd ,fill=strain)) + 
  geom_boxplot()+
  # scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  stat_compare_means()+theme_bw()
ggsave(paste0(alpha_plot_path, '/faith_pd_plot.png'), faith_pd_plot)

observed_features_plot = ggplot(data,aes(strain, observed_features ,fill=strain)) + 
  geom_boxplot()+
  # scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  stat_compare_means()+theme_bw()
ggsave(paste0(alpha_plot_path, '/observed_features_plot.png'), observed_features_plot)






