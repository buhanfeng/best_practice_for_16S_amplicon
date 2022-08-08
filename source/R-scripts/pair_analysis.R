
#
# set enviroment
#
#项目名，可以不设置。设置后，所有生成的文件前面都会加上basename
basename='batch2_EC'
base_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/batch2_EC'
data_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/data'
#command_path
command_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/microbiome_batch'
# manifest_tab的路径
manifest='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/batch2_EC/manifest.tsv'
# qiime 使用的cpu数
CPU=40
tab_tax_base_on='unoise3'
working_path=paste(base_path, tab_tax_base_on, sep = '/')
prefix = paste0(basename, '_')

#
# library
#
library(sqldf)
library(xlsx)
library(stringr)

merge_tab_path = paste0(working_path, '/post_tab_tax/', 'merge_tab.txt')
tax_path = paste0(working_path, '/post_tab_tax/', 'taxonomy.tsv')
merge_tab = read.csv(merge_tab_path, sep = '\t')
tax = read.csv(tax_path, sep = '\t')
manifest_tab=read.csv(file = manifest, sep = '\t')
Feature_ID = tax$Feature_ID

Domain = c()
Phylum = c()
Class = c()
Ord = c()
Family = c()
Genus = c()
Species = c()

text = tax[1, 2]
text_v = strsplit(text, split='; ', fixed = FALSE, perl = FALSE, useBytes = FALSE)
text_v = unlist(text_v)
print(text_v[1])

for (i in c(1:nrow(tax))) {
  text = tax[i, 2]
  text_v = strsplit(text, split='; ', fixed = FALSE, perl = FALSE, useBytes = FALSE)
  text_v = unlist(text_v)
  domain = 'unknow'
  phylum = 'unknow'
  class = 'unknow'
  order = 'unknow'
  family = 'unknow'
  genus = 'unknow'
  species = 'unknow'
  if(length(text_v)==1){
    domain = text_v[1]
  }else if(length(text_v)==2){
    domain = text_v[1]
    phylum = text_v[2]
  }else if(length(text_v)==3){
    domain = text_v[1]
    phylum = text_v[2]
    class = text_v[3]
  }else if(length(text_v)==4){
    domain = text_v[1]
    phylum = text_v[2]
    class = text_v[3]
    order = text_v[4]
  }else if(length(text_v)==5){
    domain = text_v[1]
    phylum = text_v[2]
    class = text_v[3]
    order = text_v[4]
    family = text_v[5]
  }else if(length(text_v)==6){
    domain = text_v[1]
    phylum = text_v[2]
    class = text_v[3]
    order = text_v[4]
    family = text_v[5]
    genus = text_v[6]
  }else if(length(text_v)==7){
    domain = text_v[1]
    phylum = text_v[2]
    class = text_v[3]
    order = text_v[4]
    family = text_v[5]
    genus = text_v[6]
    species=text_v[7]
  }
  Domain = append(Domain, domain)
  Phylum = append(Phylum, phylum)
  Class = append(Class, class)
  Ord = append(Ord, order)
  Family = append(Family, family)
  Genus = append(Genus, genus)
  Species = append(Species, species)
}

tax_tab = data.frame(Feature_ID,Domain,Phylum,Class,Ord,Family,Genus,Species)

id=manifest_tab$id
idd=str_c(id, collapse=',')
# pool and select need to amend here
# idd='BAT05_S3W15E,BAT05_S3W15C,BAT05_S3W16E,BAT05_S3W16C,BAT05_S3W17E,BAT05_S3W17C,BAT05_S3W18E_1,BAT05_S3W18C_1,BAT05_S3W18E_2,BAT05_S3W18C_2,BAT05_S3W19E_1,BAT05_S3W19C_1,BAT05_S3W19E_2,BAT05_S3W19C_2,BAT05_M20E,BAT05_M20C,BAT05_M21E,BAT05_M21C,BAT05_M22E,BAT05_M22C,BAT05_M23E_1,BAT05_M23C_1,BAT05_M23E_2,BAT05_M23C_2,BAT05_M24E_1,BAT05_M24C_1,BAT05_M24E_2,BAT05_M24C_2,BAT05_W26E,BAT05_W26C,BAT05_W27E,BAT05_W27C,BAT05_W28E,BAT05_W28C,BAT05_W29E_1,BAT05_W29C_1,BAT05_W29E_2,BAT05_W29C_2,BAT05_W30E_1,BAT05_W30C_1,BAT05_W30E_2,BAT05_W30C_2'
sta = paste0("select OTU_ID,", idd,",Domain,Phylum,Class,Ord,Family,Genus,Species from merge_tab, tax_tab where OTU_ID=Feature_ID order by OTU_ID")
print(sta)
init_tab = sqldf(sta)
# xlsx::write.xlsx(init_tab, file = paste(working_path, 'init_tab.xlsx', sep = '/'), row.names = F)
# write.table(init_tab, file = paste(working_path, 'init_tab.csv', sep = '/'), row.names = F, sep = '\t')


#
# Venn图
#
library(sqldf)
library(stringr)
library(VennDiagram)
library(xlsx)
# distance
library(ecodist)
library(philentropy)
# hilbertSimilarity
library(hilbertSimilarity)

# 这三个向量需要手工加内容
bi_compare_list = c("BAT02_M01E/BAT02_M01C","BAT02_M02E/BAT02_M02C","BAT02_M03E/BAT02_M03C","BAT02_W02E/BAT02_W02C","BAT02_W03E/BAT02_W03C","BAT02_W04E/BAT02_W04C")
# bi_compare_list = c("BAT03_M04E/BAT03_M04C","BAT03_M05E/BAT03_M05C","BAT03_W05E/BAT03_W05C","BAT03_W06E/BAT03_W06C","BAT03_M06E/BAT03_M06C","BAT03_M07E/BAT03_M07C","BAT03_W07E/BAT03_W07C","BAT03_W08E/BAT03_W08C")
# bi_compare_list = c("BAT05_M20E/BAT05_M20C","BAT05_M21E/BAT05_M21C","BAT05_M22E/BAT05_M22C","BAT05_M23E_1/BAT05_M23C_1","BAT05_M23E_2/BAT05_M23C_2","BAT05_M23E_1/BAT05_M23C_2","BAT05_M23E_2/BAT05_M23C_1","BAT05_M24E_1/BAT05_M24C_1","BAT05_M24E_2/BAT05_M24C_2","BAT05_M24E_1/BAT05_M24C_2","BAT05_M24E_2/BAT05_M24C_1","BAT05_W26E/BAT05_W26C","BAT05_W27E/BAT05_W27C","BAT05_W28E/BAT05_W28C","BAT05_W29E_1/BAT05_W29C_1","BAT05_W29E_2/BAT05_W29C_2","BAT05_W29E_1/BAT05_W29C_2","BAT05_W29E_2/BAT05_W29C_1","BAT05_W30E_1/BAT05_W30C_1","BAT05_W30E_2/BAT05_W30C_2","BAT05_W30E_1/BAT05_W30C_2","BAT05_W30E_2/BAT05_W30C_1")
# bi_compare_list = c("BAT05_S3W15E/BAT05_S3W15C","BAT05_S3W16E/BAT05_S3W16C","BAT05_S3W17E/BAT05_S3W17C","BAT05_S3W18E_1/BAT05_S3W18C_1","BAT05_S3W18E_2/BAT05_S3W18C_2","BAT05_S3W19E_1/BAT05_S3W19C_1","BAT05_S3W19E_2/BAT05_S3W19C_2","BAT05_M20E/BAT05_M20C","BAT05_M21E/BAT05_M21C","BAT05_M22E/BAT05_M22C","BAT05_M23E_1/BAT05_M23C_1","BAT05_M23E_2/BAT05_M23C_2","BAT05_M24E_1/BAT05_M24C_1","BAT05_M24E_2/BAT05_M24C_2","BAT05_W26E/BAT05_W26C","BAT05_W27E/BAT05_W27C","BAT05_W28E/BAT05_W28C","BAT05_W29E_1/BAT05_W29C_1","BAT05_W29E_2/BAT05_W29C_2","BAT05_W30E_1/BAT05_W30C_1","BAT05_W30E_2/BAT05_W30C_2","BAT05_S3W18E_1/BAT05_S3W18E_2","BAT05_S3W18C_1/BAT05_S3W18C_2","BAT05_S3W19E_1/BAT05_S3W19E_2","BAT05_S3W19C_1/BAT05_S3W19C_2","BAT05_M23E_1/BAT05_M23E_2","BAT05_M24E_1/BAT05_M24E_2","BAT05_M23C_1/BAT05_M23C_2","BAT05_M24C_1/BAT05_M24C_2","BAT05_W29E_1/BAT05_W29E_2","BAT05_W30E_1/BAT05_W30E_2","BAT05_W29C_1/BAT05_W29C_2","BAT05_W30C_1/BAT05_W30C_2")
tri_compare_list = c()
qua_compare_list = c()
level = c('OTU_ID',"Domain", "Phylum", "Class", "Ord", "Family", "Genus", "Species")
venn_path = paste0(feature_analysis_path, '/bi_compare_list')
if( length(bi_compare_list)!=0){
  if(dir.exists(venn_path)){
    file.remove(dir(venn_path , pattern="*"))
  }else{
    dir.create(venn_path)
  }
}
if( length(tri_compare_list)!=0 ){
  dir.create(paste0(feature_analysis_path, '/tri_compare_list'))
}
if( length(qua_compare_list)!=0 ){
  dir.create(paste0(feature_analysis_path, '/qua_compare_list'))
}
temp = "select sum(##SAM) as ##SAM, ##le from init_tab where ##SAM>0 group by ##le"
tax_level = c()
# metrics_list = c('epithelium', 'content', 'common', 'spearman', 'pearson' )
metrics_list = c('hellinger', 'jaccard', 'spearman')
level_list = c('OTU_ID',"Domain", "Phylum", "Class", "Ord", "Family", "Genus", "Species")
level = c()
metrics = c()
# BAT05_S3W15E_BAT05_S3W15C = c()
# BAT05_S3W16E_BAT05_S3W16C = c()
# BAT05_S3W17E_BAT05_S3W17C = c()
# BAT05_S3W18E_1_BAT05_S3W18C_1 = c()
# BAT05_S3W18E_2_BAT05_S3W18C_2 = c()
# BAT05_S3W19E_1_BAT05_S3W19C_1 = c()
# BAT05_S3W19E_2_BAT05_S3W19C_2 = c()
# BAT05_M20E_BAT05_M20C = c()
# BAT05_M21E_BAT05_M21C = c()
# BAT05_M22E_BAT05_M22C = c()
# BAT05_M23E_1_BAT05_M23C_1 = c()
# BAT05_M23E_2_BAT05_M23C_2 = c()
# BAT05_M24E_1_BAT05_M24C_1 = c()
# BAT05_M24E_2_BAT05_M24C_2 = c()
# BAT05_W26E_BAT05_W26C = c()
# BAT05_W27E_BAT05_W27C = c()
# BAT05_W28E_BAT05_W28C = c()
# BAT05_W29E_1_BAT05_W29C_1 = c()
# BAT05_W29E_2_BAT05_W29C_2 = c()
# BAT05_W30E_1_BAT05_W30C_1 = c()
# BAT05_W30E_2_BAT05_W30C_2 = c()
# BAT05_S3W18E_1_BAT05_S3W18E_2 = c()
# BAT05_S3W18C_1_BAT05_S3W18C_2 = c()
# BAT05_S3W19E_1_BAT05_S3W19E_2 = c()
# BAT05_S3W19C_1_BAT05_S3W19C_2 = c()
# BAT05_M23E_1_BAT05_M23E_2 = c()
# BAT05_M24E_1_BAT05_M24E_2 = c()
# BAT05_M23C_1_BAT05_M23C_2 = c()
# BAT05_M24C_1_BAT05_M24C_2 = c()
# BAT05_W29E_1_BAT05_W29E_2 = c()
# BAT05_W30E_1_BAT05_W30E_2 = c()
# BAT05_W29C_1_BAT05_W29C_2 = c()
# BAT05_W30C_1_BAT05_W30C_2 = c()

BAT02_M01E_BAT02_M01C = c()
BAT02_M02E_BAT02_M02C = c()
BAT02_M03E_BAT02_M03C = c()
BAT02_W02E_BAT02_W02C = c()
BAT02_W03E_BAT02_W03C = c()
BAT02_W04E_BAT02_W04C = c()

# BAT03_M04E_BAT03_M04C = c()
# BAT03_M05E_BAT03_M05C = c()
# BAT03_W05E_BAT03_W05C = c()
# BAT03_W06E_BAT03_W06C = c()
# BAT03_M06E_BAT03_M06C = c()
# BAT03_M07E_BAT03_M07C = c()
# BAT03_W07E_BAT03_W07C = c()
# BAT03_W08E_BAT03_W08C = c()

# BAT05_M20E_BAT05_M20C = c()
# BAT05_M21E_BAT05_M21C = c()
# BAT05_M22E_BAT05_M22C = c()
# BAT05_M23E_1_BAT05_M23C_1 = c()
# BAT05_M23E_2_BAT05_M23C_2 = c()
# BAT05_M23E_1_BAT05_M23C_2 = c()
# BAT05_M23E_2_BAT05_M23C_1 = c()
# BAT05_M24E_1_BAT05_M24C_1 = c()
# BAT05_M24E_2_BAT05_M24C_2 = c()
# BAT05_M24E_1_BAT05_M24C_2 = c()
# BAT05_M24E_2_BAT05_M24C_1 = c()
# BAT05_W26E_BAT05_W26C = c()
# BAT05_W27E_BAT05_W27C = c()
# BAT05_W28E_BAT05_W28C = c()
# BAT05_W29E_1_BAT05_W29C_1 = c()
# BAT05_W29E_2_BAT05_W29C_2 = c()
# BAT05_W29E_1_BAT05_W29C_2 = c()
# BAT05_W29E_2_BAT05_W29C_1 = c()
# BAT05_W30E_1_BAT05_W30C_1 = c()
# BAT05_W30E_2_BAT05_W30C_2 = c()
# BAT05_W30E_1_BAT05_W30C_2 = c()
# BAT05_W30E_2_BAT05_W30C_1 = c()

for (le in level_list) {
  for (me in metrics_list) {
    level = append(level, le)
    metrics = append(metrics, me)
    i = 0
    for(bi in bi_compare_list){
      sam = stringr::str_split(bi, '/')
      sam = as.vector(unlist(sam))
      sta = str_replace_all(temp, '##SAM', sam[1])
      sta = str_replace_all(sta, '##le', le)
      set1 = sqldf(sta)
      sta = str_replace_all(temp, '##SAM', sam[2])
      sta = str_replace_all(sta, '##le', le)
      set2 = sqldf(sta)
      s1t = length(set1[[le]])
      s2t = length(set2[[le]])
      ii = length(intersect(set1[[le]], set2[[le]]))
      vv = 0
      if(me == 'epithelium'){
        vv = (s1t - ii)/s1t
      }else if(me == 'content'){
        vv = (s2t - ii)/s2t
      }else if(me == 'common'){
        vv = ii/(s1t + s2t - ii)
      }else if(me == 'spearman'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        cc = cor(set[sam[1]],set[sam[2]], method = 'spearman')
        vv = cc[1,1]
      }else if(me == 'pearson'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        cc = cor(set[sam[1]],set[sam[2]], method = 'pearson')
        vv = cc[1,1]
      }else if(me == 'jaccard'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        setd = set[,c(sam[1],sam[2])]
        sett = t(setd)
        sett1 = sett[1,]/sum(sett[1,])
        sett2 = sett[2,]/sum(sett[2,])
        sett = rbind(sett1, sett2)
        vv = ecodist::distance(sett, method = "jaccard")
        vv = 1 - vv
      }else if(me == 'manhattan'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        setd = set[,c(sam[1],sam[2])]
        sett = t(setd)
        sett1 = sett[1,]/sum(sett[1,])
        sett2 = sett[2,]/sum(sett[2,])
        sett = rbind(sett1, sett2)
        vv = 1 - ecodist::distance(sett, method = "manhattan")
      }else if(me == 'hellinger'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        setd = set[,c(sam[1],sam[2])]
        sett = t(setd)
        sett1 = sett[1,]/sum(sett[1,])
        sett2 = sett[2,]/sum(sett[2,])
        sett = rbind(sett1, sett2)
        vv = 1 - philentropy::distance(sett, method = "hellinger")
      }else if(me == 'jensen-shannon'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        setd = set[,c(sam[1],sam[2])]
        sett = t(setd)
        sett1 = sett[1,]/sum(sett[1,])
        sett2 = sett[2,]/sum(sett[2,])
        sett = rbind(sett1, sett2)
        vv = 1 - ecodist::distance(sett, method = "jensen-shannon")
      }else if(me == 'bray-curtis'){
        ss = "select sum(##SAM1) as ##SAM1, sum(##SAM2) as ##SAM2, ##le from init_tab group by ##le"
        ss = str_replace_all(ss, '##SAM1', sam[1])
        ss = str_replace_all(ss, '##SAM2', sam[2])
        ss = str_replace_all(ss, '##le', le)
        set = sqldf(ss)
        setd = set[,c(sam[1],sam[2])]
        sett = t(setd)
        sett1 = sett[1,]/sum(sett[1,])
        sett2 = sett[2,]/sum(sett[2,])
        sett = rbind(sett1, sett2)
        vv = ecodist::distance(sett, method = "bray-curtis")
        vv = 1 - vv
      }
      
      i = i+1
      
      if(i==1){
        BAT02_M01E_BAT02_M01C = append(BAT02_M01E_BAT02_M01C, vv)
      }else if(i==2){
        BAT02_M02E_BAT02_M02C = append(BAT02_M02E_BAT02_M02C, vv)
      }else if(i==3){
        BAT02_M03E_BAT02_M03C = append(BAT02_M03E_BAT02_M03C, vv)
      }else if(i==4){
        BAT02_W02E_BAT02_W02C = append(BAT02_W02E_BAT02_W02C, vv)
      }else if(i==5){
        BAT02_W03E_BAT02_W03C = append(BAT02_W03E_BAT02_W03C, vv)
      }else if(i==6){
        BAT02_W04E_BAT02_W04C = append(BAT02_W04E_BAT02_W04C, vv)
      }
      
      # if(i==1){
      #   BAT03_M04E_BAT03_M04C = append(BAT03_M04E_BAT03_M04C, vv)
      # }else if(i==2){
      #   BAT03_M05E_BAT03_M05C = append(BAT03_M05E_BAT03_M05C, vv)
      # }else if(i==3){
      #   BAT03_W05E_BAT03_W05C = append(BAT03_W05E_BAT03_W05C, vv)
      # }else if(i==4){
      #   BAT03_W06E_BAT03_W06C = append(BAT03_W06E_BAT03_W06C, vv)
      # }else if(i==5){
      #   BAT03_M06E_BAT03_M06C = append(BAT03_M06E_BAT03_M06C, vv)
      # }else if(i==6){
      #   BAT03_M07E_BAT03_M07C = append(BAT03_M07E_BAT03_M07C, vv)
      # }else if(i==7){
      #   BAT03_W07E_BAT03_W07C = append(BAT03_W07E_BAT03_W07C, vv)
      # }else if(i==8){
      #   BAT03_W08E_BAT03_W08C = append(BAT03_W08E_BAT03_W08C, vv)
      # }
      
      # if(i==1){
      #   BAT05_M20E_BAT05_M20C = append(BAT05_M20E_BAT05_M20C, vv)
      # }else if(i==2){
      #   BAT05_M21E_BAT05_M21C = append(BAT05_M21E_BAT05_M21C, vv)
      # }else if(i==3){
      #   BAT05_M22E_BAT05_M22C = append(BAT05_M22E_BAT05_M22C, vv)
      # }else if(i==4){
      #   BAT05_M23E_1_BAT05_M23C_1 = append(BAT05_M23E_1_BAT05_M23C_1, vv)
      # }else if(i==5){
      #   BAT05_M23E_2_BAT05_M23C_2 = append(BAT05_M23E_2_BAT05_M23C_2, vv)
      # }else if(i==6){
      #   BAT05_M23E_1_BAT05_M23C_2 = append(BAT05_M23E_1_BAT05_M23C_2, vv)
      # }else if(i==7){
      #   BAT05_M23E_2_BAT05_M23C_1 = append(BAT05_M23E_2_BAT05_M23C_1, vv)
      # }else if(i==8){
      #   BAT05_M24E_1_BAT05_M24C_1 = append(BAT05_M24E_1_BAT05_M24C_1, vv)
      # }else if(i==9){
      #   BAT05_M24E_2_BAT05_M24C_2 = append(BAT05_M24E_2_BAT05_M24C_2, vv)
      # }else if(i==10){
      #   BAT05_M24E_1_BAT05_M24C_2 = append(BAT05_M24E_1_BAT05_M24C_2, vv)
      # }else if(i==11){
      #   BAT05_M24E_2_BAT05_M24C_1 = append(BAT05_M24E_2_BAT05_M24C_1, vv)
      # }else if(i==12){
      #   BAT05_W26E_BAT05_W26C = append(BAT05_W26E_BAT05_W26C, vv)
      # }else if(i==13){
      #   BAT05_W27E_BAT05_W27C = append(BAT05_W27E_BAT05_W27C, vv)
      # }else if(i==14){
      #   BAT05_W28E_BAT05_W28C = append(BAT05_W28E_BAT05_W28C, vv)
      # }else if(i==15){
      #   BAT05_W29E_1_BAT05_W29C_1 = append(BAT05_W29E_1_BAT05_W29C_1, vv)
      # }else if(i==16){
      #   BAT05_W29E_2_BAT05_W29C_2 = append(BAT05_W29E_2_BAT05_W29C_2, vv)
      # }else if(i==17){
      #   BAT05_W29E_1_BAT05_W29C_2 = append(BAT05_W29E_1_BAT05_W29C_2, vv)
      # }else if(i==18){
      #   BAT05_W29E_2_BAT05_W29C_1 = append(BAT05_W29E_2_BAT05_W29C_1, vv)
      # }else if(i==19){
      #   BAT05_W30E_1_BAT05_W30C_1 = append(BAT05_W30E_1_BAT05_W30C_1, vv)
      # }else if(i==20){
      #   BAT05_W30E_2_BAT05_W30C_2 = append(BAT05_W30E_2_BAT05_W30C_2, vv)
      # }else if(i==21){
      #   BAT05_W30E_1_BAT05_W30C_2 = append(BAT05_W30E_1_BAT05_W30C_2, vv)
      # }else if(i==22){
      #   BAT05_W30E_2_BAT05_W30C_1 = append(BAT05_W30E_2_BAT05_W30C_1, vv)
      # }
      
      # if(i==1){
      #   BAT05_S3W15E_BAT05_S3W15C = append(BAT05_S3W15E_BAT05_S3W15C, vv)
      # }else if(i==2){
      #   BAT05_S3W16E_BAT05_S3W16C = append(BAT05_S3W16E_BAT05_S3W16C, vv)
      # }else if(i==3){
      #   BAT05_S3W17E_BAT05_S3W17C = append(BAT05_S3W17E_BAT05_S3W17C, vv)
      # }else if(i==4){
      #   BAT05_S3W18E_1_BAT05_S3W18C_1 = append(BAT05_S3W18E_1_BAT05_S3W18C_1, vv)
      # }else if(i==5){
      #   BAT05_S3W18E_2_BAT05_S3W18C_2 = append(BAT05_S3W18E_2_BAT05_S3W18C_2, vv)
      # }else if(i==6){
      #   BAT05_S3W19E_1_BAT05_S3W19C_1 = append(BAT05_S3W19E_1_BAT05_S3W19C_1, vv)
      # }else if(i==7){
      #   BAT05_S3W19E_2_BAT05_S3W19C_2 = append(BAT05_S3W19E_2_BAT05_S3W19C_2, vv)
      # }else if(i==8){
      #   BAT05_M20E_BAT05_M20C = append(BAT05_M20E_BAT05_M20C, vv)
      # }else if(i==9){
      #   BAT05_M21E_BAT05_M21C = append(BAT05_M21E_BAT05_M21C, vv)
      # }else if(i==10){
      #   BAT05_M22E_BAT05_M22C = append(BAT05_M22E_BAT05_M22C, vv)
      # }else if(i==11){
      #   BAT05_M23E_1_BAT05_M23C_1 = append(BAT05_M23E_1_BAT05_M23C_1, vv)
      # }else if(i==12){
      #   BAT05_M23E_2_BAT05_M23C_2 = append(BAT05_M23E_2_BAT05_M23C_2, vv)
      # }else if(i==13){
      #   BAT05_M24E_1_BAT05_M24C_1 = append(BAT05_M24E_1_BAT05_M24C_1, vv)
      # }else if(i==14){
      #   BAT05_M24E_2_BAT05_M24C_2 = append(BAT05_M24E_2_BAT05_M24C_2, vv)
      # }else if(i==15){
      #   BAT05_W26E_BAT05_W26C = append(BAT05_W26E_BAT05_W26C, vv)
      # }else if(i==16){
      #   BAT05_W27E_BAT05_W27C = append(BAT05_W27E_BAT05_W27C, vv)
      # }else if(i==17){
      #   BAT05_W28E_BAT05_W28C = append(BAT05_W28E_BAT05_W28C, vv)
      # }else if(i==18){
      #   BAT05_W29E_1_BAT05_W29C_1 = append(BAT05_W29E_1_BAT05_W29C_1, vv)
      # }else if(i==19){
      #   BAT05_W29E_2_BAT05_W29C_2 = append(BAT05_W29E_2_BAT05_W29C_2, vv)
      # }else if(i==20){
      #   BAT05_W30E_1_BAT05_W30C_1 = append(BAT05_W30E_1_BAT05_W30C_1, vv)
      # }else if(i==21){
      #   BAT05_W30E_2_BAT05_W30C_2 = append(BAT05_W30E_2_BAT05_W30C_2, vv)
      # }else if(i==22){
      #   BAT05_S3W18E_1_BAT05_S3W18E_2 = append(BAT05_S3W18E_1_BAT05_S3W18E_2, vv)
      # }else if(i==23){
      #   BAT05_S3W18C_1_BAT05_S3W18C_2 = append(BAT05_S3W18C_1_BAT05_S3W18C_2, vv)
      # }else if(i==24){
      #   BAT05_S3W19E_1_BAT05_S3W19E_2 = append(BAT05_S3W19E_1_BAT05_S3W19E_2, vv)
      # }else if(i==25){
      #   BAT05_S3W19C_1_BAT05_S3W19C_2 = append(BAT05_S3W19C_1_BAT05_S3W19C_2, vv)
      # }else if(i==26){
      #   BAT05_M23E_1_BAT05_M23E_2 = append(BAT05_M23E_1_BAT05_M23E_2, vv)
      # }else if(i==27){
      #   BAT05_M24E_1_BAT05_M24E_2 = append(BAT05_M24E_1_BAT05_M24E_2, vv)
      # }else if(i==28){
      #   BAT05_M23C_1_BAT05_M23C_2 = append(BAT05_M23C_1_BAT05_M23C_2, vv)
      # }else if(i==29){
      #   BAT05_M24C_1_BAT05_M24C_2 = append(BAT05_M24C_1_BAT05_M24C_2, vv)
      # }else if(i==30){
      #   BAT05_W29E_1_BAT05_W29E_2 = append(BAT05_W29E_1_BAT05_W29E_2, vv)
      # }else if(i==31){
      #   BAT05_W30E_1_BAT05_W30E_2 = append(BAT05_W30E_1_BAT05_W30E_2, vv)
      # }else if(i==32){
      #   BAT05_W29C_1_BAT05_W29C_2 = append(BAT05_W29C_1_BAT05_W29C_2, vv)
      # }else if(i==33){
      #   BAT05_W30C_1_BAT05_W30C_2 = append(BAT05_W30C_1_BAT05_W30C_2, vv)
      # }
    }
  }
}
# res = data.frame(level, metrics, BAT05_S3W15E_BAT05_S3W15C, BAT05_S3W16E_BAT05_S3W16C, BAT05_S3W17E_BAT05_S3W17C, BAT05_S3W18E_1_BAT05_S3W18C_1, BAT05_S3W18E_2_BAT05_S3W18C_2, BAT05_S3W19E_1_BAT05_S3W19C_1, BAT05_S3W19E_2_BAT05_S3W19C_2, BAT05_M20E_BAT05_M20C, BAT05_M21E_BAT05_M21C, BAT05_M22E_BAT05_M22C, BAT05_M23E_1_BAT05_M23C_1, BAT05_M23E_2_BAT05_M23C_2, BAT05_M24E_1_BAT05_M24C_1, BAT05_M24E_2_BAT05_M24C_2, BAT05_W26E_BAT05_W26C, BAT05_W27E_BAT05_W27C, BAT05_W28E_BAT05_W28C, BAT05_W29E_1_BAT05_W29C_1, BAT05_W29E_2_BAT05_W29C_2, BAT05_W30E_1_BAT05_W30C_1, BAT05_W30E_2_BAT05_W30C_2, BAT05_S3W18E_1_BAT05_S3W18E_2, BAT05_S3W18C_1_BAT05_S3W18C_2, BAT05_S3W19E_1_BAT05_S3W19E_2, BAT05_S3W19C_1_BAT05_S3W19C_2, BAT05_M23E_1_BAT05_M23E_2, BAT05_M24E_1_BAT05_M24E_2, BAT05_M23C_1_BAT05_M23C_2, BAT05_M24C_1_BAT05_M24C_2, BAT05_W29E_1_BAT05_W29E_2, BAT05_W30E_1_BAT05_W30E_2, BAT05_W29C_1_BAT05_W29C_2, BAT05_W30C_1_BAT05_W30C_2)
res = data.frame(level, metrics, BAT02_M01E_BAT02_M01C, BAT02_M02E_BAT02_M02C, BAT02_M03E_BAT02_M03C, BAT02_W02E_BAT02_W02C, BAT02_W03E_BAT02_W03C, BAT02_W04E_BAT02_W04C)
# res = data.frame(level, metrics, BAT03_M04E_BAT03_M04C, BAT03_M05E_BAT03_M05C, BAT03_W05E_BAT03_W05C, BAT03_W06E_BAT03_W06C, BAT03_M06E_BAT03_M06C, BAT03_M07E_BAT03_M07C, BAT03_W07E_BAT03_W07C, BAT03_W08E_BAT03_W08C)
# res = data.frame(level, metrics, BAT05_M20E_BAT05_M20C, BAT05_M21E_BAT05_M21C, BAT05_M22E_BAT05_M22C, BAT05_M23E_1_BAT05_M23C_1, BAT05_M23E_2_BAT05_M23C_2, BAT05_M23E_1_BAT05_M23C_2, BAT05_M23E_2_BAT05_M23C_1, BAT05_M24E_1_BAT05_M24C_1, BAT05_M24E_2_BAT05_M24C_2, BAT05_M24E_1_BAT05_M24C_2, BAT05_M24E_2_BAT05_M24C_1, BAT05_W26E_BAT05_W26C, BAT05_W27E_BAT05_W27C, BAT05_W28E_BAT05_W28C, BAT05_W29E_1_BAT05_W29C_1, BAT05_W29E_2_BAT05_W29C_2, BAT05_W29E_1_BAT05_W29C_2, BAT05_W29E_2_BAT05_W29C_1, BAT05_W30E_1_BAT05_W30C_1, BAT05_W30E_2_BAT05_W30C_2, BAT05_W30E_1_BAT05_W30C_2, BAT05_W30E_2_BAT05_W30C_1)

res = sqldf("select * from res order by metrics")
res = sqldf("select * from res order by level")
write.table(res, file = 'D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/EC_sig/bat2.txt', col.names = T, row.names = F, quote = F, sep = '\t')

#
# plot
#
library(stringr)
library(ggplot2)
all = read.csv(file = 'D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/EC_sig/all.txt', sep = '\t')
cn = colnames(all)
sample = cn[c(-1,-2)]
# sample = factor(sample, labels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33'))
# sample = as.factor(sample)
# batch = c('batch2', 'batch2','batch2','batch2','batch2','batch2','batch3','batch3','batch3','batch3','batch3','batch3','batch3','batch3','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5','batch5')
# pheno = c('WS3','WS3','WS3','WS3','WS3','WS3','WS3','MS2','MS2','MS2','MS2','MS2','MS2','MS2','WS2','WS2','WS2','WS2','WS2','WS2','WS2','REP','REP','REP','REP','REP','REP','REP','REP','REP','REP','REP','REP')
pheno = c("M", "M", "M", "W", "W", "W", "M", "M", "W", "W", "M", "M", "W", "W", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W")
# pheno = c('WS3','WS3','WS3','WS3','WS3','WS3','WS3','MS2','MS2','MS2','MS2','MS2','MS2','MS2','WS2','WS2','WS2','WS2','WS2','WS2','WS2')
sample = apply(data.frame(sample, pheno), 1, function(x){
  str_c(x['pheno'], '_',x['sample'])
})
i = 0
for(i in 1:nrow(all)){
  iv = all[i,]
  level = iv[1]
  metrics = iv[2]
  value = iv[c(-1,-2)]
  value = as.numeric(value)
  value = round(value, 2)
  temp = data.frame(value, sample, pheno)
  p = ggplot(data=temp,mapping=aes(x=value,y=sample,fill=pheno,group=factor(1))) + 
    geom_bar(stat="identity", width = 0.5)+ 
    theme(axis.text.y = element_text(angle = 0, hjust = 0.05, vjust = 0.4)) +
    geom_text(aes(label = value, vjust = 0, hjust = 0), show.legend = TRUE)
    #xlim(0.6,1)
    # scale_x_continuous(limits=c(0.6,1))
  ggsave(filename = paste0("D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/EC_sig/draw4/",level,'_',metrics, '.png'), plot = p)
}



#
# significance
#
all = read.csv(file = 'D:\\desktop\\advance\\projects\\zebrafish_gut_microbiome\\microbiome\\EC_sig\\results\\all.txt', sep = '\t')
i = 0
kruskal_test = c()
t_test = c()
# grou1 = c("M", "M", "M", "W", "W", "W", "M", "M", "W", "W", "M", "M", "W", "W", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W")
# grou1 = c("M", "M", "M", "W", "W", "W", "M", "M", "W", "W", "M", "M", "W", "W", "M", "M", "M", "M", "M", "M", "M", "W", "W", "W", "W", "W", "W", "W")
grou1 = c("M", "M", "M", "W", "W", "W", "M", "M", "W", "W", "M", "M", "W", "W", "M", "M", "M", "M", "M", "M", "M", "W", "W", "W", "W", "W", "W", "W")

grou1 = as.factor(grou1)
for(i in 1:nrow(all)){
  iv = all[i,]
  iv = iv[c(-1,-2)]
  iv = as.numeric(iv)
  temp = data.frame(iv, grou1)
  hh = kruskal.test(iv~grou1,data=temp)
  kruskal_test = append(kruskal_test, hh$p.value)
  tt = t.test(iv~grou1,data=temp)
  t_test = append(t_test, tt$p.value)
}

all = cbind(all, kruskal_test, t_test)

write.table(all, file = 'D:\\desktop\\advance\\projects\\zebrafish_gut_microbiome\\microbiome\\EC_sig\\results\\significance.txt', col.names = T, row.names = F, quote = F, sep = '\t')

