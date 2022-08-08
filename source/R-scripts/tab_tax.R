
#
# set enviroment
#
#项目名，可以不设置。设置后，所有生成的文件前面都会加上basename
basename='batch5_EC'
base_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/batch5_EC'
data_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/data'
#command_path
command_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/microbiome_batch'
# manifest_tab的路径
manifest='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/batch5_EC/manifest.tsv'
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
sta = paste0("select OTU_ID,", idd,",Domain,Phylum,Class,Ord,Family,Genus,Species from merge_tab, tax_tab where OTU_ID=Feature_ID order by OTU_ID")
init_tab = sqldf(sta)
# xlsx::write.xlsx(init_tab, file = paste(working_path, 'init_tab.xlsx', sep = '/'), row.names = F)
write.table(init_tab, file = paste(working_path, 'init_tab.csv', sep = '/'), row.names = F, sep = '\t')

# 样本Feature数量统计：
#
#
#set enviroment
feature_analysis_path=paste(working_path, 'feature_analysis', sep = '/')
if(dir.exists(feature_analysis_path)){
  file.remove(dir(feature_analysis_path , pattern="*"))
}else{
  dir.create(feature_analysis_path)
}
id=manifest_tab$id
#
# 柱状图
#
num=c()
for (ele in id) {
  cc = sqldf(paste('select count(OTU_ID) as c from init_tab where ', ele,  '!=0 ', sep = ''))
  num = append(num, cc$c)
}
id = append(id, 'total')
total=nrow(init_tab)
num = append(num, total)
feature_bar_data = data.frame(id, num)

library(ggplot2)
p = ggplot(data=feature_bar_data,mapping=aes(x=id,y=num,fill=num,group=factor(1))) + 
  geom_bar(stat="identity")+ 
  geom_text(aes(label = num, vjust = -0.8, hjust = 0.5), show.legend = TRUE)
ggsave(filename = paste0(feature_analysis_path, '\\feature_count.png'), plot = p)
#
# Venn图
#
library(sqldf)
library(stringr)
library(VennDiagram)
library(xlsx)
# 这三个向量需要手工加内容
# bi_compare_list = c("BAT03_M04E/BAT03_M04C","BAT03_M05E/BAT03_M05C","BAT03_W05E/BAT03_W05C","BAT03_W06E/BAT03_W06C","BAT03_M06E/BAT03_M06C","BAT03_M07E/BAT03_M07C","BAT03_W07E/BAT03_W07C","BAT03_W08E/BAT03_W08C")
# bi_compare_list = c("BAT02_M01E/BAT02_M01C","BAT02_M02E/BAT02_M02C","BAT02_M03E/BAT02_M03C","BAT02_W02E/BAT02_W02C","BAT02_W03E/BAT02_W03C","BAT02_W04E/BAT02_W04C")
bi_compare_list = c("BAT05_M20E/BAT05_M20C","BAT05_M21E/BAT05_M21C","BAT05_M22E/BAT05_M22C","BAT05_M23E_1/BAT05_M23C_1","BAT05_M23E_2/BAT05_M23C_2","BAT05_M23E_1/BAT05_M23C_2","BAT05_M23E_2/BAT05_M23C_1","BAT05_M24E_1/BAT05_M24C_1","BAT05_M24E_2/BAT05_M24C_2","BAT05_M24E_1/BAT05_M24C_2","BAT05_M24E_2/BAT05_M24C_1","BAT05_W26E/BAT05_W26C","BAT05_W27E/BAT05_W27C","BAT05_W28E/BAT05_W28C","BAT05_W29E_1/BAT05_W29C_1","BAT05_W29E_2/BAT05_W29C_2","BAT05_W29E_1/BAT05_W29C_2","BAT05_W29E_2/BAT05_W29C_1","BAT05_W30E_1/BAT05_W30C_1","BAT05_W30E_2/BAT05_W30C_2","BAT05_W30E_1/BAT05_W30C_2","BAT05_W30E_2/BAT05_W30C_1")
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
SAM1 = c()
SAM2 = c()
SAM1_TOTAL = c()
SAM2_TOTAL = c()
SAM1_sub_SAM2 = c()
SAM1_I_SAM2 = c()
SAM2_sub_SAM1 = c()
group1 = c()
for(bi in bi_compare_list){
  sam = stringr::str_split(bi, '/')
  sam = as.vector(unlist(sam))
  print(sam)
  for (le in level) {
    sta = str_replace_all(temp, '##SAM', sam[1])
    sta = str_replace_all(sta, '##le', le)
    print(sta)
    set1 = sqldf(sta)
    sta = str_replace_all(temp, '##SAM', sam[2])
    sta = str_replace_all(sta, '##le', le)
    print(sta)
    set2 = sqldf(sta)
    filename = str_c(le, sam[1], sam[2], sep = '_')
    filename = paste0(feature_analysis_path, '/bi_compare_list/', filename, '.png')
    vd = venn.diagram(x=list(set1=set1[[le]], set2=set2[[le]]), filename = filename, imagetype="png", col="white", fill=c(colors()[616], colors()[38]), alpha=c(0.8, 0.8), lwd=c(1, 1), cat.dist=c(-0.07, -0.07), cat.col='black')
    tax_level = append(tax_level, le)
    SAM1 = append(SAM1, sam[1])
    SAM2 = append(SAM2, sam[2])
    s1t = length(set1[[le]])
    SAM1_TOTAL = append(SAM1_TOTAL, s1t)
    s2t = length(set2[[le]])
    SAM2_TOTAL = append(SAM2_TOTAL, s2t)
    ii = length(intersect(set1[[le]], set2[[le]]))
    SAM1_I_SAM2 = append(SAM1_I_SAM2, ii)
    SAM1_sub_SAM2 = append(SAM1_sub_SAM2, s1t - ii)
    SAM2_sub_SAM1 = append(SAM2_sub_SAM1, s2t - ii)
    type = ''
    if(stringr::str_detect(sam[1], '_M') && stringr::str_detect(sam[2], '_M') ){
      type = 'Mutant'
    }else if(stringr::str_detect(sam[1], '_W') && stringr::str_detect(sam[2], '_W') ){
      type = 'Wildtype'
    }
    group1 = append(group1, type)
  }
}
bi_summary = data.frame(tax_level,group1, SAM1, SAM2, SAM1_TOTAL, SAM2_TOTAL, SAM1_sub_SAM2, SAM1_I_SAM2, SAM2_sub_SAM1)
bi_summary['ratio_Inter'] = bi_summary['SAM1_I_SAM2']/(bi_summary['SAM1_sub_SAM2']+bi_summary['SAM1_I_SAM2']+bi_summary['SAM2_sub_SAM1'])
bi_summary['ratio_1'] = bi_summary['SAM1_sub_SAM2']/(bi_summary['SAM1_sub_SAM2']+bi_summary['SAM1_I_SAM2']+bi_summary['SAM2_sub_SAM1'])
bi_summary['ratio_2'] = bi_summary['SAM2_sub_SAM1']/(bi_summary['SAM1_sub_SAM2']+bi_summary['SAM1_I_SAM2']+bi_summary['SAM2_sub_SAM1'])
# xlsx::write.xlsx(bi_summary, file = paste0(feature_analysis_path, '/bi_compare_list/bi_summary.xlsx'), row.names = F)
write.table(bi_summary, file = paste0(feature_analysis_path, '/bi_compare_list/bi_summary.csv'), row.names = F, sep = '\t')


#
# significance analysis
#
# tax_level = c()
# fisher_p = c()
# chi_square_p = c()
ratio_inter_h_test = c()
ratio_inter_t_test = c()
ratio_1_h_test = c()
ratio_1_t_test = c()
ratio_2_h_test = c()
ratio_2_t_test = c()
temp = "select * from bi_summary where tax_level='##le'"
for (le in level) {
  sta = str_replace_all(temp, '##le', le)
  dd = sqldf(sta)
  hh = kruskal.test(ratio_Inter~group1,data=dd)
  ratio_inter_h_test = append(ratio_inter_h_test, hh$p.value)
  tt = t.test(ratio_Inter~group1,data=dd)
  ratio_inter_t_test = append(ratio_inter_t_test, tt$p.value)
  
  hh = kruskal.test(ratio_1~group1,data=dd)
  ratio_1_h_test = append(ratio_1_h_test, hh$p.value)
  tt = t.test(ratio_1~group1,data=dd)
  ratio_1_t_test = append(ratio_1_t_test, tt$p.value)
  
  hh = kruskal.test(ratio_2~group1,data=dd)
  ratio_2_h_test = append(ratio_2_h_test, hh$p.value)
  tt = t.test(ratio_2~group1,data=dd)
  ratio_2_t_test = append(ratio_2_t_test, tt$p.value)
}
bi_sig_summary = data.frame(level, ratio_inter_h_test, ratio_inter_t_test, ratio_1_h_test, ratio_1_t_test, ratio_2_h_test, ratio_2_t_test)
# xlsx::write.xlsx(bi_summary, file = paste0(feature_analysis_path, '/bi_compare_list/bi_sig_summary.xlsx'), row.names = F)
write.table(bi_sig_summary, file = paste0(feature_analysis_path, '/bi_compare_list/bi_sig_summary.csv'), row.names = F, sep = '\t')

#
# 物种分析
#
tax_analysis_path=paste(working_path, tab_tax_base_on, 'tax_analysis', sep = '/')
if(dir.exists(tax_analysis_path)){
  file.remove(dir(tax_analysis_path , pattern="*"))
}else{
  dir.create(tax_analysis_path)
}

#
# 样品各等级物种统计表
#
id=manifest_tab$id
ty = c("Domain", "Phylum", "Class", "Ord", "Family", "Genus", "Species")
ro = length(id)
co = length(ty)
dd = rep.int(0, ro*co)
xx = matrix(dd, ro, co, byrow = F) 
rownames(xx) = id
colnames(xx) = ty

for (ele in id) {
  for (t in ty) {
    print(paste0("select count(*) as c from (select count(*) as num, ",t ," from init_tab where ", ele,  '!=0 and ', t, "!='unknow' group by ", t, ")"))
    cc = sqldf(paste0("select count(*) as c from (select count(*) as num, ",t ," from init_tab where ", ele,  '!=0 and ', t, "!='unknow' group by ", t, ")"))
    print(cc$c)
    xx[ele, t] = cc$c
  }
}
tax_stats = as.data.frame(xx)
xlsx::write.xlsx(tax_stats, file = paste(tax_analysis_path, 'tax_stats.xlsx', sep = '/'))

#
# 热力图
#
library(stringr)
library(ggplot2)
id=manifest_tab$id
idd=""
for(ele in id){
  uni = ", sum(##) as ##"
  idd = paste0(idd, str_replace_all(uni, '##', ele))
}
idp=str_c(id, collapse='+')
minn='20'
maxx='2000'
condi = paste0("(", idp, ")>", minn," and (", idp, ")<",maxx,"")

png_prefix = paste0(minn, '_', maxx, '_')

# set width
width = 480
if(length(id)>12 && length(id)<18){
  width = 720
}else if(length(id)>18){
  width = 960
}

sta = paste0("select OTU_ID", idd, " from init_tab where ", condi, " group by OTU_ID")
OTU_tab = sqldf(sta)
rownames(OTU_tab) = OTU_tab$OTU_ID
draw_tab <- subset(OTU_tab, select = -OTU_ID )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "OTU", ".png"), p1)

# sta = paste0("select Domain", idd, " from init_tab where ", condi, " group by Domain")
# Domain_tab = sqldf(sta)
# rownames(Domain_tab) = Domain_tab$Domain
# draw_tab <- subset(Domain_tab, select = -Domain )
# draw_tab = log10(draw_tab+1)
# d_heatmap(draw_tab)

sta = paste0("select Phylum", idd, " from init_tab where ", condi, " group by Phylum")
Phylum_tab = sqldf(sta)
rownames(Phylum_tab) = Phylum_tab$Phylum
draw_tab <- subset(Phylum_tab, select = -Phylum )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "Phylum", ".png"), p1)

sta = paste0("select Class", idd, " from init_tab where ", condi, " group by Class")
Class_tab = sqldf(sta)
rownames(Class_tab) = Class_tab$Class
draw_tab <- subset(Class_tab, select = -Class )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "Class", ".png"), p1)

sta = paste0("select Ord", idd, " from init_tab where ", condi, " group by Ord")
Ord_tab = sqldf(sta)
rownames(Ord_tab) = Ord_tab$Ord
draw_tab <- subset(Ord_tab, select = -Ord )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "Ord", ".png"), p1)

sta = paste0("select Family", idd, " from init_tab where ", condi, " group by Family")
Family_tab = sqldf(sta)
rownames(Family_tab) = Family_tab$Family
draw_tab <- subset(Family_tab, select = -Family )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "Family", ".png"), p1)


sta = paste0("select Genus", idd, " from init_tab where ", condi, " group by Genus")
Genus_tab = sqldf(sta)
rownames(Genus_tab) = Genus_tab$Genus
draw_tab <- subset(Genus_tab, select = -Genus )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "Genus", ".png"), p1)

sta = paste0("select Species", idd, " from init_tab where ", condi, " group by Species")
Species_tab = sqldf(sta)
rownames(Species_tab) = Species_tab$Species
draw_tab <- subset(Species_tab, select = -Species )
draw_tab = log10(draw_tab+1)
p1 = d_heatmap(draw_tab, width = width)
ggsave(paste0(tax_analysis_path, "/", png_prefix, "Species", ".png"), p1)


#
# 基于 lefse 的 biomarker 分析
#





#
# plot function
#
d_heatmap = function(data, outputpictype='png', clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan", filename="XX", saveppt=F, width = 480){
  # c("euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski")
  # data = "text_upload_file/b459aac8-7fb6-4c25-9b11-69ec350f153a.txt"
  # outputprefix = "user/_visitors/result_output/pretty_heatmap2021-10-29_20:36:08/9c444841-47f4-4b71-993e-88316a872b17"
  # outputpictype = "pdf"
  renameDuplicateRowNames = TRUE
  logv = "NULL"
  log_add = 0
  scale = "none"
  annotation_row = "NULL"
  annotation_col = "NULL"
  cluster_rows = T
  cluster_cols = F
  clustering_method = "complete"
  # clustering_distance_rows = "pearson"
  # clustering_distance_cols = "pearson"
  breaks = "NA"
  breaks_mid = NULL
  breaks_digits = 2
  correlation_plot = "None"
  maximum = Inf
  minimum = -Inf
  xtics_angle = 0
  manual_color_vector = "NULL"
  fontsize = 10
  manual_annotation_colors_sidebar = "NULL"
  cutree_cols = NA
  cutree_rows = NA
  kclu = NA
  ytics = F
  xtics = TRUE
  title = ""
  # width = 480
  height = 480
  debug = FALSE
  cluster_cols_variable = "NULL"
  cluster_rows_variable = "NULL"
  remove_cluster_cols_variable_in_annocol = FALSE
  remove_cluster_rows_variable_in_annorow = FALSE
  saveppt = FALSE
  
  
  # if (data == "") {
  #   script = sub(".*=", "", commandArgs()[4])
  #   #print(script)
  #   system(paste(script, "-h"))
  #   stop("At least -f is required!")
  # }
  # For all users
  # devtools::install_github("Tong-Chen/ImageGP")
  # For chinese users
  # devtools::install_git("https://gitee.com/ct5869/ImageGP.git")
  library(ImageGP)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(vegan)
  # if (outputprefix == ""){
  #   outputprefix = data
  # }
  # filename = paste0(outputprefix,  '.pheatmap.', outputpictype)
  
  if(breaks == 'NA'){
    breaks = NA
  }else if(breaks != "quantile"){
    breaks = as.numeric(sp_string2vector(breaks))
  }
  manual_color_vector = sp_string2vector(manual_color_vector)
  
  cat(sp_current_time(), "Starting...\n")
  
  gt = sp_pheatmap(
    data = data,
    # filename = filename,
    renameDuplicateRowNames = renameDuplicateRowNames,
    logv = logv,
    log_add = log_add,
    scale = scale,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    clustering_method = clustering_method,
    clustering_distance_rows = clustering_distance_rows,
    clustering_distance_cols = clustering_distance_cols,
    breaks = breaks,
    breaks_mid = breaks_mid,
    breaks_digits = breaks_digits,
    correlation_plot = correlation_plot,
    maximum = maximum,
    minimum = minimum,
    xtics_angle = xtics_angle,
    manual_color_vector = manual_color_vector,
    fontsize = fontsize,
    manual_annotation_colors_sidebar = manual_annotation_colors_sidebar,
    cutree_cols = cutree_cols,
    cutree_rows = cutree_rows,
    kclu = kclu,
    ytics = ytics,
    xtics = xtics,
    width = width,
    height = height,
    title = title,
    cluster_cols_variable=cluster_cols_variable,
    cluster_rows_variable=cluster_rows_variable,
    remove_cluster_cols_variable_in_annocol=remove_cluster_cols_variable_in_annocol,
    remove_cluster_rows_variable_in_annorow=remove_cluster_rows_variable_in_annorow,
    debug = debug,
    saveppt = saveppt
  )
  cat(sp_current_time(), "Success.\n")
  return(gt)
}


sp_pheatmap = function (data, filename = NA, renameDuplicateRowNames = F, 
          logv = NULL, log_add = 0, scale = "none", annotation_row = NULL, 
          annotation_col = NULL, cluster_rows = FALSE, cluster_cols = FALSE, 
          cluster_cols_variable = NULL, cluster_rows_variable = NULL, 
          remove_cluster_cols_variable_in_annocol = FALSE, remove_cluster_rows_variable_in_annorow = FALSE, 
          clustering_method = "complete", clustering_distance_rows = "pearson", 
          clustering_distance_cols = "pearson", breaks = NA, breaks_mid = NULL, 
          breaks_digits = 2, correlation_plot = "None", maximum = Inf, 
          minimum = -Inf, xtics_angle = 0, manual_color_vector = NULL, 
          fontsize = 14, manual_annotation_colors_sidebar = NULL, 
          cutree_cols = NA, cutree_rows = NA, kclu = NA, ytics = TRUE, 
          xtics = TRUE, width = 0, height = 0, title = "", debug = FALSE, 
          saveppt = FALSE, ...) 
{
  if (debug) {
    argg <- c(as.list(environment()), list(...))
    print(argg)
  }
  assignInNamespace(x = "draw_colnames", value = "draw_colnames_custom", 
                    ns = asNamespace("pheatmap"))
  if (class(data) == "character") {
    data <- sp_readTable(data, row.names = 1, renameDuplicateRowNames = renameDuplicateRowNames)
  }
  else if (class(data) != "data.frame") {
    stop("Unknown input format for `data` parameter.")
  }
  numeric_check = sapply(data, is.numeric)
  non_numeric_col = names(numeric_check[numeric_check == FALSE])
  if (length(non_numeric_col) > 0) {
    stop(paste(non_numeric_col, "contains non-numeric values."))
  }
  if (!sp.is.null(logv)) {
    if (log_add == 0) {
      log_add = sp_determine_log_add(data)
    }
    data <- eval(parse(text = logv))(data + log_add)
  }
  if (!sp.is.null(manual_color_vector)) {
    manual_color_vector <- generate_color_list(manual_color_vector, 
                                               100)
  }
  else {
    manual_color_vector <- colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name = "RdYlBu")))(100)
  }
  color_length = length(manual_color_vector)
  legend_breaks = NA
  legend_labels = NA
  if (length(breaks) > 1 || !is.na(breaks)) {
    if (length(breaks) == 1 && breaks == "quantile") {
      summary_v <- summary(c(t(data)))
      summary_v[1] <- summary_v[1]
      summary_v[6] <- summary_v[6]
      if (sp.is.null(breaks_mid)) {
        breaks <- unique(c(seq(summary_v[1], summary_v[2], 
                               length = color_length/4), seq(summary_v[2], 
                                                             summary_v[3], length = color_length/4), seq(summary_v[3], 
                                                                                                         summary_v[5], length = color_length/4), seq(summary_v[5], 
                                                                                                                                                     summary_v[6], length = color_length/4 - 1)))
        legend_breaks <- summary_v
      }
      else {
        breaks_mid <- as.numeric(breaks_mid)
        breaks <- unique(c(seq(summary_v[1], breaks_mid, 
                               length = color_length/2), seq(breaks_mid, 
                                                             summary_v[6], length = color_length/2 - 1)))
        legend_breaks <- c(summary_v[1], breaks_mid, 
                           summary_v[6])
      }
    }
    else {
      legend_breaks <- breaks
      length_breaks <- length(breaks)
      if (length_breaks < color_length) {
        manual_color_vector <- generate_color_list(c(manual_color_vector[1], 
                                                     manual_color_vector[100]), length_breaks + 
                                                     1)
      }
    }
    if (breaks_digits) {
      legend_breaks <- as.numeric(prettyNum(legend_breaks, 
                                            digits = breaks_digits))
    }
    legend_labels <- legend_breaks
  }
  if (!sp.is.null(annotation_row)) {
    if (class(annotation_row) == "character") {
      annotation_row <- sp_readTable(annotation_row, row.names = 1)
      annotation_row <- annotation_row[match(rownames(data), 
                                             rownames(annotation_row)), , drop = F]
    }
    if (!sp.is.null(cluster_rows_variable)) {
      if (!cluster_rows_variable %in% colnames(annotation_row)) {
        stop(paste(cluster_rows_variable, "must be one of column names of row annotation matrix!"))
      }
    }
  }
  else {
    annotation_row <- NA
  }
  if (!sp.is.null(annotation_col)) {
    if (class(annotation_col) == "character") {
      annotation_col <- sp_readTable(annotation_col, row.names = 1)
      annotation_col <- annotation_col[match(colnames(data), 
                                             rownames(annotation_col)), , drop = F]
    }
    if (!sp.is.null(cluster_cols_variable)) {
      if (!cluster_cols_variable %in% colnames(annotation_col)) {
        stop(paste(cluster_cols_variable, "must be one of column names of column annotation matrix!"))
      }
    }
  }
  else {
    annotation_col <- NA
  }
  data[data > maximum] <- maximum
  if (minimum != -Inf) {
    data[data < minimum] <- minimum
  }
  cor_data = F
  if (scale == "row") {
    data_sd <- apply(data, 1, sd)
    data <- data[data_sd != 0, ]
  }
  if (correlation_plot == "row" || correlation_plot == "Row") {
    if (clustering_distance_rows == "pearson") {
      row_cor = cor(t(data))
    }
    else if (clustering_distance_rows == "spearman") {
      row_cor = cor(t(data), method = "spearman")
    }
    else {
      row_cor = as.data.frame(as.matrix(dist(data, method = clustering_distance_rows)))
    }
    data = round(row_cor, 2)
    annotation_col = annotation_row
    cor_data = T
  }
  else if (correlation_plot == "col" || correlation_plot == 
           "Column") {
    if (clustering_distance_cols == "pearson") {
      col_cor = cor(data)
    }
    else if (clustering_distance_cols == "spearman") {
      col_cor = cor(data, method = "spearman")
    }
    else {
      col_cor = as.data.frame(as.matrix(dist(t(data), 
                                             method = clustering_distance_cols)))
    }
    data = round(col_cor, 2)
    cor_data = T
    annotation_row = annotation_col
  }
  if (width == 0 && height == 0) {
    height = nrow(data)
    width = ncol(data) * 1.1
    if (xtics_angle == 0) {
      width = width * 1.5
    }
    if (class(annotation_row) == "data.frame") {
      width = width + ncol(annotation_row)
      width = width * 1.1
    }
    if (class(annotation_col) == "data.frame") {
      height = height + ncol(annotation_col)
      width = width * 1.1
    }
    if (cluster_rows) {
      width = width + 4
    }
    if (cluster_cols) {
      height = height + 4
    }
    if (width < 8) {
      width = 8
    }
    else if (width < 20) {
      width = 8 + (width - 8)/4
    }
    else if (width < 100) {
      width = 10 + (width - 20)/5
    }
    else {
      width = 30
    }
    if (height < 10) {
      height = 8
    }
    else if (height < 20) {
      height = 8 + (height - 8)/4
    }
    else if (height < 100) {
      height = 11 + (height - 20)/5
    }
    else {
      height = 30
    }
  }
  if (sp.is.null(manual_annotation_colors_sidebar)) {
    manual_annotation_colors_sidebar = NA
  }
  else if (class(manual_annotation_colors_sidebar) == "character") {
    manual_annotation_colors_sidebar = eval(parse(text = paste("list(", 
                                                               manual_annotation_colors_sidebar, ")")))
  }
  if (nrow(data) < 3) {
    cluster_rows = FALSE
    cluster_cols = FALSE
  }
  if (ncol(data) < 3) {
    cluster_cols = FALSE
    cluster_rows = FALSE
  }
  if (nrow(data) > 5000 & correlation_plot == "None") {
    cluster_rows = FALSE
  }
  if (ncol(data) > 5000 & correlation_plot == "None") {
    cluster_cols = FALSE
  }
  cluster_rows_results = cluster_rows
  cluster_cols_results = cluster_cols
  row_order = rownames(data)
  col_order = colnames(data)
  if (cluster_rows) {
    if (clustering_distance_rows == "pearson") {
      if (!cor_data) {
        row_cor = cor(t(data))
      }
      else {
        row_cor = data
      }
      row_dist <- as.dist(1 - row_cor)
      if (any(is.na(row_cor))) {
        row_dist = dist(data)
      }
    }
    else if (clustering_distance_rows == "spearman") {
      if (!cor_data) {
        row_cor = cor(t(data), method = "spearman")
      }
      else {
        row_cor = data
      }
      row_dist <- as.dist(1 - row_cor)
      if (any(is.na(row_cor))) {
        row_dist = dist(data)
      }
    }
    else {
      if (!cor_data) {
        dist_method = c("euclidean", "manhattan", "maximum", 
                        "canberra", "binary", "minkowski")
        if (clustering_distance_rows %in% dist_method) {
          row_dist = dist(data, method = clustering_distance_rows)
        }
        else {
          row_dist = vegdist(data, method = clustering_distance_rows)
        }
      }
      else {
        row_cor = data
        row_dist <- as.dist(1 - row_cor)
        if (any(is.na(row_cor))) {
          row_dist = dist(data)
        }
      }
    }
    cluster_rows_results = hclust(row_dist, method = clustering_method)
    if (sp.is.null(cluster_rows_variable)) {
      sv = svd(data)$v[, 1]
    }
    else {
      sv = annotation_row[[cluster_rows_variable]]
      if (remove_cluster_rows_variable_in_annorow) {
        annotation_row[[cluster_rows_variable]] <- NULL
      }
      if (length(annotation_row) == 0) {
        annotation_row = NULL
      }
    }
    dend = reorder(as.dendrogram(cluster_rows_results), 
                   wts = sv)
    cluster_rows_results <- as.hclust(dend)
    row_order = cluster_rows_results$order
  }
  if (cluster_cols) {
    if (clustering_distance_cols == "pearson") {
      if (!cor_data) {
        col_cor = cor(data)
      }
      else {
        col_cor = data
      }
      col_dist <- as.dist(1 - col_cor)
      if (any(is.na(col_cor))) {
        col_dist = dist(t(data))
      }
    }
    else if (clustering_distance_cols == "spearman") {
      if (!cor_data) {
        col_cor = cor(data, method = "spearman")
      }
      else {
        col_cor = data
      }
      col_dist <- as.dist(1 - col_cor)
      if (any(is.na(col_cor))) {
        col_dist = dist(t(data))
      }
    }
    else {
      if (!cor_data) {
        dist_method = c("euclidean", "manhattan", "maximum", 
                        "canberra", "binary", "minkowski")
        if (clustering_distance_cols %in% dist_method) {
          col_dist = dist(t(data), method = clustering_distance_cols)
        }
        else {
          col_dist = vegdist(t(data), method = clustering_distance_cols)
        }
      }
      else {
        col_cor = data
        col_dist <- as.dist(1 - col_cor)
        if (any(is.na(col_cor))) {
          col_dist = dist(t(data))
        }
      }
    }
    cluster_cols_results = hclust(col_dist, method = clustering_method)
    if (sp.is.null(cluster_cols_variable)) {
      sv = svd(data)$v[, 1]
    }
    else {
      sv = annotation_col[[cluster_cols_variable]]
      if (remove_cluster_cols_variable_in_annocol) {
        annotation_col[[cluster_cols_variable]] <- NULL
      }
      if (length(annotation_col) == 0) {
        annotation_col = NULL
      }
    }
    dend = reorder(as.dendrogram(cluster_cols_results), 
                   wts = sv)
    cluster_cols_results <- as.hclust(dend)
    col_order = cluster_cols_results$order
  }
  if (correlation_plot != "None") {
    if (cluster_rows) {
      cluster_cols_results = cluster_rows_results
      col_order = row_order
    }
    else if (cluster_cols) {
      cluster_rows_results = cluster_cols_results
      row_order = col_order
    }
  }
  if (!is.na(filename)) {
    data_order = data[row_order, col_order]
    data2 = data.frame(ID = rownames(data_order), data_order)
    write.table(data2, file = paste0(filename, ".reordered.txt"), 
                sep = "\t", quote = F, col.names = T, row.names = F)
  }
  gt <- pheatmap::pheatmap(data, kmean_k = NA, color = manual_color_vector, 
                           scale = scale, border_color = NA, cluster_rows = cluster_rows_results, 
                           cluster_cols = cluster_cols_results, cutree_rows = cutree_rows, 
                           cutree_cols = cutree_cols, kmeans_k = kclu, breaks = breaks, 
                           legend_breaks = legend_breaks, legend_labels = legend_labels, 
                           xtics_angle = xtics_angle, clustering_method = clustering_method, 
                           clustering_distance_rows = clustering_distance_rows, 
                           clustering_distance_cols = clustering_distance_cols, 
                           show_rownames = ytics, show_colnames = xtics, main = title, 
                           annotation_col = annotation_col, annotation_row = annotation_row, 
                           annotation_colors = manual_annotation_colors_sidebar, 
                           fontsize = fontsize, filename = filename, width = width, 
                           height = height, ...)
  print('here11')
  if (saveppt) {
    eoffice::topptx(gt, filename = paste0(filename, ".pptx"))
  }
  print('here22')
  return(gt)
}













