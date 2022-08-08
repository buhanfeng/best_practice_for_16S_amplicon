

#
# set enviroment
#
working_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/batch5_all'
base_on='flt_unoise3_whole_30_4'


lefse_analysis_path=paste(working_path, base_on, 'lefse_analysis', sep = '/')
if(dir.exists(lefse_analysis_path)){
  file.remove(dir(lefse_analysis_path , pattern="*"))
}else{
  dir.create(lefse_analysis_path)
}

#                                                                                                                                                                     
# library
#
library(sqldf)
library(xlsx)
library(stringr)

#
# init data
#
# if(base_on == 'unoise3'){
#   merge_tab_path = paste0(working_path, base_on, '/otutab_unoise3/merge_tab.txt')
# }else if(base_on == 'dada2'){
#   merge_tab_path = paste0(working_path, base_on, 'merge_tab.txt')
# }else if(base_on == 'uparse'){
#   merge_tab_path = paste0(working_path, base_on, 'otutab_uparse/merge_tab.txt')
# }
merge_tab_path = paste0(working_path, '/', base_on, '/post_tab_tax/merge_tab.txt')
tax_path = paste0(working_path, '/', base_on, '/post_tab_tax/taxonomy.tsv')
manifest_path = paste0(working_path, "/manifest.tsv")

merge_tab = read.csv(merge_tab_path, sep = '\t')
tax = read.csv(tax_path, sep = '\t')
manifest = read.csv(manifest_path, sep = '\t')

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
id=manifest$sampleid
id=c("BAT05_M14","BAT05_M15","BAT05_M16","BAT05_M17","BAT05_M18","BAT05_M19","BAT05_W20","BAT05_W21","BAT05_W22","BAT05_W23","BAT05_W24","BAT05_W25")
idp=str_c(id, collapse=',')
init_tab = sqldf(paste0("select OTU_ID,", idp,",Domain,Phylum,Class,Ord,Family,Genus,Species from merge_tab, tax_tab where OTU_ID=Feature_ID order by OTU_ID"))

#
# calculate relative abundance 
#
# total = 0
# for(ele in id){
#   total = total+ sum(init_tab[ele])
# }
for(ele in id){
  ss = sum(init_tab[ele])
  init_tab[ele] = init_tab[ele]/ss
}


#
# transform table (join str)
#
sample_tab = subset(init_tab, select = -c(OTU_ID,Domain,Phylum,Class,Ord,Family,Genus,Species))
tax_tab = subset(init_tab, select = c(Domain,Phylum,Class,Ord,Family,Genus,Species))
for(i in c(2:7)){
  tax_tab_s = tax_tab[ , 1:i]
  ro = nrow(tax_tab_s)
  co = ncol(tax_tab_s)
  joi = c()
  for (j in c(1:ro)) {
    roo = tax_tab_s[j, ]
    ros = str_c(roo, collapse='|')
    joi = append(joi, ros)
  }
  joi = data.frame(joi)
  if(i ==2){
    Phylum_tab = cbind(joi, sample_tab)
  }else if(i ==3){
    Class_tab = cbind(joi, sample_tab)
  }else if(i ==4){
    Ord_tab = cbind(joi, sample_tab)
  }else if(i ==5){
    Family_tab = cbind(joi, sample_tab)
  }else if(i ==6){
    Genus_tab = cbind(joi, sample_tab)
  }else if(i ==7){
    Species_tab = cbind(joi, sample_tab)
  }
}
#domain
joi = tax_tab[, 'Domain']
Domain_tab = cbind(joi, sample_tab)
# OTU_tab
OTU_tab = subset(init_tab, select = -c(Domain,Phylum,Class,Ord,Family,Genus,Species))
library(dplyr)
OTU_tab=rename(OTU_tab, joi=OTU_ID)



#
# group table
#
idd=""
for(ele in id){
  uni = ", sum(##) as ##"
  idd = paste0(idd, str_replace_all(uni, '##', ele))
}
idp=str_c(id, collapse='+')
minn='0'
maxx='1000000'
condi = paste0("(", idp, ")>", minn," and (", idp, ")<",maxx,"")

sta = paste0("select joi", idd, " from OTU_tab where ", condi, " group by joi")
OTU_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Domain_tab where ", condi, " group by joi")
Domain_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Phylum_tab where ", condi, " group by joi")
Phylum_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Class_tab where ", condi, " group by joi")
Class_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Ord_tab where ", condi, " group by joi")
Ord_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Family_tab where ", condi, " group by joi")
Family_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Genus_tab where ", condi, " group by joi")
Genus_tab = sqldf(sta)

sta = paste0("select joi", idd, " from Species_tab where ", condi, " group by joi")
Species_tab = sqldf(sta)

#
# merge table
#
group1 = c('strain', 'M', 'M', 'M', 'M', 'M', 'M', 'W', 'W', 'W', 'W', 'W', 'W')
lefse_tab = rbind(Domain_tab,Phylum_tab,Class_tab,Ord_tab,Family_tab,Genus_tab,Species_tab)
lefse_tab = rename(lefse_tab, subject_id=joi)
nn = colnames(lefse_tab)
lefse_tab = rbind(group1, nn, lefse_tab[1: nrow(lefse_tab), ])
sfile_path = paste0(lefse_analysis_path, '/lefse_tab.tsv')
write.table(lefse_tab, file = sfile_path, row.names = F, col.names = F, sep = '\t', quote = F)






