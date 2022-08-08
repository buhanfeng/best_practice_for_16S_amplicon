

#
# set enviroment
#
working_path='D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/whole_gut_all_genus_beta'
base_on='flt_unoise3_200_4'


levels_tab_path=paste(working_path, base_on, 'levels_tab', sep = '/')
if(dir.exists(levels_tab_path)){
  file.remove(dir(levels_tab_path , pattern="*"))
}else{
  dir.create(levels_tab_path)
}

#                                                                                                                                                                     
# library
#
library(sqldf)
library(xlsx)
library(stringr)


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
id=manifest$id
idp=str_c(id, collapse=',')
init_tab = sqldf(paste0("select OTU_ID,", idp,",Domain,Phylum,Class,Ord,Family,Genus,Species from merge_tab, tax_tab where OTU_ID=Feature_ID order by OTU_ID"))

type_c = c("Domain","Phylum","Class","Ord","Family","Genus","Species")

for (type in type_c) {
  sta = "select ##type as OTU_ID,##idp from init_tab where ##type!='unknow'"
  sta = str_replace_all(sta, '##idp', idp)
  sta = str_replace_all(sta, '##type', type)
  init_tab_s = sqldf(sta)
  sfile_path = paste0(levels_tab_path, '/', type,'_tab.tsv')
  write.table(init_tab_s, file = sfile_path, row.names = F, col.names = T, sep = '\t', quote = F)
}



