


path = 'D:/desktop/advance/projects/zebrafish_gut_microbiome/microbiome/deliver/Åú´Î1-5 lefse»ã×Ü/bat5_EC_S2_Not_repeat_WT_vs_MU_E_vs_C'

out = read.csv(file = paste(path, 'lefse_-c_1__-s_2__-u_3__-o_1000000_-a_0.05_-w_0.2_-s_1_--min_c_10.res.read', sep = '/'), sep = '\t')

inin = read.csv(file = paste(path, 'lefse_tab.tsv.read', sep = '/'), sep = '\t')

library(sqldf)

merge_tab = sqldf("select name,LDA,treat,LDA_adj,p_value,BAT05_M20E,BAT05_M20C,BAT05_M21E,BAT05_M21C,BAT05_M22E,BAT05_M22C,BAT05_W26E,BAT05_W26C,BAT05_W27E,BAT05_W27C,BAT05_W28E,BAT05_W28C from out,inin where inin.subject_id=out.name")

write.table(merge_tab, file = paste(path, 'merge_tab.tsv', sep =  '/'), row.names = F, quote = F)
