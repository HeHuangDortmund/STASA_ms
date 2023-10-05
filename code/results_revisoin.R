library(rprojroot)
library(data.table)
#library(imager)
root = find_root(is_git_root)
library(reticulate)
library(scanpathAggr)

#load data, evaluation with GLD
load(file.path(root, 'results_new', 'res_mit1003_minGLA_MM_Extra.Rda'))
load(file.path(root, 'results_new', 'res_osie_minGLA_MM_Extra.Rda'))

str(res_mit1003_minGLA_MM_Extra)

#load data, evaluation with scanMatch
load(file.path(root, 'results_new', 'res_mit1003_minScanMatch_MM_Extra.Rda'))
load(file.path(root, 'results_new', 'res_osie_minScanMatch_MM_Extra.Rda'))

#################################################################################
# Scasim
################################################################################

#--------------------------------------------------------------------------------
# mit1003
res_mit1003_minGLD_MM_Extra = do.call(rbind,res_mit1003_minGLA_MM_Extra)
res_mit1003_minGLD_MM_Extra = as.data.table(res_mit1003_minGLD_MM_Extra)
MultiMatch_Score_MIT1003_STASA_minGLD =  res_mit1003_minGLD_MM_Extra[,-(1:2)]

GLD_Score_MIT1003_STASA_minGLD =  res_mit1003_minGLD_MM_Extra[,c(1:2)]

MultiMatch_Score_MIT1003_STASA_minGLD$stimulusID = 1:(dim(res_mit1003_minGLD_MM_Extra)[1])
GLD_Score_MIT1003_STASA_minGLD$stimulusID = 1:(dim(res_mit1003_minGLD_MM_Extra)[1])
names(GLD_Score_MIT1003_STASA_minGLD) = c('STASA', 'minScasim', 'stimulusID')


GLD_Score_MIT1003_STASA_minGLD_long = melt(GLD_Score_MIT1003_STASA_minGLD, id.vars = 'stimulusID')


MultiMatch_Score_MIT1003_STASA = MultiMatch_Score_MIT1003_STASA_minGLD[,c(11, 1:5)]
names(MultiMatch_Score_MIT1003_STASA) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_MIT1003_STASA_long = melt(MultiMatch_Score_MIT1003_STASA, id.vars = 'stimulusID')
names(MultiMatch_Score_MIT1003_STASA_long)[2] = 'multiMatch'
MultiMatch_Score_MIT1003_STASA_long$method = 'STASA'

MultiMatch_Score_MIT1003_minGLD = MultiMatch_Score_MIT1003_STASA_minGLD[,c(11, 6:10)]
names(MultiMatch_Score_MIT1003_minGLD) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_MIT1003_minGLD_long = melt(MultiMatch_Score_MIT1003_minGLD, id.vars = 'stimulusID')
names(MultiMatch_Score_MIT1003_minGLD_long)[2] = 'multiMatch'
MultiMatch_Score_MIT1003_minGLD_long$method = 'minScasim'


MultiMatch_Score_MIT1003_STASA_minGLD_long = rbind(MultiMatch_Score_MIT1003_STASA_long, MultiMatch_Score_MIT1003_minGLD_long)





#-----------------------------------------------------------------------------
#osie
res_osie_minGLD_MM_Extra = do.call(rbind,res_osie_minGLA_MM_Extra)
res_osie_minGLD_MM_Extra = as.data.table(res_osie_minGLD_MM_Extra)
MultiMatch_Score_OSIE_STASA_minGLD =  res_osie_minGLD_MM_Extra[,-(1:2)]
GLD_Score_OSIE_STASA_minGLD =  res_osie_minGLD_MM_Extra[,c(1:2)]

MultiMatch_Score_OSIE_STASA_minGLD$stimulusID = 1:(dim(res_osie_minGLD_MM_Extra)[1])
GLD_Score_OSIE_STASA_minGLD$stimulusID = 1:(dim(res_osie_minGLD_MM_Extra)[1])
names(GLD_Score_OSIE_STASA_minGLD) = c('STASA', 'minScasim', 'stimulusID')

GLD_Score_OSIE_STASA_minGLD_long = melt(GLD_Score_OSIE_STASA_minGLD, id.vars = 'stimulusID')


MultiMatch_Score_OSIE_STASA = MultiMatch_Score_OSIE_STASA_minGLD[,c(11, 1:5)]
names(MultiMatch_Score_OSIE_STASA) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_OSIE_STASA_long = melt(MultiMatch_Score_OSIE_STASA, id.vars = 'stimulusID')
names(MultiMatch_Score_OSIE_STASA_long)[2] = 'multiMatch'
MultiMatch_Score_OSIE_STASA_long$method = 'STASA'

MultiMatch_Score_OSIE_minGLD = MultiMatch_Score_OSIE_STASA_minGLD[,c(11, 6:10)]
names(MultiMatch_Score_OSIE_minGLD) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_OSIE_minGLD_long = melt(MultiMatch_Score_OSIE_minGLD, id.vars = 'stimulusID')
names(MultiMatch_Score_OSIE_minGLD_long)[2] = 'multiMatch'
MultiMatch_Score_OSIE_minGLD_long$method = 'minScasim'


MultiMatch_Score_OSIE_STASA_minGLD_long = rbind(MultiMatch_Score_OSIE_STASA_long, MultiMatch_Score_OSIE_minGLD_long)


#################################################################################
# ScanMatch
################################################################################

#--------------------------------------------------------------------------------
# mit1003
res_mit1003_minScanMatch_MM_Extra = do.call(rbind,res_mit1003_minScanMatch_MM_Extra)


res_mit1003_minScanMatch_MM_Extra = as.data.table(res_mit1003_minScanMatch_MM_Extra)
#res_mit1003_minScanMatch_MM_Extra = res_mit1003_minScanMatch_MM_Extra[, 1:12] # fine
res_mit1003_minScanMatch_MM_Extra = res_mit1003_minScanMatch_MM_Extra[, 13:24] # grob


MultiMatch_Score_MIT1003_STASA_minScanMatch =  res_mit1003_minScanMatch_MM_Extra[,-(1:2)]
ScanMatch_Score_MIT1003_STASA_minScanMatch =  res_mit1003_minScanMatch_MM_Extra[,c(1:2)]




MultiMatch_Score_MIT1003_STASA_minScanMatch$stimulusID = 1:(dim(res_mit1003_minScanMatch_MM_Extra)[1])
ScanMatch_Score_MIT1003_STASA_minScanMatch$stimulusID = 1:(dim(res_mit1003_minScanMatch_MM_Extra)[1])
names(ScanMatch_Score_MIT1003_STASA_minScanMatch) = c('STASA', 'maxScanMatch', 'stimulusID')

ScanMatch_Score_MIT1003_STASA_minScanMatch_long = melt(ScanMatch_Score_MIT1003_STASA_minScanMatch, id.vars = 'stimulusID')


MultiMatch_Score_MIT1003_STASA = MultiMatch_Score_MIT1003_STASA_minScanMatch[,c(11, 1:5)]
names(MultiMatch_Score_MIT1003_STASA) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_MIT1003_STASA_long = melt(MultiMatch_Score_MIT1003_STASA, id.vars = 'stimulusID')
names(MultiMatch_Score_MIT1003_STASA_long)[2] = 'multiMatch'
MultiMatch_Score_MIT1003_STASA_long$method = 'STASA'

MultiMatch_Score_MIT1003_minScanMatch = MultiMatch_Score_MIT1003_STASA_minScanMatch[,c(11, 6:10)]
names(MultiMatch_Score_MIT1003_minScanMatch) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_MIT1003_minScanMatch_long = melt(MultiMatch_Score_MIT1003_minScanMatch, id.vars = 'stimulusID')
names(MultiMatch_Score_MIT1003_minScanMatch_long)[2] = 'multiMatch'
MultiMatch_Score_MIT1003_minScanMatch_long$method = 'maxScanMatch'


MultiMatch_Score_MIT1003_STASA_minScanMatch_long = rbind(MultiMatch_Score_MIT1003_STASA_long, MultiMatch_Score_MIT1003_minScanMatch_long)





#-----------------------------------------------------------------------------
#osie
res_osie_minScanMatch_MM_Extra = do.call(rbind,res_osie_minScanMatch_MM_Extra)
res_osie_minScanMatch_MM_Extra = as.data.table(res_osie_minScanMatch_MM_Extra)
#res_osie_minScanMatch_MM_Extra = res_osie_minScanMatch_MM_Extra[, 1:12] # fine
res_osie_minScanMatch_MM_Extra = res_osie_minScanMatch_MM_Extra[, 13:24] # grob



MultiMatch_Score_OSIE_STASA_minScanMatch =  res_osie_minScanMatch_MM_Extra[,-(1:2)]
ScanMatch_Score_OSIE_STASA_minScanMatch =  res_osie_minScanMatch_MM_Extra[,c(1:2)]

MultiMatch_Score_OSIE_STASA_minScanMatch$stimulusID = 1:(dim(res_osie_minScanMatch_MM_Extra)[1])
ScanMatch_Score_OSIE_STASA_minScanMatch$stimulusID = 1:(dim(res_osie_minScanMatch_MM_Extra)[1])
names(ScanMatch_Score_OSIE_STASA_minScanMatch) = c('STASA', 'maxScanMatch', 'stimulusID')

ScanMatch_Score_OSIE_STASA_minScanMatch_long = melt(ScanMatch_Score_OSIE_STASA_minScanMatch, id.vars = 'stimulusID')


MultiMatch_Score_OSIE_STASA = MultiMatch_Score_OSIE_STASA_minScanMatch[,c(11, 1:5)]
names(MultiMatch_Score_OSIE_STASA) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_OSIE_STASA_long = melt(MultiMatch_Score_OSIE_STASA, id.vars = 'stimulusID')
names(MultiMatch_Score_OSIE_STASA_long)[2] = 'multiMatch'
MultiMatch_Score_OSIE_STASA_long$method = 'STASA'

MultiMatch_Score_OSIE_minScanMatch = MultiMatch_Score_OSIE_STASA_minScanMatch[,c(11, 6:10)]
names(MultiMatch_Score_OSIE_minScanMatch) = c('stimulusID','vector',	'direction',	'length',	'position',	'duration') 
MultiMatch_Score_OSIE_minScanMatch_long = melt(MultiMatch_Score_OSIE_minScanMatch, id.vars = 'stimulusID')
names(MultiMatch_Score_OSIE_minScanMatch_long)[2] = 'multiMatch'
MultiMatch_Score_OSIE_minScanMatch_long$method = 'maxScanMatch'


MultiMatch_Score_OSIE_STASA_minScanMatch_long = rbind(MultiMatch_Score_OSIE_STASA_long, MultiMatch_Score_OSIE_minScanMatch_long)




################################################################################
# merge two data sets
################################################################################


MultiMatch_Score_MIT1003_STASA_minGLD_long$dataset = 'MIT1003'
MultiMatch_Score_OSIE_STASA_minGLD_long$dataset = 'OSIE'
MultiMatch_Score_STASA_minGLD_long = rbind(MultiMatch_Score_MIT1003_STASA_minGLD_long, MultiMatch_Score_OSIE_STASA_minGLD_long)

MultiMatch_Score_MIT1003_STASA_minScanMatch_long$dataset = 'MIT1003'
MultiMatch_Score_OSIE_STASA_minScanMatch_long$dataset = 'OSIE'
MultiMatch_Score_STASA_minScanMatch_long = rbind(MultiMatch_Score_MIT1003_STASA_minScanMatch_long, MultiMatch_Score_OSIE_STASA_minScanMatch_long)


GLD_Score_MIT1003_STASA_minGLD_long$dataset = 'MIT1003'
GLD_Score_OSIE_STASA_minGLD_long$dataset = 'OSIE'
GLD_Score_STASA_minGLD_long = rbind(GLD_Score_MIT1003_STASA_minGLD_long, GLD_Score_OSIE_STASA_minGLD_long)

ScanMatch_Score_MIT1003_STASA_minScanMatch_long$dataset = 'MIT1003'
ScanMatch_Score_OSIE_STASA_minScanMatch_long$dataset = 'OSIE'
ScanMatch_Score_STASA_minScanMatch_long = rbind(ScanMatch_Score_MIT1003_STASA_minScanMatch_long, ScanMatch_Score_OSIE_STASA_minScanMatch_long)

names(GLD_Score_STASA_minGLD_long)[2] = 'method'
names(ScanMatch_Score_STASA_minScanMatch_long)[2] = 'method'


library(ggplot2)
p1 <-
ggplot(data = MultiMatch_Score_STASA_minGLD_long,
       aes(x = multiMatch, y = value, fill = method))+
  geom_boxplot()+
  facet_wrap(~dataset)+
  ylab('score')+
  xlab('multiMatch variables')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

p2<-
ggplot(data = MultiMatch_Score_STASA_minScanMatch_long,
       aes(x = multiMatch, y = value, fill = method))+
  geom_boxplot()+
  facet_wrap(~dataset)+
  ylab('score')+
  xlab('multiMatch Variables')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




p3<-
ggplot(data = GLD_Score_STASA_minGLD_long,
       aes( y = value, x = method))+
  geom_boxplot()+
  facet_wrap(~dataset)+
  ggtitle('a)')+
  ylab('Average Scasim Distance')
  

#ScanMatch_Score_STASA_minScanMatch_long$method[ScanMatch_Score_STASA_minScanMatch_long$method == 'minScanMatch'] = 'maxScanMatch'

p4<-
ggplot(data = ScanMatch_Score_STASA_minScanMatch_long,
       aes( y = value, x = method))+
  geom_boxplot()+

  facet_wrap(~dataset)+
  ggtitle('b)')+
  ylab('Average ScanMatch Similarity')








t.test(value ~ method, data = ScanMatch_Score_STASA_minScanMatch_long[dataset == 'MIT1003'])
t.test(value ~ method, data = ScanMatch_Score_STASA_minScanMatch_long[dataset == 'OSIE'])

t.test(value ~ method, data = GLD_Score_STASA_minGLD_long[dataset == 'MIT1003'])
t.test(value ~ method, data = GLD_Score_STASA_minGLD_long[dataset == 'OSIE'])

library(lsr)
cohensD(value ~ method, data = ScanMatch_Score_STASA_minScanMatch_long[dataset == 'MIT1003'])
cohensD(value ~ method, data = ScanMatch_Score_STASA_minScanMatch_long[dataset == 'OSIE'])

cohensD(value ~ method, data = GLD_Score_STASA_minGLD_long[dataset == 'MIT1003'])
cohensD(value ~ method, data = GLD_Score_STASA_minGLD_long[dataset == 'OSIE'])


library(gridExtra)
grid.arrange(p3,p4, nrow = 1)  # Figure 9 in ms
grid.arrange(p1,p2, nrow = 1)



STASAminScanMatch = MultiMatch_Score_STASA_minScanMatch_long[method == 'STASA']
STASAminGLD = MultiMatch_Score_STASA_minGLD_long[method == 'STASA']

minScanMatch = MultiMatch_Score_STASA_minScanMatch_long[method == 'maxScanMatch']
minGLD = MultiMatch_Score_STASA_minGLD_long[method == 'minScasim']

STASAminScanMatch$method = 'STASA (to maximise ScanMatch)'
STASAminGLD$method = 'STASA (to minimise Scasim)'



#STASA_ = rbind(STASAminScanMatch, STASAminGLD, minScanMatch, minGLD)

STASA_ = rbind(STASAminScanMatch, STASAminGLD)

t.test(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'vector'])
t.test(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'direction'])
t.test(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'length'])
t.test(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'position'])
t.test(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'duration'])


cohensD(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'vector'])
cohensD(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'direction'])
cohensD(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'length'])
cohensD(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'position'])
cohensD(value ~ method, STASA_[dataset == 'MIT1003' & multiMatch == 'duration'])






t.test(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'vector'])
t.test(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'direction'])
t.test(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'length'])
t.test(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'position'])
t.test(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'duration'])


cohensD(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'vector'])
cohensD(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'direction'])
cohensD(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'length'])
cohensD(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'position'])
cohensD(value ~ method, STASA_[dataset == 'OSIE' & multiMatch == 'duration']) # 0,18 sehr gerenge unterschied





  ggplot(data = STASA_,
         aes(x = multiMatch, y = value, fill = method))+
  geom_boxplot()+
  facet_wrap(~dataset)+
  ylab('score')+
  xlab('multiMatch Variables')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="bottom")


  
  ggplot(data = STASA_,
         aes(x = multiMatch, y = value, fill = method))+
    geom_boxplot()+
    facet_wrap(~dataset)+
    ylab('score')+
    xlab('multiMatch Variables')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position="top")

  