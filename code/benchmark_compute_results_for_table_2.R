library(rprojroot)
library(data.table)
library(imager)
library(scanpathAggr)
root = find_root(is_git_root)
library(reticulate)


Sys.setenv(RETICULATE_PYTHON = file.path(root,'my_env','bin',"python")) # on linux
#Sys.setenv(RETICULATE_PYTHON = file.path(root,'python','Scripts')) # on windows

#---------1st: get Convex_Hulls

# load data Fixation_MIT1003
load(file = file.path(root, 'data','Fixation_MIT1003.Rda'))
Participants = c('ajs', 'CNG','emb','ems','ff','hp','jcw',
                 'jw','kae','krl','po','tmj','tu','ya','zb')
# -------------------------------------How many participants for each Stimulus?

# remove NA
# NA Quelle: for some Stimuli and Paticipant: there are no saccade
#            which means that the Paticipant didn't 'move' their eyes during the experiment

Fixation_MIT1003 = Fixation_MIT1003[!is.na(Duration),]
ParticipantAndStimulus = unique(Fixation_MIT1003[, c('Participant','Stimulus')])
ParticipantPerStimulus = ParticipantAndStimulus[,.N, by = Stimulus]
indexOfIncompleteStimulus = which(ParticipantPerStimulus$N < 15)
IncompleteStimulus = ParticipantPerStimulus[indexOfIncompleteStimulus,]$Stimulus
Fixation_MIT1003 = Fixation_MIT1003[!Stimulus%in%IncompleteStimulus,]


# take a single stimulus
#singleStimulus = 'i05june05_static_street_boston_p1010764'


#singleStimulus = all_stimulus[557]


aggregattionAndMultimatch = function(singleStimulus, OR_p = 0.22 ){
  Fixation_MIT1003_singleStimulus = Fixation_MIT1003[Stimulus == singleStimulus]
  
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  #OR_p = 0.22 # 0.22
  
  convex_hull = getConcaveHulls(Fixation_MIT1003_singleStimulus,
                               ExperimentDuration,
                               et)
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  
  multimatch <- import('multimatch_gaze')
  stimulusImgPath = file.path(root,'data','MIT1003','ALLSTIMULI',paste0(singleStimulus,'.jpeg'))
  img = load.image(stimulusImgPath)
  
  screensize = dim(img)[1:2]
  
  fix_vector = RepresentScanpath$representScanpath[,c('Fixation_Position_X_px','Fixation_Position_Y_px','Duration')]
  names(fix_vector) = c('start_x',  'start_y', 'duration')
  
  spLists = split(Fixation_MIT1003_singleStimulus, by = 'Participant')
  
  toFixVectorForMultimatch = function(x){
    x = x[,c('Fixation_Position_X_px','Fixation_Position_Y_px','Duration')]
    names(x) = c('start_x',  'start_y', 'duration')
    return(x)
  }
  FixVectorLists = lapply(spLists, toFixVectorForMultimatch)
  
  res = do.call(rbind,lapply(FixVectorLists, function(x) unlist(multimatch$docomparison(fix_vector, x, screensize=screensize))))
  
  return(colMeans(res, na.rm = TRUE))
  
}



ORs = seq(0.1, 0.9, by = 0.1)
all_stimulus = unique(Fixation_MIT1003$Stimulus)


res_mit1003_multimatch = list()
res_mit1003_multimatch_all = list()
for (j in 1:length(ORs)) {
  res_block_all = list()
  
  for (i in 1:length(all_stimulus)) {
    print(paste0('j: ',j,'  i:', i))
    
    TMP = try(aggregattionAndMultimatch(all_stimulus[i], ORs[j]))
    if(class(TMP) == 'try-error'){
      res_block_all[[i]] = NA
    }else{
      res_block_all[[i]] = TMP
      print(TMP)
    }
    
    
  }
  
  
  #res_j = colMeans(do.call(rbind, res_block_all), na.rm = TRUE)
  #res_mit1003_multimatch = c(res_mit1003_multimatch, res_j)
  
  res_j = do.call(rbind, res_block_all)
  res_mit1003_multimatch_all[[j]] = res_j
  
}

#save(res_mit1003_multimatch, file = file.path(root, 'results_new',' res_mit1003_multimatch.Rda'))

save(res_mit1003_multimatch_all, file = file.path(root, 'results_new',' res_mit1003_multimatch_all.Rda'))

load(file = file.path(root, 'results_new',' res_mit1003_multimatch_all.Rda'))


colMeans(res_mit1003_multimatch_all[[2]], na.rm = TRUE) #results in Table 2


#--------------------------------------------------------------------------------------



#################MIT1003SCANMATCH###############################################
rm(list = ls())

library(rprojroot)
library(data.table)
library(imager)
library(reticulate)
library(scanpathAggr)

root = find_root(is_git_root)

source(file.path(root,'code','scanpathComparision','scanMatch.R'))
source(file.path(root,'code','utilityFunctions','scanPathToSequence.R'))


load(file = file.path(root, 'data','Fixation_MIT1003.Rda'))
Participants = c('ajs', 'CNG','emb','ems','ff','hp','jcw',
                 'jw','kae','krl','po','tmj','tu','ya','zb')
# -------------------------------------How many participants for each Stimulus?

# remove stimuli with missing eye-tracking data
# NA Quelle: for some Stimuli and Paticipant: there are no saccade
#            which means that the Paticipant didn't 'move' their eyes during the experiment
Fixation_MIT1003 = Fixation_MIT1003[!is.na(Duration),]
ParticipantAndStimulus = unique(Fixation_MIT1003[, c('Participant','Stimulus')])
ParticipantPerStimulus = ParticipantAndStimulus[,.N, by = Stimulus]
indexOfIncompleteStimulus = which(ParticipantPerStimulus$N < 15)
IncompleteStimulus = ParticipantPerStimulus[indexOfIncompleteStimulus,]$Stimulus
Fixation_MIT1003 = Fixation_MIT1003[!Stimulus%in%IncompleteStimulus,]


aggregattionAndScanmatch = function(singleStimulus, OR_p = 0.5){
  
  Fixation_MIT1003_singleStimulus = Fixation_MIT1003[Stimulus == singleStimulus]
  spLists = split(Fixation_MIT1003_singleStimulus, by = 'Participant')
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  
  convex_hull = getConcaveHulls(Fixation_MIT1003_singleStimulus,
                                ExperimentDuration,
                                et)
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  
  
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  
  stimulusImgPath = file.path(root,'data','MIT1003','ALLSTIMULI',paste0(singleStimulus,'.jpeg'))
  img = load.image(stimulusImgPath)
  
  screensize = dim(img)[1:2]
  Xres = screensize[1]
  Yres = screensize[2]
  Xbin = 24
  Ybin = 18
  threshold = 3.5
  
  scanPathSequenzList = lapply(spLists, scanPathToSequence, Xbin = Xbin, Ybin = Ybin, Xres = Xres, Yres = Yres)
  maxLength = max(unlist(lapply(scanPathSequenzList, length))) # will be used for normalization of the match score
  
  RepresentScanpathSequenz = scanPathToSequence(RepresentScanpath$representScanpath, Xbin = Xbin, Ybin = Ybin, Xres = Xres, Yres = Yres)
  
  solutionMatrix = buildSolutionMatrix(Xbin, Ybin, Xres,Yres, threshold)
  
  needlemanOnlyScore = function(x){
    if(length(x) == 0){  #length(x) is zero because of the missing data, see scanPathSequenzList[[12]]
      return(NA)  
    }else{
      needleman(x, seq2 = RepresentScanpathSequenz, gap = 0, solutionMatrix = solutionMatrix)$score 
    }
  }
  
  meanScore = mean(unlist(lapply(scanPathSequenzList, needlemanOnlyScore)), na.rm = TRUE)
  
  Normalized_score = meanScore/(max(solutionMatrix))/maxLength
  return(Normalized_score)
}

#warnings are form the r package rgeos by using the method union on the 
#two objects of klass gpc.polygon : union(gpc.polygon1, gpc.polygon2)


all_stimulus = unique(Fixation_MIT1003$Stimulus)
ORs = seq(0.1, 0.9, by = 0.1)

res_mit1003_scanmatch = list()
for (j in 1:length(ORs)) {
  
  res_j = list()
  for (i in 1:length(all_stimulus)) {
    print(paste0('j:',j,'----', 'i: ',i))
    
    TMP = try(aggregattionAndScanmatch(all_stimulus[i], ORs[j]))
    if(class(TMP) == 'try-error'){
      res_j[[i]] = NA
    }else{
      print(paste0('score: ',TMP))
      res_j[[i]] = TMP
    }
    
    
  }
  res= unlist(res_j)
  res_mit1003_scanmatch[[j]] = res
}

save(res_mit1003_scanmatch, file = file.path(root, 'results_new',' res_mit1003_scanMatch.Rda'))
#lapply(res_mit1003, mean, na.rm = TRUE)
lapply(res_mit1003_scanmatch, mean, na.rm = TRUE)

load(file = file.path(root, 'results_new',' res_mit1003_scanMatch.Rda'))

getOptForEachStimulus = function(i){
  
  score_i = c()
  for (j in 1:9) {
    score_i[j] = res_mit1003_scanmatch[[j]][i]
  }
  
  return(max(score_i, na.rm = TRUE))
  
}

optscores = unlist(lapply(1:972, getOptForEachStimulus))
mean(optscores) # 0.691  Result in table 2 in ms




#################OSIE_multiMATCH###############################################
rm(list = ls())
library(rprojroot)
library(data.table)
library(imager)
root = find_root(is_git_root)
library(reticulate)
library(scanpathAggr)

#Sys.setenv(RETICULATE_PYTHON = file.path(root,'python','Scripts')) # onwindows
Sys.setenv(RETICULATE_PYTHON = file.path(root,'my_env','bin',"python")) # on linux

# load data Fixation_MIT1003
load(file = file.path(root, 'data','Fixation_OSIE.Rda'))
Participants = c('P01', 'P02','P03','P04','P05','P06','P07',
                 'P08','P09','P10','P11','P12','P13','P14','P15')
# -------------------------------------How many participants for each Stimulus?

# remove NA
# NA Quelle: for some Stimuli and Paticipant: there are no saccade
#            which means that the Paticipant didn't 'move' their eyes during the experiment

Fixation_OSIE = Fixation_OSIE[!is.na(Duration),]
ParticipantAndStimulus = unique(Fixation_OSIE[, c('Participant','Stimulus')])
ParticipantPerStimulus = ParticipantAndStimulus[,.N, by = Stimulus]


aggregattionAndMultimatch = function(singleStimulus, OR_p = 0.22 ){
  Fixation_OSIE_singleStimulus = Fixation_OSIE[Stimulus == singleStimulus]
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  
  convex_hull = getConcaveHulls(Fixation_OSIE_singleStimulus,
                               ExperimentDuration,
                               et)
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_OSIE_singleStimulus)
  
  
  #Sys.setenv(RETICULATE_PYTHON = file.path(root,'python','Scripts'))
  Sys.setenv(RETICULATE_PYTHON = file.path(root,'my_env','bin',"python"))
  
  #reticulate::py_config()
  multimatch <- import('multimatch_gaze')
  stimulusImgPath = file.path(root,'data','OSIE','master','data','stimuli',singleStimulus)
  img = load.image(stimulusImgPath)
  screensize = dim(img)[1:2]
  fix_vector = RepresentScanpath$representScanpath[,c('Fixation_Position_X_px','Fixation_Position_Y_px','Duration')]
  names(fix_vector) = c('start_x',  'start_y', 'duration')
  
  spLists = split(Fixation_OSIE_singleStimulus, by = 'Participant')
  
  toFixVectorForMultimatch = function(x){
    x = x[,c('Fixation_Position_X_px','Fixation_Position_Y_px','Duration')]
    names(x) = c('start_x',  'start_y', 'duration')
    return(x)
  }
  FixVectorLists = lapply(spLists, toFixVectorForMultimatch)
  
  res = do.call(rbind,lapply(FixVectorLists, function(x) unlist(multimatch$docomparison(fix_vector, x, screensize=screensize))))
  
  return(colMeans(res, na.rm = TRUE))
  
}



ORs = seq(0.1, 0.9, by = 0.1)
all_stimulus = unique(Fixation_OSIE$Stimulus)

res_osie_multimatch = list()
res_osie_multimatch_all = list()
for (j in 1:length(ORs)) {
  res_block_all = list()
  
  for (i in 1:length(all_stimulus)) {
    print(paste0('j: ',j,'  i:', i))
    
    TMP = try(aggregattionAndMultimatch(all_stimulus[i], ORs[j]))
    if(class(TMP) == 'try-error'){
      res_block_all[[i]] = NA
    }else{
      res_block_all[[i]] = TMP
      print(TMP)
    }
    
    
  }
  
  
  #res_j = colMeans(do.call(rbind, res_block_all), na.rm = TRUE)
  #res_osie_multimatch = c(res_osie_multimatch, res_j)
  
  res_j = do.call(rbind, res_block_all)
  res_osie_multimatch_all[[j]] = res_j
  
  
}

#res_osie_multimatch = do.call(rbind, res_osie_multimatch)

#save(res_osie_multimatch, file = file.path(root, 'results_new',' res_osie_multimatch.Rda'))
save(res_osie_multimatch_all, file = file.path(root, 'results_new',' res_osie_multimatch_all.Rda'))

load(file = file.path(root, 'results_new',' res_osie_multimatch_all.Rda'))

colMeans(res_osie_multimatch_all[[2]], na.rm = TRUE) # results in table 2 in ms



#################OSIE_scanMATCH###############################################
rm(list = ls())
library(rprojroot)
library(data.table)
library(imager)
root = find_root(is_git_root)
library(reticulate)


# load data Fixation_MIT1003
load(file = file.path(root, 'data','Fixation_OSIE.Rda'))
Participants = c('P01', 'P02','P03','P04','P05','P06','P07',
                 'P08','P09','P10','P11','P12','P13','P14','P15')

#Fixation_OSIE
# remove NA
# NA Quelle: for some Stimuli and Paticipant: there are no saccade
#            which means that the Paticipant didn't 'move' their eyes during the experiment

Fixation_OSIE = Fixation_OSIE[!is.na(Duration),]
ParticipantAndStimulus = unique(Fixation_OSIE[, c('Participant','Stimulus')])
ParticipantPerStimulus = ParticipantAndStimulus[,.N, by = Stimulus]
#indexOfIncompleteStimulus = which(ParticipantPerStimulus$N < 15)
#IncompleteStimulus = ParticipantPerStimulus[indexOfIncompleteStimulus,]$Stimulus
#Fixation_OSIE = Fixation_OSIE[!Stimulus%in%IncompleteStimulus,]


aggregattionAndScanmatch = function(singleStimulus, OR_p = 0.5){
  
  Fixation_OSIE_singleStimulus = Fixation_OSIE[Stimulus == singleStimulus]
  spLists = split(Fixation_OSIE_singleStimulus, by = 'Participant')
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  
  convex_hull = getConcaveHulls(Fixation_OSIE_singleStimulus,
                               ExperimentDuration,
                               et)
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  
  
  
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_OSIE_singleStimulus,
                                           unionAOI = TRUE)
  
  stimulusImgPath = file.path(root,'data','OSIE','master','data','stimuli',singleStimulus)
  img = load.image(stimulusImgPath)
  
  screensize = dim(img)[1:2]
  Xres = screensize[1]
  Yres = screensize[2]
  Xbin = 24
  Ybin = 18
  threshold = 3.5
  
  source(file.path(root,'code','scanpathComparision','scanMatch.R'))
  source(file.path(root,'code','utilityFunctions','scanPathToSequence.R'))
  #scanPathToSequence2 = function(x){}
  scanPathSequenzList = lapply(spLists, scanPathToSequence, Xbin = Xbin, Ybin = Ybin, Xres = Xres, Yres = Yres)
  maxLength = max(unlist(lapply(scanPathSequenzList, length))) # will be used for normalization of the match score
  
  RepresentScanpathSequenz = scanPathToSequence(RepresentScanpath$representScanpath, Xbin = Xbin, Ybin = Ybin, Xres = Xres, Yres = Yres)
  
  solutionMatrix = buildSolutionMatrix(Xbin, Ybin, Xres,Yres, threshold)
  
  needlemanOnlyScore = function(x){
    needleman(x, seq2 = RepresentScanpathSequenz, gap = 0, solutionMatrix = solutionMatrix)$score
  }
  
  meanScore = mean(unlist(lapply(scanPathSequenzList, needlemanOnlyScore)))
  Normalized_score = meanScore/(max(solutionMatrix))/maxLength
  return(Normalized_score)
  
}



all_stimulus = unique(Fixation_OSIE$Stimulus)

ORs = seq(0.1, 0.9, by = 0.1)
res_osie_scanMatch= list()
for (j in 1:length(ORs)) {
  
  res_j = list()
  for (i in 1:length(all_stimulus)) {
    print(paste0('j:',j,'----', 'i: ',i))
    
    TMP = try(aggregattionAndScanmatch(all_stimulus[i], ORs[j]))
    if(class(TMP) == 'try-error'){
      res_j[[i]] = NA
    }else{
      print(paste0('score: ',TMP))
      res_j[[i]] = TMP
    }
    
    
  }
  res= unlist(res_j)
  res_osie_scanMatch[[j]] = res
}

save(res_osie_scanMatch, file = file.path(root, 'results_new',' res_osie_scanMatch.Rda'))
load(file = file.path(root, 'results_new','res_osie_scanMatch.Rda'))


lapply(res_osie_scanMatch, mean, na.rm = TRUE)


getOptForEachStimulus = function(i){
  
  score_i = c()
  for (j in 1:9) {
    score_i[j] = res_osie_scanMatch[[j]][i]
  }
  
  return(max(score_i, na.rm = TRUE))
  
}
optscores = unlist(lapply(1:700, getOptForEachStimulus))
mean(optscores) # 0.649   Result in table 2 in ms





