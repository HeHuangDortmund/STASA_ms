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




aggregattionAndminScanmatch = function(singleStimulus, OR_p = 0.2, OR_p_grob){
  
  Fixation_OSIE_singleStimulus = Fixation_MIT1003[Stimulus == singleStimulus]
  spLists = split(Fixation_OSIE_singleStimulus, by = 'Participant')
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  
  convex_hull = getConcaveHulls(Fixation_OSIE_singleStimulus,
                                ExperimentDuration,
                                et)
  
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  
  
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  
  
  
  #stimulusImgPath = file.path(root,'data','OSIE','master','data','stimuli',singleStimulus)
  stimulusImgPath = file.path(root,'data','MIT1003','ALLSTIMULI',paste0(singleStimulus,'.jpeg'))
  img = load.image(stimulusImgPath)
  
  screensize = dim(img)[1:2]
  Xres = screensize[1]
  Yres = screensize[2]
  
  
  #  24x18 Grid
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
  
  # needlemanOnlyScore = function(x){
  #   needleman(x, seq2 = RepresentScanpathSequenz, gap = 0, solutionMatrix = solutionMatrix)$score
  # }
  # 
  
  
  
  needleman_normScore = function(x1, x2){
    needleman(x1, x2, gap = 0, solutionMatrix = solutionMatrix)$score/(max(length(x1),length(x2)))/(max(solutionMatrix))
  }
  
  needlemanOnlyScore_prenorm = function(x){
    needleman_normScore(x, x2 = RepresentScanpathSequenz)
  }
  
  
  
  #meanScore = mean(unlist(lapply(scanPathSequenzList, needlemanOnlyScore)))
  
  #Normalized_score = meanScore/(max(solutionMatrix))/maxLength
  Normalized_score_pre = mean(unlist(lapply(scanPathSequenzList, needlemanOnlyScore_prenorm)))
  # minimize scanMatch
  
  # scanMatch_distMatrix = matrix(0, nc = 15, nr = 15)
  # for (i in 1:length(scanPathSequenzList)) {
  #   for (j in 1:length(scanPathSequenzList)) {
  #     scanMatch_distMatrix[i,j] =  needleman(scanPathSequenzList[[i]], seq2 = scanPathSequenzList[[j]], gap = 0, solutionMatrix = solutionMatrix)$score
  #   }
  # }
  
  scanMatch_distMatrix_norm_pre = matrix(0, nc = 15, nr = 15)
  for (i in 1:length(scanPathSequenzList)) {
    for (j in 1:length(scanPathSequenzList)) {
      scanMatch_distMatrix_norm_pre[i,j] =  needleman_normScore(scanPathSequenzList[[i]],  scanPathSequenzList[[j]])
    }
  }
  
  ind_pre = which.max(colSums(scanMatch_distMatrix_norm_pre))
  
  scanMatch_score_max_pre = (sum(scanMatch_distMatrix_norm_pre[,ind_pre]) - 1)/14
  
  #scanMatch_distMatrix_norm = scanMatch_distMatrix/(max(solutionMatrix))/maxLength
  #ind = which.max(colSums(scanMatch_distMatrix_norm))
  #scanMatch_score_max = (sum(scanMatch_distMatrix_norm[,ind]) - 1)/14
  
  #fine = c(Normalized_score, Normalized_score_pre, scanMatch_score_max, scanMatch_score_max_pre)
  
  
  
  
  ####################################################################
  ## compute the multimatch scores for both methods
  ###################################################################
  #Sys.setenv(RETICULATE_PYTHON = file.path(root,'python','Scripts'))
  Sys.setenv(RETICULATE_PYTHON = file.path(root,'my_env','bin',"python"))
  
  #reticulate::py_config()
  multimatch <- import('multimatch_gaze')
  stimulusImgPath = file.path(root,'data','OSIE','master','data','stimuli',singleStimulus)
  
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
  
  res_minScanmatch = do.call(rbind,lapply(FixVectorLists, function(x) unlist(multimatch$docomparison(FixVectorLists[[ind_pre]], x, screensize=screensize))))
  
  multimatch_STASA = colMeans(res, na.rm = TRUE)
  multimatch_minScanMatch = colMeans(res_minScanmatch[-ind_pre,], na.rm = TRUE)
  fine = c(Normalized_score_pre, scanMatch_score_max_pre, multimatch_STASA, multimatch_minScanMatch)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##############################################################################
  #  12x8 Grid  grob
  #-----------------------------------------------------------------------------
  
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p_grob)
  
  
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  
  
  
  
  
  Xbin = 12
  Ybin = 8
  source(file.path(root,'code','scanpathComparision','scanMatch.R'))
  source(file.path(root,'code','utilityFunctions','scanPathToSequence.R'))
  #scanPathToSequence2 = function(x){}
  scanPathSequenzList = lapply(spLists, scanPathToSequence, Xbin = Xbin, Ybin = Ybin, Xres = Xres, Yres = Yres)
  #maxLength = max(unlist(lapply(scanPathSequenzList, length))) # will be used for normalization of the match score
  
  RepresentScanpathSequenz = scanPathToSequence(RepresentScanpath$representScanpath, Xbin = Xbin, Ybin = Ybin, Xres = Xres, Yres = Yres)
  
  solutionMatrix = buildSolutionMatrix(Xbin, Ybin, Xres,Yres, threshold)
  
  # needlemanOnlyScore = function(x){
  #   needleman(x, seq2 = RepresentScanpathSequenz, gap = 0, solutionMatrix = solutionMatrix)$score
  # }
  # 
  
  needleman_normScore = function(x1, x2){
    needleman(x1, x2, gap = 0, solutionMatrix = solutionMatrix)$score/(max(length(x1),length(x2)))/(max(solutionMatrix))
  }
  
  needlemanOnlyScore_prenorm = function(x){
    needleman_normScore(x, x2 = RepresentScanpathSequenz)
  }
  
  
  
  #meanScore = mean(unlist(lapply(scanPathSequenzList, needlemanOnlyScore)))
  
  #grob_Normalized_score = meanScore/(max(solutionMatrix))/maxLength
  grob_Normalized_score_pre = mean(unlist(lapply(scanPathSequenzList, needlemanOnlyScore_prenorm)))
  # minimize scanMatch
  
  # scanMatch_distMatrix = matrix(0, nc = 15, nr = 15)
  # for (i in 1:length(scanPathSequenzList)) {
  #   for (j in 1:length(scanPathSequenzList)) {
  #     scanMatch_distMatrix[i,j] =  needleman(scanPathSequenzList[[i]], seq2 = scanPathSequenzList[[j]], gap = 0, solutionMatrix = solutionMatrix)$score
  #   }
  # }
  
  scanMatch_distMatrix_norm_pre = matrix(0, nc = 15, nr = 15)
  for (i in 1:length(scanPathSequenzList)) {
    for (j in 1:length(scanPathSequenzList)) {
      scanMatch_distMatrix_norm_pre[i,j] =  needleman_normScore(scanPathSequenzList[[i]],  scanPathSequenzList[[j]])
    }
  }
  
  ind_pre_grob = which.max(colSums(scanMatch_distMatrix_norm_pre))
  
  grob_scanMatch_score_max_pre = (sum(scanMatch_distMatrix_norm_pre[,ind_pre_grob]) - 1)/14
  
  #scanMatch_distMatrix_norm = scanMatch_distMatrix/(max(solutionMatrix))/maxLength
  #ind = which.max(colSums(scanMatch_distMatrix_norm))
  #grob_scanMatch_score_max = (sum(scanMatch_distMatrix_norm[,ind]) - 1)/14
  
  #grob = c(grob_Normalized_score, grob_Normalized_score_pre, grob_scanMatch_score_max, grob_scanMatch_score_max_pre)
  
  
  
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
  
  res_minScanmatch = do.call(rbind,lapply(FixVectorLists, function(x) unlist(multimatch$docomparison(FixVectorLists[[ind_pre_grob]], x, screensize=screensize))))
  
  grob_multimatch_STASA = colMeans(res, na.rm = TRUE)
  grob_multimatch_minScanMatch = colMeans(res_minScanmatch[-ind_pre_grob,], na.rm = TRUE)
  
  
  
  grob = c(grob_Normalized_score_pre, grob_scanMatch_score_max_pre, grob_multimatch_STASA, grob_multimatch_minScanMatch)
  
  
  
  return(c(fine, grob))
  
}

all_stimulus = unique(Fixation_MIT1003$Stimulus)


ORs = matrix(c(theta_mit1003_grob,
               theta_mit1003_fine),
             nc = 2)*0.1

indNA = which(is.na(ORs[,1]) | is.na(ORs[,2]))

ORs = ORs[-indNA,]
all_stimulus = all_stimulus[-indNA]
res_mit1003_minScanMatch_MM_Extra= list()

for (i in 1:length(all_stimulus)) {
  print(paste0('i: ',i))
  
  TMP = try(aggregattionAndminScanmatch(all_stimulus[i], ORs[i,1], ORs[i,2]))
  if(class(TMP) == 'try-error'){
    res_mit1003_minScanMatch_MM_Extra[[i]] = NA
  }else{
    print(paste0('score: ',TMP))
    res_mit1003_minScanMatch_MM_Extra[[i]] = TMP
  }
}


save(res_mit1003_minScanMatch_MM_Extra, file = file.path(root, 'results','res_mit1003_minScanMatch_MM_Extra.Rda'))#TODO

