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



aggregattionAndLD = function(singleStimulus, OR_p = 0.22 ){
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
  
  
  
  
  rp_sp = data.frame(subject = "rs", trial = 1,
                     x =  RepresentScanpath$representScanpath$Fixation_Position_X_px,
                     y =  RepresentScanpath$representScanpath$Fixation_Position_Y_px,
                     duration =  RepresentScanpath$representScanpath$Duration*1000)
  
  
  scanpaths = data.frame(subject = Fixation_MIT1003_singleStimulus$Participant,
                         trial = 1,
                         x = Fixation_MIT1003_singleStimulus$Fixation_Position_X_px,
                         y = Fixation_MIT1003_singleStimulus$Fixation_Position_Y_px,
                         duration = Fixation_MIT1003_singleStimulus$Duration*1000)
  
  
  
  #stimulusImgPath = file.path(root,'data','OSIE','master','data','stimuli',singleStimulus)
  stimulusImgPath = file.path(root,'data','MIT1003','ALLSTIMULI',paste0(singleStimulus,'.jpeg'))
  img = load.image(stimulusImgPath)
  
  screensize = dim(img)[1:2]
  Xres = screensize[1]
  Yres = screensize[2]
  
  score_ld = mean(scasim(scanpaths, duration ~ x + y | subject, Xres/2, Yres/2, 60, 1/60, data2 = rp_sp))
  #score_ld = mean(scasim(scanpaths, duration ~ x + y | subject, Xres/2, Yres/2, 60, 1/60, data2 = rp_sp))
  
  d = scasim(scanpaths, duration ~ x + y | subject, Xres/2, Yres/2, 60, 1/60)
  i = which.min(colSums(d))
  score_ld_ld = colSums(d)[i]/(length(Participants) - 1)
  
  
  
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
  spLists = split(Fixation_MIT1003_singleStimulus, by = 'Participant')
  
  toFixVectorForMultimatch = function(x){
    x = x[,c('Fixation_Position_X_px','Fixation_Position_Y_px','Duration')]
    names(x) = c('start_x',  'start_y', 'duration')
    return(x)
  }
  FixVectorLists = lapply(spLists, toFixVectorForMultimatch)
  
  res = do.call(rbind,lapply(FixVectorLists, function(x) unlist(multimatch$docomparison(fix_vector, x, screensize=screensize))))
  
  res_minScanmatch = do.call(rbind,lapply(FixVectorLists, function(x) unlist(multimatch$docomparison(FixVectorLists[[i]], x, screensize=screensize))))
  
  multimatch_STASA = colMeans(res, na.rm = TRUE)
  multimatch_minScanMatch = colMeans(res_minScanmatch[-i,], na.rm = TRUE)
  
  
  return(c(score_ld, score_ld_ld, multimatch_STASA, multimatch_minScanMatch))
}
#aggregattionAndLD(singleStimulus)
ORs = min_value_theta_mit1003_gld[min_value_theta_mit1003_gld$variable == "STASA",]$V1*0.1
all_stimulus = unique(Fixation_MIT1003$Stimulus)

res_mit1003_minGLA_MM_Extra = list()
for (i in 1:length(all_stimulus)) {
  print(paste0('i: ',i))
  
  TMP = try(aggregattionAndLD(all_stimulus[i], ORs[i]))
  if(class(TMP) == 'try-error'){
    res_mit1003_minGLA_MM_Extra[[i]] = NA
  }else{
    print(paste0('score: ',TMP))
    res_mit1003_minGLA_MM_Extra[[i]] = TMP
  }
}

save(res_mit1003_minGLA_MM_Extra, file = file.path(root, 'results_new','res_mit1003_minGLA_MM_Extra.Rda'))


