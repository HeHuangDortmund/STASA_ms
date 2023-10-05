#getKConcaveHulls: function to get K ConcaveHulls by giving data with colums x and y, and giving
#                 nummer K
#           Output: a list with K elements, some elements may be NA



getKConcaveHullsForSingleCS = function(K,data_xy, eps_ini = 400){ # K >= 2
  
  if(dim(data_xy)[1] < 3){
    return(NA)
  }else{
    kmeanClusters = kmeans(data_xy, centers = K)
    
    data_clusters = list()
    for (k in 1:K) {
      data_clusters[[k]] = data_xy[kmeanClusters$cluster == k,]
    }
    source(file.path(root,'code','utilityFunctions','adaptiveDBSCAN.R'))
    lapply(data_clusters, adaptiveDBSCAN, eps_ini = eps_ini)
  }
  
  
}

getKConaveHullsForSPs = function(Fixation_singleStimulus,
                          ExperimentDuration,
                          et,
                          eps_ini = 400,
                          K){ 
  
  timepoint = seq(et, ExperimentDuration, by = et)
  n_timepoint = length(timepoint)
  concave_hull = list()
  for (i in 1:n_timepoint) {
    data_on_section = Fixation_singleStimulus[Event_Start_Trial_Time_ms <= timepoint[i] &
                                                Event_Ende_Trial_Time_ms >= timepoint[i]]
    data_on_section = data_on_section[,.(Fixation_Position_X_px, Fixation_Position_Y_px)]
    concave_hull[[i]] = getKConcaveHullsForSingleCS(data_on_section,eps_ini = eps_ini, K = K)
  }
  return(concave_hull)
}




runExampleOfGetReprastvScanpath = FALSE
if(runExampleOfGetReprastvScanpath){
  library(rprojroot)
  library(data.table)
  root = find_root(is_git_root)
  source(file.path(root,'code','utilityFunctions','getConvexHulls.R'))
  source(file.path(root,'code','utilityFunctions','getClusteringLabelsForEachSection.R'))
  
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
  singleStimulus = 'i05june05_static_street_boston_p1010764'
  
  Fixation_MIT1003_singleStimulus = Fixation_MIT1003[Stimulus == singleStimulus]
  
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  OR_p = 0.22
  

  
  res = getKConaveHullsForSPs(Fixation_singleStimulus=Fixation_MIT1003_singleStimulus,
                        ExperimentDuration,
                        et = 30/1000,
                        eps_ini = 400,
                        K=2)
  
  
  #clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  
  #RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,Fixation_MIT1003_singleStimulus)
  
  
  #prune AOI with only one or two points
  
}


