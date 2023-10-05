#DTWdistance 

# Input: Scanpath Matrix, scanpath of one of the Stimuli in Fixation_MIT1003.Rda 
# Output: ScanpathList without outlier
library(dtw)
library(rprojroot)
library(data.table)
root = find_root(is_git_root)

outlierScanpathRemoval = function(scanpaths){
  #TODO
  scanpaths = scanpaths[,c('Participant',"Fixation_Position_X_px", "Fixation_Position_Y_px")]
  Participants = unique(scanpaths$Participant)
  scanpathList = list()
  for (i in seq(along = Participants)) {
    scanpathList[[i]] = as.matrix(scanpaths[Participant == Participants[[i]],c("Fixation_Position_X_px", "Fixation_Position_Y_px")])
  }
  
  total_distance_list = c()
  for (i in seq(along = Participants)) {
    
  total_distance_list[i] = 0
    for (j in seq(along = Participants)) {
      if(j != i){
        total_distance_list[i] = total_distance_list[i] + dtw(scanpathList[[i]], scanpathList[[j]], keep = TRUE)$distance
      }
    }
  }
  
  Q1 = quantile(total_distance_list, 0.25)
  Q3 = quantile(total_distance_list, 0.75)
  ind = which(total_distance_list > Q3 + 1.5*(Q3 - Q1) | total_distance_list < Q1 - 1.5*(Q3 - Q1))
  if(length(ind) > 0){
    scanpathList = scanpathList[-ind]
    Participants = Participants[-ind]
  }

  names(scanpathList) = Participants
  #TODO return LIST
  return(scanpathList)
  
}

runExample = FALSE
if(runExample){
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
  
  scanpaths = Fixation_MIT1003_singleStimulus 
  
  scanpathsWithoutOutlier = outlierScanpathRemoval(scanpaths)
}




