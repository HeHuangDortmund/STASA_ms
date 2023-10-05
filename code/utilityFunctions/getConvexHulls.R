#---------------------------------------------------------------
# get a list of convex hull for each time section  using dbscan
#---------------------------------------------------------------

library(dbscan)
library(concaveman)
library(data.table)

# eps and minPts: The larger these two values are, 
#                 the smaller the number of clusters will be. 
#                 Choose the smallest eps and minPts so that the number of clusters is 1.
#                 -> Stat form a small eps and minPts, 
#                   increase them until number of clusters is 1
#

getConvexHulls = function(Fixation_singleStimulus,
                          ExperimentDuration,
                          et,
                          eps_ini = 400  # half screen width
                          ){ 
  
  timepoint = seq(et, ExperimentDuration, by = et)
  n_timepoint = length(timepoint)
  convex_hull = list()
  for (i in 1:n_timepoint) {
    
    eps = eps_ini
    data_on_section = Fixation_singleStimulus[Event_Start_Trial_Time_ms <= timepoint[i] &
                                                Event_Ende_Trial_Time_ms >= timepoint[i]]
    data_on_section = data_on_section[,.(Fixation_Position_X_px, Fixation_Position_Y_px)]
    
    if(dim(data_on_section)[1] < 3){ # if no or less than 3 fixation on this section, it can not form a polygon
      convex_hull[[i]] = NA
    } else{
      minPts = max(ceiling(dim(data_on_section)[1]/2), 3) # more than half points and at least 3
      db = dbscan(data_on_section, eps = eps, minPts = minPts) #
  
      n_cluster = length(unique(db$cluster[db$cluster != 0])) # why minus one ? cluster 0 are noises
      n_cluster
      
      if(n_cluster == 0){
        convex_hull[[i]] = NA
      }else{
        while (n_cluster == 1) {# decrese eps until no cluster exists
          db = dbscan(data_on_section, eps = eps, minPts = minPts)
          #hullplot(data_on_section, db, main = "DBSCAN")
          n_cluster = length(unique(db$cluster[db$cluster != 0]))
          eps = eps - 1
        }
        eps = eps + 2
        db = dbscan(data_on_section, eps = eps, minPts = minPts)
        n_cluster = length(unique(db$cluster[db$cluster != 0])) #
        
        while (n_cluster == 1) {# increse minPts until no cluster exists
          db = dbscan(data_on_section, eps = eps, minPts = minPts)
          #hullplot(data_on_section, db, main = "DBSCAN")
          n_cluster = length(unique(db$cluster[db$cluster != 0]))
          minPts = minPts + 1
        }
        minPts = minPts - 2
        
        db = dbscan(data_on_section, eps = eps, minPts = minPts)
        n_cluster = length(unique(db$cluster[db$cluster != 0]))
        
        if(n_cluster == 0){
          convex_hull[[i]] = NA
        } else{
          data_on_section[, cluster:= db$cluster]
          clusters = list()
          for (j in 1:n_cluster) {
            clusters[[j]] =  as.matrix(data_on_section[cluster == j][,1:2])
          }
          convex_hull[[i]] =  lapply(clusters, concaveman)
        }
        
      }

      
     
      # if(n_cluster == 0){
      #   convex_hull[[i]] = NA
      # } else{
      #   data_on_section[, cluster:= db$cluster]
      #   clusters = list()
      #   for (j in 1:n_cluster) {
      #     clusters[[j]] =  as.matrix(data_on_section[cluster == j][,1:2])
      #   }
      #   convex_hull[[i]] =  lapply(clusters, concaveman)
      # }
      
   
      
    }
  }
  
  return(convex_hull)
  
}


# Example:
runExample = FALSE
if(runExample == TRUE){
  
  library(rprojroot)
  library(data.table)
  root = find_root(is_git_root)
  
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
  #OR_p = 0.22
  
  convex_hull = getConvexHulls(Fixation_MIT1003_singleStimulus,
                               ExperimentDuration,
                               et)
  
}


