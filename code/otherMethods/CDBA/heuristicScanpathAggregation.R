

library(densityClust)
library(rprojroot)
library(data.table)
library(gtools)
library(dtw)
root = find_root(is_git_root)




#Input:
# scanpaths: data.table contains following columns
#             -Participant
#             -Fixation_Position_X_px
#             -Fixation_Position_Y_px
#             -Duration

# Output: 
# representSP: data.table contains following columns
#             AOIcenter_x
#             AOIcenter_y
#             Duration

heuristicScanpathAggregation = function(scanpaths){
  source(file.path(root, 'code','otherMethods','CDBA','outlierScanpathRemoval.R'))
  scanpathsWithoutOutlier = outlierScanpathRemoval(scanpaths)
  scanpathsWithoutOutlier_fixations = do.call(rbind, scanpathsWithoutOutlier)
  
  
  
  fixDist <- dist(scanpathsWithoutOutlier_fixations)
  fixClust <- densityClust(fixDist, gaussian=TRUE)
  
  xxx=sort(fixClust$rho*fixClust$delta,decreasing = TRUE)
  ind = which(fixClust$rho*fixClust$delta > mean(xxx))
  fixClust <- findClusters(fixClust, peaks = ind)
  
  #fixClust <- findClusters(fixClust, rho=2, delta=20) #TODO 2?
  
  #table(fixClust$clusters)
  
  scanpathsWithoutOutlier_fixations = cbind(scanpathsWithoutOutlier_fixations,fixClust$clusters)
  scanpathsWithoutOutlier_fixations = as.data.table(scanpathsWithoutOutlier_fixations)
  names(scanpathsWithoutOutlier_fixations)[3] = 'cluster'
  
  
  scanpathsWithoutOutlier_fixations_list = split(scanpathsWithoutOutlier_fixations, by = 'cluster')
  
  # detect the AOI center for each cluster
  source(file.path(root, 'code', 'utilityFunctions','AOIcenterDetection.R'))
  AOIcenters = lapply(scanpathsWithoutOutlier_fixations_list, AOIcenterDetection)
  AOIcenters_mat = as.data.table(cbind(names(AOIcenters),do.call(rbind,AOIcenters)))
  names(AOIcenters_mat) = c('cluster', 'AOIcenter_x','AOIcenter_y')
  AOIcenters_mat$cluster = as.numeric(AOIcenters_mat$cluster)
  AOIcenters_mat$AOIcenter_x = as.numeric(AOIcenters_mat$AOIcenter_x)
  AOIcenters_mat$AOIcenter_y = as.numeric(AOIcenters_mat$AOIcenter_y)
  
  scanpathsWithoutOutlier_fixations$Paticipant = scanpaths$Participant
  
  scanpaths_AOIs = merge(scanpathsWithoutOutlier_fixations, AOIcenters_mat, by = 'cluster', sort = FALSE)
  
  
  # to find out which combination of 2 clusters are not in data
  FindNotExistAOIpair = function(scanpaths_AOIs){
    
    allClusters = unique(scanpaths_AOIs$cluster)
    
    # all possible combination of 2 clusters
    allComb2 = apply( expand.grid(rep(list(allClusters), 2)),1, paste0,collapse = '')
    
    
    # all combination of 2 clusters in data
    scanpathAOIs_list = split(scanpaths_AOIs, by = 'Paticipant')
    
    #sp = scanpathAOIs_list[[1]]
    
    pasteComb2 = function(sp){
      combs = c()
      
      n_clu = length(sp$cluster)
      
      if(n_clu > 1){
        for (i in 1:(n_clu - 1)) {
          combs = c(combs,paste0(sp$cluster[i], sp$cluster[i+1]))
        }
        
      }
      
      return(combs)
    }
    
    all_exist_comb = unique(do.call(c,lapply(scanpathAOIs_list, pasteComb2)))
    
    noX = setdiff(allComb2,all_exist_comb)
    
    return(noX)
  }
  
  
  getCandidateSet = function(scanpaths_AOIs){
    notExisComb = FindNotExistAOIpair(scanpaths_AOIs)
    
    n_clusters = length(unique(scanpaths_AOIs$cluster))
    
    pasteComb = function(x){
      combs = c()
      
      n_clu = length(x)
      
      if(n_clu > 1){
        for (i in 1:(n_clu - 1)) {
          combs = c(combs,paste0(x[i], x[i+1]))
        }
        
      }
      
      return(combs)
    }
    
    CandidateSet = list()
    for (i in 2:n_clusters) {
      tmp = permutations(n = n_clusters, r = i, v = 1:n_clusters,repeats.allowed = TRUE)
      
      pasted = apply(tmp,1, pasteComb)
      
      if(i == 2){
        ind = which(pasted %in% notExisComb)
      }else{
        
        minifun = function(x) x %in% notExisComb
        ind = which(apply(apply(pasted, 2, minifun), 2, any))
      }
      
      CandidateSet[[i-1]] = tmp[-ind,]
    }
    
    desired_length = sum(sapply(CandidateSet, nrow))
    
    CandidateSet_list = list()
    
    for (i in 1:length(CandidateSet)) {
      n_row = dim(CandidateSet[[i]])[1]
      for (j in 1:n_row) {
        CandidateSet_list = c(CandidateSet_list,list(CandidateSet[[i]][j,]))
      }
    }
    return(CandidateSet_list)
    
    
    
  }
  
  
  candidateSet_list = getCandidateSet(scanpaths_AOIs)
  #TODO 
  
  toMatrix = function(x){
    do.call(rbind, lapply(x, function(i) AOIcenters_mat[AOIcenters_mat$cluster == i,]))[,2:3]
  }
  
  #dtw(toMatrix(candidateSet_list[[100]]),toMatrix(candidateSet_list[[120]]))$distance
  
  
  #candidateSP = candidateSet_list[[11]]
  sumDtwDistance = function(candidateSP){
    sum(sapply(scanpathsWithoutOutlier, function(x) dtw(toMatrix(candidateSP), x)$distance))
  }
  
  # res = c()
  # for (i in 1:length(candidateSet_list)) {
  #   res[i] = sumDtwDistance(candidateSet_list[[i]])
  # }
  
  
  
  
  res = lapply(candidateSet_list, sumDtwDistance)
  
  pathAOI = candidateSet_list[[which.min(res)]]
  
  
  
  representSP = toMatrix(pathAOI)
  
  
  # gaze duration analysis
  scanpaths$clusters = scanpaths_AOIs$cluster
  avgDurationOfEachAOI = scanpaths[,mean(Duration), by = 'clusters']
  names(avgDurationOfEachAOI)[2] = 'Duration'
  
  avgDurationOfEachAOI = avgDurationOfEachAOI[order(clusters)]
  
  avgDurationOfEachAOI[pathAOI,]$Duration
  
  representSP$Duration = avgDurationOfEachAOI[pathAOI,]$Duration
  
  names(representSP) = c("Fixation_Position_X_px", "Fixation_Position_Y_px", "Duration")
  
  #TODO change the column names to: "Fixation_Position_X_px", "Fixation_Position_Y_px", "Duration"
  #      representSP as list of two elements: representScanpath = representSP, AOIpolygons = NA
  representScanpath = list(representScanpath = representSP, AOIpolygons = NULL)
  return(representSP)
  
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
  
  
  # INPUT -----------------------------scanpaths
  scanpaths = Fixation_MIT1003_singleStimulus 
  
  representSP = heuristicScanpathAggregation(scanpaths)
  
  
}

