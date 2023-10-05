

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

CDBAScanpathAggregation = function(scanpaths){
  source(file.path(root, 'code','otherMethods','CDBA','outlierScanpathRemoval.R'))
  scanpathsWithoutOutlier = outlierScanpathRemoval(scanpaths)
  
  survivedlen = sapply(scanpathsWithoutOutlier, function(x) dim(x)[1])
  participants_withoutO = rep(names(scanpathsWithoutOutlier),survivedlen)
  

  
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
  
  #------------------------------------------------------------------
  scanpathsWithoutOutlier_fixations$Paticipant = participants_withoutO
  #----------------------------------------------------------------------
  scanpaths_AOIs = merge(scanpathsWithoutOutlier_fixations, AOIcenters_mat, by = 'cluster', sort = FALSE)
  #scanpaths_AOIs[order(Paticipant,cluster)]
  
  # to find out which combination of 2 clusters are not in data
  FindCandidateSet = function(scanpaths_AOIs){
    
    allClusters = sort(unique(scanpaths_AOIs$cluster))
    
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
          if(sp$cluster[i] != sp$cluster[i+1]){ #!!!!!!!!!!!!!! without repeat
            combs = c(combs,paste0(sp$cluster[i], sp$cluster[i+1]))
          }
          
        }
        
      }
      
      return(combs)
    }
    
    all_exist_comb = unique(do.call(c,lapply(scanpathAOIs_list, pasteComb2)))
    
    noX = setdiff(allComb2,all_exist_comb)
    
    condidata_set = list()
    for (i in 1:length(allClusters)) {
      cb = paste0(allClusters[i], 1:length(allClusters))
      condidata_set[[i]] = which(!cb %in% noX)
    }
    
    
    return(list(condidataSet = condidata_set, noXcomb = noX))
  }
  
  
  
  CandidateSetAndNoX = FindCandidateSet(scanpaths_AOIs)
  
 
  
  
  #getCandidateSet = function(scanpaths_AOIs){
  #notExisComb = FindNotExistAOIpair(scanpaths_AOIs)
    
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
    
    

  toMatrix = function(x){
      do.call(rbind, lapply(x, function(i) AOIcenters_mat[AOIcenters_mat$cluster == i,]))[,2:3]
    }
    
    # compute the nearst AOI for each AOI in the condidat set
  NNB = c()
  for (i in 1:n_clusters) {
      
    if(sum(CandidateSetAndNoX$condidataSet[[i]] == i) > 0){
        
        ind = which(CandidateSetAndNoX$condidataSet[[i]] == i)
        candi = CandidateSetAndNoX$condidataSet[[i]][-ind]
    }else{
        candi = CandidateSetAndNoX$condidataSet[[i]]
    }
      
      
    NNB[i] = candi[which.min(dist(toMatrix(i),toMatrix(candi)))]
     
      
  }
    
    
    
  optimal_sp = list()
    for (i in 2:n_clusters) {
      
      #-----------------------------------------------------------------------------------
      # find an initial sp
      aRandomSP = sample(n_clusters, i, replace = TRUE)
      while (!all(!pasteComb(aRandomSP) %in% CandidateSetAndNoX$noXcomb)) {
        aRandomSP = sample(n_clusters, i, replace = TRUE)
      }
      Dis  = sum(sapply(scanpathsWithoutOutlier, function(x) dtw(toMatrix(aRandomSP), x)$distance))
      repeat{
        
        aRandomSP2 = NNB[aRandomSP] # update to the neasrt neighor
        newDis = sum(sapply(scanpathsWithoutOutlier, function(x) dtw(toMatrix(aRandomSP), x)$distance))
        
        verbesserung = newDis - Dis
        #print(verbesserung)
        if(verbesserung > 0) break
        Dis = newDis
        aRandomSP = aRandomSP2
      
      }
      
      optimal_sp[[i-1]] = list(aRandomSP, Dis)
    }
    
    
  alldist = c()
  for (i in 1:(n_clusters -1)) {
     alldist[i] = optimal_sp[[i]][[2]]
  }
    
    
  pathAOI = optimal_sp[[which.min(alldist)]][[1]]
    
    
  
  representSP = toMatrix(pathAOI)
  
  
  # gaze duration analysis
  
  

  scanpaths_noOutlier = scanpaths[Participant %in%  names(scanpathsWithoutOutlier)]
  #----------------------------------------------------------------
  scanpaths_noOutlier$clusters = scanpaths_AOIs$cluster
  #-------------------------------------------------------------
  
  avgDurationOfEachAOI = scanpaths_noOutlier[,mean(Duration), by = 'clusters']
  names(avgDurationOfEachAOI)[2] = 'Duration'
  
  avgDurationOfEachAOI = avgDurationOfEachAOI[order(clusters)]
  
  avgDurationOfEachAOI[pathAOI,]$Duration
  
  representSP$Duration = avgDurationOfEachAOI[pathAOI,]$Duration
  
  names(representSP) = c("Fixation_Position_X_px", "Fixation_Position_Y_px", "Duration")
  
  #TODO change the column names to: "Fixation_Position_X_px", "Fixation_Position_Y_px", "Duration"
  #      representSP as list of two elements: representScanpath = representSP, AOIpolygons = NA
  representScanpath = list(representScanpath = representSP, AOIpolygons = NULL)
  return(representScanpath)
  
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
  
  unique(Fixation_MIT1003$Stimulus)
  # take a single stimulus
  singleStimulus = 'i05june05_static_street_boston_p1010764'
  singleStimulus = unique(Fixation_MIT1003$Stimulus)[200]
  singleStimulus = 'i1182314083'
  
  Fixation_MIT1003_singleStimulus = Fixation_MIT1003[Stimulus == singleStimulus]
  
  
  # INPUT -----------------------------scanpaths
  scanpaths = Fixation_MIT1003_singleStimulus 
  
  representSP = CDBAScanpathAggregation(scanpaths)
  
  
}

