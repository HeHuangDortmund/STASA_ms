#Visual attention identification using random walks based eye tracking protocols

#https://ieeexplore.ieee.org/document/7416925/authors#authors
# X. Chen and Z. Chen, "Visual attention identification using random walks based eye tracking protocols," 
# 2015 IEEE Global Conference on Signal and Information Processing (GlobalSIP), 2015, pp. 6-9, doi: 10.1109/GlobalSIP.2015.7416925.

# AOI center detection 

#input: data.frame of 2D coordinates of fixtions of a single cluster (AOI)
#outout: center of the input coordinates (vector of 2 elements)



AOIcenterDetection = function(singleCluster, sigma = 1, r_def = 100, alpha = 0.2, thresholdT = 0.0005){
  
  
  if(dim(singleCluster)[1] == 1){
    return(c(singleCluster[[1]],singleCluster[[2]]))
  }
  
  #fixDist <- matrix(dist(singleCluster, diag = TRUE, upper = TRUE), nc = dim(singleCluster)[1])
  fixDist = dist(singleCluster, diag = TRUE, upper = TRUE)
  fixDist = as.matrix(fixDist)
  
   n = dim(fixDist)[1]
  #The transition probability from fixation i to j :
  q_ij = matrix(0, nr = n, nc = n)
  for (i in 1:n) {
    nenner = sum(exp(-sigma*fixDist[i,]))
    for (j in 1:n) {
      q_ij[i,j] = exp(-sigma*fixDist[i,j])/nenner
    }
  }
  
  
  Density_r = function(r = r_def){
    rowSums(fixDist <= r)
  }
  w_i = Density_r()/sum(Density_r())
  
  #TODO
  #initialization of l_i  
  l_i = rep(1/n, n)
  
  
  repeat{
    eta = 0
    for (i in 1:n) {
      sum_1_in_7 = 0
      for (j in 1:n) {
        sum_1_in_7 = sum_1_in_7 + (1 - (1 - alpha)*l_i[i])*l_i[j]*q_ij[j,i]
      }
      eta = eta + sum_1_in_7 + (1 - alpha)*l_i[i]*w_i[i]
    }
    
    l_iPLUS1 = numeric(n)
    for (i in 1:n) {
      sum_1_in_6 = 0
      for (j in 1:n) {
        sum_1_in_6 = sum_1_in_6 + (1 - (1 - alpha)*l_i[i])*l_i[j]*q_ij[j,i]
      }
      l_iPLUS1[i] = (sum_1_in_6 + (1- alpha)*l_i[i]*w_i[i])/eta
    }
    
    
    if(all(abs(l_iPLUS1 - l_i) < thresholdT)){
      break
    } else{
      l_i = l_iPLUS1
    }
    
  }
  
  AOIcenter = c(sum(l_i*singleCluster[,1]), sum(l_i*singleCluster[,2]))

  return(AOIcenter)
  
}



runExample = FALSE

if(runExample){
  
  
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
  
  source(file.path(root, 'code','otherMethods','CDBA','outlierScanpathRemoval.R'))
  scanpathsWithoutOutlier = outlierScanpathRemoval(scanpaths)
  
  
  scanpathsWithoutOutlier_fixations = do.call(rbind, scanpathsWithoutOutlier)
  
  
  
  fixDist = dist(scanpathsWithoutOutlier_fixations)
  fixClust = densityClust(fixDist, gaussian=TRUE)
  fixClust = findClusters(fixClust, rho=2, delta=2)
  table(fixClust$clusters)
  scanpathsWithoutOutlier_fixations = cbind(scanpathsWithoutOutlier_fixations,fixClust$clusters)
  scanpathsWithoutOutlier_fixations = as.data.table(scanpathsWithoutOutlier_fixations)
  names(scanpathsWithoutOutlier_fixations)[3] = 'cluster'
  
  
  singleCluster = scanpathsWithoutOutlier_fixations[cluster == 3][,-3]
  
  AOIcenterDetection(singleCluster)
}





