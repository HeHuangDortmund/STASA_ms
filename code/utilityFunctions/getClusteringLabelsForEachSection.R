#--------------------------------------------------------------------------
# get Clustering labels for each time section
#        -i)   convex hulls to gpcPolygon
#        -ii)  compute the similarities of two adjecent polygons
#        -iii) similar polygons together
#-------------------------------------------------------------------------

# Function -getClusteringLabelsForEachSection
# Argument: 
#       - convex_hull: output of getConvexHulls()
#       - OR_p : threshold of the similarity
# Value:
#       vector of the cluster labels of each time cross section, the same length
#       with convex_hull

getClusteringLabelsForEachSection = function(convex_hull,OR_p, compareFirstAndRest = FALSE){
  library(geometry)
  library(rgeos)
  gpcPoly = list()
  for (i in 1:length(convex_hull)) {
    if(is.na(convex_hull[[i]])){
      gpcPoly[[i]] = NA
    } else{
      n_ch = length(convex_hull[[i]])
      if(n_ch == 1){
        tmp = as.data.frame(convex_hull[[i]])
        names(tmp) = c("x", "y")
        gpcPoly[[i]] = as(tmp, "gpc.poly")
      } else{
        chtmp = list()
        for (k in 1:n_ch) {
          chtmp[[k]] = as.data.frame(convex_hull[[i]][[k]])
          names(chtmp[[k]]) = c("x", "y")
          chtmp[[k]] = as(chtmp[[k]], "gpc.poly")
        }
        
        CHtmp = chtmp[[1]]
        for (h in 2:n_ch) {
          CHtmp = union(CHtmp, chtmp[[h]])
        }
        gpcPoly[[i]] = CHtmp
        
      }
    }
  }
  
  
  
  
  
  n_plg = length(gpcPoly)
  
  ## 2.1 similarity of each 2 neighbors
  similarities = c()
  NAorNOt = is.na(gpcPoly)
  for (i in 1:(n_plg-1)) {
    #if(is.na(gpcPoly[[i]]) | is.na(gpcPoly[[i+1]])){
    if(NAorNOt[i] | NAorNOt[i+1]){
      similarities[i] = 0
    }else{
      area_union = area.poly(union(gpcPoly[[i]], gpcPoly[[i+1]]))
      area_intersec = area.poly(intersect(gpcPoly[[i]], gpcPoly[[i+1]]))
      similarity = area_intersec/area_union
      similarities[i] = similarity
    }
  }
  
  
  
  #------------------------rewrite--------------------------------------------
  clusterings = c()
  firstCluster = min(which(similarities > 0))
  clusterings[firstCluster] = 1
  #clusterings[firstCluster+1] = 1
  if(firstCluster > 1){
    clusterings[1:(firstCluster-1)] = 0
  }
  
  
  clus = 1
  for (i in (firstCluster):(n_plg-1)) {
    if(similarities[i]>= OR_p){
      clusterings[i+1] = clus
    }else if(similarities[i] == 0 & is.na(gpcPoly[[i+1]])){
      clusterings[i+1] = 0
    }else if(similarities[i] < OR_p & !is.na(gpcPoly[[i+1]])){
      clus = clus + 1
      clusterings[i+1] = clus
    }
  }
  #---------------------------------------------------------------------
  

  #-------------------next step -----------------------------------------------
  # compare the similarity between fist section and the rest sections within the cluster,
  # 
  # instead of neighbors
  if(compareFirstAndRest == TRUE){
    clusering_uniq = unique(clusterings)
    clusering_uniq = clusering_uniq[clusering_uniq!=0]
    
    for (i in seq(along = clusering_uniq)) {
      
      ind = which(clusterings == clusering_uniq[i])
      clusterBegin = min(ind)
      
      similarity2 = c()
      if(length(ind) > 2){
        for (j in 2:length(ind)) {
          area_union = area.poly(union(gpcPoly[[clusterBegin]], gpcPoly[[ind[j]]]))
          area_intersec = area.poly(intersect(gpcPoly[[clusterBegin]], gpcPoly[[ind[j]]]))
          similarity2[j-1] = area_intersec/area_union
          
        }
      }
      
      
      newclusterBegins = c()
      
      while (length(which(similarity2 < OR_p)) != 0) {
        xxx= min(which(similarity2 < OR_p)) + 1
        newclusterBegin = ind[xxx]
        newclusterBegins = c(newclusterBegins, newclusterBegin)
        ind = ind[-(1:(xxx-1))]
        
        similarity2 = c()
        if(length(ind) > 2){
          for (j in 2:length(ind)) {
            area_union = area.poly(union(gpcPoly[[newclusterBegin]], gpcPoly[[ind[j]]]))
            area_intersec = area.poly(intersect(gpcPoly[[newclusterBegin]], gpcPoly[[ind[j]]]))
            similarity2[j-1] = area_intersec/area_union
            
          }
        }
        
      }
      
      #newBreaks = c(newBreaks,newclusterBegins)
      if(!is.null(newclusterBegins)){
        
        ind2 = which(clusterings == clusering_uniq[i])
        FirstSec = min(ind2)
        LastSec = max(ind2)
        
        breaks = c(FirstSec,newclusterBegins,LastSec)
        
        for (i in 1:(length(breaks)-1)) {
          clusterings[breaks[i]:breaks[i+1]] = breaks[i]*100
        }
        
      }
      
    }
    
    
    
    
    newClusters = unique(clusterings)[unique(clusterings) != 0]
    
    
    clusterings2 = clusterings
    for (i in 1:length(newClusters)) {
      clusterings2[clusterings == newClusters[i]] = i
    }
    
    
    
    
    
    
    #--------------------------------------------------------------------------
    #similarities >= OR_p
    #plot(similarities, type = "l")
    # clusterings = c()
    # firstCluster = min(which(similarities > 0))
    # clusterings[firstCluster] = 1
    # clusterings[1:(firstCluster-1)] = 0
    # clus = 1
    # for (i in firstCluster:(n_plg-1)){
    #   
    #   if(similarities[i]>= OR_p){
    #     clusterings[i+1] = clus
    #   }else if(similarities[i] == 0){
    #     clusterings[i+1] = 0
    #   }else{
    #     clus = clus +1
    #     clusterings[i+1] = clus
    #   }
    #   
    # }
    return(clusterings2)
  }else{
    return(clusterings)
  }
}


# Example

runExampleOfgetClusteringLabel = FALSE

if(runExampleOfgetClusteringLabel){
  
  
  library(rprojroot)
  library(data.table)
  root = find_root(is_git_root)
  source(file.path(root,'code','utilityFunctions','getConvexHulls.R'))
  
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
  
  convex_hull = getConvexHulls(Fixation_MIT1003_singleStimulus,
                               ExperimentDuration,
                               et)
  
  # 2nd: getClusteringLabelsForEachSection
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  
}



getClusteringLabelsForEachSection2 = function(convex_hull,OR_p, compareFirstAndRest = FALSE){
  library(geometry)
  library(rgeos)
  gpcPoly = list()
  for (i in 1:length(convex_hull)) {
    if(is.na(convex_hull[[i]])){
      gpcPoly[[i]] = NA
    } else{
      n_ch = length(convex_hull[[i]])
      if(n_ch == 1){
        tmp = as.data.frame(convex_hull[[i]])
        names(tmp) = c("x", "y")
        gpcPoly[[i]] = as(tmp, "gpc.poly")
      } else{
        chtmp = list()
        for (k in 1:n_ch) {
          chtmp[[k]] = as.data.frame(convex_hull[[i]][[k]])
          names(chtmp[[k]]) = c("x", "y")
          chtmp[[k]] = as(chtmp[[k]], "gpc.poly")
        }
        
        CHtmp = chtmp[[1]]
        for (h in 2:n_ch) {
          CHtmp = union(CHtmp, chtmp[[h]])
        }
        gpcPoly[[i]] = CHtmp
        
      }
    }
  }
  
  
  
  
  
  n_plg = length(gpcPoly)
  
  ## 2.1 similarity of each 2 neighbors
  similarities = c()
  NAorNOt = is.na(gpcPoly)
  for (i in 1:(n_plg-1)) {
    #if(is.na(gpcPoly[[i]]) | is.na(gpcPoly[[i+1]])){
    if(NAorNOt[i] | NAorNOt[i+1]){
      similarities[i] = 0
    }else{
      area_union = area.poly(union(gpcPoly[[i]], gpcPoly[[i+1]]))
      area_intersec = area.poly(intersect(gpcPoly[[i]], gpcPoly[[i+1]]))
      similarity = area_intersec/area_union
      similarities[i] = similarity
    }
  }
  
  
  #------------------------rewrite--------------------------------------------
  clusterings = c()
  firstCluster = min(which(similarities > 0))
  clusterings[firstCluster] = 1
  #clusterings[firstCluster+1] = 1
  if(firstCluster > 1){
    clusterings[1:(firstCluster-1)] = 0
  }
  
  
  clus = 1
  for (i in (firstCluster):(n_plg-1)) {
    if(similarities[i]>= OR_p){
      clusterings[i+1] = clus
    }else if(similarities[i] == 0 & is.na(gpcPoly[[i+1]])){
      clusterings[i+1] = 0
    }else if(similarities[i] < OR_p & !is.na(gpcPoly[[i+1]])){
      clus = clus + 1
      clusterings[i+1] = clus
    }
  }
  #---------------------------------------------------------------------
  
  
  #-------------------next step -----------------------------------------------
  # compare the similarity between fist section and the rest sections within the cluster,
  # 
  # instead of neighbors
  if(compareFirstAndRest == TRUE){
    clusering_uniq = unique(clusterings)
    clusering_uniq = clusering_uniq[clusering_uniq!=0]
    
    for (i in seq(along = clusering_uniq)) {
      
      ind = which(clusterings == clusering_uniq[i])
      clusterBegin = min(ind)
      
      similarity2 = c()
      if(length(ind) > 2){
        for (j in 2:length(ind)) {
          area_union = area.poly(union(gpcPoly[[clusterBegin]], gpcPoly[[ind[j]]]))
          area_intersec = area.poly(intersect(gpcPoly[[clusterBegin]], gpcPoly[[ind[j]]]))
          similarity2[j-1] = area_intersec/area_union
          
        }
      }
      
      
      newclusterBegins = c()
      
      while (length(which(similarity2 < OR_p)) != 0) {
        xxx= min(which(similarity2 < OR_p)) + 1
        newclusterBegin = ind[xxx]
        newclusterBegins = c(newclusterBegins, newclusterBegin)
        ind = ind[-(1:(xxx-1))]
        
        similarity2 = c()
        if(length(ind) > 2){
          for (j in 2:length(ind)) {
            area_union = area.poly(union(gpcPoly[[newclusterBegin]], gpcPoly[[ind[j]]]))
            area_intersec = area.poly(intersect(gpcPoly[[newclusterBegin]], gpcPoly[[ind[j]]]))
            similarity2[j-1] = area_intersec/area_union
            
          }
        }
        
      }
      
      #newBreaks = c(newBreaks,newclusterBegins)
      if(!is.null(newclusterBegins)){
        
        ind2 = which(clusterings == clusering_uniq[i])
        FirstSec = min(ind2)
        LastSec = max(ind2)
        
        breaks = c(FirstSec,newclusterBegins,LastSec)
        
        for (i in 1:(length(breaks)-1)) {
          clusterings[breaks[i]:breaks[i+1]] = breaks[i]*100
        }
        
      }
      
    }
    
    
    
    
    newClusters = unique(clusterings)[unique(clusterings) != 0]
    
    
    clusterings2 = clusterings
    for (i in 1:length(newClusters)) {
      clusterings2[clusterings == newClusters[i]] = i
    }
    
    
    
    
    
    
    #--------------------------------------------------------------------------
    #similarities >= OR_p
    #plot(similarities, type = "l")
    # clusterings = c()
    # firstCluster = min(which(similarities > 0))
    # clusterings[firstCluster] = 1
    # clusterings[1:(firstCluster-1)] = 0
    # clus = 1
    # for (i in firstCluster:(n_plg-1)){
    #   
    #   if(similarities[i]>= OR_p){
    #     clusterings[i+1] = clus
    #   }else if(similarities[i] == 0){
    #     clusterings[i+1] = 0
    #   }else{
    #     clus = clus +1
    #     clusterings[i+1] = clus
    #   }
    #   
    # }
    return(clusterings2)
  }else{
    return(list(clusterings = clusterings, similarities = similarities))
  }
}

