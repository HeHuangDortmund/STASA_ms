# getRepresentScanpath: Funktion to compute the respresentative scanpath
# Arguments: 
#          - convex_hull: list of convex_hulls
#          - clusterings: vector of the labels of the clustering with the same length with convex_hull
#          - et : time interval of two cross sections
# Value:
#          - list with two elements:
#            -$representScanpath
#            -$AOIpolygons: list of class'SpatialPolygons'

library(prevR)

getRepresentScanpath = function(convex_hull, clusterings, et,
                                Fixation_singleStimulus, 
                                unionAOI = TRUE,
                                durationrec =FALSE){
  

    coordsToSpatialPolygons = function(x){
    if(is.na(x)){
      sps = NA
    }else{
      p = Polygon(x[[1]])
      ps = Polygons(list(p),1)
      sps = SpatialPolygons(list(ps))
    }
    return(sps)
  }
  SpatialPolygonsofHulls = lapply(convex_hull, coordsToSpatialPolygons)
  #table(clusterings)
  #labels = names(table(clusterings))[-1] # remove cluster 0
  labels = names(table(clusterings))[names(table(clusterings))!=0]

  if(unionAOI == FALSE){
    intersectPolygon_list = list()
    for (i in seq(along = labels)) {
      ind = which(clusterings == labels[i])
      #get intersections of the convex hulls with index 'ind'
      if(length(ind) == 1){
        intersectPolygon = SpatialPolygonsofHulls[[ind]]
      } else{
        intersectPolygon = SpatialPolygonsofHulls[[ind[1]]]
        for (j in ind[-1]) {
          intersectPolygon = gIntersection(intersectPolygon, 
                                           SpatialPolygonsofHulls[[j]])
          
          if(class(intersectPolygon) == "SpatialCollections"){
            # some times the results is not a SpatialPolygons
            #  e.g. collection of polygons and lines
            intersectPolygon = intersectPolygon@polyobj
          }
          
        }
        
        
      }
      intersectPolygon_list = c(intersectPolygon_list, intersectPolygon)
    }
  }else{
    overlapPolygon_list = list()
    for (i in seq(along = labels)) {
      ind = which(clusterings == labels[i])
      #get intersections of the convex hulls with index 'ind'
      if(length(ind) == 1){
        overlapPolygon = SpatialPolygonsofHulls[[ind]]
      } else{
        overlapPolygon = SpatialPolygonsofHulls[[ind[1]]]
        for (j in ind[-1]) {
          overlapPolygon = gUnion(overlapPolygon, 
                                           SpatialPolygonsofHulls[[j]])
          
          if(class(overlapPolygon) == "SpatialCollections"){
            # some times the results is not a SpatialPolygons
            #  e.g. collection of polygons and lines
            overlapPolygon = overlapPolygon@polyobj
          }
          
        }
        
        
      }
      overlapPolygon_list = c(overlapPolygon_list, overlapPolygon)
    
    }
  }
  

  #clusterCentroid = intersectPolygon_list[[1]]@polygons[[1]]@Polygons[[1]]@labpt
  #clusterDuration = table(clusterings)[-1]*et
  clusterDuration = table(clusterings)[names(table(clusterings))!=0]*et
  #length(intersectPolygon_list)
  
  tmpFun = function(x){
    
    if(class(x) == "SpatialPoints"){
      
      if(!is.null(dim(x@coords))){
        return(colMeans(x@coords))
      }else{
        return(x@coords)
      }
      
      
    }else if(class(x) == "SpatialLines"){
      return(colMeans(x@lines[[1]]@Lines[[1]]@coords))
    }else{
      return(x@polygons[[1]]@Polygons[[1]]@labpt)
    }
    
  }
  
  if(unionAOI == FALSE){
    Fixation_Position = do.call(rbind,lapply(intersectPolygon_list, tmpFun))
  }else{
    Fixation_Position = do.call(rbind,lapply(overlapPolygon_list, tmpFun))
  }


  indLogi = clusterDuration >= 0.1  #  fixation should not be short as 0.1 second?
 
  
  start_end = do.call(rbind,lapply(seq(along = labels), function(i) range(which(clusterings == labels[i])*et)))
 
  if(length(labels) == 1 | sum(indLogi) == 1){
    representScanpath = c(Fixation_Position[indLogi,],start_end[indLogi,],clusterDuration[indLogi])
    representScanpath = data.frame(Fixation_Position_X_px = representScanpath[1], 
                                   Fixation_Position_Y_px = representScanpath[2],
                                   Event_Start_Trial_Time_ms= representScanpath[3], 
                                   Event_Ende_Trial_Time_ms= representScanpath[4],
                                   Duration = representScanpath[5])
  }else{
    representScanpath = cbind(Fixation_Position[indLogi,],start_end[indLogi,],clusterDuration[indLogi])
    representScanpath = as.data.frame(representScanpath)
    names(representScanpath) = c("Fixation_Position_X_px", "Fixation_Position_Y_px","Event_Start_Trial_Time_ms", "Event_Ende_Trial_Time_ms","Duration")
  }
  
  
  
  #return(list(representScanpath = representScanpath, AOIpolygons = intersectPolygon_list))  
  #------------------------------------------------------------------------------------
  # 
  #----------------------------------------------------------------------------------------------HHHHH
  
  if(durationrec == TRUE){
    ind = list()
    if(unionAOI == FALSE){
      for (i in 1:length(intersectPolygon_list)) {
        ind[[i]] = point.in.SpatialPolygons(Fixation_singleStimulus$Fixation_Position_X_px,
                                            Fixation_singleStimulus$Fixation_Position_Y_px,
                                            intersectPolygon_list[[i]])
      }
    }else{
      for (i in 1:length(overlapPolygon_list)) {
        ind[[i]] = point.in.SpatialPolygons(Fixation_singleStimulus$Fixation_Position_X_px,
                                            Fixation_singleStimulus$Fixation_Position_Y_px,
                                            overlapPolygon_list[[i]])
      }
    }
    
    
    
    fix_duration = c()
    
    for (i in 1:length(ind)) {
      fix_duration[i] = mean(Fixation_singleStimulus$Duration[ind[[i]]], na.rm = TRUE)
    }
    
    
    fix_duration = fix_duration[indLogi]
    
    tosub = !is.na(fix_duration)
    representScanpath$Duration[tosub] = fix_duration[tosub]
  }
  
 

  if(unionAOI == TRUE){
    AOIpolygons = overlapPolygon_list[indLogi]
  }else{
    AOIpolygons = intersectPolygon_list[indLogi]
    }
  
  #AOIpolygons = ifelse(unionAOI, overlapPolygon_list[indLogi], intersectPolygon_list[indLogi])
  return(list(representScanpath = representScanpath, AOIpolygons = AOIpolygons))  
  # 
}


# Example

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
  
  convex_hull = getConvexHulls(Fixation_MIT1003_singleStimulus,
                               ExperimentDuration,
                               et)
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,Fixation_MIT1003_singleStimulus)
  RepresentScanpath$representScanpath
  length(RepresentScanpath$AOIpolygons)
  #prune AOI with only one or two points
  
}


#str(RepresentScanpath)
