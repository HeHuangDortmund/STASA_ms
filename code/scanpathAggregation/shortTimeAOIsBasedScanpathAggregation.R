library(dbscan)
library(concaveman)
library(rgl)
library(rgeos)
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


#preparing dara for plot, duplicate each row, 1st row as start time, snd row as end time of the fixation
data_for_plot = Fixation_MIT1003_singleStimulus[rep(seq_len(nrow(Fixation_MIT1003_singleStimulus)), each = 2), ]
# even rows 
index = seq(2, nrow(data_for_plot), 2)
data_for_plot$Event_Start_Trial_Time_ms[index] = data_for_plot$Event_Ende_Trial_Time_ms[index]
data_for_plot[,Event_Ende_Trial_Time_ms:=NULL]

# clustering using DBSCAN for each time point

ExperimentDuration = 3  # 3 seconds
et = 30/1000 # 30 ms
OR_p = 0.22


#---------------------------------------------------------------
# get a list of convex hull for each time section  using dbscan
#---------------------------------------------------------------

getConvexHulls = function(Fixation_singleStimulus,
                          ExperimentDuration,
                          et,
                          OR_p){
  
  #ExperimentDuration = 3  # 3 seconds
  #et = 30/1000 # 30 ms
  #OR_p = 0.22
  timepoint = seq(et, ExperimentDuration, by = et)
  n_timepoint = length(timepoint)
  
  convex_hull = list()
  
  for (i in 1:n_timepoint) {
    
    data_on_section = Fixation_singleStimulus[Event_Start_Trial_Time_ms <= timepoint[i] &
                                                        Event_Ende_Trial_Time_ms >= timepoint[i]]
    data_on_section = data_on_section[,.(Fixation_Position_X_px, Fixation_Position_Y_px)]
    
    if(dim(data_on_section)[1] == 0){
      convex_hull[[i]] = NA
    } else{
      db = dbscan(data_on_section, eps = 30, minPts = 3)
      #hullplot(data_on_section, db, main = "DBSCAN")
      n_cluster = length(unique(db$cluster)) - 1 # why minus one ? cluster 0 are noises
      
      
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
  }
  
  return(convex_hull)
  
}


convex_hull = getConvexHulls(Fixation_MIT1003_singleStimulus,
                             ExperimentDuration,
                             et,
                             OR_p)

#--------------------------------------------------------------------------
# get Clustering labels for each time section
#        -i)   convex hulls to gpcPolygon
#        -ii)  compute the similarities of two adjecent polygons
#        -iii) similar polygons together
#-------------------------------------------------------------------------

getClusteringLabelsForEachSection = function(convex_hull){
  library(geometry)
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
  for (i in 1:(n_plg-1)) {
    if(is.na(gpcPoly[[i]]) | is.na(gpcPoly[[i+1]])){
      similarities[i] = 0
    }else{
      area_union = area.poly(union(gpcPoly[[i]], gpcPoly[[i+1]]))
      area_intersec = area.poly(intersect(gpcPoly[[i]], gpcPoly[[i+1]]))
      similarity = area_intersec/area_union
      similarities[i] = similarity
    }
  }
  
  
  
  #similarities >= OR_p
  #plot(similarities, type = "l")
  clusterings = c()
  firstCluster = min(which(similarities > 0))
  clusterings[firstCluster] = 1
  clusterings[1:(firstCluster-1)] = 0
  clus = 1
  for (i in firstCluster:(n_plg-1)){
    
    if(similarities[i]>= OR_p){
      clusterings[i+1] = clus
    }else if(similarities[i] == 0){
      clusterings[i+1] = 0
    }else{
      clus = clus +1
      clusterings[i+1] = clus
    }
    
  }
  return(clusterings)
  
}

clusterings = getClusteringLabelsForEachSection(convex_hull)


# getRepresentScanpath: Funktion to compute the respresentative scanpath
#   -note- need et.
# Arguments: 
#          - convex_hull: list of convex_hulls
#          - clusterings: vector of the labels of the clustering with the same length with convex_hull
# Value:
#          - list with two elements:
#            -$representScanpath
#            -$AOIpolygons: list of class'SpatialPolygons'

getRepresentScanpath = function(convex_hull, clusterings){
  
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
  labels = names(table(clusterings))[-1] # remove cluster 0

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
  
  #clusterCentroid = intersectPolygon_list[[1]]@polygons[[1]]@Polygons[[1]]@labpt
  clusterDuration = table(clusterings)[-1]*et
  

  Fixation_Position = do.call(rbind,lapply(intersectPolygon_list, function(x) x@polygons[[1]]@Polygons[[1]]@labpt))
  start_end = do.call(rbind,lapply(seq(along = labels), function(i) range(which(clusterings == labels[i])*et)))
  
  representScanpath = cbind(Fixation_Position,start_end,clusterDuration)
  representScanpath = as.data.frame(representScanpath)
  names(representScanpath) = c(
    "Fixation_Position_X_px", 
    "Fixation_Position_Y_px",
    "Event_Start_Trial_Time_ms", 
    "Event_Ende_Trial_Time_ms",
    "Duration")
  
  return(list(representScanpath = representScanpath, AOIpolygons = intersectPolygon_list))  
}


res = getRepresentScanpath(convex_hull, clusterings)





#clusterings
#length(clusterings)

library(RColorBrewer)
n = length(unique(clusterings))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)
col_clusterings = col[clusterings+1]

#---------------------------------------------------------------------------------------
# visualisation
point_list = list()
n_list = nrow(data_for_plot)/2
for (i in 1:n_list) {
  point_list[[i]] = data_for_plot[(2*i-1):(2*i),]
}



ROIS = list()
for (i in 1:length(convex_hull)) {
  
  if(is.na(convex_hull[[i]])){
    ROIS[[i]] = NA
  }else{
    n_cl = length(convex_hull[[i]])
    ROIS[[i]] = list()
    for (j in 1:n_cl) {
      poly =  convex_hull[[i]][[j]]
      poly3d = cbind(poly, timepoint[i])
      colnames(poly3d) = c("x", "y", "z")
      poly3d = as.data.frame(poly3d)
      ROIS[[i]][[j]] = poly3d
    }
  }
}



#------------------------------------------------------------ 
##################

library(imager)
img = load.image(file.path(root,'data','MIT1003','ALLSTIMULI',paste0(singleStimulus,'.jpeg')))

#draw lines and polygons on image
#plot(img)
#abline(h = 600, col = 'red')
#circle <- (Xc(img) - 200)^2 + (Yc(img) - 350)^2 < 150^2 #250, 350 center
#highlight(circle)]
#polygon(convex_hull[[2]][[1]][,1],convex_hull[[2]][[1]][,2], col = 'red')



img_grey = grayscale(img)[,,1,1]
rownames(img_grey) = 1:1024
colnames(img_grey) = 1:768
px = rep(1:1024, 768)
py = rep(1:1024, each = 768)
pz = rep(10, 1024*768)
col = as.vector(img_grey)
plot3d(x = px, y = py, z = pz, xlim = c(0,1024), ylim = c(0,768), zlim = c(0,3),col = grey(col, 0.5))


for (i in 1:length(point_list)) {
  plot3d(
    x = point_list[[i]]$Fixation_Position_X_px, 
    y = point_list[[i]]$Fixation_Position_Y_px, 
    z = point_list[[i]]$Event_Start_Trial_Time_ms, 
    type="l", add = TRUE, col = "blue")
}


for (i in 1:length(ROIS)) {
  
  if(!is.na(ROIS[[i]])){
    n_cls = length(ROIS[[i]])
    for (j in 1:n_cls) {
      plot3d(x=ROIS[[i]][[j]]$x, y=ROIS[[i]][[j]]$y,z=ROIS[[i]][[j]]$z, type="l", col = col_clusterings[i], add = TRUE)
    }
  }
  
}




