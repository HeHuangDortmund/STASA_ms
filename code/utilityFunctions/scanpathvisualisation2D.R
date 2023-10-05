library(ggplot2)
library(imager)
library(rprojroot)
library(data.table)
root = find_root(is_git_root)

#background_image: A function used in ggplot for ploting the background image
background_image <- function(raster.img){
  annotation_raster(raster.img,
                    xmin = -Inf, xmax = Inf,
                    ymin = -Inf, ymax = Inf)
}



# scanpathVisualization2D: Function for visualisation the scanpath on stimulus
# Arguments:
#     -scanPath:the a single scanpath or the scanpaths of a group of participants of a stimulus.
#               It is a data.table contains the following columns:
#               - Fixation_Position_X_px
#               - Fixation_Position_Y_px
#               - Duration
#               - Participant, if it is the scanpaths for a group
#     -representScanpath: output of the function getRepresentScanpath()
#     -stimulusImgPath: file path of the stimulus
# Value: a ggplot object

scanpathVisualization2D = function(scanPath=NULL,
                                   representScanpath=NULL,
                                   stimulusImgPath){
  
  img = load.image(stimulusImgPath)
  
  if(!is.null(scanPath)){
    scanPath$Fixation_Position_Y_px = height(img) - scanPath$Fixation_Position_Y_px
  }
  
  if(!is.null(representScanpath)){
    representScanpath$representScanpath$Fixation_Position_Y_px = height(img) - representScanpath$representScanpath$Fixation_Position_Y_px
    #TODO
    # representScanpath$AOIpolygons ---Y coordinates must  also be changed 
  }
  


  #if(representative) cat('draw a representative scanpath')
  if(is.null(representScanpath)){ # only scanpath or scanpaths
    #check scanPath is from a single Paticipant or a group
    if('Participant' %in% names(scanPath)){
      if(length(unique(scanPath$Participant)) > 1){
        singleSP = FALSE
        cat('draw scanpath form a group of participants')
      }else{
        singleSP = TRUE
        cat('draw scanpath form a single participants')
      }
    }else{
      singleSP = TRUE
      cat('draw scanpath form a single participants')
    }
    
    if(singleSP == FALSE){
      # add a column named Order
      order = do.call(c,sapply(scanPath[,.N, by = Participant]$N, function(n) 1:n))
      scanPath$Order = order
      img = load.image(stimulusImgPath)
      p <- ggplot(scanPath, aes(x = Fixation_Position_X_px, 
                                y = Fixation_Position_Y_px, 
                                col = Participant)) +
        background_image(img) +
        geom_line()+
        geom_point(aes(size = Duration)) +
        ylim(0, height(img)) +
        xlim(0, width(img)) + 
        theme(legend.position = "none")
      #geom_text(aes(label = Order, y = Fixation_Position_Y_px + 10))
      
    }else{
      # add a column named order
      scanPath$Order = 1:dim(scanPath)[1]
      img = load.image(stimulusImgPath)
      randomColor = rainbow(26, s=.6, v=.9)[sample(1:26,1)]
      p <-ggplot(scanPath, aes(x = Fixation_Position_X_px, y = Fixation_Position_Y_px)) +
        background_image(img) +
        geom_line(col = randomColor)+
        geom_point(aes(size = Duration), col = randomColor)+
        ylim(0, height(img)) +
        xlim(0, width(img))+
        theme(legend.position = "none")
      
    }
    
    
  }else if(is.null(scanPath)){ # only representative scanpath ???????????????????
    img = load.image(stimulusImgPath)
    scanPath_rep = representScanpath$representScanpath
    scanPath_rep$Order = 1:(dim(scanPath_rep)[1])
    
    
    if(!is.null(representScanpath$AOIpolygons)){
      AOIs_coords = list()
      for(i in 1:length(representScanpath$AOIpolygons)){
        AOIs_coords[[i]] = representScanpath$AOIpolygons[[i]]@polygons[[1]]@Polygons[[1]]@coords
      }
      #AOIs = RepresentScanpath$AOIpolygons
      nrows = sapply(AOIs_coords, function(x) dim(x)[1])
      AOIs = rep(paste0('AOI',1:length(nrows)), nrows)
      AOIcoord = do.call(rbind,AOIs_coords)
      AOIcoord = as.data.frame(cbind(AOIcoord, AOIs))
      names(AOIcoord) = c('x', 'y', 'AOIs')
      AOIcoord$x = as.numeric(AOIcoord$x)
      AOIcoord$y = as.numeric(AOIcoord$y)
      
      #str(img)
      
      
      p <- ggplot()+
        background_image(img) +
        #geom_line(data=scanPath, aes(x = Fixation_Position_X_px, 
        #                             y = Fixation_Position_Y_px, 
        #                             col = Participant))+
        #geom_point(data=scanPath,aes(x = Fixation_Position_X_px, 
        #                             y = Fixation_Position_Y_px,
        #                             size = Duration,
        #                             col = Participant)) +
        geom_path(data = scanPath_rep, aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px), col = 'red', size = 2)+
        geom_point(data = scanPath_rep,aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px,
                                           size = Duration),col = 'red') +
        geom_label(data = scanPath_rep,aes(label = Order, 
                                          x = Fixation_Position_X_px + 30, 
                                          y = Fixation_Position_Y_px + 30),col = 'red') + 
        #geom_polygon(data = AOIcoord, aes(x = x, y = y,col = AOIs, alpha = 0)) +
        ylim(0, height(img)) +
        xlim(0, width(img))+
        theme(legend.position = "none")
      
      
    }else{ # if no AOIs in represent scanpath, than do not draw AOI polygons
      p <- ggplot()+
        background_image(img) +
        geom_path(data = scanPath_rep, aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px), col = 'red', size = 2)+
        geom_point(data = scanPath_rep,aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px,
                                           size = Duration),col = 'red') +
        geom_label(data = scanPath_rep,aes(label = Order, 
                                          x = Fixation_Position_X_px + 30, 
                                          y = Fixation_Position_Y_px + 30),col = 'red') + 
        #geom_polygon(data = AOIcoord, aes(x = x, y = y,col = AOIs, alpha = 0.2)) + 
        ylim(0, height(img)) +
        xlim(0, width(img))+
        theme(legend.position = "none")
    }
    
    
    
    # p <- ggplot(scanPath_rep, aes(x = Fixation_Position_X_px, 
    #                      y = Fixation_Position_Y_px)) +
    #   background_image(img) +
    #   geom_line(col = 'blue')+
    #   geom_point(aes(size = Duration),col = 'blue') +
    #   geom_text(aes(label = Order, y = Fixation_Position_Y_px + 1), col = 'red')+
    #   geom_polygon(data = AOIcoord, aes(x = x, y = y,col = AOIs, alpha = 0.2))
  } else{ # scanpaths representative scanpath together
    
    # add a column named order
    #order = do.call(c,sapply(scanPath[,.N, by = Participant]$N, function(n) 1:n))
    #scanPath$Order = order
    scanPath_rep = representScanpath$representScanpath
    scanPath_rep$Order = 1:(dim(scanPath_rep)[1])
    
    if(!is.null(representScanpath$AOIpolygons)){
      AOIs_coords = list()
      for(i in 1:length(representScanpath$AOIpolygons)){
        AOIs_coords[[i]] = representScanpath$AOIpolygons[[i]]@polygons[[1]]@Polygons[[1]]@coords
      }
      #AOIs = RepresentScanpath$AOIpolygons
      nrows = sapply(AOIs_coords, function(x) dim(x)[1])
      AOIs = rep(paste0('AOI',1:length(nrows)), nrows)
      AOIcoord = do.call(rbind,AOIs_coords)
      AOIcoord = as.data.frame(cbind(AOIcoord, AOIs))
      names(AOIcoord) = c('x', 'y', 'AOIs')
      AOIcoord$x = as.numeric(AOIcoord$x)
      AOIcoord$y = as.numeric(AOIcoord$y)
      scanPath_rep$Participant = 'represent'
      img = load.image(stimulusImgPath)
      p <- ggplot()+
        background_image(img) +
        geom_path(data=scanPath, aes(x = Fixation_Position_X_px, 
                                     y = Fixation_Position_Y_px, 
                                     col = Participant))+
        geom_point(data=scanPath,aes(x = Fixation_Position_X_px, 
                                     y = Fixation_Position_Y_px,
                                     size = Duration,
                                     col = Participant)) +
        geom_path(data = scanPath_rep, aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px), col = 'red', size = 2)+
        geom_point(data = scanPath_rep,aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px,
                                           size = Duration),col = 'red') +
        geom_label(data = scanPath_rep,aes(label = Order, 
                                          x = Fixation_Position_X_px + 30, 
                                          y = Fixation_Position_Y_px + 30), col = 'red', size = 2.5) + 
        #geom_polygon(data = AOIcoord, aes(x = x, y = y,col = AOIs, alpha = 0))+
        ylim(0, height(img)) +
        xlim(0, width(img))+
        theme(legend.position = "none")
    }else{
      scanPath_rep$Participant = 'represent'
      img = load.image(stimulusImgPath)
      p <- ggplot()+
        background_image(img) +
        geom_path(data=scanPath, aes(x = Fixation_Position_X_px, 
                                     y = Fixation_Position_Y_px, 
                                     col = Participant))+
        geom_point(data=scanPath,aes(x = Fixation_Position_X_px, 
                                     y = Fixation_Position_Y_px,
                                     size = Duration,
                                     col = Participant)) +
        geom_path(data = scanPath_rep, aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px), col = 'red', size = 2)+
        geom_point(data = scanPath_rep,aes(x = Fixation_Position_X_px, 
                                           y = Fixation_Position_Y_px,
                                           size = Duration),col = 'red') +
        geom_label(data = scanPath_rep,aes(label = Order, 
                                          x = Fixation_Position_X_px + 30, 
                                          y = Fixation_Position_Y_px + 30), col = 'red', size = 2.5) + 
        #geom_polygon(data = AOIcoord, aes(x = x, y = y,col = AOIs, alpha = 0.2))+
        ylim(0, height(img)) +
        xlim(0, width(img))+
        theme(legend.position = "none")
    }
    
   
  }
  return(p)
}
  


# Example

runExampleVisual = FALSE

if(runExampleVisual){
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
 
  singleStimulus = 'i1182314083' #!!!!!!!!!!!!!!!!!!
  
  Fixation_MIT1003_singleStimulus = Fixation_MIT1003[Stimulus == singleStimulus]
  stimulusImgPath = file.path(root,'data','MIT1003','ALLSTIMULI',paste0(singleStimulus,'.jpeg'))
  
  
  
  ExperimentDuration = 3  # 3 seconds
  et = 30/1000 # 30 ms
  OR_p = 0.22 # 0.22
  OR_p = 0.72
  
  source(file.path(root,'code','utilityFunctions','getConvexHulls.R'))
  source(file.path(root,'code','utilityFunctions','getClusteringLabelsForEachSection.R'))
  source(file.path(root,'code','utilityFunctions','getRepresentScanpath.R'))
  
  
  
  convex_hull = getConvexHulls(Fixation_MIT1003_singleStimulus,
                               ExperimentDuration,
                               et)
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  
  
  
  scanpathVisualization2D(scanPath = Fixation_MIT1003_singleStimulus, 
                          stimulusImgPath = stimulusImgPath)  
  
  
  
  scanpathVisualization2D(scanPath = Fixation_MIT1003_singleStimulus,
                          representScanpath = RepresentScanpath,
                          stimulusImgPath = stimulusImgPath)    
  
  ###############################Multiplots###############################################
  library(cowplot)
  
  
  
  source(file.path(root,'code','utilityFunctions','getConvexHulls.R'))
  source(file.path(root,'code','utilityFunctions','getClusteringLabelsForEachSection.R'))
  source(file.path(root,'code','utilityFunctions','getRepresentScanpath.R'))
  
  
  
  convex_hull = getConvexHulls(Fixation_MIT1003_singleStimulus,
                               ExperimentDuration,
                               et)

  
  OR_p = 0.42
  
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  similarities4 = getClusteringLabelsForEachSection2(convex_hull, OR_p)[[2]]
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  rep4 = RepresentScanpath$representScanpath
  
  p4 <-
    scanpathVisualization2D(scanPath = Fixation_MIT1003_singleStimulus,
                            representScanpath = RepresentScanpath,
                            stimulusImgPath = stimulusImgPath) 
  
  OR_p = 0.62
  
  
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  similarities6 = getClusteringLabelsForEachSection2(convex_hull, OR_p)[[2]]
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  rep6 = RepresentScanpath$representScanpath
  
  
  p6 <-
    scanpathVisualization2D(scanPath = Fixation_MIT1003_singleStimulus,
                            representScanpath = RepresentScanpath,
                            stimulusImgPath = stimulusImgPath) 
  
  
  OR_p = 0.82
  
  
  clusterings = getClusteringLabelsForEachSection(convex_hull, OR_p)
  similarities8 = getClusteringLabelsForEachSection2(convex_hull, OR_p)[[2]]
  RepresentScanpath = getRepresentScanpath(convex_hull, clusterings, et,
                                           Fixation_MIT1003_singleStimulus)
  rep8 = RepresentScanpath$representScanpath
  
  
  p8 <-
    scanpathVisualization2D(scanPath = Fixation_MIT1003_singleStimulus,
                            representScanpath = RepresentScanpath,
                            stimulusImgPath = stimulusImgPath) 
 
  library(latex2exp)

  p4 <- p4 + 
    #labs(title = expression(theta==0.4)) + 
    xlab('') + ylab('') +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  p6 <- p6 + 
    #labs(title = expression(theta==0.6)) + 
    xlab('') + ylab('') +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  p8 <- p8 + 
    #labs(title = expression(theta==0.8)) + 
    xlab('') + ylab('') +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())

  
  #plot_grid(p2, p4, p6, p8, nrow = 2)
  #ggsave("C:/Users/huang/Documents/He/Clustering and Permutationtest for trajectories/DDBSCAN_papaja/results_theta/thetas_label.png",
  #       width = 12,
  #       height = 6)
  
  
  #plot_grid(p4, p6, p8, nrow = 1)
  #plot_grid(p4, p6, p8, nrow = 3)
  
  similarities = data.frame(similarities = similarities4, time = ((1:length(similarities4))-1)*0.03)
  
  breaks4 = unique(sort(c(0, rep4$Event_Start_Trial_Time_ms - 0.03, 3,
                   rep4$Event_Ende_Trial_Time_ms - 0.06)))
  labels4 = c(expression(b[1]),expression(e[1]),expression(b[2]),expression(e[2]),
             expression(b[3]),expression(e[3]),expression(b[4]),expression(e[4]),
             expression(b[5]),expression(e[5]),expression(b[6]),expression(e[6]),
             expression(b[7]),expression(e[7]),expression(b[8]),expression(e[8]),
             '')
  
  p4_simi <- ggplot(similarities,aes(x = time, y = similarities))+
    geom_line()+
    geom_point()+
    geom_hline(yintercept = 0.42, col = 'red')+
    #geom_vline(xintercept = rep4$Event_Start_Trial_Time_ms - 0.03, col = 'red')+
    #geom_vline(xintercept = rep4$Event_Ende_Trial_Time_ms - 0.06, col = 'blue')+
    annotate("rect", xmin=rep4$Event_Start_Trial_Time_ms - 0.03, xmax=rep4$Event_Ende_Trial_Time_ms - 0.06, ymin=0, ymax=Inf, alpha=0.2, fill="red")+
    xlab('time [s]')+
    scale_x_continuous(breaks= breaks4, labels = labels4) +
    scale_y_continuous(breaks= c(0, 0.42, 1), labels = c(0, expression(paste(theta,'=0.4')),1))+
    theme(axis.text.y  = element_text(angle=90, colour = c('black','red','black')))
  

  
  breaks6 = unique(sort(c(0, rep6$Event_Start_Trial_Time_ms - 0.03, 3,
                          rep6$Event_Ende_Trial_Time_ms - 0.06)))
  labels6 = c('',
              expression(b[1]),expression(e[1]),expression(b[2]),expression(e[2]),
              expression(b[3]),expression(e[3]),expression(b[4]),expression(e[4]),
              expression(b[5]),expression(e[5]),
              '')
  
  p6_simi <- ggplot(similarities,aes(x = time, y = similarities))+
    geom_line()+
    geom_point()+
    geom_hline(yintercept = 0.62, col = 'red')+
    #geom_vline(xintercept = rep6$Event_Start_Trial_Time_ms - 0.03, col = 'red')+
    #geom_vline(xintercept = rep6$Event_Ende_Trial_Time_ms - 0.06, col = 'blue')+
    annotate("rect", xmin=rep6$Event_Start_Trial_Time_ms - 0.03, xmax=rep6$Event_Ende_Trial_Time_ms - 0.06, ymin=0, ymax=Inf, alpha=0.2, fill="red")+
    xlab('time [s]')+
    scale_x_continuous(breaks= breaks6, labels = labels6)+
    scale_y_continuous(breaks= c(0, 0.62, 1), labels = c(0, expression(paste(theta,'=0.6')),1))+
    theme(axis.text.y  = element_text(angle=90, colour = c('black','red','black')))
  
  
  breaks8 = unique(sort(c(0, rep8$Event_Start_Trial_Time_ms - 0.03, 3,
                          rep8$Event_Ende_Trial_Time_ms - 0.06)))
  labels8 = c('',
              expression(b[1]),expression(e[1]),expression(b[2]),expression(e[2]),
              expression(b[3]),expression(e[3]),
              '')
  
  p8_simi <- ggplot(similarities,aes(x = time, y = similarities))+
    geom_line()+
    geom_point()+
    geom_hline(yintercept = 0.82, col = 'red')+
    #geom_vline(xintercept = rep8$Event_Start_Trial_Time_ms - 0.03, col = 'red')+
    #geom_vline(xintercept = rep8$Event_Ende_Trial_Time_ms - 0.06, col = 'blue')+
    annotate("rect", xmin=rep8$Event_Start_Trial_Time_ms - 0.03, xmax=rep8$Event_Ende_Trial_Time_ms - 0.06, ymin=0, ymax=Inf, alpha=0.2, fill="red")+
    xlab('time [s]')+
    scale_x_continuous(breaks= breaks8, labels = labels8) +
    scale_y_continuous(breaks= c(0, 0.82, 1), labels = c(0, expression(paste(theta,'=0.8')),1))+
    theme(axis.text.y  = element_text(angle=90, colour = c('black','red','black')))
  
  
  plot_grid(p4,  p4_simi,p6, p6_simi, p8, p8_simi, nrow = 3)
  
  ggsave("C:/Users/huang/Documents/He/Clustering and Permutationtest for trajectories/DDBSCAN_papaja/results_theta/thetas_simi.png",
         width = 12,
         height = 9)
    
  
  


}






#ggsave(file = file.path(root, 'results','some_rsps',"istatic_boston_street_april_p1010271.png"), width = 16, height = 9)


