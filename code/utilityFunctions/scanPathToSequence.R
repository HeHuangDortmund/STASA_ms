

scanPathToSequence = function(x, Xbin, Ybin, Xres, Yres){
  
  NumberOfBins = c(Xbin,Ybin) # 8 rows and 12 columns
  GridAois = matrix(nr = NumberOfBins[2], nc = NumberOfBins[1])
  for (i in 1:NumberOfBins[2]) {
    for (j in 1:NumberOfBins[1]) {
      GridAois[i,j] = paste0(letters[i],LETTERS[j])
    }
  }
  
  
  X_in_bin = ceiling(x$Fixation_Position_X_px/ (Xres/Xbin))
  X_in_bin[X_in_bin > Xbin] = Xbin
  Y_in_bin = ceiling(x$Fixation_Position_Y_px/ (Yres/Ybin))
  Y_in_bin[Y_in_bin > Ybin] = Ybin
  
  AOIs = c()
  for (i in 1:length(X_in_bin)) {
    AOIs = c(AOIs, GridAois[Y_in_bin[i],X_in_bin[i]])
  }
  
  
  # 100ms (0.1s) bins
  #number_of_Aois = floor(x$Event_Ende_Trial_Time_ms/0.1 - x$Event_Start_Trial_Time_ms/0.1)
  number_of_Aois = floor(x$Duration/0.1)
  
  AOI_sequence = rep(AOIs, number_of_Aois)
  return(AOI_sequence)
}


runExample = FALSE
if(runExample){
  library(rprojroot)
  library(data.table)
  library(imager)
  root = find_root(is_git_root)
  library(reticulate)
  
  source(file.path(root,'code','utilityFunctions','getConvexHulls.R'))
  source(file.path(root,'code','utilityFunctions','getClusteringLabelsForEachSection.R'))
  

  
  # load data Fixation_OSIE
  load(file = file.path(root, 'data','Fixation_OSIE.Rda'))
  Participants = c('P01', 'P02','P03','P04','P05','P06','P07',
                   'P08','P09','P10','P11','P12','P13','P14','P15')
  

  # remove NA
  # NA Quelle: for some Stimuli and Paticipant: there are no saccade
  #            which means that the Paticipant didn't 'move' their eyes during the experiment
  
  Fixation_OSIE = Fixation_OSIE[!is.na(Duration),]

  
  # take a single stimulus
  singleStimulus = '1001.jpg'
  
  Fixation_OSIE_singleStimulus = Fixation_OSIE[Stimulus == singleStimulus]
  spLists = split(Fixation_OSIE_singleStimulus, by = 'Participant')
  Xres = 1024
  Yres = 768
  Xbin = 12
  Ybin = 8
  scanPathToSequence(x = spLists[[1]], Xbin, Ybin, Xres, Yres)
  
}
