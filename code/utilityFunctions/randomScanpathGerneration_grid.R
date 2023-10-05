# ramdom scanpath generation.

# -  We generated firstly a random ground truth scanpath with 4 fixations on a stimulus with size 800x800, 
#    which is divide into 4x4 grids. 
# -  We selected out 4 of the 16 grid centers randomly as the 
#    positions of the fixations for the ground truth scanpath. 
# -  The Duration of each fixation is set to 0.5 second, and
#    every two fixations are separated by one third of a second. 
# -  Thus the total time for the ground truth scanpath is 3 seconds,
#    and the total fixation time is 2 seconds.

# Then we add gauss noise to the ground truth scanpath, in which
# the fixation positions were perturbed with noise drawn from a Gaussian distribution, 
# with a standard deviation of one third of the width of the grid.
# In the same way, the start time and the end time
# of the fixations were perturbed with noise drawn from a Gaussian distribution with a standard deviation
# of 1/9 seconds. 
# We randomly generate 50 such scanpaths with noises and then apply our algorithm to these scanpaths 
# to observe how much the resulting representative scanpath is similar to the ground truth scanpath.




library(mvtnorm)
library(data.table)

#randomScanpathGenerator: function to generate a scanpath randomly with giving experiment design
#Input:  ....
#      -noiseStandErrorCoeffiziet = 3, 3 sigma...
#Output: list of 2 Elements
#        - a random generated groundtruch scanpath
#        - scanpaths with gauss noise

randomScanpathGenerator = function(mySeed = 100,
                                   Xres = 800,
                                   Yres = 800 ,
                                   Xbin = 4,
                                   Ybin = 4,
                                   FixationDuration = 1/2,
                                   SaccadeDuation = 1/3,
                                   viewLast = 3,
                                   numberOfFixations = 4,
                                   numberOfscanpaths = 50,
                                   noiseStandErrorCoeffiziet = 3){ 
  
  SIGMA_spatial = min(Xres/Xbin, Yres/Ybin)/noiseStandErrorCoeffiziet  # SD of the spatial noise 
  SIGMA_temporal = SaccadeDuation/noiseStandErrorCoeffiziet            # SD of the temporal noise
  
  BinCenters_x = Xres/Xbin/2 + (0:(Xbin-1))*Xres/Xbin
  BinCenters_y = Yres/Ybin/2 + (0:(Ybin-1))*Yres/Ybin
  centerOftheBins = expand.grid(BinCenters_x, BinCenters_y)
  set.seed(mySeed)
  fixation_x_y = centerOftheBins[sample(1:nrow(centerOftheBins), numberOfFixations),]
  
  fixation_b = (0:(numberOfFixations - 1))*(FixationDuration + SaccadeDuation)
  fixation_e = (0:(numberOfFixations - 1))*(FixationDuration + SaccadeDuation) + FixationDuration
  
  Fixation_groundTruth = data.frame(cbind(fixation_x_y, fixation_b, fixation_e))
  # Aim Data structure
  # "Participant"               "Stimulus"                  "Fixation_Position_X_px"    "Fixation_Position_Y_px"   
  # "Event_Start_Trial_Time_ms" "Event_Ende_Trial_Time_ms"  "Duration"
  names(Fixation_groundTruth) = c("Fixation_Position_X_px",    
                                  "Fixation_Position_Y_px", 
                                  "Event_Start_Trial_Time_ms",
                                  "Event_Ende_Trial_Time_ms")
  
  
  Fixation_groundTruth$Duration = Fixation_groundTruth$Event_Ende_Trial_Time_ms - Fixation_groundTruth$Event_Start_Trial_Time_ms
  
  # spatial Gauss noise in 2D
  GaussNoise2D = rmvnorm(numberOfFixations, mean = c(0, 0), sigma = matrix(c(SIGMA_spatial, 0, 0, SIGMA_spatial), nrow = 2))
  # Temporal Gauss noise 
  GaussNoise1D = rnorm(numberOfFixations*2, mean = 0, sd = SIGMA_temporal)
  
  
  randomScanpaths = list()
  for (i in 1:numberOfscanpaths) {
    GaussNoise2D = rmvnorm(numberOfFixations, mean = c(0, 0), sigma = matrix(c(SIGMA_spatial, 0, 0, SIGMA_spatial), nrow = 2))
    GaussNoise1D = rnorm(numberOfFixations*2, mean = 0, sd = SIGMA_temporal)
    GaussNoise_be = matrix(GaussNoise1D, nc = 2)
    NOISE = data.frame(cbind(GaussNoise2D, GaussNoise_be))
    names(NOISE) = c("Fixation_Position_X_px",    
                     "Fixation_Position_Y_px", 
                     "Event_Start_Trial_Time_ms",
                     "Event_Ende_Trial_Time_ms") 
    
    NOISE$Duration = NOISE$Event_Ende_Trial_Time_ms - NOISE$Event_Start_Trial_Time_ms
    
    randomScanpath = Fixation_groundTruth + NOISE
    randomScanpath$Event_Start_Trial_Time_ms[randomScanpath$Event_Start_Trial_Time_ms < 0] = 0
    randomScanpath$Event_Ende_Trial_Time_ms[randomScanpath$Event_Ende_Trial_Time_ms > viewLast] = viewLast
    randomScanpath$Duration = randomScanpath$Event_Ende_Trial_Time_ms - randomScanpath$Event_Start_Trial_Time_ms
    
    randomScanpath$Participant = ifelse(i<10, paste0('P0',i), paste0('P',i))
    randomScanpaths[[i]] = randomScanpath
  }
  
  randomScanpathsInDataframe = do.call(rbind, randomScanpaths)
  randomScanpathsInDataframe = randomScanpathsInDataframe[randomScanpathsInDataframe$Duration > 0,]
  
  return(list(groundTruthScanpath = Fixation_groundTruth, randomScanpathsInDataframe = data.table(randomScanpathsInDataframe)))
}







