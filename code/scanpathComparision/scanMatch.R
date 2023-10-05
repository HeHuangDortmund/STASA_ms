# Tow Steps:
# bulid solutionmatrix and apply Needleman-Wunsch algorithmus



#thresholdCalculate 

thresholdCalculate = function( Xres,
                               Yres, 
                               Xbin,
                               Ybin,
                               Fixation_of_singleStimulus){
  spLists = split(Fixation_of_singleStimulus, by = 'Participant')
  saccade_amplitude_in_pixels_calcu = function(x){
    #x = spLists[[1]]
    max(max(abs(x$Fixation_Position_X_px - Xres/2)),
        max(abs(x$Fixation_Position_Y_px - Yres/2)))
  }
  #Threshold = 2 * std(saccade amplitude in pixels) / (Xres / Xbin)
  threshold = 2*sd(unlist(lapply(spLists, saccade_amplitude_in_pixels_calcu)))/ (Xres / Xbin)
  return(threshold)
}



#buildSolutionMatrix: function to build solution matrix related to the Euclidean
#                     distance between bins
# Par: 
#    - Xbin, Ybin: number of Xbins and Ybins
#    - Xres,Yres: resolution of the stimulus, like 1024 768
#    - threshold: Threshold = 2 * std(saccade amplitude in pixels) / (Xres / Xbin)
# Res:
#    - the SolutionMatrix

buildSolutionMatrix = function(Xbin=24, Ybin=18, Xres,Yres, threshold=3.5){
  NumberOfBins = c(Xbin,Ybin) # 8 rows and 12 columns
  GridAois = matrix(nr = NumberOfBins[2], nc = NumberOfBins[1])
  for (i in 1:NumberOfBins[2]) {
    for (j in 1:NumberOfBins[1]) {
      GridAois[i,j] = paste0(letters[i],LETTERS[j])
    }
  }
  
  solutionMatrix = matrix(0, nr = prod(NumberOfBins), nc = prod(NumberOfBins))
  rownames(solutionMatrix) = GridAois
  colnames(solutionMatrix) = GridAois
  
  for (i in 1:prod(NumberOfBins)) {
    for (j in 1:prod(NumberOfBins)) {
      
      ind_i = which(GridAois == rownames(solutionMatrix)[i],arr.ind = TRUE)
      ind_j = which(GridAois == colnames(solutionMatrix)[j],arr.ind = TRUE)
      
      
      solutionMatrix[i,j] = sqrt((ind_i[1] - ind_j[1])^2 + (ind_i[2] - ind_j[2])^2)
    }
  }
  
  # A substitution matrix based on the inverse distance between RoIs 
  #will have value ranging from 0 (between RoIs at the diagonal opposites) 
  #to sqrt(Xbin^2 + Ybin^2) ( between the two same RoIs). As a rule of thumb, 
  #it is suggest that you leave the Gap Value to 0 and instead adjust the Threshold 
  #value (the values in the matrix are equal to the 
  #old values - sqrt(Xbin^2 + Ybin^2)  - Threshold).  (or ?)
  
  #The Threshold value can be set by relating it to the variability of
  #the saccade amplitude of your data set and can be calculated as follow:
  
  solutionMatrix =  max(solutionMatrix) - solutionMatrix - threshold
  return(solutionMatrix)
}


#needleman: Function, implementation of Needleman-Wunsch Algorithmus
#Par:
#   -seq1, seq2: e.g. seq1 = c('aA','dD','cD','hE') #seq2 = c('bB','eE','dE','hE','dD')
#   -gap: gap penalty, default is zero
#   -solutionMatrix: results of buildSolutionMatrix()
#Res:
# A list of:
#    - aligned_seqs
#    - score of the match (unnormalized)
#    - score_matrix  
#    - movement_matrix 

needleman = function(seq1, seq2, gap = 0, solutionMatrix){
  
  
  
  #seq1 = c('aA','dD','cD','hE')
  #seq2 = c('bB','eE','dE','hE','dD')
  len1 = length(seq1); len2 = length(seq2)
  
  # Initialize matrix M (for scores)
  M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(M) = c("-", seq1) # assign seq chars to matrix names
  colnames(M) = c("-", seq2) # assign seq chars to matrix names
  M[1, ] = cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
  M[, 1] = cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites
  
  # Initialize matrix D (for directions)
  D = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(D) = c("-", seq1) # assign seq chars to matrix names
  colnames(D) = c("-", seq2) # assign seq chars to matrix names
  D[1, ] = rep("hor") # Fill 1st row with "hor" for horizontal moves
  D[, 1] = rep("ver") # Fill 1st col with "ver" for vertical moves
  type = c("dia", "hor", "ver") # Lookup vector
  
  # Compute scores and save moves
  for (i in 2:(len1 + 1)){# for every (initially zero) row
    for (j in 2:(len2 + 1)){# for every (initially zero) col
      
      
      hor = M[i, j - 1] + gap # horizontal move = gap for seq1
      ver = M[i - 1, j] + gap # vertical move = gap for seq2 
      dia = M[i - 1, j - 1] + solutionMatrix[rownames(M)[i], colnames(M)[j]]
      

      #dia = ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch) 
      #             M[i - 1, j - 1] + match, 
      #             M[i - 1, j - 1] + mismatch)
      M[i, j] = max(dia, hor, ver) # Save current (best) score in M
      D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
    }
  } 
  
  # Backtracing
  align1 = c(); align2 = c() # Note: length of final alignments is unknown at this point
  
  while(i > 1 && j > 1){
    
    if(D[i, j] == "dia") {
      align1 = c(rownames(M)[i], align1)
      align2 = c(colnames(M)[j], align2)
      j = j - 1; i = i - 1  # update indices
    } else if (D[i, j] == "ver") {
      align1 = c(rownames(M)[i], align1)
      align2 = c("-", align2) # vertical movement = gap for seq2
      i = i - 1 # update indices
    } else if (D[i, j] == "hor") {
      align1 = c("-", align1) # horizontal movement = gap for seq1
      align2 = c(colnames(M)[j], align2) 
      j = j - 1 # update indices
    } 
    
  }
  
  # Prepare output
  return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2),
              score = M[nrow(M), ncol(M)], score_matrix = M, movement_matrix = D))
  
}

runExample = FALSE
if(runExample){
  threshold = 3.5
  Xres = 1024
  Yres = 768
  Xbin = 24
  Ybin = 18
  solutionMatrix = buildSolutionMatrix(Xbin, Ybin, Xres,Yres, threshold)
  seq1 = c('bD','bG','bG','bG','cD','cD','hB','hF','eH','eH','hJ','dL','dL','eK','cI','cH','bL','aL','cI')
  
  seq2 = c('dF','cC','bG','aG','cD','cC','hB','hB','hF','eH','eH','hJ','eK','aL','cH','cH')
  res = needleman(seq1, seq2, 0, solutionMatrix)
  # a normalization of the score is performed
  maxLength = 20 # length of the longest sequence
  Normalized_score = res$score/(max(solutionMatrix))/maxLength
  Normalized_score
  
}






