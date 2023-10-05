
# MIT1003 250Hz
# Fixation data preprocessing 
# into the tabular format with the following columns:
# c("Participant", 
#"Stimulus"
#"Event_Start_Trial_Time_ms", 
#"Event_Ende_Trial_Time_ms", 
#"Fixation_Position_X_px", 
#"Fixation_Position_Y_px",
#"Duration")

library(R.matlab)
library(saccades)
library(rprojroot)
library(data.table)
root = find_root(is_git_root)

# MIT1003 **************************************************************!!

Hz = 250
# get Participant and stimulus name
Participants = c('ajs', 'CNG','emb','ems','ff','hp','jcw',
                 'jw','kae','krl','po','tmj','tu','ya','zb')
Stimuli.mat = list.files(file.path(root, 'data', 'MIT1003','DATA','DATA','ajs'))
Stimuli = sapply(Stimuli.mat, function(x) unlist(strsplit(x,split = '[.]'))[1])

prtcpt_stml_cmb = data.table(Participants = rep(Participants, each = 1003),
                             Stimuli = rep(Stimuli, 15))

Participant = 'emb'
Stimulus = 'i1032358'

Participant = 'ajs'
Stimulus = 'i1032358'

getFixation <- function(Participant, Stimulus){
  MIT1003_path = file.path(root, 'data', 'MIT1003','DATA','DATA',Participant,
                           paste0(Stimulus,'.mat'))
  MIT1003 = readMat(MIT1003_path)
  #fixationData  = MIT1003[[1]][,,1]$DATA[,,1]$fixationData[,,1]
  
  if(TRUE){
    #get fixation from eyedata 
    eyeData  = MIT1003[[1]][,,1]$DATA[,,1]$eyeData
    n = dim(eyeData)[1]
    Hz = 250
    range(eyeData)
    eyeDataSamples = data.frame(time = (0:(n-1))/Hz, trial = 'trial', x = eyeData[,1], y = eyeData[,2])
    #fix = detect.fixations(eyeDataSamples) #ERROR if no saccades.
    try.res = try(fix <- detect.fixations(eyeDataSamples))
    if(class(try.res)=="try-error"){
      res = data.table(Participant = Participant,
                       Stimulus = Stimulus,
                       Fixation_Position_X_px = NA,
                       Fixation_Position_Y_px = NA,
                       Event_Start_Trial_Time_ms = NA, 
                       Event_Ende_Trial_Time_ms = NA,  
                       Duration = NA)
    } else{
      fix = fix[!(fix$x < 0 | fix$y < 0),]
      res = data.table(Participant = Participant,
                       Stimulus = Stimulus,
                       Fixation_Position_X_px = fix$x,
                       Fixation_Position_Y_px = fix$y,
                       Event_Start_Trial_Time_ms = fix$start, 
                       Event_Ende_Trial_Time_ms = fix$end,  
                       Duration = fix$dur)
      }
  
  }
  return(res)
}

#getFixation('ajs','i1032358')
#getFixation('emb','i1032358')






Fixation_MIT1003 = do.call(rbind,apply(prtcpt_stml_cmb, 1, function(x) getFixation(x[1], x[2])))

save(Fixation_MIT1003, file = file.path(root, 'data','Fixation_MIT1003.Rda'))



