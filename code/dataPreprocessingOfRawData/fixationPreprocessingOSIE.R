library(saccades)
library(rprojroot)
library(data.table)
library(R.matlab)
root = find_root(is_git_root)



# Aim Data structure
# "Participant"               "Stimulus"                  "Fixation_Position_X_px"    "Fixation_Position_Y_px"   
# "Event_Start_Trial_Time_ms" "Event_Ende_Trial_Time_ms"  "Duration"
#READ .mat

# OSIE

OSIE_path = file.path(root ,'data','OSIE','master','data','eye',
                      'fixations.mat')



# library to read matlab data formats into R


fixations = readMat(OSIE_path)

getFixationForStimulusi = function(i){
  Stimulus_stmls_i = fixations$fixations[[i]][[1]][[1]][1,1]
  fixations_stmls_i = fixations$fixations[[i]][[1]][[2]]
  
  FIXIATIONS = list()
  for (j in 1:15) {
 
    P_j = as.data.frame(t(do.call(rbind, fixations_stmls_i[[j]][[1]])))
    names(P_j) = c("Fixation_Position_X_px", "Fixation_Position_Y_px", 'Duration')
    if(j < 10){
      Participant_j = paste0('P0',j)
    }else{
      Participant_j = paste0('P',j)
    }
    #TODO: add fake b_i and e_i
    n_fix = dim(P_j)[1]
    saccad_dur = (3000 - sum(P_j$Duration))/n_fix
    Event_Start_Trial_Time_ms = c(0,cumsum(P_j$Duration)[1:(n_fix-1)]) + (0:(n_fix-1))*saccad_dur
    Event_Ende_Trial_Time_ms = Event_Start_Trial_Time_ms + P_j$Duration
    
    Event_Start_Trial_Time_ms = round(Event_Start_Trial_Time_ms/1000,3)
    Event_Ende_Trial_Time_ms = round(Event_Ende_Trial_Time_ms/1000, 3)
    Participant = rep(Participant_j,n_fix)
    Stimulus = rep(Stimulus_stmls_i, n_fix)
    FIXIATIONS[[j]] = data.frame(Participant = Participant,
                                 Stimulus = Stimulus,
                                 Fixation_Position_X_px = P_j$Fixation_Position_X_px,
                                 Fixation_Position_Y_px = P_j$Fixation_Position_Y_px,
                                 Event_Start_Trial_Time_ms = Event_Start_Trial_Time_ms,
                                 Event_Ende_Trial_Time_ms = Event_Ende_Trial_Time_ms,
                                 Duration = P_j$Duration/1000)

    
  }
  FIXIATIONS_mat = do.call(rbind,FIXIATIONS)
  
  return(FIXIATIONS_mat)

}

n_stimuli = length(fixations$fixations)
Fixation_OSIE = do.call(rbind, lapply(1:n_stimuli, getFixationForStimulusi))
Fixation_OSIE = as.data.table(Fixation_OSIE)
save(Fixation_OSIE, file = file.path(root, 'data','Fixation_OSIE.Rda'))
