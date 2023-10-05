#install.packages("reticulate")
library(reticulate)
library(rprojroot)
library(data.table)
root = find_root(is_git_root)
Sys.setenv(RETICULATE_PYTHON = file.path(root,'python','Scripts'))

#reticulate::py_config()


multimatch <- import('multimatch_gaze')

# start_x, start_y, duration
s1 = data.frame(x_start, y_start, duration)
s2 = data.frame(x_start, y_start, duration)

load(file = file.path(root, 'data','Fixation_MIT1003.Rda'))


Fixation_MIT1003 = Fixation_MIT1003[!is.na(Duration),]
#Example:
s1 = Fixation_MIT1003[Participant == 'ajs' & Stimulus == 'i05june05_static_street_boston_p1010764']
s2 = Fixation_MIT1003[Participant == 'CNG' & Stimulus == 'i05june05_static_street_boston_p1010764']

fix_vector1 = data.frame(start_x = s1$Fixation_Position_X_px, start_y = s1$Fixation_Position_Y_px, duration = s1$Duration)
fix_vector2 = data.frame(start_x = s2$Fixation_Position_X_px, start_y = s2$Fixation_Position_Y_px, duration = s2$Duration)


multimatch$docomparison(fix_vector1, fix_vector2, screensize=c(1280, 720))

                        
                        