#################
## PREPARING R ##
#################

# Removes all objects currently stored (equivalent to 'clear' in Stata)
rm(list=ls())

# Sets seed to assure full replicability
set.seed(1) 

# Installs packages necessary to calculate IOp, tres and forests (equivalent to 'ssc install xxx' in Stata)
# You only need to run the next three lines once, after that, you can put a '#' in front of the line
# install.packages("party") # To calculate trees and forests
# install.packages("IC22") # To calculate MLD
# install.packages("haven") # To load .dta files

# Loads the packages. You need to run those lines once every time you open R
library("party")
library("IC2")
library("haven")

# Set working directory (use "/" not "\")
setwd("C:/Users/bqp576/Desktop/ML/Sebastian")



################################
## LOADING & RESTRICTING DATA ##
################################

# Load the data
sample = read_dta("PISA_clean_v1.dta")

# Restrict sample to only 2012. You can remove this or make further restrictions based on treatment and control variables as you wish.
subsample = sample[sample$year==2012,]

# Outcome variable
outcome = "stdpvmath3"

# Circumstances. Edit as you please
circumstances = "gender+age+migr_background+lang_home+high_isced+nbook2_pisasum+fam_Hiseini_all+famstruc+f_labsit4+m_labsit4+privtut+
  first_language+migr_1+migr_2+migr_3+f_isced+m_isced+nbook_stand2+f_isei_all+m_isei_all+f_abi+m_abi+
  wealth+escs+hedres+cultpos+
  st_wl_ownrm+st_wl_stdypla+st_wl_internt+st_wl_txtbook+st_wl_classlit+st_wl_nhandy+st_wl_ntv+st_wl_ncomp+st_wl_ncar+st_wl_nmusins"

# Counts number of circumstances. Will be useful later.
(numberofcircumstances = sum(charToRaw(circumstances) == charToRaw('+'))+1) 

####################################
## CREATING AN 'OPPORTUNITY TREE' ##
####################################

# Converting categorical circumstances such that R understands they are categorical. 
# These were the only two I could find. Binary variables or ordinal variables are fine.
subsample$f_labsit4=as.factor(subsample$f_labsit4)
subsample$m_labsit4=as.factor(subsample$m_labsit4)


# Type '?ctree' and '?ctree_control' for documentation (equivalent to 'help ctree' in Stata)
opportunitytree = ctree(as.formula(paste('subsample[[outcome]] ~', circumstances)), # Outcome variables and circumstances, defined above
                        data=subsample, # Name of dataset
                        control=ctree_control(mincriterion=0.99, # 1 minus p-value, governs when to make splits
                                              testtype="Bonferroni", # Corrects for multiple hypothesis testing
                                              maxsurrogate=3, # Number of surrogate variables to be used to allocate observations into buckets when they are missing in the splitting variable
                                              maxdepth=10)) # The maximum number of steps from the top to the bottom of the tree
                                                            # The trees below will look prettier if you restrict this to, say, 4, 
                                                            # but I would always use a very high number when calculating IOp.


# Store type means from the tree
subsample$typemeans_tree <- predict(opportunitytree)[,1]

# Calculate R-squared
RSS_tree = sum(subsample$w_fstuwt*(subsample[[outcome]]-subsample$typemeans_tree)^2) # Residual sum of squares
TSS_tree = sum(subsample$w_fstuwt*(subsample[[outcome]]-mean(subsample[[outcome]]))^2) # Total sum of squares
(IOp_tree = 1-RSS_tree/TSS_tree) # IOp = R-squared

#####################################
## PLOTTING THE 'OPPORTUNITY TREE' ##
#####################################

plot(opportunitytree, type="simple")
# Again, this will look nicer if you rerun the tree with maxdepth=4


######################################
## CREATING AN 'OPPORTUNITY FOREST' ##
######################################

# Here I show how to produce an opportunity forest based on certain parameter choices.
# I have selected the parameters based on minimizing the out-of-bag-error. I show how to do this later in this R script.
# NB: TAKES ABOUT 1 MINUTE TO RUN ON MY COMPUTER. It will take longer if you use larger dataset.

opportunityforest = cforest(as.formula(paste('subsample[[outcome]] ~', circumstances)), # Outcome variables and circumstances, defined above,
                            data=subsample, # Name of dataset
                            control = cforest_control(mincriterion = 0, # 1 minus p-value, governs when to make splits. With forests, it often makes sense to continue splitting regardless of the p-value
                                              mtry=0.9*numberofcircumstances, # Number of circumstances to consider at each splitting point. Good choices are usually approximiately equivalent to the squareroot of the number of circumstances.
                                              fraction=0.5, # Fraction of observations to be used for each tree
                                              replace=FALSE, # Whether you sample with replacement when choosing the observations to be used for each tree.
                                              trace=TRUE, # Whether a status bar showing how many of the trees have been calculated is shown in the ouput.
                                              ntree=200, # Number of trees in the forest. This greatly influences computing time.
                                              testtype="Bonferroni", # Corrects for multiple hypothesis testing.
                                              maxsurrogate=1, # Number of surrogate variables to be used to allocate observations into buckets when they are missing in the splitting variable.
                                              minsplit=5)) # The minimum number of observations there has to be in a type before a further split is allowed. Governs the depth of the tree together with alpha.


# Store type means from the tree. Take about 10 seconds on my computer
subsample$typemeans_forest <- predict(opportunityforest)[,1]

# Calculate R-squared
RSS_forest = sum(subsample$w_fstuwt*(subsample[[outcome]]-subsample$typemeans_forest)^2) # Residual sum of squares
TSS_forest = sum(subsample$w_fstuwt*(subsample[[outcome]]-mean(subsample[[outcome]]))^2) # Total sum of squares
(IOp_forest = 1-RSS_forest/TSS_forest) # IOp = R-squared


##############################
## VARIABLE IMPORTANCE PLOT ##
##############################

# Calculates the importance of each circumstance in predicting the outcome 
# NB: TAKES ABouT 10 seconds on my computer.
varimp_absolute = varimp(opportunityforest)

# Adjust the variable importante measure, such that the circumstance with the highest importantce gets a value of 1
varimp_relative = round(varimp_absolute/max(varimp_absolute),3)

# Plot the results
dotchart(varimp_relative[order(varimp_relative)],
         main = "Circumstance Importance Plot", 
         xlab = "Circumstance importance (1 = max)")


###########################################################
## FINDING OPTIMAL FOREST BY MINIMIZING OUT-OF-BAG ERROR ##
###########################################################

# There are at least five different parameters that can be tweaked to find the optimal forests:
# 1. alpha (p-value for a split to be considered),
# 2. ntree (tree size), 
# 3. maxsurrogate (number of surrogate splits to allow for), 
# 4. mtry (circumstances considered at each split),
# 5. minsplit (minimum number of observations for a split to be allowed)
# Tuning over a 5-dimensional grid of all these variables would be incredibly time consuming (it would probably take several days).

# Based on some preliminary experimentation, below I restrict the tuning parameters as follows: 
# 1. Fix alpha to 1 (i.e. such that it is not-binding and minsplit is used to govern the depth of the trees instead.)
# 2. Fix ntree to 200. Beyond that, there is little improvments in the fit
# 3. Fix masurrogotae at 1. This seemed to give the best fit.

# I tune over a two-dimensional grid of values of mtry={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1}*numberofcircumstances 
# and values of minssplit={5,10,20,30,40,50}
tuning_ntree = c(200)
tuning_maxsurrogate = c(1)
tuning_mtry = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)*numberofcircumstances
tuning_alpha = c(1)
tuning_minsplit = c(5,10,20,30,40,50)

# Tells you how many forests you are computing 
(length = length(tuning_ntree)*length(tuning_maxsurrogate)*length(tuning_mtry)*length(tuning_alpha)*length(tuning_minsplit))

# Vector which will store the final results. We want to chose the row that minimizes the mse (mean squared error).
RESULTS = data.frame(matrix(NA,length,6))
colnames(RESULTS) <- c("ntree","maxsurrogate","mtry","alpha","minsplit","mse")

# Starts looping over all tuning parameters
# NB: THIS MAY TAKE SEVERAL HOURS TO RUN DEPENDING ON THE SPEED OF YOUR COMPUTER
row = 1
for (ntr in tuning_ntree) {
for (sur in tuning_maxsurrogate) {
for (mtr in tuning_mtry) {
for (alp in tuning_alpha) {
for (spl in tuning_minsplit) {

RESULTS[row,1]=ntr
RESULTS[row,2]=sur
RESULTS[row,3]=mtr
RESULTS[row,4]=alp
RESULTS[row,5]=spl
print(row)
print(c(ntr,sur,mtr,alp,spl))

opportunityforest = cforest(as.formula(paste('subsample[[outcome]] ~', circumstances)), # Outcome variables and circumstances, defined above,
                            data=subsample, # Name of dataset
                            control = cforest_control(mincriterion = 1-alp, # 1 minus p-value governing when to make splits. With forests, it often makes sense to continue splitting regardless of the p-value
                                                      mtry=mtr, # Number of circumstances to consider at each splitting point. Good choices are usually approximiately equivalent to the squareroot of the number of circumstances.
                                                      fraction=0.5, # Fraction of observations to be used for each tree
                                                      replace=FALSE, # Does not sample with replacement when choosing the observations to include
                                                      trace=TRUE, # Allows you to trace how many trees have been calculated
                                                      ntree=ntr, # Number of trees in the forest. This greatly influences computing time
                                                      testtype="Bonferroni", # Corrects for multiple hypothesis testing
                                                      maxsurrogate=sur, # Allows for surrogate variables to be used to allocate observations into buckets when they are missing in the splitting variable
                                                      minsplit=spl)) # The minimum number of observations there has to be in a type before a further split is allowed. Governs the depth of the tree together with alpha.

subsample$predictedvalue=predict(opportunityforest,OOB=TRUE)
subsample$squareerror <- (subsample[[outcome]]-subsample$predictedvalue)^2
# STORES RMSE IN THE RELEVANT CELL
RESULTS[row,6] <- mean(subsample$squareerror^(1/2), na.rm=TRUE)
save(RESULTS,file="TuningForest.Rda")

row = row + 1

}
}
}
}
}

# Optimal parameters
# Optimal mtry: 0.9*numberofcircumstances
RESULTS[RESULTS$mse==min(RESULTS$mse),3]
# Optimal minsplit: 5
RESULTS[RESULTS$mse==min(RESULTS$mse),5]