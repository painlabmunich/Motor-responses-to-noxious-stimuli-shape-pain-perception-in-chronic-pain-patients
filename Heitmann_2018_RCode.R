#################################################################################################
# R code needed to perform a moderated mediation analysis as outlined in Heitmann et al., 2018. #
#################################################################################################

# Software requirements: R, including the *lme4* and *mediation* packages 
# Input: nx1 comma separated (.csv) data vector containing preprocessed single trial data (see Heitmann et al., 2018: trials only used for analysis if stimulus modality (pain/touch) correctly identified, rating > 0,  reaction time within two standard deviations of the individual mean)
# - first row contains variable headings (separated by commas): SubjectNr, Pain/Touch, Intensity, ReactionTime, Rating 
# - subsequent rows contain single trial data with (again separated by commas) the following columns: subject number, stimulus modality (1 = pain, 2 = touch), level of stimulus intensity (1-3), reaction time, rating 
# see "example_datafile.csv" for a template file


# load necessary R packages
library(lme4)
library(mediation)

# set working directory
setwd('path')

# load data
dat = read.csv("example_datafile.csv")


######################
# data preprocessing #
######################

# define variable *Pain.Touch* as a factor
dat$Pain.Touch = factor(dat$Pain.Touch, levels = 1:2, labels = c("Pain", "Touch"))

# z-transform response variables *Rating* and *ReactionTime* for each modality separately while keeping original data for later analysis
dat.Orig = dat
dat$Rating[dat$Pain.Touch == "Pain"] = scale(dat$Rating[dat$Pain.Touch == "Pain"])
dat$Rating[dat$Pain.Touch == "Touch"] = scale(dat$Rating[dat$Pain.Touch == "Touch"])
dat$ReactionTime[dat$Pain.Touch == "Pain"] = scale(dat$ReactionTime[dat$Pain.Touch == "Pain"])
dat$ReactionTime[dat$Pain.Touch == "Touch"] = scale(dat$ReactionTime[dat$Pain.Touch == "Touch"])

# center variable *Intensity* around 0 to treat as numeric variable
dat$Intensity = dat$Intensity - 2


################################
# moderated mediation analysis #
################################

# manually redefine random effects (in the later *mediation* function, a random effect cannot have the same name as the moderator)
dat$Touch = 1 * (dat$Pain.Touch == "Touch")
dat$IntTouch = dat$Intensity * dat$Touch
dat$RatTouch = dat$Rating * dat$Touch
dat$RTTouch = dat$ReactionTime * dat$Touch

# fit models needed for the moderated mediation analyses using full models including all random effects
# mediation models ( A1 = *perception-behavior-model*, A2 = *behavior-perception-model*)
A1.med.fit = lmer(Rating ~ Intensity * Pain.Touch + 
                    (1 + Intensity + Touch + IntTouch | SubjectNr), data = dat)
A2.med.fit = lmer(ReactionTime ~ Intensity * Pain.Touch + 
                    (1 + Intensity + Touch + IntTouch | SubjectNr), data = dat)
## observation models
A1.obs.fit = lmer(ReactionTime ~ Intensity * Pain.Touch + Rating * Pain.Touch + 
                    (1 + Rating + Touch + RatTouch + Intensity + IntTouch | SubjectNr), data = dat)
A2.obs.fit = lmer(Rating ~ Intensity * Pain.Touch + ReactionTime * Pain.Touch + 
                    (1 + Intensity + Touch + IntTouch + ReactionTime + RTTouch | SubjectNr)
                  , data = dat)

# perform moderated mediation analysis (since a moderator is used, separately for *Pain* and *Touch*)
# mediation analysis for A1 (the *perception-behavior-model*)
A1.ma.Pain = mediate(A1.med.fit, A1.obs.fit, treat = "Intensity", 
                     mediator = "Rating", covariates = list(Pain.Touch = "Pain"), sims = 1000)
A1.ma.Touch = mediate(A1.med.fit, A1.obs.fit, treat = "Intensity", 
                      mediator = "Rating", covariates = list(Pain.Touch = "Touch"), sims = 1000)
## Mediation analysis for A2 (the *behavior-perception-model*)
A2.ma.Pain = mediate(A2.med.fit, A2.obs.fit, treat = "Intensity", 
                     mediator = "ReactionTime", covariates = list(Pain.Touch = "Pain"), sims = 1000)
A2.ma.Touch = mediate(A2.med.fit, A2.obs.fit, treat = "Intensity", 
                      mediator = "ReactionTime", covariates = list(Pain.Touch = "Touch"), sims = 1000)


###################################
# read out results and statistics #
###################################

# show summary of all models (ACME = mediation effect, ADE = direct effect)
summary(A1.ma.Pain)
summary(A1.ma.Touch)
summary(A2.ma.Pain)
summary(A2.ma.Touch)

# calculate p-values for the differences of direct effects, mediation effects and proportions mediated between modalities (*Pain* and *Touch*)
# perception-behavior-model, mediation effect
sim_diff_A1_ACME = A1.ma.Pain$d0.sims - A1.ma.Touch$d0.sims
mean(sim_diff_A1_ACME > 0)
# perception-behavior-model, direct effect
sim_diff_A1_ADE = A1.ma.Pain$z0.sims - A1.ma.Touch$z0.sims
mean(sim_diff_A1_ADE > 0)
# perception-behavior-model, proportion mediated
sim_diff_A1_PropMed = A1.ma.Pain$n0.sims - A1.ma.Touch$n0.sims
mean(sim_diff_A1_PropMed > 0)
# behavior-perception-model, mediation effect
sim_diff_A2_ACME = A2.ma.Pain$d0.sims - A2.ma.Touch$d0.sims
mean(sim_diff_A2_ACME < 0)
# behavior-perception-model, direct effect
sim_diff_A2_ADE = A2.ma.Pain$z0.sims - A2.ma.Touch$z0.sims
mean(sim_diff_A2_ADE > 0)
# behavior-perception-model, proportion mediated
sim_diff_A2_PropMed = A2.ma.Pain$n0.sims - A2.ma.Touch$n0.sims
mean(sim_diff_A2_PropMed < 0)

# calculate standard deviation of ratings and reaction times for both modalities to transform coefficients and confidence intervals back into original units
rating.SD.pain = sd(dat.Orig$Rating[dat.Orig$Pain.Touch == "Pain"])
rating.SD.touch = sd(dat.Orig$Rating[dat.Orig$Pain.Touch == "Touch"])
RT.SD.pain = sd(dat.Orig$ReactionTime[dat.Orig$Pain.Touch == "Pain"])
RT.SD.touch = sd(dat.Orig$ReactionTime[dat.Orig$Pain.Touch == "Touch"])

# transform all coefficients back in original units, transfer proportion mediated to percentage and collect results in *results_orig_units*
results_orig_units = mat.or.vec(12,3)
# perception-behavior-model, mediation effect, *Pain*
results_orig_units[1,1] = A1.ma.Pain$d0*RT.SD.pain
results_orig_units[1,2] = A1.ma.Pain$d0.ci["2.5%"]*RT.SD.pain
results_orig_units[1,3] = A1.ma.Pain$d0.ci["97.5%"]*RT.SD.pain
# perception-behavior-model, mediation effect, *Touch*
results_orig_units[2,1] = A1.ma.Touch$d0*RT.SD.touch
results_orig_units[2,2] = A1.ma.Touch$d0.ci["2.5%"]*RT.SD.touch
results_orig_units[2,3] = A1.ma.Touch$d0.ci["97.5%"]*RT.SD.touch
cat('ACME A1 Touch: ', results_orig_units[2,1], ' ( CI = [', results_orig_units[2,2], ', ', results_orig_units[2,3], '])')
# perception-behavior-model, direct effect, *Pain*
results_orig_units[3,1] = A1.ma.Pain$z0*RT.SD.pain
results_orig_units[3,2] = A1.ma.Pain$z0.ci["2.5%"]*RT.SD.pain
results_orig_units[3,3] = A1.ma.Pain$z0.ci["97.5%"]*RT.SD.pain
# perception-behavior-model, direct effect, *Touch*
results_orig_units[4,1] = A1.ma.Touch$z0*RT.SD.touch
results_orig_units[4,2] = A1.ma.Touch$z0.ci["2.5%"]*RT.SD.touch
results_orig_units[4,3] = A1.ma.Touch$z0.ci["97.5%"]*RT.SD.touch
# perception-behavior-model, proportion mediated, *Pain*
results_orig_units[5,1] = A1.ma.Pain$n0*100
results_orig_units[5,2] = A1.ma.Pain$n0.ci["2.5%"]*100
results_orig_units[5,3] = A1.ma.Pain$n0.ci["97.5%"]*100
# perception-behavior-model, proportion mediated, *Touch*
results_orig_units[6,1] = A1.ma.Touch$n0*100
results_orig_units[6,2] = A1.ma.Touch$n0.ci["2.5%"]*100
results_orig_units[6,3] = A1.ma.Touch$n0.ci["97.5%"]*100
# behavior-perception-model, mediation effect, *Pain*
results_orig_units[7,1] = A2.ma.Pain$d0*rating.SD.pain
results_orig_units[7,2] = A2.ma.Pain$d0.ci["2.5%"]*rating.SD.pain
results_orig_units[7,3] = A2.ma.Pain$d0.ci["97.5%"]*rating.SD.pain
# behavior-perception-model, mediation effect, *Touch*
results_orig_units[8,1] = A2.ma.Touch$d0*rating.SD.touch
results_orig_units[8,2] = A2.ma.Touch$d0.ci["2.5%"]*rating.SD.touch
results_orig_units[8,3] = A2.ma.Touch$d0.ci["97.5%"]*rating.SD.touch
# behavior-perception-model, direct effect, *Pain*
results_orig_units[9,1] = A2.ma.Pain$z0*rating.SD.pain
results_orig_units[9,2] = A2.ma.Pain$z0.ci["2.5%"]*rating.SD.pain
results_orig_units[9,3] = A2.ma.Pain$z0.ci["97.5%"]*rating.SD.pain
# behavior-perception-model, direct effect, *Touch*
results_orig_units[10,1] = A2.ma.Touch$z0*rating.SD.touch
results_orig_units[10,2] = A2.ma.Touch$z0.ci["2.5%"]*rating.SD.touch
results_orig_units[10,3] = A2.ma.Touch$z0.ci["97.5%"]*rating.SD.touch
# behavior-perception-model, proportion mediated, *Pain*
results_orig_units[11,1] = A2.ma.Pain$n0*100
results_orig_units[11,2] = A2.ma.Pain$n0.ci["2.5%"]*100
results_orig_units[11,3] = A2.ma.Pain$n0.ci["97.5%"]*100
# behavior-perception-model, proportion mediated, *Touch*
results_orig_units[12,1] = A2.ma.Touch$n0*100
results_orig_units[12,2] = A2.ma.Touch$n0.ci["2.5%"]*100
results_orig_units[12,3] = A2.ma.Touch$n0.ci["97.5%"]*100
