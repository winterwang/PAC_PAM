
# Reanalysis of using the data from Mabe ----------------------------------

library(R2OpenBUGS)
library(BRugs)
library(coda)

#------------------------------------------------------------------------------
# THE MODEL

# Specify the model in BUGS language, but save it as a string in R:
modelString = "
  # BUGS model specification begins ...
  model {
      # data model for the Mabe RCT data 
      # For PAC arm 
      Y_pac ~ dbin(r_pac, X_pac) 
            # Y_pac is the number of success in treatment PAC    
            # X_pac is the number of patients in treatment PAC
            # r_pac is the probability of success rate in treatment PAC
      logit(r_pac) <- lr_pac
            # lr_pac is the log(odds) of success in treatment PAC

      # For PAM arm
      Y_pam ~ dbin(r_pam, X_pam)
            # Y_pam is the number of success in treatment PAM    
            # X_pam is the number of patients in treatment PAM
            # r_pam is the probability of success rate in treatment PAM
      logit(r_pam) <- lr_pam
            # lr_pam is the log(odds) of success in treatment PAM

      # calculate odds ratio for PAM vs. PAC
      lr_pam <- lr_pac + lor
            # lor is the log(odds ratio, OR) comparing PAM against PAC
      OR <- exp(lor)
            # convert lor back to OR
      # critical values calculation 
      P.pam.better <- step(OR - 1) 
            # 1 if or >= 1, which favors PAM treatment
            # 0 if or < 1, which favors PAC treatment

      # Priors 
      lr_pac ~ dnorm(0, 0.3)  # vague priors for log success(%) in PAC arm equal 
      lor    ~ dnorm(0, 0.33) # vague priors for log(OR)
                              # this is a neutral normal prior centred on no effect; 
                              # where 95% of the prior OR values lays between 1/30, and 30
  }
  # ... BUGS model specification ends.
" # close quote to end modelString


# data from the mabe RCT 
data.list <- list(
  X_pac = 137,       # there were 137 patients recruited for PAC treatment 
  Y_pac = 89,        #             89 patients eradication succeeded
  X_pam = 169,       #            169 patients recruited for PAM treatment
  Y_pam = 163        #            163 patients eradication succeeded
)

# initial values 

initi1 <- list(lr_pac = 2.197225,  # means log(0.9/(1-0.9)) success (%) = 0.9 in PAC
               lor = 0)      # means OR = 1, no difference between PAM and PAC
initi2 <- list(lr_pac = -6.906755,  # means log(0.001/(1-0.001)) success (%) = 0.001 in PAC
               lor = 5)      # means OR = 148, strong effect difference favors PAM


# Write the modelString to a file, using R commands:
writeLines(modelString,con="model.txt")


# Check the model
modelCheck( "model.txt" )

# Load the data
modelData(bugsData(data.list))

# Complie the data with 2 chains
modelCompile(numChains = 2) 

# set initial values
modelInits(bugsData(initi1))
modelInits(bugsData(initi2))

# set monitors on nodes of interest
pars <- c("r_pac", "r_pam", "OR", "lor")#, "P.pam.better")
samplesSet(pars)


# Generate 1000 iterations
modelUpdate(1000)


samplesHistory("*", mfrow = c(4,1), ask = FALSE)
samplesHistory("*", mfrow = c(4,1), ask = FALSE, beg = 51)

# Gelman-Rubin statistics
postsamples <- buildMCMC("*")
gelman.diag(postsamples)
gelman.plot(postsamples)

# Generate another 100000 iterations 
modelUpdate(50000)
sample.statistics <- samplesStats("*", beg = 1001)
print(sample.statistics)

# #         mean      sd  MC_error val2.5pc  median val97.5pc start sample
# OR    14.6400 7.00600 0.0248700   6.0510 13.0900   32.5700  1001 100000
# lor    2.5890 0.42830 0.0015680   1.8000  2.5720    3.4830  1001 100000
# r_pac  0.6547 0.04012 0.0001444   0.5743  0.6555    0.7310  1001 100000
# r_pam  0.9595 0.01476 0.0000470   0.9261  0.9613    0.9831  1001 100000


par(mar = c(3,2,3,1))
samplesAutoC("*", mfrow = c(2, 2), chain = 1, beg = 1001,
             ask = FALSE, lag.max = 50)
samplesAutoC("*", mfrow = c(2, 2), chain = 2, beg = 1001,
             ask = FALSE, lag.max = 50)

#### PUT THE SAMPLED VALUES IN ONE R DATA FRAME:
chain <- data.frame(r_pac  = samplesSample("r_pac"),
                    r_pam  = samplesSample("r_pam"),
                    OR = samplesSample("OR"),
                    lor = samplesSample("lor"))#, 
#P.pam.better = samplesSample("P.pam.better"))


#### PLOT THE DENSITY OF THE SAMPLED VALUES

plot(density(chain$OR[1001:102000]), main = "OR chains 1:2 sample: 100000", 
     ylab = "", xlab = "", col = "red", ylim = c(0, 0.4), xlim = c(6, 200))

# call R2Openbugs ---------------------------------------------------------


# dire.dir <- "/home/takeshi/ドキュメント/githubprojects/PAC_PAM/"
# output.dir <- "bugsoutput/"
# drug.odir <- paste(dire.dir, output.dir, "mabeRCT/", sep = "")
# pars <- c("r_pac", "r_pam", "OR", "lor", "P.pam.better")
# 
# mabeRCT <- bugs(data = data.list, parameters.to.save = pars,
#                       model.file = paste(dire.dir, "model.txt", sep = ""),
#                 　　　inits = list(initi1, initi2), 
#                       n.chains = 2, n.iter = 100000,
#                       n.burnin = 1000, DIC = T, working.directory = drug.odir,
#                       codaPkg = TRUE)
