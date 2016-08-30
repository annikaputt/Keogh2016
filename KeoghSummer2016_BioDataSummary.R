############################
# KeoghSummer2016_BioDataSummary.R
# Keogh2016.Rproj

# Calculates mean length, mean weight, mean condition factor for all species per unique site
# Also calculates 95% CI for above *Assumes a normal distribution

# Files sourced:
# Files produced:

# Created August 26, 2016
# A Putt
#############################
library(plyr)

# Source the data
source("KeoghSummer2016_DataUpload.R")

# Calculate condition factor K and scale to approximately 1
biodata$ConditionFactor_K <- 100000*( (biodata$Weight_g)/(biodata$ForkLength_mm^3) )

# Summarize and calculate mean length, weight, condition, and 95% CI
biodatasummary <- ddply(biodata,c("uniquesite","Species"),summarize,
                        N = length(Species),
                        AvgLength = mean(ForkLength_mm,na.rm=TRUE),
                        AvgLengthSD = sd(ForkLength_mm,na.rm=TRUE),
                        AvgWeight = mean(Weight_g,na.rm=TRUE),
                        AvgWeightSD = sd(Weight_g,na.rm=TRUE),
                        AvgCondition = mean(ConditionFactor_K,na.rm=TRUE),
                        AvgConditionSD = sd(ConditionFactor_K,na.rm=TRUE)
)


# Add in the upper and lower 95% CI values for the summary data
CIFunc <- function(sd,N) { error <- qnorm(0.975)*sd/sqrt(N) }

biodatasummary$AvgLengthUpperCI <- biodatasummary$AvgLength + CIFunc(biodatasummary$AvgLengthSD,biodatasummary$N)
biodatasummary$AvgLengthLowerCI <- biodatasummary$AvgLength - CIFunc(biodatasummary$AvgLengthSD,biodatasummary$N)
biodatasummary$AvgWeightUpperCI <- biodatasummary$AvgWeight + CIFunc(biodatasummary$AvgWeightSD,biodatasummary$N)
biodatasummary$AvgWeightLowerCI <- biodatasummary$AvgWeight - CIFunc(biodatasummary$AvgWeightSD,biodatasummary$N)
biodatasummary$AvgCondUpperCI   <- biodatasummary$AvgCondition + CIFunc(biodatasummary$AvgConditionSD,biodatasummary$N)
biodatasummary$AvgCondLowerCI   <- biodatasummary$AvgCondition - CIFunc(biodatasummary$AvgConditionSD,biodatasummary$N)

write.csv(biodatasummary,"R_ResultOutput/BioDataSummary.csv",row.names=FALSE)


