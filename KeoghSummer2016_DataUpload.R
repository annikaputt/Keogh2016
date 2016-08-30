############################
# KeoghSummer2016_DataUpload.R
# Keogh2016.Rproj

# Uploads bio data and count summaries of the 2016 Keogh 2 and 3-pass depletion electroshocking

# Files sourced:
# Files produced:

# Created August 26, 2016
# A Putt
#############################

# Upload the bio data and check to make sure that R has selected good default classes
biodata <- read.csv("Data_CSV/BioData.csv",head=TRUE)
#str(biodata)
# Add unique column with site name
biodata$uniquesite <- as.factor(paste(biodata$Location,biodata$SiteID))

# Upload count summaries
counts <- read.csv("Data_CSV/CountSummary.csv",head=TRUE)
#str(counts)
# Add unique column with site name
counts$uniquesite <- as.factor(paste(counts$Location,counts$SiteID))

# Upload Sample Location info for areas
locations <- read.csv("Data_CSV/SamplingLocations.csv",head=TRUE)
locations$uniquesite <- as.factor(paste(locations$Location,locations$SiteID))

