# Preparing the different data sets for use of mem

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

############################
# France:

# get data:
dat_fr0 <- read.csv("../Data/raw/ili_fr.csv")

# bring in format for mem:
dat_fr <- data.frame(value = dat_fr0$inc)
dat_fr$location <- dat_fr0$geo_name
dat_fr$week <- dat_fr0$week %% 100
dat_fr$year <- (dat_fr0$week - dat_fr$week)/100
dat_fr$season[dat_fr$week <= 40] <- paste0(dat_fr$year[dat_fr$week <= 40] - 1, "/", dat_fr$year[dat_fr$week <= 40])
dat_fr$season[dat_fr$week > 40] <- paste0(dat_fr$year[dat_fr$week > 40], "/", dat_fr$year[dat_fr$week > 40] + 1)
dat_fr <- dat_fr[order(dat_fr$location, dat_fr$year, dat_fr$week), ]

# remove half-seasons at the ends and pandemic season:
dat_fr <- subset(dat_fr, year <= 2019)
dat_fr <- subset(dat_fr, !(season %in% c("1984/1985", "2009/2010", "2019/2020")))
# exclude Corse:
dat_fr <- subset(dat_fr, location != "CORSE")

# split by region to do some additional stuff:
dat_fr_scaled <- NULL
regions <- unique(dat_fr$location)

# a place to store peak values:
peaks <- matrix(NA, ncol = 33, nrow = length(regions))
rownames(peaks) <- regions

# run through regions:
for (reg in regions) {
  # subset to region
  subs <- subset(dat_fr, location == reg)
  # create variable week_in_season
  for(i in 1:nrow(subs)) subs$week_in_season[i] <- sum(subs$season[1:i] == subs$season[i])
  # re-scale:
  # compute scaling factors
  peaks_reg <- aggregate(subs$value, by = list(subs$season), FUN = max)
  factor <- 100/mean(peaks_reg$x)
  # apply scaling
  subs$value <- factor*subs$value
  peaks[reg, ] <- factor*peaks_reg$x
  
  # add to data set to be written out:
  if(is.null(dat_fr_scaled)){
    dat_fr_scaled <- subs
  }else{
    dat_fr_scaled <- rbind(dat_fr_scaled, subs)
  }
}

# check means of peaks are indeed 100 in each region:
rowMeans(peaks)
apply(peaks, 1, sd)

# introduce variable combining region and season (needed for re-shaping)
dat_fr_scaled$location_season <- paste(dat_fr_scaled$location, dat_fr_scaled$season, sep = "_")
# select relevant columns
dat_fr_scaled <- dat_fr_scaled[, c("location_season", "week_in_season", "value")]
# restrict to first 30 weeks of season:
dat_fr_scaled <- subset(dat_fr_scaled, week_in_season <= 30)
# adapt column names:
colnames(dat_fr_scaled)[colnames(dat_fr_scaled) == "value"] <- "inc"

# bring to wide format:
dat_fr_wide <- reshape(dat_fr_scaled, direction = "wide", idvar = "week_in_season", timevar = "location_season", v.names = "inc")

# select relevant columns:
dat_for_mem_fr <- (dat_fr_wide[, grepl(pattern = "/", colnames(dat_fr_wide), fixed = TRUE)])
# re-name columns:
colnames(dat_for_mem_fr) <- gsub("/", ".", gsub("inc.", "", colnames(dat_for_mem_fr)))

# write out:
write.csv(dat_for_mem_fr, file = "../Data/for_mem/ili_mem_fr.csv", row.names = FALSE)


############################
# US:

# get data:
dat_us_states <- read.csv("../Data/raw/ili_us_states.csv")
dat_us_national <- read.csv("../Data/raw/ili_us_national.csv")
dat_us_national$REGION <- "National"
dat_us0 <- rbind(dat_us_states, dat_us_national)

head(dat_us0)
# bring in format for mem:
dat_us <- dat_us0[, c("REGION", "YEAR", "WEEK", "X..WEIGHTED.ILI")]
colnames(dat_us) <- c("location", "year", "week", "value")
dat_us$season <- round(dat_us$year + dat_us$week/52) - 1
dat_us$season <- paste(dat_us$season, dat_us$season + 1, sep = "/")
dat_us <- subset(dat_us, !season %in% c("1997/1998", "2009/2010", "2018/2019"))

# split by region to do some additional stuff:
dat_us_scaled <- NULL
regions <- unique(dat_us$location)

# a place to store peak values:
peaks <- matrix(NA, ncol = length(unique(dat_us$season)), nrow = length(regions))
rownames(peaks) <- regions

# run through regions:
for (reg in regions) {
  # subset to region
  subs <- subset(dat_us, location == reg)
  # create variable week_in_season
  for(i in 1:nrow(subs)) subs$week_in_season[i] <- sum(subs$season[1:i] == subs$season[i])
  # re-scale:
  # compute scaling factors
  peaks_reg <- aggregate(subs$value, by = list(subs$season), FUN = max)
  factor <- 100/mean(peaks_reg$x)
  # apply scaling
  subs$value <- factor*subs$value
  peaks[reg, ] <- factor*peaks_reg$x
  
  # add to data set to be written out:
  if(is.null(dat_us_scaled)){
    dat_us_scaled <- subs
  }else{
    dat_us_scaled <- rbind(dat_us_scaled, subs)
  }
}

# check means of peaks are indeed 100 in each region:
rowMeans(peaks)
apply(peaks, 1, sd)

# introduce variable combining region and season (needed for re-shaping)
dat_us_scaled$location_season <- paste(dat_us_scaled$location, dat_us_scaled$season, sep = "_")
# remove spaces:
dat_us_scaled$location_season <- gsub(" ", "", dat_us_scaled$location_season, fixed = TRUE)

# select relevant columns
dat_us_scaled <- dat_us_scaled[, c("location_season", "week_in_season", "value")]
# restrict to first 30 weeks of season:
# dat_us_scaled <- subset(dat_us_scaled, week_in_season <= 30)
# adapt column names:
colnames(dat_us_scaled)[colnames(dat_us_scaled) == "value"] <- "inc"

# bring to wide format:
dat_us_wide <- reshape(dat_us_scaled, direction = "wide", idvar = "week_in_season", timevar = "location_season", v.names = "inc")

# select relevant columns:
dat_for_mem_us <- (dat_us_wide[, grepl(pattern = "/", colnames(dat_us_wide), fixed = TRUE)])
# re-name columns:
colnames(dat_for_mem_us) <- gsub("/", ".", gsub("inc.", "", colnames(dat_for_mem_us)))

# fill in NAs with zeros:
dat_for_mem_us[is.na(dat_for_mem_us)] <- 0

# write out:
write.csv(dat_for_mem_us, file = "../Data/for_mem/ili_mem_us.csv", row.names = FALSE)

