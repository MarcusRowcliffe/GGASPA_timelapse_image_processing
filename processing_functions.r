# Read in a csv file of exif data
# - sort to ensure correct order
# - Reformat date-time data to allow numeric actions
# - add full filepath
# - add stepwise time differences (timestep)
# - add time of day
# - add date
read.exifcsv <- function(file, ...){
  file %>%
    read.csv(...) %>%
    arrange(Directory, FileName) %>%
    mutate(timestamp=as.POSIXct(timestamp,
                                tz="UTC"),
           FilePath = paste(Directory, FileName, sep="/"),
           time = hms::as_hms(timestamp),
           date = as.Date(timestamp))
}

# Visualise operational times on a single panel (takes a while)
plot_schedule <- function(dat){
  par(mfrow=c(1,1), mar=c(2,4,2,1))
  locs <- sort(unique(dat$location))
  nlocs <- length(locs)
  plot(as.factor(dat$location), dat$timestamp, xaxt="n")
  axis(1, 1:nlocs, locs, las=3, cex.axis=0.7)
  axis(3, 1:nlocs, locs, las=3, cex.axis=0.7)
  for(i in 1:nlocs) lines(rep(i,2), c(0,1e12), col="grey80")
}

# Extract the nth subdirectory from each of a vector of full directory names
extract_nth_subdir <- function(dirs, n){
  splitdirs <- strsplit(dirs, .Platform$file.sep)
  unlist(lapply(splitdirs, function(y) y[n]))
}

# Plot location schedules, panel per location
plot_schedule_by_loc <- function(dat){
  par(mfrow=c(4, 6), mar=c(1,3,1,1))
  locs <- sort(unique(dat$location))
  for(l in locs){
    i <- dat$location==l
    plot(dat$timestamp[i], main=l, xaxt="n")
  }
}

# Plot time difference frequencies, panel per location
hist_timestep <- function(dat){
  par(mfrow=c(4, 6), mar=c(3,1,1,1))
  locs <- sort(unique(dat$location))
  for(l in locs){
    i <- dat$location==l
    x <- dat %>%
      filter(i & timestep>0) %>%
      pull(timestep)
    if(length(x)>1) hist(log10(x), breaks=100, main=l) else
      plot(0, type="n", main=l)
  }
}

# Difference by category - like diff but recognises a category indicator,
# but padding category switches with missing values (including last value,
# so that output length is equal to length x)
diffbycat <- function(x, cat){
  td <- as.numeric(diff(x))
  td[head(cat, -1) != tail(cat, -1)] <- NA
  c(td, NA)
}

# Modal average by frequency
modefrq <- function(x){
  tab <- sort(table(x), decr=TRUE)
  as.numeric(names(tab)[1])
}

# Timlapse index - the proportion of records falling within the n most frequent
# time difference categories
tl_index <- function(x, n=5){
  tdif_counts <- sort(table(x), decreasing = TRUE)
  if(length(x)<100)
    return(0) else{
      if(length(tdif_counts)<n) n <- length(tdif_counts)
      return(sum(tdif_counts[1:n]) / sum(tdif_counts))
    }
}

# Calculates location-specific timelapse indices, returns a dataframe with
# columns location and tli (see tl_index)
test_timelapse <- function(dat, plot=FALSE){
  tltest <- dat %>%
    group_by(location) %>%
    summarise(tli=tl_index(timestep, 5))
  if(plot==TRUE){
    locs <- sort(unique(dat$location))
    nlocs <- length(locs)
    par(mfrow=c(1,1), mar=c(2,4,2,1))
    plot(as.factor(tltest$location), tltest$tli, cex.axis=0.7, las=3)
    axis(3, 1:nlocs, locs, las=3, cex.axis=0.7)
    for(i in 1:nlocs) lines(rep(i,2), c(0,2), col="grey80")
  }
  tltest
}

# Correct timestamp in dat given location and time increment data taken from
# the ith row of resets.
# dat must contain columns:
#   location (character)
#   timestamp (posix)
# resets must contain columns:
#   location: character
#   startrow: numeric giving the row within dat after which the reset occurs
#   endrow: numeric row at which to end the correction; if NA defaults to last
#       row of that location
#   inc1, inc2 (numeric seconds by which to adjust timestamps for respectively
#     the first and subsequent timestamps relative to the last timestamp before
#     the reset).
# Returns dat with timestamp corrected.
correct_timestamp <- function(dat, resets){

  correct <- function(i){
    j <- resets$startrow[i]
    dat$timestamp[j+1] <- dat$timestamp[j] + resets$inc1[i]
    k <- resets$endrow[i]
    if(is.na(k))
      k <- tail(which(dat$location == resets$location[i]), 1)
    kk <- (j+2):k
    dat$timestamp[kk] <- dat$timestamp[kk] +
      difftime(dat$timestamp[j], dat$timestamp[j+2], tz="UTC", units="secs") +
      resets$inc2[i]
    dat
  }

  for(i in 1:nrow(resets)) dat <- correct(i)
  mutate(dat,
         timestep = diffbycat(timestamp, location) / 60,
         date = as.Date(timestamp),
         time = hms::as_hms(timestamp))
}

# Create dataframe of location-wise day and night time lapse intervals
get_intervals <- function(dat){
  intervals <- dat %>%
    group_by(location) %>%
    summarise(day = 60*modefrq(timestep))
  thrshld <- intervals$day[match(dat$location, intervals$location)] + 60
  nightintervals <- with(subset(dat, timestep*60 > thrshld),
                         tapply(timestep, location, modefrq) * 60)
  intervals$night <- nightintervals[match(intervals$location, names(nightintervals))]
  intervals
}

# Returns a subset of image records in dat that look like timelapse.
# Judges timelapse by setting out a regular sequence of times for each
# location/date combination (based on the location-specific day and night
# timelapse intervals contained in intervals), and finding the images with
# times that most closely match these within a tolerance threshold (tol).
# dat: A dataframe of image metadata with (at least) columns: location, date,
#   timestamp and timestep.
# intervals: a dataframe with columns location, day and night, giving
#   respectively the location ID, daytime time lapse interval and night time
#   interval
# tol: tolerance - seconds from a sequence time at which to accept an image
#   as timelapse
# verbose: if TRUE, prints location_date combination being processed
filter_timelapse <- function(dat, intervals, tol=180, verbose=FALSE){

  is.tl <- function(pd){
    if(verbose) print(pd)
    times <- subset(dat, loc_date==pd)$timestamp
    loc <- unlist(strsplit(pd, "_"))[1]
    start <- min(times, na.rm=T)
    interval <- intervals %>%
      filter(location==loc) %>%
      pull(day)
    timeseq <- seq(start, start+24*60^2, interval)
    i <- expand.grid(as.numeric(times), as.numeric(timeseq)) %>%
      apply(1, diff) %>%
      abs() %>%
      matrix(nrow=length(times)) %>%
      apply(2, function(x) which(x==min(x) & x<tol)) %>%
      unlist()
    1:length(times) %in% i
  }

  dat$loc_date <- with(dat, paste(location, date, sep="_"))
  ld <- unique(dat$loc_date)
  index <- pbapply::pbsapply(ld, is.tl) %>%
    unlist()
  dat %>%
    filter(index) %>%
    mutate(timestep = diffbycat(timestamp, location) / 60)
}

get_species_frequencies <- function(dat){
  res <- subdat4 %>%
    group_by(location, species) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = species,
                values_from = n)
  res[is.na(res)] <- 0
}
