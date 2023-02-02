# prep data for analysis

library(openxlsx)

# load databases

# obs & ppl come from Bret's thlhp script

obs <- read.csv("data/observations.csv")
ppl <- read.csv("data/people.csv")

# id_tab comes from Ben

id_tab <- read.xlsx("data/Data_export_7Sept21.xlsx")

# part 1: age-specific fertility, using "etic" data (ppl/obs)
# part 2: reproducing prox analyses, using "emic" data (id_tab)

# make an id-level dataset

d <- data.frame(id = ppl[!is.na(ppl$apoe_genotype), ]$pid) 

# merge in the apoe data

d$apoe <- ppl[match(d$id, ppl$pid), ]$apoe

# add demographic data, check everyone is "female"

d$sex <- ppl[match(d$id, ppl$pid), ]$male
d$sex <- ifelse(d$sex == 0, "female", "male")
all(d$sex == "female")

# add demographic data, birth year and death year

d$birth_year <- ppl[match(d$id, ppl$pid), ]$date_of_birth
d$birth_year <- substr(d$birth_year, 1, 4)
d$birth_year <- as.numeric(d$birth_year)

# drop the one ind without a birth year
d <- d[!is.na(d$birth_year), ]

d$death_year <- ppl[match(d$id, ppl$pid), ]$date_of_death
d$death_year[!is.na(d$death_year)] <- substr(d$death_year[!is.na(d$death_year)], 1, 4)
d$death_year[which(d$death_year == 9999)] <- NA

# overwrite the death years with NA if death is "disputed"
d$death_disputed <- ppl[match(d$id, ppl$pid), ]$death_disputed
d$death_year[d$death_disputed == 1] <- NA
d$death_disputed <- NULL
d$death_year <- as.numeric(d$death_year)

# add the most recent year an individual was observed

d$last_obs <- ppl[match(d$id, ppl$pid), ]$date_last_obs
d$last_obs <- substr(d$last_obs, 1, 4)
d$last_obs <- as.numeric(d$last_obs)

# compute a max_age variable, i.e. last obs - birth year

for (i in 1:nrow(d)) {
  
  if (is.na(d$death_year[i])) {
    d$max_age[i] <- d$last_obs[i] - d$birth_year[i]
  }
  
  if (!is.na(d$death_year[i])) {
    d$max_age[i] <- d$death_year[i] - d$birth_year[i]
  }
  
}

# remove the individual that appears -2?
table(d$max_age)

# for the individual that appears -2, correct death year to last obs

d[d$id == 855, ]$death_year <- NA
d[d$id == 855, ]$max_age <- d$last_obs[d$id == 855] - d$birth_year[d$id == 855]

nrow(d) # we have 843 individuals (772 are "in sample")
# we don't subset, because we dont need "completed" fertility

d$in_sample <- ppl[match(d$id, ppl$pid), ]$in_sample

# create the fertility data

# add "birth_year" to ppl table for calculating conceptions

ppl$birth_year <- substr(ppl$date_of_birth, 1, 4)
table(ppl$birth_year)
ppl$birth_year <- as.numeric(ppl$birth_year) # NA warning which is ok

# now get the demographic lifetime data
# for n individuals over y years of life

n <- nrow(d)
y <- max(d$max_age + 1)

# conceptions

# first check no female appears as father
intersect(d[d$sex == "female", ]$id, ppl$father_pid[!is.na(ppl$father_pid)])

con <- array(dim = c(n, y))
d$n_babies <- NA

# if female, check if they appear as "mother_id"
# conceive in the year that corresponds to a baby's birth year
# do not conceive in all the other years they are alive

# also the proximate determinants data
# afr is the min baby birth year
# alr is the max baby birth year
# ibi = alr - afr / n_babies

d$afr <- NA
d$alr <- NA
d$ibi <- NA

for (i in 1:n) {
  
  babies <- ppl[which(ppl$mother_pid == d$id[i]), ]  
  
  # we have to remove babies without a year of birth
  babies <- babies[!is.na(babies$date_of_birth), ] 
  
  # we have to remove babies with a year of birth BEFORE mother yob
  babies <- babies[babies$birth_year > d$birth_year[i], ]
  
  age_at_birth <- babies$birth_year - d$birth_year[i]
  # we want to filter women who "give birth" before age 12, or after age 59
  # these are likely errors
  age_at_birth <- age_at_birth[age_at_birth > 12 & age_at_birth < 59] 
  
  # fill in 0s when alive
  alive_till <- d$max_age[i] # same as max_age
  con[i, 1:alive_till] <- 0
  con[i, age_at_birth] <- 1
  
  # now n babies is how many years they gave birth
  d$n_babies[i] <- length(age_at_birth)
  
  if (d$n_babies[i] > 0) {
    
    d$afr[i] <- min(age_at_birth)
    d$alr[i] <- max(age_at_birth)
    d$ibi[i] <- (d$alr[i] - d$afr[i]) / d$n_babies[i]
    
  }
  
  
}

# we decided to use the data from the id_tab/ ben_tab for the prox vars

# add the matching variable

d$id_not_anon <- ppl[match(d$id, ppl$pid), ]$old_pid

d$ibi_b <- id_tab[match(d$id_not_anon, id_tab$pid), ]$avgibi
d$loss_b <- id_tab[match(d$id_not_anon, id_tab$pid), ]$abortos

# we don't have the AFR from this tab -.- 

# get the height data
d$height <- NA

for (i in 1:nrow(d)) {
  
  heights <- obs[obs$pid == d$id[i], ]$height
  d$height[i] <- median(heights, na.rm = TRUE)
  
}

# compute BMI and take the median

d$bmi <- NA
obs$bmi <- obs$weight / (obs$height/100)^2

for (i in 1:nrow(d)) {
  
  bmis <- obs[obs$pid == d$id[i], ]$bmi
  d$bmi[i] <- median(bmis, na.rm = TRUE)
  
}



# add the ben source too

d$height_b <- id_tab[match(d$id_not_anon, id_tab$pid), ]$talladepie

# make a new d, as we need for stan

d <- list(id = d$id, 
          allele = d$apoe,
          con = con, 
          n_babies = d$n_babies, 
          afr = d$afr, 
          alr = d$alr, 
          ibi = d$ibi, 
          bmi = d$bmi,
          ibi_b = d$ibi_b, 
          loss_b = d$loss_b, 
          height_b = d$height_b)

save(d, file = "d.robj")



