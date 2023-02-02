
library(longformer)

#################
# load database #
#################

db <- load_database("thlhp_db/data")
# commit 8cfd9ee, 3 Aug 2021

print("loaded thlhp_db")

reg <- db$people
cens <- db$visits
coms <- db$communities

# germline changes
reg$population[is.na(reg$population)] <- "Tsimane"

# label miscarriages / stillbirths
reg$is_miscarriage_stillbirth <- as.numeric(reg$first_name == "Aborto")
drop <- which(reg$first_name == "Aborto")
expect_true(length(drop) == 591)

############################
## start with medical data #
############################

m <- db$medical_cheqseg

o <- order(m$Fecha)
m <- m[o,]

m$Edad[m$Edad > 200] <- NA
m$date_of_birth <- reg$date_of_birth[match(m$pid, reg$pid)]
m$age_reg <- as.numeric(as.Date(m$Fecha) - as.Date(m$date_of_birth))/365

# clean height
m$TallaDePie <- as.numeric(m$TallaDePie)
m$TallaDePie[which(m$TallaDePie > 250 | m$TallaDePie < 110)] <- NA

m$TallaDePie[which(m$TallaDePie < 130 & m$age_reg > 18)] <- NA
m$TallaDePie[which(m$TallaDePie > 180 & m$age_reg < 18)] <- NA

tar <- which(m$TallaDePie > 110 + 5 * m$age_reg)
m$TallaDePie[tar] <- NA

check <- sort(unique(m$pid[which(m$TallaDePie > 178)]))
# bad: 7GMW7R, TLYP3M, BTJPCA

drop <- which(m$vid %in% c("7GMW7R", "TLYP3M", "BTJPCA"))
m$TallaDePie[drop] <- NA

# clean weight
m$Peso <- as.numeric(m$Peso)
m$Peso[which(m$Peso > 130 | m$Peso < 30)] <- NA

# child weights are not in kilograms, they must be in grams or something

tar <- which(m$Peso > 100 & !m$pid %in% c("NPTW", "K3PL", "TAN3", "JCWD", "VTPD", "4HTQ"))
m$Peso[tar] <- NA
# nearly all of these are BS, but a few stand out as legit-looking
# and a few i cant' tell since they involve individual measurments...

tar <- which(m$pid == "DAND")
m$Peso[tar] <- NA
m$Grasa[tar] <- NA
m$TallaDePie[tar] <- NA

#### temporary renamings

m$name <- paste(m$Nombre, m$Apellido1, m$Apellido2)
m$date <- m$Fecha
m$height <- m$TallaDePie
m$weight <- m$Peso
m$doctor <- m$Medico
m$community_id <- m$ComID
m$community <- m$Comunidad
m$family_id <- m$FamID
m$systolic_bp <- m$PAenHGSys
m$diastolic_bp <- m$PAenHGDias
m$body_temp <- m$Temperatura

######### create observations table from medical #########

obs <- select(m,
  pid, vid, name, date, height, weight
)

obs <- rename(obs,
  height_med = height,
  weight_med = weight
)

o <- order(obs$date)
obs <- obs[o,]


##########################
# add anthropometry data #
##########################

a <- db$anthropometry

a$date_of_birth <- reg$date_of_birth[match(a$pid, reg$pid)]
a$age_reg <- as.numeric(as.Date(a$date) - as.Date(a$date_of_birth))/365

# unbelievable values
a$weight[a$weight > 300] <- NA

tar <- which(a$age_reg < 8 & a$weight > 50)
a$weight[tar] <- NA

tar <- which(a$weight > 100 & !a$pid %in% !m$pid %in% c("NPTW", "K3PL", "TAN3", "JCWD", "VTPD", "4HTQ"))
# none of these look credible
check <- sort(unique(a$pid[tar]))
a$weight[tar] <- NA

a$name <- paste(a$first_name, a$last_name_1, a$last_name_2)

######### now add observations from anthropometry #########

# so <1% of vids are duplicated here, again can treat as one row per vid

obs$height_a <- NA
obs$weight_a <- NA

# fill out existing observations

tar <- which(obs$vid %in% a$vid)
link <- match(obs$vid[tar], a$vid)

obs$height_a[tar] <- a$height[link]
obs$weight_a[tar] <- a$weight[link]

# add new observations

tar <- which(!is.na(a$vid) & !a$vid %in% obs$vid)

new_obs <- a[tar,]

new_obs <- select(new_obs, date, vid, pid, name, height, weight)

new_obs <- rename(new_obs,
  height_a = height,
  weight_a = weight
)

obs <- bind_rows(obs, new_obs)

o <- order(obs$date)
obs <- obs[o,]

##########################################
# add medical seguimiento from 2002-2006 #
##########################################

mold <- db$medical_seguimiento # medical from 2002 to 2006

mold$name <- paste(mold$nombre, mold$apellido1, mold$apellido2)
mold$date <- mold$fecha

mold$height <- as.numeric(mold$Altura) # warning
mold$weight <- as.numeric(mold$Peso) # warning

######### now add observations from old medical #########

# about 4% of observations have duplicated vids, meh

obs$height_mold <- NA
obs$weight_mold <- NA

# fill out existing observations

tar <- which(!is.na(obs$vid) & obs$vid %in% mold$vid)
link <- match(obs$vid[tar], mold$vid)

obs$height_mold[tar] <- mold$height[link]
obs$weight_mold[tar] <- mold$weight[link]

# add new observations

tar <- which(!mold$vid %in% obs$vid & !is.na(mold$vid))

new_obs <- mold[tar,]

new_obs <- select(new_obs, date, vid, pid, name, height, weight)

new_obs <- rename(
  new_obs,
  weight_mold = weight,
  height_mold = height
)

obs <- bind_rows(obs, new_obs)

o <- order(obs$date)
obs <- obs[o,]


# ############################
# # add apoe4 data to people #
# ############################

library(readxl)

apoe4 <- read.csv("dna_fertility_db.csv")
apoe4_new <- read_xlsx("Data_export_7Sept21.xlsx")

reg$apoe_genotype <- apoe4$apoe_e234[match(reg$pid, apoe4$pid)]
reg$in_sample <- as.numeric(reg$pid %in% apoe4_new$pid)

#####################
# define last visit #
#####################

reg$date_last_obs <- NA

for (i in 1:nrow(reg)) {
  my_obs_dates <- cens$date[which(cens$pid == reg$pid[i] & cens$has_data)]
  if (length(my_obs_dates) > 0) {
    reg$date_last_obs[i] <- max(my_obs_dates)
  }
  if (i %% 100 == 0) print(i)
}


#######################
# transform variables #
#######################

# clean up dates of observation
obs$date <- gsub("\\s*", "", obs$date)

drop <- which(is.na(obs$date))
length(drop)
expect_true(length(drop) == 1)
obs <- obs[-drop,]

drop <- which(is.na(obs$pid))
length(drop)
expect_true(length(drop) == 318)
obs <- obs[-drop,]

expect_true(sum(is.na(as.Date(obs$date))) == 0)
expect_true(all(nchar(obs$date) == 10))
expect_true(sum(is.na(obs$pid)) == 0)
expect_true(sum(is.na(obs$vid)) == 0)

obs <- select(obs, vid, date, pid, name, everything())

# combine height
# three sources of height:
# "height_med"  "height_a"    "height_mold"
# medical looks the most up-to-date in this case
obs$height <- obs$height_med
tar <- which(is.na(obs$height))
obs$height[tar] <- obs$height_a[tar]
tar <- which(is.na(obs$height))
obs$height[tar] <- obs$height_mold[tar]
drop_cols <- c("height_med", "height_a", "height_mold")
obs <- select(obs, -all_of(drop_cols))

# combine weight
obs$weight <- obs$weight_med
tar <- which(is.na(obs$weight))
obs$weight[tar] <- obs$weight_a[tar]
tar <- which(is.na(obs$weight))
obs$weight[tar] <- obs$weight_mold[tar]
tar <- which(is.na(obs$weight))
drop_cols <- c("weight_med", "weight_a", "weight_mold")
obs <- select(obs, -all_of(drop_cols))

####################
# final subsetting #
####################

people <- reg[which(reg$population == "Tsimane"),]
# subset to only Tsimane
drop <- which(!obs$pid %in% people$pid)
obs <- obs[-drop,]

obs$date_of_birth <- people$date_of_birth[match(obs$pid, people$pid)]
obs$age <- as.numeric(as.Date(obs$date) - as.Date(obs$date_of_birth))/365
obs$age <- floor(obs$age)

######################
# assess missingness #
######################

# expect_true(nrow(obs) == 64694)
# expect_true(nrow(people) == 25956)
# expect_true(mean(is.na(obs$family_id)) < 0.22)
# expect_true(mean(is.na(obs$community_id)) < 0.01)
# expect_true(mean(is.na(obs$spanish)) < 0.66)
# expect_true(mean(is.na(obs$school)) < 0.71)
# expect_true(mean(is.na(obs$body_temp)) < 0.36)
# expect_true(mean(is.na(obs$resting_hr)) < 0.77)
# expect_true(mean(is.na(obs$height)) < 0.33)
# expect_true(mean(is.na(obs$weight)) < 0.35)
# expect_true(mean(is.na(obs$doctor)) < 0.37)
# expect_true(mean(is.na(obs$age)) < 0.01)


###############
# save output #
###############

write.csv(people, "./people.csv", row.names = FALSE)
write.csv(obs, "./observations.csv", row.names = FALSE)

print("analysis tables saved")
