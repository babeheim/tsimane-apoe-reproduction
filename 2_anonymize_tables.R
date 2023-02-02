
rm(list = ls())

library(testthat)
library(dplyr)

obs <- read.csv("observations.csv", stringsAsFactors = FALSE)
people <- read.csv("people.csv", stringsAsFactors = FALSE)

# consistency checks

expect_true(all(obs$pid %in% people$pid))

# anonymize individual pids
people$old_pid <- people$pid
people$pid <- sample(1:nrow(people))

################

# cascade new ids out
people$father_pid <- people$pid[match(people$father_pid, people$old_pid)]
people$mother_pid <- people$pid[match(people$mother_pid, people$old_pid)]
obs$pid <- people$pid[match(obs$pid, people$old_pid)]

# consistency checks

check <- sample(1:nrow(obs), 1)
obs$name[check]
people$first_name[people$pid == obs$pid[check]]

# drop unnecessary/identifying variables
people <- select(people, pid, old_pid, male, date_of_birth,
  date_of_death, death_disputed, father_pid, mother_pid, has_data,
  date_last_obs, apoe_genotype, in_sample, is_miscarriage_stillbirth)
obs <- select(obs, -vid, -name, -date_of_birth)

write.csv(people, "people.csv", row.names = FALSE)
write.csv(obs, "observations.csv", row.names = FALSE)

