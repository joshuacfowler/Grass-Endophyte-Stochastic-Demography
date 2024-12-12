## Authors: Josh and Tom	## Grass endophyte population model
## Purpose: Create a script that imports Endodemog data, perform all raw data manipulation to set up data lists for Survival and Growth models and for the flowering tiller and seed production models,	
## and create an .RData object that can be loaded for analysis	
## Last Update: Oct 31, 2019
######################################################
library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)
#library(rnoaa)
library(SPEI)


##############################################################################
####### Reading in endo_demog_long ------------------------------
##############################################################################
## tom path for shared dropbox
# setwd("C:/Users/tm9")

# This is the main compiled data sheet that we will merge with the flowering and seed data
LTREB_endodemog <- 
  read.csv(file = "~Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
year_factor_key <- c('2007' = 1, '2008' = 2, '2009' = 3, 
                       '2010' = 4, '2011' = 5, '2012' = 6, 
                       '2013' = 7, '2014' = 8, '2015' = 9, 
                       '2016' = 10, '2017' = 11, "2018" = 12,
                       "2019" = 13, "2020" = 14, "2021" = 15)
species_factor_key <- c("AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                        "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                        "POSY" = 7)
plot_factor_key <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5,"6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10,
                     "11" = 11, "12" = 12, "13" = 13, "14" = 14, "15" = 15, "16" = 16, "17" = 17, "18" = 18, "19" = 19, "20" = 20,
                     "31" = 21, " 32" = 22, "33" = 23, "34" = 24, "35" = 25, "36" = 26, "37" = 27, "38" = 28, "39" = 29, "40" = 30,
                     "91" = 31, " 92" = 32, " 93" = 33, " 94" = 34, " 95" = 35, " 96" = 36, " 97" = 37, "98" = 38, " 99" = 39, "100" = 40,
                     "101" = 41, "102" = 42, "103" = 43, "104" = 44, "105" = 45, "106" = 46, "107" = 47, "108" = 48, "109" = 49, "110" = 50,
                     "111" = 51, "112" = 52, "113" = 53, "114" = 54, "115" = 55, "116" = 56, "117" = 57, "118" = 58, "119" = 59, "120" = 60,
                     "121" = 61, "122" = 62, "123" = 63, "124" = 64, "125" = 65, "126" = 66, "127" = 67, "128" = 68, "129" = 69, "130" = 70,
                     "141" = 71, "142" = 72, "143" = 73, "144" = 74, "145" = 75, "146" = 76, "147" = 77, "148" = 78, "149" = 79, "150" = 80,
                     "151" = 81, "152" = 82, "153" = 83, "154" = 84, "155" = 85, "156" = 86, "157" = 87, "158" = 88)


LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t = na_if(size_t, 0)) %>% 
  mutate(size_t1 = na_if(size_t1, 0)) %>%  
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>%
  mutate(endo_01 = as.integer(case_when(species == "ELRI" & grepl("2178R", id) ~ 0, #This is for a typo where two 2017 ELRI recruits (labeled 2178R and 209_2178R) were  were labeled as eplus but should be eminus. (this still leaves the original endo column as labelled but corrects the endo_01 column
                                        endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
  mutate(species = case_when(species == "ELVI" & plot == 101 ~ "ELRI", species == "ELVI" & plot != 101 ~ "ELVI",  # This is for the compiled data where a ELRI data point in plot 101, tag 2004 is labelled as ELVI
                                   species == "ELRI" ~ "ELRI",
                                   species == "FESU" ~ "FESU",
                                   species == "AGPE" ~ "AGPE",
                                   species == "POAL" ~ "POAL",
                                   species == "POSY" ~ "POSY",
                                   species == "LOAR" ~ "LOAR"))  %>%    
  mutate(species_index = as.integer(recode(species, !!!species_factor_key))) %>% 
  mutate(year_t_index = as.integer(recode(year_t, !!!year_factor_key))) %>%             
  mutate(year_t1_index = as.integer(recode(year_t1, !!!year_factor_key))) %>%               
  mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, 
                                          origin == "R" ~ 1, 
                                          origin != "R" | origin != "O" ~ 1))) %>%
  mutate(plot_fixed = as.integer(case_when(species == "LOAR" & id == "39_1B" ~ "39", # This is for a copy error with LOAR individual 39_1B, which assigned it to plots 41-44, as well as 
                                 plot != "R" ~ as.character(plot), 
                                 plot == "R" ~ as.character(origin)))) %>% 
  mutate(plot_index = as.integer(recode(plot_fixed, !!! plot_factor_key))) %>% 
  mutate(surv_t1 = as.integer(case_when(surv_t1 == 1 ~ 1,
                                   surv_t1 == 0 ~ 0,
                                   is.na(surv_t1) & birth == year_t1 ~ 1,
                                   is.na(surv_t1) & size_t1 > 0 ~ 1))) %>% 
  mutate(size_t1 = case_when(birth == year_t1 & surv_t1 == 1 & is.na(size_t1) ~ 1, 
                             TRUE ~ size_t1),
         size_t = case_when(birth == year_t1-1 & !is.na(seed_t) & is.na(size_t) ~ 1,
                            TRUE ~ size_t)) %>%  # This is for a few copy errors where seedlings have no size_t1 data entered, or a 0 entered.
  filter(duplicated(.) == FALSE)
# dim(LTREB_data)

# #This is a list of all the original plot endophyte status
# # There are a few plants, particularly in LOAR that have NA's or potentially have the incorrect endo status for their plots. leaving this for now.
# LTREB_endophyte_plot_numbers <- LTREB_data %>% 
#   group_by(species, plot_fixed) %>% 
#   summarize(endophyte_status = as.integer(mean(endo_01, na.rm = T)))
#   filter(!is.nan(endophyte_status))


########################################################################################################################
###### Reading in raw excel files which have the detailed Reproduction data.---------------------------
########################################################################################################################
# read in raw data from POAL
POAL_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POALnew complete with 2016 data.xlsx", sheet = "POAL", .name_repair = "universal")
POAL_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POALnew complete with 2016 data.xlsx", sheet = "POAL (NEW) recruits", .name_repair = "universal")
POAL_data_old <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POALold complete with 2016 data.xlsx", sheet = "POAL (OLD)", .name_repair = "universal")
POAL_data_old_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POALold complete with 2016 data.xlsx", sheet = "POAL (OLD) recruits", .name_repair = "universal")

# read in data from POSY
POSY_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POSYnew complete with 2016 data.xlsx", sheet = "POSY", .name_repair = "universal")
POSY_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POSYnew complete with 2016 data.xlsx", sheet = "POSY (NEW) recruits", .name_repair = "universal")
POSY_data_old <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POSYold complete with 2016 data.xlsx", sheet = "POSY(Old)", .name_repair = "universal")
POSY_data_old_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/POSYold complete with 2016 data.xlsx", sheet = "POSY(Old)recruits", .name_repair = "universal")

# read in data from LOAR
LOAR_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR", .name_repair = "universal")
LOAR_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR recruits", .name_repair = "universal")
LOAR_data_seed2008 <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR seeds 2008", skip=1, .name_repair = "universal")
LOAR_data_seed2009 <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/LOAR_for_demog_long.xlsx", sheet = "LOAR seeds 2009", .name_repair = "universal")

# Read in data from FESU
FESU_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU", .name_repair = "universal")
FESU_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/FESU Datasheet complete 7 13 16.xlsx", sheet = "FESU recruits", .name_repair = "universal")

# Read in data from ELVI
ELVI_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI(IN) originals up to 2016.xlsx", sheet = "ELVI", .name_repair = "universal")
ELVI_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVI recruits", .name_repair = "universal")
ELVI_data_seed <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/ELVI (IN) FINAL 3 10 16 updated and checked.xlsx", sheet = "ELVISeeds", .name_repair = "universal")
ELVI_data_r <- ELVI_data_r %>% 
  mutate(tag = paste(Plot, RecruitNo, sep = "_")) 
  

# Read in data from ELRI
ELRI_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI", .name_repair = "universal")
ELRI_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI recruits", .name_repair = "universal")
ELRI_data_seed2009 <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/ELRI data up to 2016.xlsx", sheet = "ELRI infl 09", .name_repair = "universal")
ELRI_data_r <- ELRI_data_r %>% 
  mutate(tag = paste(PLOT, RecruitNo, sep = "_")) 
# Read in data from AGPE
AGPE_data <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/AGPE2016_final.xlsx", sheet = "AGPE", .name_repair = "universal")
AGPE_data_r <- read_xlsx("Dropbox/EndodemogData/rawdatafilesbyspecies/AGPE2016_final.xlsx", sheet = "AGPE recruits", .name_repair = "universal")

# A list of the columns that are redundant for now, or will be used in our seed estimates
 # pmain <- POAL_data %>%  
# dplyr::select(-survive4may08,-`Surv4/11`, -Aphids1, -aphids2, -aphid3, -contains("SB"),-notes2,
 #         -notes2__1, -notes3, -notes4, -notes5,-`notes5/4/2008`,-notes6, -notes7, 
 #         -notes8, -notes9, -othernotes3, -data, -`Planting Notes`, -TAG, -Endocheck, 
 #         -EndoDateCheck, -EndoDateCheck_Day, -EndoDateCheck_Month, -EndoDateCheck_Year,
 #         -`Coll Date`, -TotTillers11sep08, -tilleradjust, -VisualEST1, -endoyr, 
 #         -seed2surv, -seed3surv, -seed4surv, -actseed2, -contains("Est"), -contains("Hbv"), 
 #         -contains("Lvs"), -contains("Infl"), -contains("CB"), -contains("Prop"))



#################################################################################################################################################################
# Pulling out the seed production estimates. These are not measured for all plants, and so will go into a separate dataframe------------------------------
#################################################################################################################################################################
# Pulling out the seed production estimates for the "New" POAL data --------
pseed <- POAL_data %>% 
  mutate(Birth = year(as.character(Planted.Date))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA2", "seeds_InflB2", 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))

pseed$year<- ifelse(pseed$variable == "seed2007", 2007, ifelse(pseed$variable == "seed2008", 2008,ifelse(pseed$variable == "seeds_InflA2", 2009, ifelse(pseed$variable  == "seeds_InflB2", 2009, ifelse(pseed$variable  == "seed2010", 2010, ifelse(pseed$variable  == "seed2011", 2011, ifelse(pseed$variable  == "seed2012", 2012, ifelse(pseed$variable  == "seed2013", 2013,ifelse(pseed$variable == "seed2014", 2014,ifelse(pseed$variable == "seed2015", 2015,ifelse(pseed$variable  == "seed2016", 2016, NA)))))))))))

pseed1 <- pseed %>% 
  filter(!is.na(seed))
 # View(pseed1)

pspike <- POAL_data %>% 
  mutate("Birth" = year(as.character(Planted.Date))) %>% 
  mutate(spike2007 = NA, spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spikelets_InflA2", "spikelets_InflB2", "spikelets_inflA3", 
                       "spikelets_inflB3", "spikelets_inflC3", "spikelets_inflA4",
                       "spikelets_inflA5", "spikelets_inflA6", "spikelets_InflB6", 
                       "spikelets_inflC6", "spikelets_inflA7", "spikelets_InflB7", 
                       "spikelets_inflC7","spikelets_inflA8", "spikelets_InflB8", 
                       "spikelets_inflC8","spikelets_inflA9", "spikelets_InflB9", 
                       "spikelets_inflC9"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))
pspike$year<- ifelse(pspike$variable == "spike2007", 2007, ifelse(pspike$variable == "spike2008", 2008,ifelse(pspike$variable == "spikelets_InflA2", 2009, ifelse(pspike$variable  == "spikelets_InflB2", 2009, ifelse(pspike$variable  == "spikelets_inflA3", 2010, ifelse(pspike$variable  == "spikelets_inflB3", 2010, ifelse(pspike$variable  == "spikelets_inflC3", 2010, ifelse(pspike$variable  == "spikelets_inflA4", 2011,ifelse(pspike$variable == "spikelets_inflA5", 2012,ifelse(pspike$variable == "spikelets_inflA6", 2013, ifelse(pspike$variable == "spikelets_InflB6", 2013, ifelse(pspike$variable == "spikelets_inflC6", 2013, ifelse(pspike$variable  == "spikelets_inflA7", 2014, ifelse(pspike$variable == "spikelets_InflB7", 2014, ifelse(pspike$variable == "spikelets_inflC7", 2014, ifelse(pspike$variable == "spikelets_inflA8", 2015, ifelse(pspike$variable == "spikelets_InflB8", 2015, ifelse(pspike$variable ==   "spikelets_inflC8", 2015, ifelse(pspike$variable == "spikelets_inflA9", 2016, ifelse(pspike$variable == "spikelets_InflB9", 2016, ifelse(pspike$variable == "spikelets_inflC9", 2016, NA)))))))))))))))))))))
# View(pspike)
pspike1 <- pspike %>%  
  filter(!is.na(spikelets))
# View(pseed1)

pflw <- POAL_data %>% 
  mutate("Birth" = year(as.character(Planted.Date))) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
pflw$year<- ifelse(pflw$variable == "flw2007", 2007, ifelse(pflw$variable == "Flwtillers1", 2008, ifelse(pflw$variable  == "FlwTillers2", 2009, ifelse(pflw$variable  == "FlwTillers3", 2010, ifelse(pflw$variable  == "FlwTillers4", 2011, ifelse(pflw$variable  == "FlwTillers5", 2012, ifelse(pflw$variable  == "FlwTillers6", 2013,ifelse(pflw$variable == "FlwTillers7", 2014,ifelse(pflw$variable == "FlwTillers8", 2015,ifelse(pflw$variable  == "FlwTillers9", 2016, NA))))))))))

# View(pflw)

pseedmerge_ss <- full_join(pspike1, pseed1, by = c( "plot", "pos", "tag", "Endo", 
                                         "Birth","year", "tillerid"))
# View(pseedmerge_ss)

pseedmerge_ssf <- merge(pseedmerge_ss, pflw, by = c("plot", "pos", "tag", "Endo", 
                                                    "Birth","year"), all = TRUE)
# View(pseedmerge_ssf)



# Pulling out the seed, spikelet and flw tiller info from the POAL New Recruits
rseed <-POAL_data_r %>% 
  mutate("Loc'n" = NA, "TRT" = NA) %>% 
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010","seed2011", "seed2012", "seed2013", "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B"))
rseed$year<- ifelse(rseed$variable == "seed2007", 2007, ifelse(rseed$variable == "seed2008", 2008, ifelse(rseed$variable == "seed2009", 2009, ifelse(rseed$variable == "seed2010", 2010, ifelse(rseed$variable  == "seed2011", 2011, ifelse(rseed$variable  == "seed2012", 2012, ifelse(rseed$variable  == "seed2013", 2013, ifelse(rseed$variable  == "seed2014", 2014, ifelse(rseed$variable  == "seed2015", 2015, ifelse(rseed$variable  == "seed2016", 2016, NA))))))))))
rseed1 <- rseed %>% 
  filter(!is.na(seed))
# View(rseed1)

rspike <- POAL_data_r %>% 
  mutate("Loc'n" = NA, "TRT" = NA) %>% 
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", "spikelets_inflA5...20", "spikelets_inflA5...26", "spikelets_infl14", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",))
rspike$year<- ifelse(rspike$variable == "spike2007", 2007, ifelse(rspike$variable == "spike2008", 2008, ifelse(rspike$variable == "spike2009", 2009, ifelse(rspike$variable == "spike2010", 2010, ifelse(rspike$variable == "spike2011", 2011, ifelse(rspike$variable == "spikelets_inflA5...20", 2012, ifelse(rspike$variable  == "spikelets_inflA5...26", 2013, ifelse(rspike$variable  == "spikelets_infl14", 2014, ifelse(rspike$variable  == "spikelets_infl15", 2015, ifelse(rspike$variable  == "spikelets_infl16", 2016, NA))))))))))
rspike1 <- rspike %>% 
  filter(!is.na(spikelets))
# View(rspike1)

# We already have the FlwTiller data within rflw dataframe
rflw <- POAL_data_r %>%
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
rflw$year<- ifelse(rflw$variable == "flw2007", 2007, ifelse(rflw$variable == "flw2008", 2008, ifelse(rflw$variable == "flw2009", 2009, ifelse(rflw$variable == "FLWtiller10", 2010, ifelse(rflw$variable == "FLWtiller11", 2011, ifelse(rflw$variable  == "FLWtiller12", 2012, ifelse(rflw$variable  == "FLWtiller13", 2013, ifelse(rflw$variable  == "FLWtiller14", 2014, ifelse(rflw$variable  == "FLWtiller15", 2015, ifelse(rflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(rflw)


rseedmerge_ss <- full_join(rspike1, rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                              "Birth","year", "tillerid"))
# View(rseedmerge_ss)

rseedmerge_ssf <- merge(rseedmerge_ss, rflw, by = c("plot", "pos", "tag", "Endo", 
                                                    "Birth","year"), all = TRUE)
# View(rseedmerge_ssf)


# Pulling out the seed production estimates for the "Old" POAL data --------
pold_seed <- POAL_data_old %>%
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA3", 'seeds_InflB3', 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

pold_seed$year<- ifelse(pold_seed$variable == "seed2007", 2007, ifelse(pold_seed$variable == "seed2008", 2008, ifelse(pold_seed$variable  == "seeds_InflA3", 2009, ifelse(pold_seed$variable == "seeds_InflB3", 2009, ifelse(pold_seed$variable == "seed2010", 2010, ifelse(pold_seed$variable  == "seed2011", 2011,ifelse(pold_seed$variable  == "seed2012", 2012, ifelse(pold_seed$variable  == "seed2013", 2013,ifelse(pold_seed$variable == "seed2014", 2014,ifelse(pold_seed$variable == "seed2015", 2015,ifelse(pold_seed$variable  == "seed2016", 2016, NA)))))))))))

pold_seed1 <- pold_seed %>% 
  filter(!is.na(seed))

# View(pold_seed1)

pold_spike <- POAL_data_old%>% 
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate(spike2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "Spikelets_InflA1", "Spikelets_InflB1","spikelets_InflA3",
                       "spikelets_InflB3", "spikelets_InflA4","spikelets_InflB4",
                       "spikelets_InflC4","spikelets_InflA5", "spikelets_inflB5", 
                       "spikelets_inflC5", "spikelets_inflA6", 
                       "spikelets_inflA7", "spikelets_InflB7",  "spikelets_inflC7",
                       "spikelets_inflA8", "spikelets_InflB8", "spikelets_inflC8",
                       "spikelets_inflA9", "spikelets_InflB9","spikelets_inflC9",
                       "spikelets_inflA10", "spikelets_InflB10","spikelets_inflC10"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("InflA", variable) ~ "A", 
                              grepl("InflB", variable) ~ "B",
                              grepl("InflC", variable) ~ "C",
                              grepl("inflA", variable) ~ "A", 
                              grepl("inflB", variable) ~ "B",
                              grepl("inflC", variable) ~ "C"))
pold_spike$year<- ifelse(pold_spike$variable == "spike2007", 2007, ifelse(pold_spike$variable == "Spikelets_InflA1", 2008, ifelse(pold_spike$variable == "Spikelets_InflB1", 2008, ifelse(pold_spike$variable == "spikelets_InflA3", 2009, ifelse(pold_spike$variable  == "spikelets_InflB3", 2009, ifelse(pold_spike$variable  == "spikelets_InflA4", 2010, ifelse(pold_spike$variable  == "spikelets_InflB4", 2010, ifelse(pold_spike$variable  == "spikelets_InflC4", 2010, ifelse(pold_spike$variable  == "spikelets_InflA5", 2011, ifelse(pold_spike$variable  == "spikelets_inflB5", 2011, ifelse(pold_spike$variable  == "spikelets_inflC5", 2011,ifelse(pold_spike$variable == "spikelets_inflA6", 2012, ifelse(pold_spike$variable == "spikelets_inflA7", 2013, ifelse(pold_spike$variable == "spikelets_InflB7", 2013, ifelse(pold_spike$variable == "spikelets_inflC7", 2013, ifelse(pold_spike$variable  == "spikelets_inflA8", 2014, ifelse(pold_spike$variable == "spikelets_InflB8", 2014, ifelse(pold_spike$variable == "spikelets_inflC8", 2014, ifelse(pold_spike$variable == "spikelets_inflA9", 2015, ifelse(pold_spike$variable == "spikelets_InflB9", 2015, ifelse(pold_spike$variable ==   "spikelets_inflC9", 2015, ifelse(pold_spike$variable == "spikelets_inflA10", 2016, ifelse(pold_spike$variable == "spikelets_InflB10", 2016, ifelse(pold_spike$variable == "spikelets_inflC10", 2016, NA))))))))))))))))))))))))

pold_spike1 <- pold_spike %>% 
  filter(!is.na(spikelets))
# View(pold_spike1)

poldflw <- POAL_data_old %>% 
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
poldflw$year<- ifelse(poldflw$variable == "flw2007", 2007, ifelse(poldflw$variable == "Flwtillers1", 2008, ifelse(poldflw$variable  == "FlwTillers3", 2009, ifelse(poldflw$variable  == "FlwTillers4", 2010, ifelse(poldflw$variable  == "FlwTillers5", 2011, ifelse(poldflw$variable  == "FlwTillers6", 2012, ifelse(poldflw$variable  == "FlwTillers7", 2013,ifelse(poldflw$variable == "FlwTillers8", 2014,ifelse(poldflw$variable == "FlwTillers9", 2015,ifelse(poldflw$variable  == "FlwTillers10", 2016, NA))))))))))
# View(poldflw)

pold_seedmerge_ss <- full_join(pold_spike1, pold_seed1, by = c( "plot", "pos", "tag", "Endo", 
                                                          "Birth","year", "tillerid"))
# View(pold_seedmerge_ss)

pold_seedmerge_ssf <- merge(pold_seedmerge_ss, poldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth","year"), all = TRUE)
# View(pold_seedmerge_ssf)



# Pulling out the seed, spikelet and flw tiller info from the POAL New Recruits
rold_seed <-POAL_data_old_r %>% 
  rename("tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate("Birth" = case_when(Date == "Q2010" ~ "2010",
                                  Date != "Q2010" ~ Date)) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", 
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B")) 
rold_seed$year<- ifelse(rold_seed$variable == "seed2007", 2007, ifelse(rold_seed$variable == "seed2008", 2008, ifelse(rold_seed$variable == "seed2008", 2008, ifelse(rold_seed$variable == "seed2009", 2009,ifelse(rold_seed$variable == "seed2010", 2010, ifelse(rold_seed$variable  == "seed2011", 2011, ifelse(rold_seed$variable  == "seed2012", 2012, ifelse(rold_seed$variable  == "seed2013", 2013, ifelse(rold_seed$variable  == "seed2014", 2014, ifelse(rold_seed$variable  == "seed2015", 2015, ifelse(rold_seed$variable  == "seed2016", 2016, NA)))))))))))

rold_seed1 <- rold_seed %>% 
  filter(!is.na(seed))
# View(rold_seed1)

rold_spike <- POAL_data_old_r %>% 
  rename("tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate("Birth" = case_when(Date == "Q2010" ~ "2010",
                                  Date != "Q2010" ~ Date))  %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike1_ = spike1, spike2_ = spike2) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", "spikelets_inflA5", "spikelets_infl13","spike1_","spike2_", "spikelets_infl15", 
                       "spikelets_infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("spike1_",variable)~"A",
                              grepl("spike2_", variable)~"B"))
rold_spike$year<- ifelse(rold_spike$variable == "spike2007", 2007, ifelse(rold_spike$variable == "spike2008", 2008, ifelse(rold_spike$variable == "spike2009", 2009, ifelse(rold_spike$variable == "spike2010", 2010, ifelse(rold_spike$variable == "spike2011", 2011, ifelse(rold_spike$variable == "spikelets_inflA5", 2012, ifelse(rold_spike$variable  == "spikelets_infl13", 2013, ifelse(rold_spike$variable  == "spike1_", 2014, ifelse(rold_spike$variable == "spike2_", 2014, ifelse(rold_spike$variable  == "spikelets_infl15", 2015, ifelse(rold_spike$variable  == "spikelets_infl16", 2016, NA)))))))))))
rold_spike1 <- rold_spike %>% 
  filter(!is.na(spikelets))
# View(rold_spike1)

roldflw <- POAL_data_old_r %>%
  rename("tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate("Birth" = case_when(Date == "Q2010" ~ "2010",
                                  Date != "Q2010" ~ Date)) %>%   
  mutate(flw2007 = NA, flw2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "FLWTiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
roldflw$year<- ifelse(roldflw$variable == "flw2007", 2007, ifelse(roldflw$variable == "flw2008", 2008, ifelse(roldflw$variable == "FLWTiller09", 2009, ifelse(roldflw$variable == "FLWtiller10", 2010, ifelse(roldflw$variable == "FLWtiller11", 2011, ifelse(roldflw$variable  == "FLWtiller12", 2012, ifelse(roldflw$variable  == "FLWtiller13", 2013, ifelse(roldflw$variable  == "FLWtiller14", 2014, ifelse(roldflw$variable  == "FLWtiller15", 2015, ifelse(roldflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(roldflw)

rold_seedmerge_ss <- full_join(rold_spike1, rold_seed1, by = c( "plot", "pos", "tag", "Endo", 
                                                          "Birth","year", "tillerid"))
# View(rold_seedmerge_ss)

rold_seedmerge_ssf <- merge(rold_seedmerge_ss, roldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth","year"), all = TRUE)
# View(rold_seedmerge_ssf)
# POAL(Old) recruits data includes tag # 10_19, which doesn't have any data collected for the reproductive data and is not present in the endo_demog_long file. This will be filtered out later when we merge the reproductive output with endo_demog_long 

# View(rold_seedmerge_ssf)

pseedmerge_ssf <- pseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                  "tillerid", "flw","spikelets", "seed")]
pold_seedmerge_ssf <- pold_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                           "tillerid", "flw","spikelets", "seed")]
rseedmerge_ssf <- rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                   "tillerid", "flw","spikelets", "seed")]
rold_seedmerge_ssf <- rold_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                           "tillerid", "flw","spikelets", "seed")]


str(pseedmerge_ssf)
str(rseedmerge_ssf)
str(pold_seedmerge_ssf)
str(rold_seedmerge_ssf)

POALrepro <- pseedmerge_ssf %>% 
  rbind(rseedmerge_ssf) %>% 
  rbind(pold_seedmerge_ssf) %>% 
  rbind(rold_seedmerge_ssf) %>% 
  mutate(species = "POAL")
# POALrepro <- POALrepro[!(is.na(POALrepro$flw)),]





# Combining seed productions measurements across years for the “New” POSY data -------------

## recoding for the year of measurement
## merging these measurements into one dataframe
po_seed <- POSY_data %>% 
  mutate("Birth" = year(as.character(Planted.Date))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA2", "seeds_InflB2", 
                       "seed2010", "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
         
po_seed$year<- ifelse(po_seed$variable == "seed2007", 2007, ifelse(po_seed$variable == "seed2008", 2008,ifelse(po_seed$variable == "seeds_InflA2", 2009, ifelse(po_seed$variable  == "seeds_InflB2", 2009, ifelse(po_seed$variable  == "seed2010", 2010, ifelse(po_seed$variable  == "seed2011", 2011, ifelse(po_seed$variable  == "seed2012", 2012, ifelse(po_seed$variable  == "seed2013", 2013,ifelse(po_seed$variable == "seed2014", 2014,ifelse(po_seed$variable == "seed2015", 2015,ifelse(po_seed$variable  == "seed2016", 2016, NA)))))))))))
po_seed1 <- po_seed %>% 
  filter(!is.na(seed))
# View(po_seed1)

po_spike <- POSY_data %>% 
  mutate("Birth" = year(as.character(Planted.Date))) %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spikelets_InflA2", "spikelets_InflB2", "spikelets_inflA3", 
                       "spikelets_inflB3", "spikelets_inflA4","spikelets_inflB4",
                       "spikelets_inflA5", "spikelets_inflA6", "spikelets_InflB6", 
                       "spikelets_inflC6", "spikelets_inflA7", "spikelets_InflB7", 
                       "spikelets_inflC7","spikelets_inflA8", "spikelets_InflB8", 
                       "spikelets_inflC8", "spike2016"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B", variable) ~ "B",
                              grepl("C", variable) ~ "C"))

po_spike$year<- ifelse(po_spike$variable == "spike2007", 2007, ifelse(po_spike$variable == "spike2008", 2008,ifelse(po_spike$variable == "spikelets_InflA2", 2009, ifelse(po_spike$variable  == "spikelets_InflB2", 2009, ifelse(po_spike$variable  == "spikelets_inflA3", 2010, ifelse(po_spike$variable  == "spikelets_inflB3", 2010, ifelse(po_spike$variable  == "spikelets_inflB3", 2010, ifelse(po_spike$variable  == "spikelets_inflA4", 2011, ifelse(po_spike$variable == "spikelets_inflB4", 2011, ifelse(po_spike$variable == "spikelets_inflA5", 2012,ifelse(po_spike$variable == "spikelets_inflA6", 2013, ifelse(po_spike$variable == "spikelets_InflB6", 2013, ifelse(po_spike$variable == "spikelets_inflC6", 2013, ifelse(po_spike$variable  == "spikelets_inflA7", 2014, ifelse(po_spike$variable == "spikelets_InflB7", 2014, ifelse(po_spike$variable == "spikelets_inflC7", 2014, ifelse(po_spike$variable == "spikelets_inflA8", 2015, ifelse(po_spike$variable == "spikelets_InflB8", 2015, ifelse(po_spike$variable ==   "spikelets_inflC8", 2015, ifelse(po_spike$variable == "spike2016", 2016, NA))))))))))))))))))))
po_spike1 <- po_spike %>% 
  filter(!is.na(spikelets))
# View(po_spike1)


po_flw <- POSY_data %>% 
  mutate("Birth" = year(as.character(Planted.Date))) %>% 
  mutate(flw2007 = 0, flw2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "flw2016"), 
       value.name = "flw") %>% 
  filter(!is.na(plot))
po_flw$year<- ifelse(po_flw$variable == "flw2007", 2007, ifelse(po_flw$variable == "Flwtillers1", 2008, ifelse(po_flw$variable  == "FlwTillers2", 2009, ifelse(po_flw$variable  == "FlwTillers3", 2010, ifelse(po_flw$variable  == "FlwTillers4", 2011, ifelse(po_flw$variable  == "FlwTillers5", 2012, ifelse(po_flw$variable  == "FlwTillers6", 2013,ifelse(po_flw$variable == "FlwTillers7", 2014,ifelse(po_flw$variable == "FlwTillers8", 2015,ifelse(po_flw$variable  == "flw2016", 2016, NA))))))))))
# View(po_flw)

po_seedmerge_ss <- full_join(po_spike1, po_seed1, by = c("plot", "pos", "tag", "Endo", 
                                                   "Birth","year","tillerid"))
# View(po_seedmerge_ss)


po_seedmerge_ssf <- merge(po_seedmerge_ss,po_flw, by = c("plot", "pos", "tag", "Endo", 
                                                          "Birth","year"), all = TRUE)

# View(po_seedmerge_ssf)




# Combining repro measurements across years for the "New" POSY recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
po_rseed <- POSY_data_r %>% 
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2009", "seed2010", "seed2011",
                       "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
po_rseed$year<- ifelse(po_rseed$variable == "seed2007", 2007, ifelse(po_rseed$variable == "seed2008", 2008, ifelse(po_rseed$variable == "seed2009", 2009, ifelse(po_rseed$variable  == "seed2010", 2010, ifelse(po_rseed$variable  == "seed2011", 2011, ifelse(po_rseed$variable  == "seed2012", 2012, ifelse(po_rseed$variable  == "seed2013", 2013,ifelse(po_rseed$variable == "seed2014", 2014,ifelse(po_rseed$variable == "seed2015", 2015,ifelse(po_rseed$variable  == "seed2016", 2016, NA))))))))))

po_rseed1 <- po_rseed %>% 
  filter(!is.na(seed))
# View(po_rseed1)

po_rspike <- POSY_data_r %>% 
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA,) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", "spike2011", 
                       "spikelets_inflA5...26", "spikelets_inflA5...32", 
                       "spikelets_infl14","spike1...45", "spike2...46", 
                       "spike3...47", "spike1...54", "spike2...55", "spike3...56"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                              grepl("spike2", variable) ~ "B",
                              grepl("spike3", variable) ~ "C"))
po_rspike$year<- ifelse(po_rspike$variable  == "spike2007", 2007, ifelse(po_rspike$variable  == "spike2008", 2008, ifelse(po_rspike$variable  == "spike2009", 2009, ifelse(po_rspike$variable  == "spike2010", 2010, ifelse(po_rspike$variable == "spike2011", 2011, ifelse(po_rspike$variable  == "spikelets_inflA5...26", 2012, ifelse(po_rspike$variable == "spikelets_inflA5...32", 2013, ifelse(po_rspike$variable == "spikelets_InflB6", 2013, ifelse(po_rspike$variable == "spikelets_inflC6", 2013, ifelse(po_rspike$variable  == "spikelets_infl14", 2014, ifelse(po_rspike$variable == "spike1...45", 2015, ifelse(po_rspike$variable == "spike2...46", 2015, ifelse(po_rspike$variable ==   "spike3...47", 2015, ifelse(po_rspike$variable == "spike1...54", 2016, ifelse(po_rspike$variable == "spike2...55", 2016, ifelse(po_rspike$variable == "spike3...56", 2016, NA))))))))))))))))

po_rspike1 <- po_rspike %>% 
  filter(!is.na(spikelets))# View(po_rspike)


po_rflw <- POSY_data_r %>%
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = 0, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
po_rflw$year<- ifelse(po_rflw$variable == "flw2007", 2007, ifelse(po_rflw$variable == "flw2008", 2008, ifelse(po_rflw$variable == "flw2009", 2009, ifelse(po_rflw$variable == "FLWtiller10", 2010, ifelse(po_rflw$variable == "FLWtiller11", 2011, ifelse(po_rflw$variable  == "FLWtiller12", 2012, ifelse(po_rflw$variable  == "FLWtiller13", 2013, ifelse(po_rflw$variable  == "FLWtiller14", 2014, ifelse(po_rflw$variable  == "FLWtiller15", 2015, ifelse(po_rflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(po_rflw)


po_rseedmerge_ss <- full_join(po_rspike1, po_rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                           "Birth","year","tillerid"))
# View(po_rmerge_ss)

po_rseedmerge_ssf <- merge(po_rseedmerge_ss, po_rflw, by = c("plot", "pos", "tag", "Endo", 
                                                         "Birth","year"), all = TRUE)
# View(po_rseedmerge_ssf)






# Combining repro measurements across years for the “Old” POSY data ------------------


## Combining data across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe

po_oldseed <- POSY_data_old %>%
  mutate(seed2007 = NA, seed2008 = NA,  seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                 is.na(Date) ~ 2007)) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", 
                  "Birth", "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA3", "seeds_InflB3", "seed2010", 
                       "seed2011", "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

po_oldseed$year<- ifelse(po_oldseed$variable == "seed2007", 2007, ifelse(po_oldseed$variable == "seed2008", 2008, ifelse(po_oldseed$variable  == "seeds_InflA3", 2009, ifelse(po_oldseed$variable  == "seeds_InflB3",2009, ifelse(po_oldseed$variable == "seed2010", 2010, ifelse(po_oldseed$variable  == "seed2011", 2011, ifelse(po_oldseed$variable  == "seed2012", 2012, ifelse(po_oldseed$variable  == "seed2013", 2013,ifelse(po_oldseed$variable == "seed2014", 2014,ifelse(po_oldseed$variable == "seed2015", 2015,ifelse(po_oldseed$variable  == "seed2016", 2016, NA)))))))))))
po_oldseed1 <- po_oldseed %>% 
  filter(!is.na(seed))
# View(po_oldseed1)

po_oldspike <- POSY_data_old %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(spike2007 = NA, ) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", 
                  "Birth", "TRT", "Plant"), 
       measure.var = c("spike2007", "Spikelets_InflA1", "Spikelets_InflB1", 
                       "spikelets_InflA3","spikelets_InflB3", 
                       "spikelets_InflA4","spikelets_inflB4", "spikelets_InflA5", 
                       "spikelets_inflB5", "spikelets_inflA6", "spikelets_inflA7",
                       "spikelets_InflB7","spikelets_inflC7","spikelets_inflA8",
                       "spikelets_InflB8","spikelets_inflC8","spikelets_inflA9",
                       "spikelets_InflB9","spikelets_inflC9","spikelets_inflA10",
                       "spikelets_InflB10","spikelets_inflC10"), 
       value.name = "spikelets") %>% 
    mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                                grepl("B",variable) ~ "B",
                                grepl("C", variable) ~ "C"))

po_oldspike$year<- ifelse(po_oldspike$variable == "spike2007", 2007, ifelse(po_oldspike$variable == "Spikelets_InflA1", 2008,  ifelse(po_oldspike$variable == "Spikelets_InflB1", 2008, ifelse(po_oldspike$variable  == "spikelets_InflA3", 2009, ifelse(po_oldspike$variable  == "spikelets_InflB3", 2009, ifelse(po_oldspike$variable  == "spikelets_InflA4", 2010, ifelse(po_oldspike$variable  == "spikelets_inflB4", 2010, ifelse(po_oldspike$variable  == "spikelets_InflA5", 2011, ifelse(po_oldspike$variable  == "spikelets_inflB5", 2011, ifelse(po_oldspike$variable  == "spikelets_inflA6", 2012, ifelse(po_oldspike$variable  == "spikelets_inflA7", 2013, ifelse(po_oldspike$variable  == "spikelets_InflB7", 2013, ifelse(po_oldspike$variable  == "spikelets_inflC7", 2013, ifelse(po_oldspike$variable == "spikelets_inflA8", 2014, ifelse(po_oldspike$variable == "spikelets_InflB8", 2014, ifelse(po_oldspike$variable == "spikelets_inflC8", 2014, ifelse(po_oldspike$variable == "spikelets_inflA9", 2015,ifelse(po_oldspike$variable == "spikelets_InflB9", 2015, ifelse(po_oldspike$variable == "spikelets_inflC9", 2015, ifelse(po_oldspike$variable  == "spikelets_inflA10", 2016, ifelse(po_oldspike$variable  == "spikelets_InflB10", 2016, ifelse(po_oldspike$variable  == "spikelets_inflC10", 2016,NA))))))))))))))))))))))
po_oldspike1 <- po_oldspike %>% 
  filter(!is.na(spikelets), spikelets != 0)
# View(po_oldspike1)

po_oldflw <- POSY_data_old %>% 
  rename("plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n",
                  "Birth", "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
po_oldflw$year<- ifelse(po_oldflw$variable == "flw2007", 2007, ifelse(po_oldflw$variable == "FlwTillers1", 2008, ifelse(po_oldflw$variable  == "FlwTillers3", 2009, ifelse(po_oldflw$variable  == "FlwTillers4", 2010, ifelse(po_oldflw$variable  == "FlwTillers5", 2011, ifelse(po_oldflw$variable  == "FlwTillers6", 2012, ifelse(po_oldflw$variable  == "FlwTillers7", 2013,ifelse(po_oldflw$variable == "FlwTillers8", 2014,ifelse(po_oldflw$variable == "FlwTillers9", 2015,ifelse(po_oldflw$variable  == "FlwTillers10", 2016, NA))))))))))
# View(po_oldflw)


po_oldseedmerge_ss <- full_join(po_oldseed1, po_oldspike1, by = c("plot", "pos", "tag", "Endo", 
                                                                "Birth","year","tillerid"))
# View(po_oldmerge_ss)

po_oldseedmerge_ssf <- merge(po_oldseedmerge_ss, po_oldflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth","year"), all = TRUE)
# View(po_oldseedmerge_ssf)



# Combining  measurements across years for the “Old” POSY recruits data --------

## Combining data for recruits across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
po_roldseed <- POSY_data_old_r %>%
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("seed1", variable) ~ "A", 
                              grepl("seed2",variable) ~ "B",
                              grepl("seed3", variable) ~ "C"))
po_roldseed$year<- ifelse(po_roldseed$variable == "seed2007", 2007, ifelse(po_roldseed$variable == "seed2008", 2008, ifelse(po_roldseed$variable == "seed2009", 2009, ifelse(po_roldseed$variable == "seed2010", 2010, ifelse(po_roldseed$variable == "seed2011", 2011, ifelse(po_roldseed$variable  == "seed2012", 2012, ifelse(po_roldseed$variable  == "seed2013", 2013, ifelse(po_roldseed$variable  == "seed2014", 2014, ifelse(po_roldseed$variable  == "seed2015", 2015, ifelse(po_roldseed$variable  == "seed2016", 2016, NA))))))))))
po_roldseed1 <- po_roldseed %>% 
  filter(!is.na(seed))
# View(po_roldseed1)


po_roldspike <- POSY_data_old_r %>%
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 =NA, spike2010 = NA, spike2011 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010", 
                       "spike2011", "spikelets_inflA12","spikelets_inflA13",
                       "spikelets_inflb13","spikelets_inflc13",
                       "spike1...39", "spike2...40", "spike3...41",
                       "spike1...48", "spike2...49", "spike3...50",
                       "spike1...57", "spike2...58", "spike3...59"),
       value.name = "spikelets") %>% 
    mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                                grepl("spike2",variable) ~ "B",
                                grepl("spike3", variable) ~ "C",
                                grepl("inflA", variable) ~ "A",
                                grepl("inflb", variable) ~ "B",
                                grepl("inflc", variable) ~ "C"))
  
po_roldspike$year<- ifelse(po_roldspike$variable == "spike2007", 2007, ifelse(po_roldspike$variable == "spike2008", 2008, ifelse(po_roldspike$variable == "spike2009", 2009, ifelse(po_roldspike$variable == "spike2010", 2010, ifelse(po_roldspike$variable == "spike2011", 2011, ifelse(po_roldspike$variable  == "spikelets_inflA12", 2012, ifelse(po_roldspike$variable  == "spikelets_inflA13", 2013, ifelse(po_roldspike$variable  == "spikelets_inflb13", 2013, ifelse(po_roldspike$variable  == "spikelets_inflc13", 2013, ifelse(po_roldspike$variable  == "spike1...39", 2014, ifelse(po_roldspike$variable  == "spike2...40", 2014, ifelse(po_roldspike$variable  == "spike3...41", 2014, ifelse(po_roldspike$variable  == "spike1...48", 2015, ifelse(po_roldspike$variable  == "spike2...49", 2015, ifelse(po_roldspike$variable  == "spike3...50", 2015, ifelse(po_roldspike$variable  == "spike1...57", 2016, ifelse(po_roldspike$variable  == "spike2...58", 2016,ifelse(po_roldspike$variable  == "spike3...59", 2016,NA))))))))))))))))))
po_roldspike1 <- po_roldspike %>% 
  filter(!is.na(spikelets))
#View(po_roldspike1)

po_roldflw <- POSY_data_old_r %>%
  rename("Birth" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLWtiller10","FLWtiller11", "FLWTiller12", 
                       "FLWTiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") %>% 
  distinct() # there are a couple of tag ids that are duplicated
po_roldflw$year<- ifelse(po_roldflw$variable == "flw2007", 2007, ifelse(po_roldflw$variable == "flw2008", 2008, ifelse(po_roldflw$variable == "flw2009", 2009, ifelse(po_roldflw$variable == "FLWtiller10", 2010, ifelse(po_roldflw$variable == "FLWtiller11", 2011, ifelse(po_roldflw$variable  == "FLWTiller12", 2012, ifelse(po_roldflw$variable  == "FLWTiller13", 2013, ifelse(po_roldflw$variable  == "FLWtiller14", 2014, ifelse(po_roldflw$variable  == "FLWtiller15", 2015, ifelse(po_roldflw$variable  == "FLWtiller16", 2016, NA))))))))))
# View(po_roldflw)

po_roldseedmerge_ss <-  full_join(po_roldspike1, po_roldseed1, by = c("plot", "pos", "tag", "Endo", 
                                                                   "Birth","year","tillerid"))
# View(po_roldseedmerge_ss)

po_roldseedmerge_ssf <- merge(po_roldseedmerge_ss, po_roldflw, by = c("plot", "pos", "tag", "Endo", 
                                                                     "Birth","year"), all = TRUE)
# View(po_roldseedmerge_ssf)






# Combining the old and new and original and recruit POSY repro dataframes ---------
po_seedmerge_ssf <- po_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                   "tillerid", "flw","spikelets", "seed")]
po_oldseedmerge_ssf <- po_oldseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                   "tillerid", "flw","spikelets", "seed")]
po_rseedmerge_ssf <- po_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                   "tillerid", "flw","spikelets", "seed")]
po_roldseedmerge_ssf <- po_roldseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                           "tillerid", "flw","spikelets", "seed")]


str(po_seedmerge_ssf)
str(po_rseedmerge_ssf)
str(po_oldseedmerge_ssf)
str(po_roldseedmerge_ssf)

POSYrepro <- po_seedmerge_ssf %>% 
  rbind(po_rseedmerge_ssf) %>% 
  rbind(po_oldseedmerge_ssf) %>% 
  rbind(po_roldseedmerge_ssf) %>% 
  mutate(species = "POSY")

# POSYrepro <- POSYrepro[!(is.na(POSYrepro$flw)),]
# View(POSYrepro)





# Combining repro measurements across years for the LOAR data ------------------

## Combining measurements across years for reproduction measurements

# there is a separate sheet with the main raw data for seeds from 2008. It is pretty messy, but I am removing the columns that contain duplicated spikelet information as well as a few rows that have questionable data
LOAR_data_seed2008_long <- LOAR_data_seed2008 %>% 
  dplyr::select(-contains("__1"), -contains("Avg"), -PropClaviceps1, -contains("Bug"), -tag2, -Notes,-contains("SeedNotes1"), -contains("totunfilled")) %>% 
  rename(Tag = Tag...2) %>% 
  filter(!is.na(Tag), Spikelets...3 != "-") %>% 
  group_by(Tag) %>% mutate(tillerid2008 = as.character(LETTERS[row_number()]),  year = 2008) %>% 
  mutate(seed2008 = as.numeric(Unfilled.Green...7) + as.numeric(Unfilled.Brown...6) + Filled...4) %>% 
  rename(spikelets2008 = Spikelets...3)

LOAR_data_seed2008_seed <- LOAR_data_seed2008_long %>% 
  dplyr::select(Tag, tillerid2008, year, seed2008)
LOAR_data_seed2008_spike <- LOAR_data_seed2008_long %>% 
  dplyr::select(Tag, tillerid2008, year, spikelets2008)
# there is a separate sheet with the main raw data for seeds from 2009. It is pretty messy, but I am removing the columns that contain duplicated spikelet information as well as a few rows that have questionable data
LOAR_data_seed2009_long <-LOAR_data_seed2009 %>% 
  dplyr::select(-contains("__1"), -contains("Avg"), -contains("__2"), -contains("Bug"), -Notes, -contains("TOTunfilled"), -TOTUnfilled) %>% 
  rename(Tag = Tag...2) %>% 
  filter(!is.na(Tag)) %>% 
  group_by(Tag) %>% mutate(tillerid2009 = as.character(LETTERS[row_number()]), year = 2009) %>% 
  mutate(seed2009 = Unfilled.Green...7 + Unfilled.Brown...6 + Filled...4) %>% 
  rename(spikelets2009 = Spikelets...3)

LOAR_data_seed2009_seed <- LOAR_data_seed2009_long %>% 
  dplyr::select(Tag, tillerid2009, year, seed2009)
LOAR_data_seed2009_spike <- LOAR_data_seed2009_long %>% 
  dplyr::select(Tag, tillerid2009, year, spikelets2009)

# I'm going to merge these two years into the spikelet and seed data respectively after producing the "long" version of the rest of the years.
# There is really no raw seed data in the LOAR excel sheet, just calculations from the spikelet counts, and some averages (some of which I think are taking the wrong cells, or omiting the *.12 rate of seeds per spikelet)
lseed <- LOAR_data%>% 
  rename("plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG" ) %>% 
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(seed2007 = NA,
         seed2008 = NA,
         seed2009 = NA,
         seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009",
                       "seed2010", "seed2011", "seed2012", "seed2013","seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
    mutate(tillerid = case_when(grepl("a", variable) ~ "A", 
                                grepl("b",variable) ~ "B",
                                grepl("c", variable) ~ "C"))
lseed$year <- ifelse(lseed$variable == "seed2007", 2007, ifelse(lseed$variable == "seed2008", 2008, ifelse(lseed$variable  == "seed2009", 2009, ifelse(lseed$variable  == "seed2010", 2010, ifelse(lseed$variable  == "seed2011", 2011, ifelse(lseed$variable  == "seed2012", 2012, ifelse(lseed$variable  == "seed2013", 2013,ifelse(lseed$variable == "seed2014", 2014,ifelse(lseed$variable == "seed2015", 2015,ifelse(lseed$variable  == "seed2016", 2016, NA))))))))))
# View(lseed)

# merging in the 2008 and 2009 seed data
lseed_1 <- full_join(lseed, LOAR_data_seed2008_seed, by = c("tag" = "Tag", "year" = "year"))


lseed_2 <- full_join(lseed_1, LOAR_data_seed2009_seed, by = c("tag" = "Tag", "year" = "year"))

lseed_2$seed_new <- ifelse(!is.na(lseed_2$seed), lseed_2$seed, ifelse(!is.na(lseed_2$seed2008), lseed_2$seed2008, ifelse(!is.na(lseed_2$seed2009), lseed_2$seed2009, NA)))
lseed_2$tillerid_new <- ifelse(!is.na(lseed_2$tillerid), lseed_2$tillerid, ifelse(!is.na(lseed_2$tillerid2008), lseed_2$tillerid2008, ifelse(!is.na(lseed_2$tillerid2009), lseed_2$tillerid2009, NA)))
lseed_3 <- lseed_2 %>% 
  dplyr::select(plot, pos, tag, Endo, Birth, variable, year, seed_new, tillerid_new) %>% 
  rename(seed = seed_new, tillerid = tillerid_new) %>% 
  filter(!is.na(seed))


# 2011 doesn't seem to have clear spikelet/infl data. In the LOAR_for_demog_long.xlsx file, there is a column for good_seeds/inf, and a column for number of infl counted, so presumably it is a mean. There are some spikelet counts in the LOAR_to2016_complete.xlsx for what I believe is 2011. One example: plant 601, has listed 42 spike/infl in the later DB, but has 0 good seeds/infl listed in the former.
lspike <- LOAR_data %>% 
  rename("plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(spike2007 = NA,
         spike2008 = NA,
         spike2009 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"), 
       measure.var = c("spike2007", "spike2008", "spike2009",
                       "spikelets_inflA3", "spikelets_inflB3",
                       "spikelets_inflA5", "spikelets_inflB5",
                       "spikelets_inflA6", "spikelets_inflB6",
                       "spikelets_inflA7", "spikelets_inflB7",
                       "spikelets_inflA8", "spikelets_inflB8",
                       "spikelets_inflA9", "spikelets_inflB9"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
lspike$year<- ifelse(lspike$variable == "spike2007", 2007, ifelse(lspike$variable == "spike2008", 2008,  ifelse(lspike$variable == "spike2009", 2009, ifelse(lspike$variable  == "spikelets_inflA3", 2010, ifelse(lspike$variable  == "spikelets_inflB3", 2010, ifelse(lspike$variable  == "spikelets_inflA5", 2012, ifelse(lspike$variable  == "spikelets_inflB5", 2012, ifelse(lspike$variable  == "spikelets_inflA6", 2013, ifelse(lspike$variable  == "spikelets_inflB6", 2013, ifelse(lspike$variable  == "spikelets_inflA7", 2014, ifelse(lspike$variable  == "spikelets_inflB7", 2014, ifelse(lspike$variable == "spikelets_inflA8", 2015, ifelse(lspike$variable == "spikelets_inflB8", 2015, ifelse(lspike$variable == "spikelets_inflA9", 2016, ifelse(lspike$variable == "spikelets_inflB9", 2016, NA)))))))))))))))

# View(lspike)

lspike_1 <- full_join(lspike, LOAR_data_seed2008_spike, by = c("tag" = "Tag", "year" = "year"))
lspike_2 <- full_join(lspike_1, LOAR_data_seed2009_spike, by = c("tag" = "Tag", "year" = "year"))

lspike_2$spikelets_new <- ifelse(!is.na(lspike_2$spikelets), lspike_2$spikelets, ifelse(!is.na(lspike_2$spikelets2008), lspike_2$spikelets2008, ifelse(!is.na(lspike_2$spikelets2009), lspike_2$spikelets2009, NA)))
lspike_2$tillerid_new <- ifelse(!is.na(lspike_2$tillerid), lspike_2$tillerid, ifelse(!is.na(lspike_2$tillerid2008), lspike_2$tillerid2008, ifelse(!is.na(lspike_2$tillerid2009), lspike_2$tillerid2009, NA)))
lspike_3 <- lspike_2 %>% 
  dplyr::select(plot, pos, tag, Endo, Birth, variable, year, spikelets_new, tillerid_new) %>% 
  rename(spikelets = spikelets_new, tillerid = tillerid_new) %>% 
  filter(!is.na(spikelets))


lflw <- LOAR_data %>% 
  rename("plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate("Birth" = case_when(!is.na(Date) ~ year(as.character(Date)),
                                  is.na(Date) ~ 2007)) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Birth"), 
       measure.var = c("flw2007", "FlwTillers1", "FLWTiller2", "FLWTillers3", "FLWTillers4",
                       "FLWTillers5", "FLWTillers6", "FLWTillers7", 
                       "FLWTillers8", "FLWTillers9"), 
       value.name = "flw") 
lflw$year<- ifelse(lflw$variable == "flw2007", 2007, ifelse(lflw$variable == "FlwTillers1", 2008, ifelse(lflw$variable  == "FLWTiller2", 2009, ifelse(lflw$variable  == "FLWTillers3", 2010, ifelse(lflw$variable  == "FLWTillers4", 2011, ifelse(lflw$variable  == "FLWTillers5", 2012, ifelse(lflw$variable  == "FLWTillers6", 2013,ifelse(lflw$variable == "FLWTillers7", 2014,ifelse(lflw$variable == "FLWTillers8", 2015,ifelse(lflw$variable  == "FLWTillers9", 2016, NA))))))))))
# View(lflw)

l_seedmerge_ss <- full_join(lspike_3, lseed_3, by = c("plot", "pos", "tag", "Endo", 
                                                      "Birth", "year", "tillerid"))
# View(l_seedmerge_ss)


l_seedmerge_ssf <- merge(l_seedmerge_ss, lflw, by = c("plot", "pos", "tag", "Endo", 
                                                      "Birth", "year"), all = TRUE)

# View(l_seedmerge_ssf)





# Combining repro measurements across years for the LOAR recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
l_rseed <- LOAR_data_r %>%
  rename("Birth" = "Date", "plot" = "Plot", "pos" = "Recruit.ID", "Endo" = "endo") %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"), 
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed")   %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))
l_rseed$year<- ifelse(l_rseed$variable == "seed2007", 2007, ifelse(l_rseed$variable == "seed2008", 2008, ifelse(l_rseed$variable == "seed2009", 2009,ifelse(l_rseed$variable == "seed2010", 2010, ifelse(l_rseed$variable == "seed2011", 2011, ifelse(l_rseed$variable  == "seed2012", 2012, ifelse(l_rseed$variable  == "seed2013", 2013, ifelse(l_rseed$variable  == "seed2014", 2014, ifelse(l_rseed$variable  == "seed2015", 2015, ifelse(l_rseed$variable  == "seed2016", 2016, NA))))))))))
l_rseed_1 <- l_rseed %>% 
  filter(!is.na(seed))
# View(l_rseed_1)


l_rspike <- LOAR_data_r %>%
  rename("Birth" = "Date", "plot" = "Plot", "pos" = "Recruit.ID", "Endo" = "endo") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", 
                       "spike2010","spike2011", "spike2012",
                       "spike1...30", "spike2...31", "spike1...37", "spike2...38",
                       "spike1_15", "spike2_15", "spike1_16", "spike2_16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("spike1", variable) ~ "A", 
                              grepl("spike2",variable) ~ "B"))
l_rspike$year<- ifelse(l_rspike$variable == "spike2007", 2007, ifelse(l_rspike$variable == "spike2008", 2008, ifelse(l_rspike$variable == "spike2009", 2009, ifelse(l_rspike$variable == "spike2010", 2010, ifelse(l_rspike$variable == "spike2011", 2011, ifelse(l_rspike$variable  == "spike2012", 2012, ifelse(l_rspike$variable  == "spike1...30", 2013, ifelse(l_rspike$variable  == "spike2...31", 2013, ifelse(l_rspike$variable  == "spike1...37", 2014, ifelse(l_rspike$variable  == "spike2...38", 2014, ifelse(l_rspike$variable  == "spike1_15", 2015,ifelse(l_rspike$variable  == "spike2_15", 2015, ifelse(l_rspike$variable  == "spike1_16", 2016,ifelse(l_rspike$variable  == "spike2_16", 2016, NA))))))))))))))
l_rspike_1 <- l_rspike %>% 
  filter(!is.na(spikelets))
# View(l_rspike_1)

l_rflw <- LOAR_data_r %>%
  rename("Birth" = "Date", "plot" = "Plot", "pos" = "Recruit.ID", "Endo" = "endo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FlwTillers10","Flw11", "Flw12", 
                       "Flw13", "Flw14", "Flw15",
                       "Flw16"),
       value.name = "flw") 
l_rflw$year<- ifelse(l_rflw$variable == "flw2007", 2007, ifelse(l_rflw$variable == "flw2008", 2008, ifelse(l_rflw$variable == "flw2009", 2009, ifelse(l_rflw$variable == "FlwTillers10", 2010, ifelse(l_rflw$variable == "Flw11", 2011, ifelse(l_rflw$variable  == "Flw12", 2012, ifelse(l_rflw$variable  == "Flw13", 2013, ifelse(l_rflw$variable  == "Flw14", 2014, ifelse(l_rflw$variable  == "Flw15", 2015, ifelse(l_rflw$variable  == "Flw16", 2016, NA))))))))))
# View(l_rflw)

l_rseedmerge_ss <- full_join(l_rspike_1, l_rseed_1, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Birth", "year", "tillerid"))
# View(l_rseedmerge_ss)

l_rseedmerge_ssf <- merge(l_rseedmerge_ss, l_rflw, by = c( "plot", "pos", "tag", "Endo", 
                                                    "Birth", "year"),all = TRUE)
# View(l_rseedmerge_ssf)

# Combining recruit and original plant data
l_seedmerge_ssf <- l_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                       "tillerid", "flw","spikelets", "seed")]
l_rseedmerge_ssf <- l_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                         "tillerid", "flw","spikelets", "seed")]


LOARrepro <- l_seedmerge_ssf %>% 
  rbind(l_rseedmerge_ssf) %>% 
  mutate(species = "LOAR")

# LOARrepro <- LOARrepro[!(is.na(LOARrepro$flw)),]
# View(LOARrepro)







# Combining repro measurements across years for the FESU data ------------------

## Combining measurements across years for Seed, Spike, and Flowering using melt
## Recoding those measurements for the year they are taken

fseed <- FESU_data %>%
  rename("Birth" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(seed2007 = NA, seed2008 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", "seeds_InflA2", "seeds_InflB2",  
                       "seed2010", "seed2011","seed2012",
                       "seed2013", "seed2014", "seed2015", 
                       "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

fseed$year<- ifelse(fseed$variable == "seed2007", 2007, ifelse(fseed$variable == "seed2008", 2008, ifelse(fseed$variable  == "seeds_InflA2", 2009, ifelse(fseed$variable  == "seeds_InflB2", 2009, ifelse(fseed$variable  == "seed2010", 2010, ifelse(fseed$variable  == "seed2011", 2011, ifelse(fseed$variable  == "seed2012", 2012, ifelse(fseed$variable  == "seed2013", 2013,ifelse(fseed$variable == "seed2014", 2014,ifelse(fseed$variable == "seed2015", 2015,ifelse(fseed$variable  == "seed2016", 2016, NA)))))))))))
fseed1 <- fseed %>% 
  filter(!is.na(seed))
# View(fseed1)

fspike <- FESU_data %>% 
  rename("Birth" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(spike2007 = NA, spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "spike2008", "spikelets_inflA2", "spikelets_InflB2",
                       "spikelets_inflA3", "spikelets_inflB3",
                       "spikelets_inflA4", "spikelets_inflB4",
                       "spikelets_inflA5", "spikelets_inflB5",
                       "spikelets_inflA6", "spikelets_inflB6",
                       "spikelets_inflA7", "spikelets_inflB7",
                       "spikelets_inflA8", "spikelets_inflB8",
                       "spikelets_inflA9", "spikelets_inflB9"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

fspike$year<- ifelse(fspike$variable == "spike2007", 2007, ifelse(fspike$variable == "spike2008", 2008, ifelse(fspike$variable  == "spikelets_inflA2", 2009, ifelse(fspike$variable  == "spikelets_InflB2", 2009, ifelse(fspike$variable  == "spikelets_inflA3", 2010, ifelse(fspike$variable  == "spikelets_inflB3", 2010, ifelse(fspike$variable  == "spikelets_inflA4", 2011, ifelse(fspike$variable  == "spikelets_inflB4", 2011, ifelse(fspike$variable  == "spikelets_inflA5", 2012, ifelse(fspike$variable  == "spikelets_inflB5", 2012, ifelse(fspike$variable  == "spikelets_inflA6", 2013, ifelse(fspike$variable  == "spikelets_inflB6", 2013, ifelse(fspike$variable == "spikelets_inflA7", 2014, ifelse(fspike$variable == "spikelets_inflB7", 2014, ifelse(fspike$variable == "spikelets_inflA8", 2015, ifelse(fspike$variable == "spikelets_inflB8", 2015, ifelse(fspike$variable  == "spikelets_inflA9", 2016, ifelse(fspike$variable  == "spikelets_inflB9", 2016, NA))))))))))))))))))
fspike1 <- fspike %>% 
  filter(!is.na(spikelets))
# View(fspike1)


fflw <- FESU_data %>% 
  rename("Birth" = "Planteddate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(flw2007 = 0, flw2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "flw2008", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
fflw$year<- ifelse(fflw$variable == "flw2007", 2007, ifelse(fflw$variable == "flw2008", 2008, ifelse(fflw$variable  == "FlwTillers2", 2009, ifelse(fflw$variable  == "FlwTillers3", 2010, ifelse(fflw$variable  == "FlwTillers4", 2011, ifelse(fflw$variable  == "FlwTillers5", 2012, ifelse(fflw$variable  == "FlwTillers6", 2013,ifelse(fflw$variable == "FlwTillers7", 2014,ifelse(fflw$variable == "FlwTillers8", 2015,ifelse(fflw$variable  == "FlwTillers9", 2016, NA))))))))))
# View(fflw)

f_seedmerge_ss <- full_join(fspike1, fseed1, by = c("plot", "pos", "tag", "Endo", 
                                                  "Birth","year","tillerid"))
# View(f_seedmerge_ss)

f_seedmerge_ssf <- merge(f_seedmerge_ss, fflw, by = c("plot", "pos", "tag", "Endo", 
                                                      "Birth","year"), all = TRUE)
# View(f_seedmerge_ssf)


# Combining repro measurements across years for the FESU recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
f_rseed <- FESU_data_r %>%
  rename("Birth" = "Date", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C"))

f_rseed$year<- ifelse(f_rseed$variable == "seed2007", 2007, ifelse(f_rseed$variable == "seed2008", 2008, ifelse(f_rseed$variable == "seed2009", 2009, ifelse(f_rseed$variable == "seed2010", 2010, ifelse(f_rseed$variable == "seed2011", 2011, ifelse(f_rseed$variable  == "seed2012", 2012, ifelse(f_rseed$variable  == "seed2013", 2013, ifelse(f_rseed$variable  == "seed2014", 2014, ifelse(f_rseed$variable  == "seed2015", 2015, ifelse(f_rseed$variable  == "seed2016", 2016, NA))))))))))
f_rseed1 <- rseed %>% 
  filter(!is.na(seed))
# View(f_rseed1)

f_rspike <- FESU_data_r %>%
  rename("Birth" = "Date", "Endo" = "endo") %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012",
                       "spikelets_infla13", "spikelets_inflb13",
                       "spikelets_infla14", "spikelets_inflb14",
                       "spikelets_infla15", "spikelets_inflb15",
                       "spikelets_infla16", "spikelets_inflb16",
                       "spikelets_inflc16"
                       ),
       value.name = "spikelets")   %>% 
  mutate(tillerid = case_when(grepl("a", variable) ~ "A", 
                              grepl("b",variable) ~ "B",
                              grepl("c", variable) ~ "C"))

f_rspike$year<- ifelse(f_rspike$variable == "spike2007", 2007, ifelse(f_rspike$variable == "spike2008", 2008, ifelse(f_rspike$variable == "spike2009", 2009, ifelse(f_rspike$variable == "spike2010", 2010, ifelse(f_rspike$variable == "spike2011", 2011, ifelse(f_rspike$variable  == "spike2012", 2012, ifelse(f_rspike$variable  == "spikelets_infla13", 2013, ifelse(f_rspike$variable  == "spikelets_inflb13", 2013, ifelse(f_rspike$variable  == "spikelets_infla14", 2014, ifelse(f_rspike$variable  == "spikelets_inflb14", 2014, ifelse(f_rspike$variable  == "spikelets_infla15", 2015, ifelse(f_rspike$variable  == "spikelets_inflb15", 2015, ifelse(f_rspike$variable  == "spikelets_infla16", 2016, ifelse(f_rspike$variable  == "spikelets_inflb16", 2016, ifelse(f_rspike$variable  == "spikelets_inflc16", 2016,NA)))))))))))))))
f_rspike1 <- f_rspike %>% 
  filter(!is.na(spikelets))
# View(f_rspike1)

f_rflw <- FESU_data_r %>%
  rename("Birth" = "Date", "Endo" = "endo") %>% 
  mutate(flw2007 = NA, flw2008 = NA, flw2009 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "flw2009", "FLW10", "FLW11", "FLW12", 
                       "FLW13", "FLW14", "FLW15",
                       "FLW16"),
       value.name = "flw") %>% 
  distinct() # There are a few duplicated rows of tag id's
f_rflw$year<- ifelse(f_rflw$variable == "flw2007", 2007, ifelse(f_rflw$variable == "flw2008", 2008,ifelse(f_rflw$variable == "flw2009", 2009,ifelse(f_rflw$variable == "FLW10", 2010, ifelse(f_rflw$variable == "FLW11", 2011, ifelse(f_rflw$variable  == "FLW12", 2012, ifelse(f_rflw$variable  == "FLW13", 2013, ifelse(f_rflw$variable  == "FLW14", 2014, ifelse(f_rflw$variable  == "FLW15", 2015, ifelse(f_rflw$variable  == "FLW16", 2016, NA))))))))))
# View(f_rflw)

f_rseedmerge_ss <- full_join(f_rspike1, f_rseed1, by = c("plot", "pos", "tag", "Endo", 
                                                       "Birth","year","tillerid"))
# View(f_rseedmerge_ss)

f_rseedmerge_ssf <- merge(f_rseedmerge_ss, f_rflw, by = c("plot", "pos", "tag", "Endo", 
                                                          "Birth","year"), all = TRUE)
# View(f_rseedmerge_ssf)








# Combining the original and recruit FESU repro dataframes ----------------------
f_seedmerge_ssf <- f_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                             "tillerid", "flw","spikelets", "seed")]

f_rseedmerge_ssf <- f_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                               "tillerid", "flw","spikelets", "seed")]

FESUrepro <- f_seedmerge_ssf %>% 
  rbind(f_rseedmerge_ssf) %>% 
  mutate(species = "FESU")

# FESUrepro <- FESUrepro[!(is.na(FESUrepro$flw)),]

# View(FESUrepro)






# Combining repro measurements across years for the ELVI data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken
ELVI_seed_tiller <- merge(ELVI_data, ELVI_data_seed, by = "TAG") #there is a separate sheet with the main raw data for seeds. These are recorded as seeds/infl

elviseed <- ELVI_seed_tiller %>%
  rename("Birth" = "Planted.Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>% 
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(seed2007 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "Seeds_Infl08", 
                       "INfl1Seeds09","Infl2Seeds09","Infl3Seeds09","Infl4Seeds09", 
                       "Infl1Seeds10","Infl2Seeds10",
                       "Infl1Seeds11","Infl2Seeds11",
                       "Seeds1_Infl12", "Seeds2_Infl12",
                       "Seeds1_Infl13", "Seeds2_Infl13",
                       "Seeds1_Infl14", "Seeds2_Infl14",
                       "Seeds1_Infl15", "Seeds2_Infl15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "seed")    %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1S", variable) ~ "A", 
                              grepl("fl2S",variable) ~ "B",
                              grepl("fl3S", variable) ~ "C",
                              grepl("fl4S", variable) ~ "D",
                              grepl("Seeds_", variable) ~ "A",
                              grepl("Seeds1", variable) ~ "A",
                              grepl("Seeds2", variable) ~ "B"))

elviseed$year<- ifelse(elviseed$variable == "seed2007", 2007, ifelse(elviseed$variable == "Seeds_Infl08", 2008, ifelse(elviseed$variable  == "INfl1Seeds09", 2009, ifelse(elviseed$variable  == "Infl2Seeds09", 2009, ifelse(elviseed$variable  == "Infl3Seeds09", 2009, ifelse(elviseed$variable  == "Infl4Seeds09", 2009, ifelse(elviseed$variable  == "Infl1Seeds10", 2010, ifelse(elviseed$variable  == "Infl2Seeds10", 2010,ifelse(elviseed$variable  == "Infl1Seeds11", 2011,ifelse(elviseed$variable  == "Infl2Seeds11", 2011, ifelse(elviseed$variable  ==  "Seeds1_Infl12", 2012, ifelse(elviseed$variable  ==  "Seeds2_Infl12", 2012, ifelse(elviseed$variable  ==  "Seeds1_Infl13", 2013, ifelse(elviseed$variable  ==  "Seeds2_Infl13", 2013, ifelse(elviseed$variable =="Seeds1_Infl14", 2014, ifelse(elviseed$variable =="Seeds2_Infl14", 2014, ifelse(elviseed$variable == "Seeds1_Infl15", 2015, ifelse(elviseed$variable == "Seeds2_Infl15", 2015, ifelse(elviseed$variable  =="Seeds1_Infl16", 2016, ifelse(elviseed$variable  =="Seeds2_Infl16", 2016, NA))))))))))))))))))))
elviseed1 <- elviseed %>% 
  filter(!is.na(seed))
# View(elviseed1)

elvispike <- ELVI_seed_tiller %>% 
  rename("Birth" = "Planted.Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012","spike2013", 
                       "spike2014", "spike2015", "spike2016"), 
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))

elvispike$year<- ifelse(elvispike$variable == "spike2007", 2007, ifelse(elvispike$variable == "spike2008", 2008, ifelse(elvispike$variable  == "spike2009", 2009, ifelse(elvispike$variable  == "spike2010", 2010, ifelse(elvispike$variable  == "spike2011", 2011, ifelse(elvispike$variable  == "spike2012", 2012, ifelse(elvispike$variable  == "spike2013", 2013, ifelse(elvispike$variable == "spike2014", 2014, ifelse(elvispike$variable == "spike2015", 2015, ifelse(elvispike$variable  == "spike2016", 2016, NA))))))))))
elvispike1 <- elvispike %>% 
  filter(!is.na(spikelets))
# View(elvispike1)

elviflw <- ELVI_seed_tiller %>% 
  rename("Birth" = "Planted.Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG") %>%
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(flw2007 = 0) %>% 
 melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers08", "FlwTillers09", "FlwTillers10", 
                       "FlwTillers11", "FlwTillers12", "FlwTillers13", 
                       "FlwTillers14", "FlwTillers15", "FlwTillers16"), 
       value.name = "flw") 
elviflw$year<- ifelse(elviflw$variable == "flw2007", 2007,ifelse(elviflw$variable == "FlwTillers08", 2008, ifelse(elviflw$variable  == "FlwTillers09", 2009, ifelse(elviflw$variable  == "FlwTillers10", 2010, ifelse(elviflw$variable  == "FlwTillers11", 2011, ifelse(elviflw$variable  == "FlwTillers12", 2012, ifelse(elviflw$variable  == "FlwTillers13", 2013,ifelse(elviflw$variable == "FlwTillers14", 2014,ifelse(elviflw$variable == "FlwTillers15", 2015,ifelse(elviflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elviflw)

elvi_seedmerge_ss <- full_join(elvispike1, elviseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                        "Birth","year","tillerid"))
# View(elvi_seedmerge_ss)

elvi_seedmerge_ssf <- merge(elvi_seedmerge_ss, elviflw, by = c("plot", "pos", "tag", "Endo", 
                                                               "Birth","year"), all = TRUE)
# View(elvi_seedmerge_ssf)



# Combining repro measurements across years for the ELVI recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
elvi_rseed <- ELVI_data_r %>%
  rename("Birth" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", 
                       "seed2010", "seed2011", "seed2012", 
                       "seed2013", "seed2014", 
                       "seeds.infl1.15","seeds.infl2.15", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B"))

elvi_rseed$year<- ifelse(elvi_rseed$variable == "seed2007", 2007,ifelse(elvi_rseed$variable == "seed2008", 2008,ifelse(elvi_rseed$variable == "seed2009", 2009, ifelse(elvi_rseed$variable == "seed2010", 2010, ifelse(elvi_rseed$variable == "seed2011", 2011, ifelse(elvi_rseed$variable  == "seed2012", 2012, ifelse(elvi_rseed$variable  == "seed2013", 2013, ifelse(elvi_rseed$variable  == "seed2014", 2014, ifelse(elvi_rseed$variable  == "seeds.infl1.15", 2015, ifelse(elvi_rseed$variable  == "seeds.infl2.15", 2015, ifelse(elvi_rseed$variable == "seed2016", 2016, NA)))))))))))
elvi_rseed1 <- elvi_rseed %>% 
  filter(!is.na(seed))
# View(elvi_rseed1)

elvi_rspike <- ELVI_data_r %>%
  rename("Birth" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA, spike2014 = NA, spike2015 = NA, spike2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010","spike2011", "spike2012", 
                       "spike2013", "spike2014", "spike2015", "spike2016"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("C", variable) ~ "C",
                              grepl("fl1", variable) ~ "A", 
                              grepl("fl2",variable) ~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))
elvi_rspike$year<- ifelse(elvi_rspike$variable == "spike2007", 2007, ifelse(elvi_rspike$variable == "spike2008", 2008, ifelse(elvi_rspike$variable == "spike2009", 2009, ifelse(elvi_rspike$variable == "spike2010", 2010, ifelse(elvi_rspike$variable == "spike2011", 2011, ifelse(elvi_rspike$variable  == "spike2012", 2012, ifelse(elvi_rspike$variable  == "spike2013", 2013, ifelse(elvi_rspike$variable  == "spike2014", 2014, ifelse(elvi_rspike$variable  == "spike2015", 2015, NA)))))))))
elvi_rspike1 <- elvi_rspike %>% 
  filter(!is.na(spikelets))
# View(elvi_rspike1)

elvi_rflw <- ELVI_data_r %>%
  rename("Birth" = "Year", "plot" = "Plot", "pos" = "RecruitNo") %>%
  mutate(flw2007 = NA, flw2008 =NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "FlwTillers09", "FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15"),
       value.name = "flw") 
elvi_rflw$year<- ifelse(elvi_rflw$variable == "flw2007", 2007, ifelse(elvi_rflw$variable == "flw2008", 2008, ifelse(elvi_rflw$variable == "FlwTillers09", 2009, ifelse(elvi_rflw$variable == "FlwTillers10", 2010, ifelse(elvi_rflw$variable == "FlwTillers11", 2011, ifelse(elvi_rflw$variable  == "FlwTillers12", 2012, ifelse(elvi_rflw$variable  == "FlwTillers13", 2013, ifelse(elvi_rflw$variable  == "FlwTillers14", 2014, ifelse(elvi_rflw$variable  == "FlwTillers15", 2015, ifelse(elvi_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elvi_rflw)

elvi_rseedmerge_ss <- full_join(elvi_rspike1, elvi_rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                                 "Birth","year","tillerid"))
# View(elvi_rseedmerge_ss)

elvi_rseedmerge_ssf <- merge(elvi_rseedmerge_ss, elvi_rflw, by = c("plot", "pos", "tag", "Endo", 
                                                                   "Birth","year"), all = TRUE)
# View(elvi_rseedmerge_ssf)




# Combining the  original and recruit ELVI repro dataframes ---------
elvi_seedmerge_ssf <- elvi_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                           "tillerid", "flw","spikelets", "seed")]

elvi_rseedmerge_ssf <- elvi_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                             "tillerid", "flw","spikelets", "seed")]


ELVIrepro <- elvi_seedmerge_ssf %>% 
  rbind(elvi_rseedmerge_ssf) %>% 
  mutate(species = "ELVI") %>% 
  mutate(spikelets = seed)

# ELVIrepro <- ELVIrepro[!(is.na(ELVIrepro$flw)),]

# View(ELVIrepro)








# Combining repro measurements across years for the ELRI data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken
ELRI_seed_tiller <- merge(ELRI_data, ELRI_data_seed2009, by.x = c("PLOT", "POS", "TAG", "ENDO", "Planted.Date","TRT", "Loc.n", "Plant"), by.y = c("plot", "pos", "tag", "Endo", "Planted.Date","TRT", "Loc.n", "Plant"), all.x = TRUE) #there is a separate sheet with the main raw data for seeds.


elriseed <- ELRI_seed_tiller %>%
  rename("Birth" = "Planted.Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>% 
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(seed2007 = NA, seed2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seed2008", 
                       "Infl1Seeds2", "Infl2Seeds2", "Infl3Seeds2", "Infl4Seeds2", 
                       "Seeds_Infl10","seeds_1_11","seeds_2_11",
                       "seeds_1_12","seeds_2_12",
                       "seeds_1_13","seeds_2_13",
                       "seeds_1_14","seeds_2_14",
                       "seeds_1_15","seeds_2_15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("_1_", variable) ~ "A",
                              grepl("_2_", variable) ~ "B", 
                              grepl("fl1S",variable) ~ "A",
                              grepl("fl2S",variable)~ "B",
                              grepl("fl3S", variable) ~ "C",
                              grepl("fl4S", variable) ~ "D",
                              grepl("Infl10", variable) ~ "A"))

elriseed$year<- ifelse(elriseed$variable == "seed2007", 2007, ifelse(elriseed$variable == "seed2008", 2008, ifelse(elriseed$variable  == "Infl1Seeds2", 2009, ifelse(elriseed$variable  == "Infl2Seeds2", 2009, ifelse(elriseed$variable  == "Infl3Seeds2", 2009, ifelse(elriseed$variable  == "Infl4Seeds2", 2009, ifelse(elriseed$variable  == "Seeds_Infl10", 2010, ifelse(elriseed$variable  == "seeds_1_11", 2011, ifelse(elriseed$variable  == "seeds_2_11", 2011, ifelse(elriseed$variable  == "seeds_1_12", 2012, ifelse(elriseed$variable  == "seeds_2_12", 2012,ifelse(elriseed$variable  == "seeds_1_13", 2013,ifelse(elriseed$variable  == "seeds_2_13", 2013, ifelse(elriseed$variable == "seeds_1_14", 2014,ifelse(elriseed$variable == "seeds_2_14", 2014, ifelse(elriseed$variable == "seeds_1_15", 2015, ifelse(elriseed$variable == "seeds_2_15", 2015,ifelse(elriseed$variable  == "Seeds1_Infl16", 2016,ifelse(elriseed$variable  == "Seeds2_Infl16", 2016, NA)))))))))))))))))))
elriseed1 <- elriseed %>% 
  filter(!is.na(seed), seed != ".")
# View(elriseed1)

elrispike <- ELRI_seed_tiller %>% 
  rename("Birth" = "Planted.Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(spike2007 = NA, spike2008 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT"), 
       measure.var = c("spike2007", "spike2008", "Infl1Seeds2", "Infl2Seeds2", "Infl3Seeds2", "Infl4Seeds2", 
                       "Seeds_Infl10","seeds_1_11","seeds_2_11",
                       "seeds_1_12","seeds_2_12",
                       "seeds_1_13","seeds_2_13",
                       "seeds_1_14","seeds_2_14",
                       "seeds_1_15","seeds_2_15",
                       "Seeds1_Infl16", "Seeds2_Infl16"), 
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("_1_", variable) ~ "A",
                              grepl("_2_", variable) ~ "B", 
                              grepl("fl1",variable) ~ "A",
                              grepl("fl2",variable)~ "B",
                              grepl("fl3", variable) ~ "C",
                              grepl("fl4", variable) ~ "D"))
elrispike$year<- ifelse(elrispike$variable == "spike2007", 2007, ifelse(elrispike$variable == "spike2008", 2008, ifelse(elrispike$variable  == "Infl1Seeds2", 2009, ifelse(elrispike$variable  == "Infl2Seeds2", 2009, ifelse(elrispike$variable  == "Infl3Seeds2", 2009, ifelse(elrispike$variable  == "Infl4Seeds2", 2009, ifelse(elrispike$variable  == "Seeds_Infl10", 2010, ifelse(elrispike$variable  == "seeds_1_11", 2011, ifelse(elrispike$variable  == "seeds_2_11", 2011, ifelse(elrispike$variable  == "seeds_1_12", 2012, ifelse(elrispike$variable  == "seeds_2_12", 2012,ifelse(elrispike$variable  == "seeds_1_13", 2013,ifelse(elrispike$variable  == "seeds_2_13", 2013, ifelse(elrispike$variable == "seeds_1_14", 2014,ifelse(elrispike$variable == "seeds_2_14", 2014, ifelse(elrispike$variable == "seeds_1_15", 2015, ifelse(elrispike$variable == "seeds_2_15", 2015,ifelse(elrispike$variable  == "Seeds1_Infl16", 2016,ifelse(elrispike$variable  == "Seeds2_Infl16", 2016, NA)))))))))))))))))))
elrispike1 <- elrispike %>% 
  filter(!is.na(spikelets), spikelets != ".")
# View(elrispike1)

elriflw <- ELRI_seed_tiller %>% 
  rename("Birth" = "Planted.Date", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG", "Endo" = "ENDO") %>%
  mutate(Birth = year(as.character(Birth))) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "FlwTillers08", "FlwTillers09.x", "FLwTillers10", 
                       "FlwTiller11", "FlwTiller12", "FlwTiller13", 
                       "FlwTiller14", "FlwTiller15", "FlwTillers16"), 
       value.name = "flw") 
elriflw$year<- ifelse(elriflw$variable == "flw2007", 2007, ifelse(elriflw$variable == "FlwTillers08", 2008, ifelse(elriflw$variable  == "FlwTillers09.x", 2009, ifelse(elriflw$variable  == "FLwTillers10", 2010, ifelse(elriflw$variable  == "FlwTiller11", 2011, ifelse(elriflw$variable  == "FlwTiller12", 2012, ifelse(elriflw$variable  == "FlwTiller13", 2013,ifelse(elriflw$variable == "FlwTiller14", 2014,ifelse(elriflw$variable == "FlwTiller15", 2015,ifelse(elriflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elriflw)

elri_seedmerge_ss <- full_join(elrispike1, elriseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                            "Birth","year","tillerid"))
# View(elri_seedmerge_ss)

elri_seedmerge_ssf <- merge(elri_seedmerge_ss, elriflw, by = c( "plot", "pos", "tag", "Endo", 
                                                                "Birth","year"), all = TRUE)
# View(elri_merge_ssf)


# Combining repro measurements across years for the ELRI recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
elri_rseed <- ELRI_data_r %>%
  rename("Birth" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", 
                       "seed2011", "seed2012", "seed2013",
                       "seeds1_14","seeds2_14", "seeds3_14",
                       "seeds1_15","seeds2_15", "seeds3_15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "seed") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("seeds1_", variable) ~ "A",
                              grepl("seeds2_", variable) ~ "B"))
elri_rseed$year<- ifelse(elri_rseed$variable == "seed2007", 2007,ifelse(elri_rseed$variable == "seed2008", 2008,ifelse(elri_rseed$variable == "seed2009", 2009,ifelse(elri_rseed$variable == "seed2010", 2010, ifelse(elri_rseed$variable == "seed2011", 2011, ifelse(elri_rseed$variable  == "seed2012", 2012, ifelse(elri_rseed$variable  == "seed2013", 2013, ifelse(elri_rseed$variable  == "seeds1_14", 2014, ifelse(elri_rseed$variable  == "seeds2_14", 2014, ifelse(elri_rseed$variable  == "seeds3_14", 2014, ifelse(elri_rseed$variable  == "seeds1_15", 2015, ifelse(elri_rseed$variable  == "seeds2_15", 2015, ifelse(elri_rseed$variable  == "seeds3_15", 2015, ifelse(elri_rseed$variable == "Seeds1_Infl16", 2016, ifelse(elri_rseed$variable == "Seeds2_Infl16", 2016,NA)))))))))))))))
elri_rseed1 <- elri_rseed %>% 
  filter(!is.na(seed))
# View(elri_rseed1)

elri_rspike <- ELRI_data_r %>%
  rename("Birth" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA, spike2011 = NA, spike2012 = NA, spike2013 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", "spike2010",
                       "spike2011", "spike2012","spikelets_infla13","spikelets_inflb13", 
                       "seeds1_14","seeds2_14", "seeds3_14",
                       "seeds1_15","seeds2_15", "seeds3_15",
                       "Seeds1_Infl16", "Seeds2_Infl16"),
       value.name = "spikelets") %>% 
  mutate(tillerid = case_when(grepl("Seeds1_", variable) ~ "A", 
                              grepl("Seeds2_",variable) ~ "B",
                              grepl("seeds1_", variable) ~ "A",
                              grepl("seeds2_", variable) ~ "B",
                              grepl("infla", variable) ~ "A",
                              grepl("inflb", variable) ~ "B"))
elri_rspike$year<- ifelse(elri_rspike$variable == "spike2007", 2007, ifelse(elri_rspike$variable == "spike2008", 2008, ifelse(elri_rspike$variable == "spike2009", 2009, ifelse(elri_rspike$variable == "spike2010", 2010, ifelse(elri_rspike$variable == "spike2011", 2011, ifelse(elri_rspike$variable  == "spike2012", 2012, ifelse(elri_rspike$variable  == "spikelets_infla13", 2013, ifelse(elri_rspike$variable  == "spikelets_inflb13", 2013, ifelse(elri_rspike$variable  == "seeds1_14", 2014, ifelse(elri_rspike$variable  == "seeds2_14", 2014, ifelse(elri_rspike$variable  == "seeds3_14", 2014, ifelse(elri_rspike$variable  == "seeds1_2015", 2015,ifelse(elri_rspike$variable  == "seeds2_2015", 2015,ifelse(elri_rspike$variable  == "seeds3_2015", 2015, ifelse(elri_rspike$variable == "Seeds1_Infl16", 2016,ifelse(elri_rspike$variable == "Seeds2_Infl16", 2016, NA))))))))))))))))
elri_rspike1 <- elri_rspike %>% 
  filter(!is.na(spikelets))
# View(elri_rspike1)

elri_rflw <- ELRI_data_r %>%
  rename("Birth" = "Date", "plot" = "PLOT", "pos" = "RecruitNo", "Endo" = "endo") %>%
  mutate(seed2007 = NA, seed2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "FLWtiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLW13", "FLW14", "FLW15", "FlwTillers16"),
       value.name = "flw") 
elri_rflw$year<- ifelse(elri_rflw$variable == "seed2007", 2007, ifelse(elri_rflw$variable == "seed2008", 2008, ifelse(elri_rflw$variable == "FLWtiller09", 2009, ifelse(elri_rflw$variable == "FLWtiller10", 2010, ifelse(elri_rflw$variable == "FLWtiller11", 2011, ifelse(elri_rflw$variable  == "FLWtiller12", 2012, ifelse(elri_rflw$variable  == "FLW13", 2013, ifelse(elri_rflw$variable  == "FLW14", 2014, ifelse(elri_rflw$variable  == "FLW15", 2015, ifelse(elri_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(elri_rflw)

elri_rseedmerge_ss <- full_join(elri_rspike1, elri_rseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                                 "Birth","year","tillerid"))
# View(elri_rseedmerge_ss)

elri_rseedmerge_ssf <- merge(elri_rseedmerge_ss, elri_rflw, by = c( "plot", "pos", "tag", "Endo", 
                                                                    "Birth","year"), all = TRUE)
# View(elri_rseedmerge_ssf)






# Combining the  original and recruit ELRI repro dataframes ---------
elri_seedmerge_ssf <- elri_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                           "tillerid", "flw","spikelets", "seed")]

elri_rseedmerge_ssf <- elri_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                             "tillerid", "flw","spikelets", "seed")]


ELRIrepro <- elri_seedmerge_ssf %>% 
  rbind(elri_rseedmerge_ssf) %>% 
  mutate(species = "ELRI")

# ELRIrepro <- ELRIrepro[!(is.na(ELRIrepro$flw)),]

# View(ELRIrepro)







# Combining repro measurements across years for the AGPE data ------------------

## Combining measurements across years for Seed, Spikelet, and Flowering using melt
## Recoding those measurements for the year they are taken
### The repro data for AGPE is pretty funky. These are recorded as averages from multiple tillers as opposed to as counts with tiller id's for the other species.
### The 2016 spikelet values are total spikelets for the plants, but all but one of the plants have only one tiller, meaning that for most of the plants the data is essentially spikelets per inflorescence. Currently I left this in, but the plant is plot 120, Pos 17, tag 2397. 
### spikelet data is also recorded as tot spikes for a few years of recruits data.
agpeseed <- AGPE_data %>%
  rename("Birth" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG...4") %>% 
  mutate("Birth" = year(as.character(Birth))) %>% 
  mutate(seed2007 = NA,  seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"),
       measure.var = c("seed2007", "seeds_spikelet1", "seeds_spikelet2", 
                       "seeds_spikelet3", "seed2011",  "seed2012", "seed2013", 
                       "seed2014", "seed2015", "seed2016"),
       value.name = "seed") %>% 
  mutate(tillerid = "multitillermean")
agpeseed$year<-  ifelse(agpeseed$variable == "seed2007", 2007, ifelse(agpeseed$variable == "seeds_spikelet1", 2008, ifelse(agpeseed$variable  == "seeds_spikelet2", 2009, ifelse(agpeseed$variable  == "seeds_spikelet3", 2010, ifelse(agpeseed$variable  == "seed2011", 2011, ifelse(agpeseed$variable  == "seed2012", 2012, ifelse(agpeseed$variable  == "seed2013", 2013,ifelse(agpeseed$variable == "seed2014", 2014,ifelse(agpeseed$variable == "seed2015", 2015,ifelse(agpeseed$variable  == "seed2016", 2016, NA))))))))))
agpeseed1 <- agpeseed %>% 
  filter(!is.na(seed))
# View(agpeseed1)

agpespike <- AGPE_data %>% 
  rename("Birth" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG...4") %>%
  mutate("Birth" = year(as.character(Birth))) %>% 
  mutate(spike2007 = NA,) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("spike2007", "Spikelets_tiller1", "Spikelets_Infl2", "no_total_spikelets_infl3",
                       "avg_spikelets4", "avg_spikelets5", "avg_spikelets6", 
                       "avg_spikelets7", "spike_infl8", "TotSpikelets9"), 
       value.name = "spikelets")  %>% 
  mutate(tillerid = "multitillermean")
agpespike$year<- ifelse(agpespike$variable == "spike2007", 2007, ifelse(agpespike$variable == "Spikelets_tiller1", 2008, ifelse(agpespike$variable  == "Spikelets_Infl2", 2009, ifelse(agpespike$variable  == "no_total_spikelets_infl3", 2010, ifelse(agpespike$variable  == "avg_spikelets4", 2011, ifelse(agpespike$variable  == "avg_spikelets5", 2012, ifelse(agpespike$variable  == "avg_spikelets6", 2013, ifelse(agpespike$variable == "avg_spikelets7", 2014, ifelse(agpespike$variable == "spike_infl8", 2015, ifelse(agpespike$variable  == "TotSpikelets9", 2016, NA))))))))))
agpespike1 <- agpespike %>% 
  filter(!is.na(spikelets))
# View(agpespike1)
 

agpeflw <- AGPE_data %>% 
  rename("Birth" = "PlantedDate", "plot" = "PLOT", "pos" = "POS", 
         "tag" = "TAG...4") %>%
  mutate("Birth" = year(as.character(Birth))) %>% 
  mutate(flw2007 = 0) %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc.n", "Birth", 
                  "TRT", "Plant"), 
       measure.var = c("flw2007", "Flwtillers1", "FlwTillers2", "FLWTiller3...78", 
                       "FLWTiller4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
agpeflw$year<- ifelse(agpeflw$variable == "flw2007", 2007, ifelse(agpeflw$variable == "Flwtillers1", 2008, ifelse(agpeflw$variable  == "FlwTillers2", 2009, ifelse(agpeflw$variable  == "FLWTiller3", 2010, ifelse(agpeflw$variable  == "FLWTiller4", 2011, ifelse(agpeflw$variable  == "FlwTillers5", 2012, ifelse(agpeflw$variable  == "FlwTillers6", 2013,ifelse(agpeflw$variable == "FlwTillers7", 2014,ifelse(agpeflw$variable == "FlwTillers8", 2015,ifelse(agpeflw$variable  == "FlwTillers9", 2016, NA))))))))))
# View(agpeflw)

agpe_seedmerge_ss <- full_join(agpespike1, agpeseed1, by = c( "plot", "pos", "tag", "Endo", 
                                                    "Birth","year","tillerid"))
# View(agpe_seedmerge_ss)

agpe_seedmerge_ssf <- merge(agpe_seedmerge_ss, agpeflw, by = c("plot", "pos", "tag", "Endo", 
                                                       "Birth","year"), all = TRUE)
# View(agpe_seedmerge_ssf)


# Combining repro measurements across years for the AGPE recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
### The spikelets data is also reported in Total spikelets for several years. The first two years there aren't any flwing tillers, and then in 2011 there are only single tillers, so this could essentially be spikelets/infl
### 2013 has lots of spikelet data, but this it is recorded as totals for the plant within an equatio, so I am calculating avgs from it with the number of flw tillers
### 2014 has data recorded as avg spikelets
agpe_rseed <- AGPE_data_r %>%
  rename("Birth" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(seed2007 = NA, seed2008 = NA, seed2009 = NA, seed2010 = NA, seed2011 = NA, seed2012 = NA, seed2013 = NA, seed2014 = NA, seed2015 = NA, seed2016 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("seed2007", "seed2008", "seed2009", "seed2010", "seed2011", "seed2012", "seed2013", "seed2014", 
                       "seed2015", "seed2016"),
       value.name = "seed")  %>% 
  mutate(tillerid = "multitillermean")
agpe_rseed$year<- ifelse(agpe_rseed$variable == "seed2007", 2007, ifelse(agpe_rseed$variable == "seed2008", 2008, ifelse(agpe_rseed$variable == "seed2009", 2009, ifelse(agpe_rseed$variable == "seed2010", 2010, ifelse(agpe_rseed$variable == "seed2011", 2011, ifelse(agpe_rseed$variable  == "seed2012", 2012, ifelse(agpe_rseed$variable  == "seed2013", 2013, ifelse(agpe_rseed$variable  == "seed2014", 2014, ifelse(agpe_rseed$variable  == "seed2015", 2015, ifelse(agpe_rseed$variable == "seed2016", 2016, NA))))))))))
agpe_rseed1 <- agpe_rseed %>% 
  filter(!is.na(seed))
# View(agpe_rseed1)

agpe_rspike <- AGPE_data_r %>%
  rename("Birth" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(spike2007 = NA, spike2008 = NA, spike2009 = NA, spike2010 = NA) %>% 
  mutate(avgspike13 = TotSpikelets13/FlwTillers13) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("spike2007", "spike2008", "spike2009", 
                       "spike2010","TotSpikelets11", "SpikeletsA12", 
                       "avgspike13", "avgspikepertiller14", "spikepertillerA15","spikepertillerB15",
                       "spikepertillerA16","spikepertillerB16"),
       value.name = "spikelets")  %>% 
  mutate(tillerid = case_when(grepl("A", variable) ~ "A", 
                              grepl("B",variable) ~ "B",
                              grepl("avg", variable) ~ "multitillermean",
                              grepl("Spikelets", variable) ~ "multitillermean"))
agpe_rspike$year<-  ifelse(agpe_rspike$variable == "spike2007", 2007, ifelse(agpe_rspike$variable == "spike2008", 2008, ifelse(agpe_rspike$variable == "spike2009", 2009, ifelse(agpe_rspike$variable == "spike2010", 2010, ifelse(agpe_rspike$variable == "TotSpikelets11", 2011, ifelse(agpe_rspike$variable  == "SpikeletsA12", 2012, ifelse(agpe_rspike$variable  == "avgspike13", 2013, ifelse(agpe_rspike$variable  == "avgspikepertiller14", 2014, ifelse(agpe_rspike$variable  == "spikepertillerA15", 2015, ifelse(agpe_rspike$variable  == "spikepertillerB15", 2015,ifelse(agpe_rspike$variable == "spikepertillerA16", 2016, ifelse(agpe_rspike$variable == "spikepertillerB16", 2016, NA))))))))))))
agpe_rspike1 <- agpe_rspike %>% 
  filter(!is.na(spikelets))
# View(agpe_rspike1)

agpe_rflw <- AGPE_data_r %>%
  rename("Birth" = "birth", "plot" = "Plot", 
         "pos" = "RecruitNo", "Endo" = "endo", "tag" = "Tag") %>%
  mutate(flw2007 = NA, flw2008 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth"),
       measure.var = c("flw2007", "flw2008", "FlwTillers09", "FlwTillers10","FlwTillers11", "FlwTillers12", 
                       "FlwTillers13", "FlwTillers14", "FlwTillers15", "FlwTillers16"),
       value.name = "flw") 
agpe_rflw$year<- ifelse(agpe_rflw$variable == "flw2007", 2007, ifelse(agpe_rflw$variable == "flw2008", 2008, ifelse(agpe_rflw$variable == "FlwTillers09", 2009, ifelse(agpe_rflw$variable == "FlwTillers10", 2010, ifelse(agpe_rflw$variable == "FlwTillers11", 2011, ifelse(agpe_rflw$variable  == "FlwTillers12", 2012, ifelse(agpe_rflw$variable  == "FlwTillers13", 2013, ifelse(agpe_rflw$variable  == "FlwTillers14", 2014, ifelse(agpe_rflw$variable  == "FlwTillers15", 2015, ifelse(agpe_rflw$variable  == "FlwTillers16", 2016, NA))))))))))
# View(agpe_rflw)

agpe_rseedmerge_ss <- full_join(agpe_rspike1, agpe_rseed1, by = c( "plot", "pos", "tag", "Endo", "Birth", "year", "tillerid"))
# View(agpe_rseedmerge_ss)

agpe_rseedmerge_ssf <- merge(agpe_rseedmerge_ss, agpe_rflw, by = c("plot", "pos", "tag", "Endo", "Birth", "year"), all = TRUE)
# View(agpe_rseedmerge_ssf)






# Combining the  original and recruit AGPE repro dataframes ---------
agpe_seedmerge_ssf <- agpe_seedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                           "tillerid", "flw","spikelets", "seed")]

agpe_rseedmerge_ssf <- agpe_rseedmerge_ssf[c("plot", "pos", "tag", "Endo", "Birth","year",
                                             "tillerid", "flw","spikelets", "seed")]


AGPErepro <- agpe_seedmerge_ssf %>% 
  rbind(agpe_rseedmerge_ssf) %>% 
  mutate(species = "AGPE")



# View(AGPErepro)

###### Bind together the different species repro datasets and merge with endo_demog_long
LTREB_repro <- AGPErepro %>% 
  merge(ELRIrepro, by = c("plot", "pos","tag", "Endo", "Birth", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(ELVIrepro, by = c("plot", "pos","tag", "Endo", "Birth", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(FESUrepro, by = c("plot", "pos","tag", "Endo", "Birth", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(LOARrepro, by = c("plot", "pos","tag", "Endo", "Birth", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(POALrepro, by = c("plot", "pos","tag", "Endo", "Birth", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  merge(POSYrepro, by = c("plot", "pos","tag", "Endo", "Birth", "year", "species", "flw", "seed", "spikelets", "tillerid"), all = TRUE) %>% 
  mutate(flw = as.numeric(flw), seed = as.numeric(seed), spikelets = as.numeric(spikelets))
  
# There are a few instances within the reproductive data where there is a 0 flw measurement with spikelet data. Often these are recorded in the raw data with a note about how seeds were possibly mislabeled.
# There are more instances where there is a flw measurement but spikelet data is not recorded. In some instances, this is likely due to forgetting to collect the data, and there are sometimes in the raw data 
# files where I tried to extract the spikelet count and there were only opaque seed estimates and no specific spikelet counts recorded (although they possibly were at some point to generate the seed estimate). 
# An example of the latter would be POAL plot 3, tag 41, in year 2013, where the note marks that this tiller was skipped.
# I'm going to filter out cases where seeds or spikelet > 0 and flw = 0
LTREB_repro_flw_spike_mismatches <- LTREB_repro %>% 
  filter(flw==0 & spikelets>0 | flw > 0 & is.na(spikelets))

###########################################################################################################################################################################
###### Cleaning up seed and spikelet information and merging with LTREB_data that will be used in the reproductive kernels -----------------
###########################################################################################################################################################################
LTREB_repro1 <- LTREB_repro %>% 
  rename(endo = Endo) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>%
  mutate(year = as.integer(year)) %>% 
  mutate(`birth` = as.integer(Birth)) %>% 
  mutate(plot_fixed = as.integer(plot)) %>% 
  mutate(spikelets_fixed = as.numeric(case_when(species == "ELVI" & flw > 0 ~ as.character(seed), # ELRI and ELVI has data stored as seed/infl
                                                flw > 0 & !is.na(spikelets) ~ as.character(spikelets),
                                     flw == 0 & spikelets == 0 ~ NA_character_,
                                     flw == 0 & spikelets != 0 ~ NA_character_,
                                     flw == 0 & is.na(spikelets) ~ NA_character_))) %>% 
  mutate(tillerid_fixed = case_when(!is.na(tillerid) & spikelets >= 0 ~ tillerid,
                                    !is.na(tillerid) & seed >= 0 ~ tillerid,
                               is.na(tillerid) & !is.na(spikelets) ~ "A",
                               is.na(tillerid) & is.na(spikelets) ~ NA_character_)) %>% 
  distinct()
write_csv(LTREB_repro1,"~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_repro1.csv")

# spreading out the spikelet info by tiller to create single row per year per individual.
LTREB_repro_wide <- LTREB_repro1 %>%
  dplyr::select(plot_fixed, pos, tag, endo_01, birth, year, species, flw, spikelets_fixed, tillerid_fixed) %>% 
  spread(key = tillerid_fixed, value = spikelets_fixed) %>% 
  rename(spike_a_t1 = A, spike_b_t1 = B, spike_c_t1 = C, spike_d_t1 = D,  spikelets_AGPE_mean = multitillermean) %>% 
  dplyr::select(-'<NA>')
  
  
# View(LTREB_repro2)
dim(LTREB_repro_wide)
table(LTREB_repro_wide$species, LTREB_repro_wide$year, !is.na(LTREB_repro_wide$flw))
table(is.na(LTREB_repro1$flw), LTREB_repro1$flw>0, !is.na(LTREB_repro1$spikelets))

# I am going to try to merge and then add the lagged repro variables at the end
LTREB_repro_t1 <- LTREB_repro_wide%>% 
  rename(flw_t1 = flw,
         year_t1 = year)


# merge the reproductive data with LTREB long data file for recent data and for size information for the reproductive model
# This is endodemoglong which is stored in LTREB_data
head(LTREB_data)

LTREB_data1 <- LTREB_data %>% 
  dplyr::select(-contains("seed"), -plot, -endo)



# now we can merge the two datasets together.
LTREB_repro_combo <- LTREB_data1 %>% 
  left_join(LTREB_repro_t1,
            by = c("plot_fixed" = "plot_fixed", "pos" = "pos",
                   "id" = "tag", "species" = "species",
                   "endo_01" = "endo_01",
                   "year_t1" = "year_t1")) %>% 
  rename("birth" = "birth.x", "birth_fromrepro" = "birth.y",
         "spike_a_t1_long" = "spike_a_t1.x", "spike_a_t1_fromrepro" = "spike_a_t1.y",
         "spike_b_t1_long" = "spike_b_t1.x", "spike_b_t1_fromrepro" = "spike_b_t1.y",
         "spike_c_t1_long" = "spike_c_t1.x", "spike_c_t1_fromrepro" = "spike_c_t1.y",
         "spike_d_t1_fromrepro" = "spike_d_t1", 
         "flw_t1_long" = "flw_t1.x", "flw_t1_fromrepro" = "flw_t1.y")
# View(LTREB_repro_combo)



# Now dplyr::select the correct repro data from long or from the raw files into new master columns
LTREB_full_to2018 <- LTREB_repro_combo %>% 
  mutate(FLW_COUNT_T1 = as.integer(case_when(surv_t1 == 0 ~ NA_integer_, !is.na(flw_t1_fromrepro) & !is.na(flw_t1_long) ~ as.integer(flw_t1_fromrepro),
                                       is.na(flw_t1_fromrepro) & !is.na(flw_t1_long) ~ as.integer(flw_t1_long),
                                       !is.na(flw_t1_fromrepro) & is.na(flw_t1_long) ~ as.integer(flw_t1_fromrepro))),   # In this case, where both datasets had data entered, spot checking showed that they have they were identical
         FLW_STAT_T1 = as.integer(case_when(FLW_COUNT_T1 > 0 & surv_t1 == 1 ~ 1, FLW_COUNT_T1 == 0 & surv_t1 == 1 ~ 0, is.na(FLW_COUNT_T1) & surv_t1 == 1 ~ 0)),
         SPIKE_A_T1 =  case_when(is.na(spike_a_t1_fromrepro) & !is.na(spike_a_t1_long) ~ as.numeric(spike_a_t1_long),
                                 !is.na(spike_a_t1_fromrepro) & is.na(spike_a_t1_long) ~ as.numeric(spike_a_t1_fromrepro),                                                                          
                                 !is.na(spike_a_t1_fromrepro) & !is.na(spike_a_t1_long) ~ as.numeric(spike_a_t1_fromrepro)),
         SPIKE_B_T1 =  case_when(is.na(spike_b_t1_fromrepro) & !is.na(spike_b_t1_long) ~ as.numeric(spike_b_t1_long),
                                 !is.na(spike_b_t1_fromrepro) & is.na(spike_b_t1_long) ~ as.numeric(spike_b_t1_fromrepro),                                                                          
                                 !is.na(spike_b_t1_fromrepro) & !is.na(spike_b_t1_long) ~ as.numeric(spike_b_t1_fromrepro)),
         SPIKE_C_T1 =  case_when(is.na(spike_c_t1_fromrepro) & !is.na(spike_c_t1_long) ~ as.numeric(spike_c_t1_long),
                                 !is.na(spike_c_t1_fromrepro) & is.na(spike_c_t1_long) ~ as.numeric(spike_c_t1_fromrepro),                                                                          
                                 !is.na(spike_c_t1_fromrepro) & !is.na(spike_c_t1_long) ~ as.numeric(spike_c_t1_fromrepro)),
         SPIKE_D_T1 = spike_d_t1_fromrepro,
         SPIKE_AGPE_MEAN_T1 = spikelets_AGPE_mean) # <- This last column is because AGPE has some years of spikelet data collected as a mean, so I am keeping it separate from the other data.
LTREB_full_to2018_lag <- LTREB_full_to2018 %>%         
group_by(id) %>% 
  mutate(FLW_COUNT_T = as.integer(dplyr::lag(FLW_COUNT_T1, n = 1, default = NA)),
         FLW_STAT_T = as.integer(dplyr::lag(FLW_STAT_T1, n = 1, default = NA)),
         SPIKE_A_T = dplyr::lag(SPIKE_A_T1, n = 1, default = NA),
         SPIKE_B_T = dplyr::lag(SPIKE_B_T1, n = 1, default = NA),
         SPIKE_C_T = dplyr::lag(SPIKE_C_T1, n = 1, default = NA),
         SPIKE_D_T = dplyr::lag(SPIKE_D_T1, n = 1, default = NA),
         SPIKE_AGPE_MEAN_T = dplyr::lag(SPIKE_AGPE_MEAN_T1, n = 1, default = NA)) %>% 
  dplyr::select(plot_fixed, plot_index, pos, id, species, species_index, 
                endo_01, endo_index, origin_01, birth,
                year_t1, year_t1_index,
                surv_t1, size_t1, logsize_t1,
                FLW_COUNT_T1, FLW_STAT_T1,
                SPIKE_A_T1, SPIKE_B_T1, SPIKE_C_T1, SPIKE_D_T1, SPIKE_AGPE_MEAN_T1,
                year_t, year_t_index, size_t, logsize_t, 
                FLW_COUNT_T, FLW_STAT_T,
                SPIKE_A_T, SPIKE_B_T, SPIKE_C_T, SPIKE_D_T, SPIKE_AGPE_MEAN_T)



##############################################################################
####### Here we will merge in 2019 and 2020 data ------------------------------
##############################################################################
# 2019 census data
AGPE_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "AGPE")
ELRI_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "ELRI")
ELVI_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "ELVI")
FESU_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "FESU")
LOAR_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "LOAR")
POAL_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "POAL")
POSY_2019_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2019/LTREB_data_2019.xlsx", sheet = "POSY")
# # Now we can merge all the different species together.
LTREB_update_data <- AGPE_2019_data %>% 
  merge(AGPE_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(ELRI_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(ELVI_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(FESU_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(LOAR_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(POAL_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(POSY_2019_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE)

# 2020 census data
AGPE_2020_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "AGPE")
ELRI_2020_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "ELRI")
ELVI_2020_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "ELVI")
FESU_2020_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "FESU")
POAL_2020_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "POAL")
POSY_2020_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2020/LTREB_data_2020.xlsx", sheet = "POSY")
# # Now we can merge all the different species together.
LTREB_update_data <- LTREB_update_data %>%
  merge(AGPE_2020_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(ELRI_2020_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(ELVI_2020_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(FESU_2020_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(POAL_2020_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(POSY_2020_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE)

# 2021 census data
AGPE_2021_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "AGPE")
ELRI_2021_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "ELRI")
ELVI_2021_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "ELVI")
FESU_2021_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "FESU")
POAL_2021_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "POAL")
POSY_2021_data <- read_xlsx("Dropbox/EndodemogData/Field Data/2021/LTREB_data_2021.xlsx", sheet = "POSY")
# # Now we can merge all the different species together.
LTREB_update_data <- LTREB_update_data %>%
  merge(AGPE_2021_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
  merge(ELRI_2021_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
  merge(ELVI_2021_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
  merge(FESU_2021_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>% 
  merge(POAL_2021_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE) %>%
  merge(POSY_2021_data, by = c("species", "origin", "plot", "pos","id",  "birth_year", "observation_year", "species", "distance_A", "distance_B", "survival", "size_tillers", "flowering_tillers", "spikelets_A", "spikelets_B", "spikelets_C", "notes"), all = TRUE)



# # We need to do a little bit of cleaning up for some of the missing tags and changing variable names.
# 
LTREB_update_cleaned <- LTREB_update_data %>%
  rename(year_t1 = observation_year, surv_t1 = survival,
         size_t1 = size_tillers, FLW_COUNT_T1 = flowering_tillers,
         SPIKE_A_T1 = spikelets_A, SPIKE_B_T1 = spikelets_B, SPIKE_C_T1 = spikelets_C,
         dist_a = distance_A, dist_b = distance_B, birth = birth_year) %>% 
    mutate(birth = as.integer(birth)) %>% #There are NA's for a few birth years where we have found new but older plants
    mutate(plot_fixed = as.integer(plot)) %>% 
    mutate(plot_index = as.integer(recode(plot_fixed, !!!plot_factor_key))) %>% 
    mutate(size_t1 = na_if(size_t1, 0)) %>%
    mutate(size_t1, logsize_t1 = log(size_t1)) %>%
    mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1))) %>%
    mutate(FLW_STAT_T1 = case_when(FLW_COUNT_T1 == 0 ~ 0,
                                   FLW_COUNT_T1 > 0 ~1)) %>% 
    mutate(species_index = as.integer(recode(species, !!!species_factor_key))) %>% 
    mutate(year_t1_index = as.integer(recode(year_t1, !!!year_factor_key))) %>%
    mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, # original transplanted plants from 2007
                                            origin == "O21" ~ 0, # original transplanted plants from 2021. In our model, we have an origin effect for recruit vs transplant, and we are treating the transplants the same.
                                            origin == "R" ~ 1, # recruit
                                            origin != "R" | origin != "O" | origin != "O21" ~ 1))) %>% 
    mutate(surv_t1 = as.integer(case_when(surv_t1 == 1 ~ 1,
                                          surv_t1 == 0 ~ 0,
                                          is.na(surv_t1) & birth == year_t1 ~ 1))) %>% 
    filter(!is.na(size_t1) & surv_t1 == 1 | surv_t1 == 0) # filtering out mismatches where there is no size data entry but the survival was recorded; this is sometimes TNF or new recruits where I think it was an oversight in data entry where they are likely a small size like 1 tiller. For now, I am just removing them

# There are a few plants that were not found in 2019 or 2020, but found in the following year. Here I am updating the survival of those plants to be 1, although we do not have size data for them.
LTREB_update_tnfs <- LTREB_update_cleaned %>% 
  group_by(id) %>% 
  filter(year_t1[1] <= year_t1[2] & surv_t1[1] == 0) %>% 
  mutate(surv_t1 = case_when(year_t1 == min(year_t1) ~ 1,
                             TRUE ~ as.numeric(surv_t1)))
dim(LTREB_update_tnfs)
# now rebind the these updated rows to the rest of the data
LTREB_update_cleaned_tnf<- filter(LTREB_update_cleaned, !(id %in% LTREB_update_tnfs$id)) %>% 
  bind_rows(LTREB_update_tnfs)


LTREB_update <- LTREB_update_cleaned_tnf %>% 
  mutate(year_t = year_t1 - 1,
         year_t_index = year_t1_index -1) %>% 
  dplyr::select(plot_fixed, plot_index, pos, id, species, species_index, 
                origin_01, birth,
                year_t1, year_t1_index,
                surv_t1, size_t1, logsize_t1,
                FLW_COUNT_T1, FLW_STAT_T1,
                SPIKE_A_T1, SPIKE_B_T1, SPIKE_C_T1,
                year_t, year_t_index, dist_a, dist_b)
# Assigning plot endo status to the 2019 data
LTREB_plot_endo_status <- LTREB_full_to2018 %>% 
  group_by(plot_fixed, plot_index) %>% 
  summarize(endo_01 = round(mean(endo_01, na.rm = T)),
            endo_index = round(mean(endo_index, na.rm = T)))

LTREB_update <- LTREB_update %>% 
  left_join(LTREB_plot_endo_status)

# Now we can merge our post-2018 data with our full dataframe
LTREB_full_update <- LTREB_full_to2018_lag %>%
  full_join(LTREB_update)


##############################################################################
####### Then we will add lagged variables to have the measurements in time t 
##############################################################################

# I'm creating the lagged variable from the combined data because the 2019 data by itself doesn't have 2018 data
# This leaves some NA's in the year t variables, so I am creating from the merged data frame and then selecting the values without NA's from the new column and from the endo_demog_long's column.
# There is probably a smoother way to do some of this, but I think this works okay.
LTREB_full_update_lag <- LTREB_full_update %>% 
  group_by(id) %>% 
  mutate(SIZE_T_NEW = dplyr::lag(size_t1, n = 1, default = NA),
         LOGSIZE_T_NEW = dplyr::lag(logsize_t1, n = 1, default = NA),
         FLW_COUNT_T_NEW = as.integer(dplyr::lag(FLW_COUNT_T1, n = 1, default = NA)),
         FLW_STAT_T_NEW = as.integer(dplyr::lag(FLW_STAT_T1, n = 1, default = NA)),
         SPIKE_A_T_NEW  = dplyr::lag(SPIKE_A_T1, n = 1, default = NA),
         SPIKE_B_T_NEW  = dplyr::lag(SPIKE_B_T1, n = 1, default = NA),
         SPIKE_C_T_NEW  = dplyr::lag(SPIKE_C_T1, n = 1, default = NA),
         SPIKE_D_T_NEW  = dplyr::lag(SPIKE_D_T1, n = 1, default = NA),
         SPIKE_AGPE_MEAN_T_NEW  = dplyr::lag(SPIKE_AGPE_MEAN_T1, n = 1, default = NA)) %>% 
  mutate(size_t = case_when(!is.na(size_t) ~ size_t,
                           !is.na(SIZE_T_NEW) ~ SIZE_T_NEW),
         logsize_t = case_when(!is.na(logsize_t) ~ logsize_t,
                            !is.na(LOGSIZE_T_NEW) ~ LOGSIZE_T_NEW),
         FLW_COUNT_T = case_when(!is.na(FLW_COUNT_T) ~ FLW_COUNT_T,
                                 !is.na(FLW_COUNT_T_NEW) ~ FLW_COUNT_T_NEW),
         FLW_STAT_T = case_when(!is.na(FLW_STAT_T) ~ FLW_STAT_T,
                                 !is.na(FLW_STAT_T_NEW) ~ FLW_STAT_T_NEW),
         SPIKE_A_T = case_when(!is.na(SPIKE_A_T) ~ SPIKE_A_T,
                                !is.na(SPIKE_A_T_NEW) ~ SPIKE_A_T_NEW),
         SPIKE_B_T = case_when(!is.na(SPIKE_B_T) ~ SPIKE_B_T,
                               !is.na(SPIKE_B_T_NEW) ~ SPIKE_B_T_NEW),
         SPIKE_C_T = case_when(!is.na(SPIKE_C_T) ~ SPIKE_C_T,
                               !is.na(SPIKE_C_T_NEW) ~ SPIKE_C_T_NEW),
         SPIKE_D_T = case_when(!is.na(SPIKE_D_T) ~ SPIKE_D_T,
                               !is.na(SPIKE_D_T_NEW) ~ SPIKE_D_T_NEW),
         SPIKE_AGPE_MEAN_T = case_when(!is.na(SPIKE_AGPE_MEAN_T) ~ SPIKE_AGPE_MEAN_T,
                                       !is.na(SPIKE_AGPE_MEAN_T_NEW) ~ SPIKE_AGPE_MEAN_T_NEW),
         ) %>% 
  dplyr::select(plot_fixed, plot_index, pos, id, species, species_index, 
                endo_01, endo_index, origin_01, birth,
                year_t1, year_t1_index,
                surv_t1, size_t1, logsize_t1,
                FLW_COUNT_T1, FLW_STAT_T1,
                SPIKE_A_T1, SPIKE_B_T1, SPIKE_C_T1, SPIKE_D_T1, SPIKE_AGPE_MEAN_T1,
                year_t, year_t_index, size_t, logsize_t, 
                FLW_COUNT_T, FLW_STAT_T,
                SPIKE_A_T, SPIKE_B_T, SPIKE_C_T, SPIKE_D_T, SPIKE_AGPE_MEAN_T,
                dist_a, dist_b)

dim(LTREB_full_update_lag)

##############################################################################
####### Merging in the endophyte checks ------------------------------
##############################################################################

LTREB_endo_check <- read_csv(file = "Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_status.csv") %>%  
  dplyr::select(-recno,-check, -...11) %>% 
  rename("origin_from_check" = "origin", "endo_status_from_check" = "status", "plot_endo_for_check" = "endo") %>% 
  mutate(origin_01 = as.integer(case_when(origin_from_check == "O" ~ 0, 
                                          origin_from_check == "R" ~ 1))) %>% 
  mutate(plot_endo_for_check = as.integer(recode(plot_endo_for_check, "plus" = 1, "minus" = 0))) %>% 
  mutate(endo_mismatch = plot_endo_for_check - endo_status_from_check) # 0 = no change, >0 = loss of endophyte from positive plot, <0 = gain of endophyte in negative plot
  # Metadata from Jenn
#   recno:	unique record number
#   species:	four letter species code
#   origin:	"O" = original plant from greenhouse, "R" = recruit
#   plot:	plot number
#   pos:	for "O" = position of plant in plot, for "R" = recruit number on tag
#   id:	unique identifier for each plant, "O" = single unique number, "R" = concatenation of plot number and tag number separated by "_"
#   status:	0= no endophyte found via microscopy on leaf peel at 200X, 1= endophyte detected
#   date_status:	year leaf peel was taken
#   endo:	plot level endophyte status: "minus" = no endophyte in original planting, "plus" = endophyte present in original planting
#   check:	"x" indicates the incorrect status was detected given original plot level treatment

# There are two plants that are checked but are not present in the endo_demog_long dataset
setdiff(LTREB_endo_check$id,LTREB_full_update_lag$id)

LTREB_full_2 <- LTREB_full_update_lag %>% 
  left_join(LTREB_endo_check, by = c("species" = "species", "plot_fixed" = "plot", "pos" = "pos", "origin_01" = "origin_01", "id" = "id"))

# here are some summaries of the amount of changes in endophyte status
LTREB_status_changes <- LTREB_full_2 %>% 
  group_by(species, plot_fixed, endo_01) %>% 
  summarize(same = sum(endo_mismatch == 0, na.rm = TRUE),
            lose_endo = sum(endo_mismatch > 0, na.rm = TRUE),
            gain_endo = sum(endo_mismatch <0, na.rm = TRUE)) %>% 
  mutate(percent_gain = (gain_endo/same)*100, percent_lose = (lose_endo/same)*100)
  
LTREB_status_changes_species <- LTREB_full_2 %>% 
  distinct(species, plot_fixed, id, endo_01, endo_mismatch) %>% 
  group_by(species, endo_01) %>% 
  summarize(same = sum(endo_mismatch == 0, na.rm = TRUE),
            lose_endo = sum(endo_mismatch > 0, na.rm = TRUE),
            gain_endo = sum(endo_mismatch <0, na.rm = TRUE)) %>% 
  mutate(percent_gain = (gain_endo/same)*100, percent_lose = (lose_endo/same)*100)

LTREB_status_changes_endo <- LTREB_full_2 %>% 
  distinct(species, plot_fixed, id, endo_01, endo_mismatch) %>% 
  group_by() %>% 
  summarize(same = sum(endo_mismatch == 0, na.rm = TRUE),
            lose_endo = sum(endo_mismatch > 0, na.rm = TRUE),
            gain_endo = sum(endo_mismatch <0, na.rm = TRUE)) %>% 
  mutate(percent_gain = (gain_endo/same)*100, percent_lose = (lose_endo/same)*100, faithfulness = 100-((percent_gain + percent_lose)/2))

# Here is a list of all the original plot endophyte statuses

LTREB_endophyte_plot_numbers <- LTREB_full_2 %>% 
  group_by(species, plot_fixed) %>% 
  summarize(endophyte_status = mean(plot_endo_for_check, na.rm = T)) %>% 
  filter(!is.nan(endophyte_status))
# write_csv(LTREB_endophyte_plot_numbers, path = "LTREB_endophyte_plot_numbers.csv")

##############################################################################
####### Merging in the location data ------------------------------
##############################################################################

LTREB_distances <- read_csv(file = "Dropbox/EndodemogData/Fulldataplusmetadata/endo_distance_tubeid.csv",
                            col_types = cols( species = col_character(),
                            origin = col_character(),
                            plot = col_integer(),
                            pos = col_character(),
                            id = col_character(),
                            dist_a = col_double(),
                            dist_b = col_double(),
                            tubeid = col_character(),
                            notes = col_character(),
                            date_dist = col_character())) %>% 
  dplyr::select(species, origin, plot, pos, id, dist_a, dist_b, date_dist) %>% 
  rename("origin_from_distance" = "origin") %>% 
  mutate(origin_01 = as.integer(case_when(origin_from_distance == "O" ~ 0, 
                                          origin_from_distance == "R" ~ 1))) %>% 
  filter(!is.na(dist_a), !is.na(dist_b)) %>% 
  mutate(duplicate = duplicated(id)) %>% #There are three LOAR id's that have two measurements, one from may 2018 and one from sept 2018: 40_F5, 33_4B, 33_12
  filter(!(species == "LOAR" & date_dist == "may_18" & id %in% c("40_F5", "33_4B", "33_12"))) # I am dplyr::selecting the september measurements for these id's. The distances are different but similar

# Here are the plant id's that are in the distance file but not the long file
setdiff(LTREB_distances$id, LTREB_full_2$id)

##############################################################################
####### Getting climate data and calculating SPEI  ------------------------------
##############################################################################

# Using rnoaa to get weather station data

# From the user manual:
# • PRCP: Precipitation, in tenths of millimeters
# • TAVG: Average temperature, in tenths of degrees Celsius
# • TMAX: Maximum temperature, in tenths of degrees Celsius
# • TMIN: Minimum temperature, in tenths of degrees Celsius

#The lat lon of Lilly-Dickey Woods
lat_lon_df <- data.frame(id = "ldw",
                         lat = 39.2359000,
                         lon = -86.2181000)
# Finding the nearest weather stations within 60 km to Lilly-Dickey Woods
mon_near_ldw <-  meteo_nearby_stations(
    lat_lon_df = lat_lon_df,
    lat_colname = "lat",
    lon_colname = "lon",
    var = "PRCP",
    year_min = 1880,
    year_max = 2020,
    radius = 30, #kilometers
  ) 

# Looking at the temporal coverage of these stations

# mon_near_ldw_coverage <- mon_near_ldw$ldw %>% full_join(meteo_coverage(meteo_pull_monitors(mon_near_ldw$ldw$id, var = "PRCP"))$detail)

# plotting the time coverage of the different stations. Ordered by distance from ldw

# mon_near_ldw_coverage$id <- reorder(mon_near_ldw_coverage$id, mon_near_ldw_coverage$distance)

# mon_near_ldw_coverage %>% 
#   filter(prcp >0) %>% 
# ggplot()+
#   geom_point(aes(y = id, x = date)) + theme_minimal()

# USC00126056 is the closest in Nashville (~1 km away) but has gap in the history and doesn't have coverage more recent than 2019
# USC00120784 is the station near Bloomington (27 km away), with good coverage with the best coverage into the past (although is seems there is a gap in the middle) 
# USC00121747 is the station near Columbus (26 km away in other direction) with seemingly even better coverage into the past
# Both of these seem to be continuing to be update (have current year data)

# Here is a better plot of all the variables from these stations which makes it clear the coverage between Columbus and Bloomington is very similar
monitors <- c("USC00126056", "USC00120784", "USC00121747")
# autoplot(meteo_coverage(meteo_pull_monitors(monitors)))

# Getting lat, long and elevation for each station
station_data <- ghcnd_stations() %>% 
  filter(id %in% monitors) %>% 
  dplyr::select(id, name, latitude, longitude, elevation) %>% 
  distinct()
  

# extract precipitation data for the three stations
climate <- 
  meteo_pull_monitors(
    monitors = monitors,
    date_min = "1895-01-01",
    date_max = Sys.Date()
  ) %>% 
  left_join(station_data) # merging the station name, the distance to Lilly-Dickey Woods, and the lat long for each weather station


ppt_r2 <- cor(climate$prcp, climate$prcp, use = "complete.obs")^2
#Some graphs to look at how the different station's values compare 
# climate_2 <- climate %>%
#   dplyr::select(id, date, prcp, tobs) %>%
#   distinct()  %>%
#   pivot_wider(names_from = id, values_from = c(prcp, tobs)) %>% 
#   rename(ppt_nashville = prcp_USC00126056,ppt_bloomington = prcp_USC00120784,
#          temp_nashville = tobs_USC00126056,temp_bloomington = tobs_USC00120784) %>% 
#   mutate(year = year(date)) %>%
#   filter(year>=2007)
#   
# ppt_r2 <- cor(climate_2$ppt_nashville, climate_2$ppt_bloomington, use = "complete.obs")^2
# temp_r2 <- cor(climate_2$temp_nashville, climate_2$temp_bloomington, use = "complete.obs")^2
# 
# climate_summary <- lm(climate_2$ppt_nashville ~ 1+climate_2$ppt_bloomington)
# summary(climate_summary)
# # 
# # 
# climate %>%
#   dplyr::select(id, date, prcp, tmax) %>%
#   distinct() %>%
#   pivot_wider(names_from = c(id), values_from = c(prcp, tmax)) %>%
#   ggplot()+
#   geom_point(aes(x = prcp_USC00126056, y = prcp_USC00120784), col = "red")+ #nashville vs bloomington
#   # geom_point(aes(x = prcp_USC00126056, y = prcp_USC00121747), col = "blue")+ #nashville vs columbus
#   # geom_point(aes(x = prcp_USC00121747, y = prcp_USC00120784), col = "green")+ #columbus vs bloominton
#   geom_abline(slope = 1, intercept = 0)
# # 
# # climate %>% 
#   dplyr::select(id, date, prcp, tmax) %>% 
#   distinct() %>% 
#   pivot_wider(id_cols = c(id,date), names_from = c(id), values_from = c(prcp, tmax)) %>% 
#   ggplot()+
#   geom_point(aes(x = tmax_USC00126056, y = tmax_USC00120784), col = "red")+ #nashville vs bloomington
#   geom_point(aes(x = tmax_USC00126056, y = tmax_USC00121747), col = "blue")+ #nashville vs columbus
#   geom_point(aes(x = tmax_USC00121747, y = tmax_USC00120784), col = "green")+#columbus vs bloominton
#   geom_abline(slope = 1, intercept = 0)
#   


# Looking at how many years for each of the stations have the full 12 months of data
# This is not very elegant, but Bloomington has the best recent data, while Columbus has the best historic data, but they are fairly similar
# month_count <- climate %>%
#   mutate(year = year(date), month = month(date)) %>%
#   group_by(id, year, month) %>%
#   summarize(prcp = sum(prcp)) %>%
#   group_by(id, year) %>%
#   summarize(n())
# 
# ggplot(data = month_count)+
#   geom_histogram(aes(x = as.integer(`n()`))) + facet_wrap(~id)
# 
# ggplot(data = filter(month_count, year >2005))+
#   geom_histogram(aes(x = as.integer(`n()`))) + facet_wrap(~id)
# 

# census months for each species to define climate year
census_months <- data.frame(species = c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY"),
                            census_month = c(9,7,7,5,7,5,5)) # 2021 climate data is missing months 6,7

LTREB_full_3 <- LTREB_full_2 %>% 
  left_join(census_months)


# Getting climate data for species with census month (just bloomington)

climate_census_month <- climate %>% 
  rename(station_id = id) %>% #renaming this so it's not confusing with plant id's
  filter(station_id == "USC00120784") %>% 
  mutate(year = year(date), month = month(date), day = day(date)) %>% 
  filter(year >1895) %>% 
  mutate(tmax = tmax/10, tmin = tmin/10, prcp = prcp/10) %>%  # Converting temp. to Celsius, and precip. to mm (from 10ths of units)
  mutate(tmean = case_when(!is.na(tmax) & !is.na(tmin) ~ (tmax + tmin)/2,
                           !is.na(tmax) & is.na(tmin) ~ tmax,
                           is.na(tmax) & !is.na(tmin) ~ tmin)) %>%# there are some NA's in the Max or Min temp at various times, and so that gives us NA's in the mean temp too, but I filter those out
  mutate(tmean_forthorthwaite = case_when(tmean <= 0 ~ 0,
                                          TRUE ~ tmean)) %>% 
  # mutate(PET = thornthwaite(tmean, unique(latitude)), #PET is Potential Evapotranspiration
  #        BAL = prcp - PET) %>% # BAL is Climatic Water Balance
  group_by(latitude, longitude, elevation, station_id, name, year, month) %>% 
  dplyr::summarize(monthly_ppt = sum(prcp, na.rm = TRUE),
                   monthly_tmean = mean(tmean, na.rm = TRUE),
                   monthly_PET = thornthwaite(mean(tmean_forthorthwaite), mean(latitude), na.rm = TRUE), # I have to give the formula only temperatures >0, so I transformed it to that, and there are still some NA's where this is just no data at all.
                   monthly_BAL = monthly_ppt - monthly_PET) 
  # filter(!is.na(monthly_tmean), !is.na(monthly_PET), !is.na(monthly_BAL)) 

# calulate SPEI, we'll use the 12 month spei which is calculated as a 12 month lag from each month
# and the 3 month spei
  # I am removing NA's which means that spei is calculated without those months for some of the time intervals, but I think thiis is better than if we were to drop the months from the database. Dropping the months leads to the column having some skipped months, and the spei's would be calculated with the shifted set of months.
climate_census_month$spei12 <- spei(climate_census_month$monthly_BAL, 12, na.rm = TRUE)$fitted
climate_census_month$spei3 <- spei(climate_census_month$monthly_BAL, 3, na.rm = TRUE)$fitted
climate_census_month$spei24 <- spei(climate_census_month$monthly_BAL, 24, na.rm = TRUE)$fitted
climate_census_month$spei1 <- spei(climate_census_month$monthly_BAL, 1, na.rm = TRUE)$fitted


# Making a dataframe with columns for each species and their census month
spei_census_month_spp <- climate_census_month %>% 
  crossing(census_months) %>% 
  mutate(climate_year = as.numeric(ifelse(month > census_month, year+1, year))) %>%  
  # filter(month == census_month) %>% 
  filter(case_when(year == 2021 & species == "ELVI" | year == 2021 & species  == "ELRI" | year == 2021 & species == "LOAR" ~ month == census_month+1,
                   TRUE ~ month == census_month)) %>% # Here we are taking just the spei starting from the census month, which should cover the climate  year  preceding the census (For 2021, there is no climate data from june, july, so for ELVI, ELRI, and LOAR, I am taking august data for this year)
  mutate(climate_year = case_when(year == 2021 & species == "ELVI" | year == 2021 & species  == "ELRI" | year == 2021 & species == "LOAR" ~ 2021,
                                  TRUE ~ climate_year)) %>% # this is fixinig the label for the climate year for that month fix
  dplyr::select(species, year, climate_year, census_month, spei1, spei3, spei12, spei24)

# Making a dataframe with annual precipitation and temperature for the census year
ppt_temp_census_annual_spp <- climate_census_month %>% 
  crossing(census_months) %>% 
  mutate(climate_year = as.numeric(ifelse(month > census_month, year+1, year))) %>% 
  group_by(species,census_month, climate_year) %>% 
  dplyr::summarize(annual_temp = mean(monthly_tmean, na.rm = T),
                   annual_precip = sum(monthly_ppt, na.rm = T))

# Now merging the climate data to the year for each species
LTREB_full_climate <- LTREB_full_3 %>% 
  left_join(spei_census_month_spp, by = c("species" = "species",  "year_t1" = "climate_year", "census_month" = "census_month")) %>% 
  left_join(ppt_temp_census_annual_spp, by = c("species" = "species",  "year_t1" = "climate_year", "census_month" = "census_month")) 


##############################################################################
####### This is the main dataframe that is used to fit vital rate models  ------------------------------
##############################################################################

LTREB_full <- LTREB_full_climate %>% 
  left_join(LTREB_distances, by = c("species" = "species","pos" = "pos", "plot_fixed" = "plot", "origin_01" = "origin_01", "id" = "id")) %>% 
  mutate(dist_a = case_when(!is.na(dist_a.x) ~ dist_a.x,
                            TRUE ~ dist_a.y),
         dist_b = case_when(!is.na(dist_b.x) ~ dist_b.x,
                            TRUE ~ dist_b.y)) %>% 
  mutate(spei1 = as.numeric(spei1), spei3 = as.numeric(spei3), spei12 = as.numeric(spei12), spei24 = as.numeric(spei24)) %>% # I don't know why but this was giving an error when trying to write the file cause it was saving the column as a list
  dplyr::select(-duplicate, -origin_from_check, -origin_from_distance, -date_status, -date_dist, -contains(".x"), -contains(".y")) # I'm removing some of the extraneous variable. We also have distance data in the new field data that needs to be merged in.
# write_csv(LTREB_full,file = "Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_full.csv")

## Tom is loading this in, bypassing above code
# tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
# LTREB_full <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_full.csv"))
# 
# LTREB_findtypo <- LTREB_full %>% 
#   filter(species == "ELRI", plot_fixed == 109)

# LTREB_plot_density <- LTREB_full %>% 
#   group_by(species, year_t, plot_fixed) %>% 
#   summarize(plants = n(),
#             density = plants/2.517) %>%  # 2.517 is the average distance between measured posts 
#   group_by() %>% 
#   summarize(avg_density = mean(density))

  


