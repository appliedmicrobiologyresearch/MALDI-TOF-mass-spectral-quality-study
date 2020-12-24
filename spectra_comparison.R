# Aline Cuénod 2020
# This script evaluates the spectra acquired for this study. 
# Spectra features were extracted using the following scripts: 
# 'eval_ascii_axima.R' outputting 'read_out_sum_shimadzu_spectra_test1.csv', 
# 'eval_ascii_microflex.R' outputting 'read_out_sum_microflex_spectra_test1.csv'
# 'eval_ascii_microflex_calibration.R' outputting one script per run

# Species identification was performed three databases, and summarised using the following scripts, each outputting a .csv file per analysis, all stored in the same directory
# 'readout_microflex_Biotyper_output.R' for species identification by the MALDI Biotyper database
# 'readout_VitekMS_ID.R' for species identification by the VitekMS database
# 'summarise_marker_based_args.R' for species identification by the PAPMID database

# Besides these files, this script requires the following metadata file
# 'Strains_groups_assigned.csv' indicating which strain belongs to which group

setwd('')

# load packages
library('ggplot2')
library('ggpubr')
library('tidyr')
library('dplyr')
library('tidyverse')
library('stringr')
library('gridExtra')
library('cowplot')

#import which strain belonges to which group
groups<-read.csv2('./04_Strains/Strains_groups_assigned.csv', sep=',')
groups['Numbering_Shipment_1']<-ifelse(nchar(groups$Numbering_Shipment_1)==1, paste0('0', groups$Numbering_Shipment_1), groups$Numbering_Shipment_1)

#import Vitek MS evaluation
eval_shimadzu<-read.csv2('./08_Data/01_Spectra/04_outputs/spectra_test/read_out_sum_shimadzu_spectra_test1.csv')
eval_shimadzu['MALDI']<-'Axima Confidence'

# correct for typos
eval_shimadzu['Method']<-gsub('j2 ', 'j2', eval_shimadzu$Method)

# Remove spectra from the run "test_qty" these was found to be contaminated and repeated
eval_shimadzu<-eval_shimadzu[!duplicated(eval_shimadzu),]
eval_shimadzu<-eval_shimadzu[!eval_shimadzu$Method=="test_qty_conta",]
eval_shimadzu<-eval_shimadzu[!duplicated(eval_shimadzu$spectra_name),]
# Harmonise columnnames
colnames(eval_shimadzu)[colnames(eval_shimadzu) == "strain_number"] <- "strainnumber"

# Import Microflex evaluation
eval_bruker<-read.csv2('./08_Data/01_Spectra/04_outputs/spectra_test/read_out_sum_microflex_spectra_test1.csv')
eval_bruker['MALDI']<-'microflex Biotyper'

# Harmonise columnnames
colnames(eval_bruker)[colnames(eval_bruker) == "strain_number"] <- "strainnumber"

# Read out the method from the file name. Three different samplepreparations were assessed 'smear', 'as25' and 'brukerextraction' and with different quantities. all experiments were repeated three times (biological replicates). 
# Additionaly, each strain was measured at day 1,2,3,4,5,6 days of incubation time. 
eval_bruker['Method']<-gsub('\\.\\d{8}\\-\\d{4}.*', '', eval_bruker$spectra)
eval_bruker['Method']<-ifelse(grepl('qty1', eval_bruker$Method), "quantity_1", 
                              ifelse(grepl('qty-2', eval_bruker$Method), "quantity_2", 
                                     ifelse(grepl('qty-3', eval_bruker$Method), "quantity_3", 
                                            ifelse(grepl('smear$|smear\\-rep', eval_bruker$Method), "smear_1", 
                                                   ifelse(grepl('smear.2', eval_bruker$Method), "smear_2", 
                                                          ifelse(grepl('smear.3', eval_bruker$Method), "smear_3", 
                                                                 ifelse(grepl('^j1.*as25(\\-rep)*$', eval_bruker$Method), "as25_1", 
                                                                        ifelse(grepl('as25.2', eval_bruker$Method), "as25_2", 
                                                                               ifelse(grepl('as25.3', eval_bruker$Method), "as25_3", 
                                                                                      ifelse(grepl('extraction$', eval_bruker$Method), "brukerextraction_1",  
                                                                                             ifelse(grepl('extraction(\\_)*2', eval_bruker$Method), "brukerextraction_2", 
                                                                                                    ifelse(grepl('extraction3', eval_bruker$Method), "brukerextraction_3", 
                                                                                                           ifelse(grepl('j2', eval_bruker$Method), "j2",
                                                                                                                  ifelse(grepl('j3', eval_bruker$Method), "j3", 
                                                                                                                         ifelse(grepl('j4', eval_bruker$Method), "j4",
                                                                                                                                ifelse(grepl('j5', eval_bruker$Method), "j5",
                                                                                                                                       ifelse(grepl('j6', eval_bruker$Method), "j6",NA)))))))))))))))))


# Remove duplicated spectra. These have been exported twice from fid files, once as 'as25' and once as 'day 1'.  
eval_bruker<-eval_bruker[!grepl('j1\\-as25\\.1\\.j1\\-aline\\-as25', eval_bruker$spectra),] 

# Remove '.1' tag from the spectra
eval_bruker['spectra']<- gsub('\\.1$', '', eval_bruker$spectra)

## Remove spectra from the run "test_qty" these was found to be contaminated and repeated
#eval_bruker<-eval_bruker[!eval_bruker$Method=="test_qty_conta",]

# extract the infomation on the run and the position from the file name
eval_bruker['run']<-gsub('(.*)(\\d{8}\\-\\d{4})(.*)', '\\2', eval_bruker$spectra)
eval_bruker['position']<- gsub('(.*)(0\\_[[:alnum:]]{2,3})', '\\2', eval_bruker$spectra)
eval_bruker['position']<- gsub('\\.1$', '', eval_bruker$position)

# Add 'spectra short' column to merge
eval_bruker['spectra_short']<- gsub('(.*\\.)(BRU.*)(\\.0\\_[[:alnum:]]{2,3})*', '\\2', eval_bruker$spectra)

# Remove all empty and duplicated rows
eval_bruker<-eval_bruker[!is.na(eval_bruker$spectra),]
eval_bruker<-eval_bruker[!duplicated(eval_bruker[,c('run', 'position','spectra_short')]),]

# Remove spectra which have no position indicated
eval_bruker<-eval_bruker[!nchar(eval_bruker$position)>5,]

# On day 6, strain 32 is missing, what is labelled as 'day 6 - strain 32' are measurements from day 4 (identical spectra / run / position)
eval_bruker<-eval_bruker[!(eval_bruker$Method == 'j6' & eval_bruker$strainnumber == 32),]

# Add the species identification by the MALDI Biotyper database. The outputs were saved as '.html' files and convrted into .csv files with the script 'readout_microflex_Biotyper_output.R'
# define path to were the outputs are stored
path_brukerrep = './08_Data/01_Spectra/04_outputs/spectra_test/brukerreport_eval/'

# read in all output files
files<-list.files(path_brukerrep, pattern ='.csv')
brukerID_eval<-lapply(paste0(path_brukerrep,files), function(i){
  read.csv(i,stringsAsFactors = FALSE, sep=';')
  })

# merge all brukerreports to one file
brukerID_eval <- brukerID_eval %>% reduce(full_join, by = colnames(brukerID_eval[[1]]))

# add 'bruker' tag to all columnnnames to make them unique and to be able to distinguish these from other database outputs
colnames(brukerID_eval)[colnames(brukerID_eval) %in% c("Score", "Species", "spectra", "species_identified", "genus_identified", "n_species", "n_genera", "n_species_over_2", "n_genera_over_2", "diff", "correct_genus_identified", "correct_species_identified")] <- paste0('brukerDB.', colnames(brukerID_eval)[colnames(brukerID_eval) %in% c("Score", "Species", "spectra", "species_identified", "genus_identified", "n_species", "n_genera", "n_species_over_2", "n_genera_over_2", "diff", "correct_genus_identified", "correct_species_identified")])

# Merge the species ID to the spectra evaluation
eval_bruker<-merge(eval_bruker, brukerID_eval, by=intersect(colnames(eval_bruker), colnames(brukerID_eval)), all.x=T)

# Merge the Bruker and the Shimadzu dataset
eval<-merge(eval_bruker, eval_shimadzu, by=intersect(colnames(eval_bruker), colnames(eval_shimadzu)), all = T)

# Extract and harmonise the format of the 'day' describing how many overnight cultures the strain has been grown before the measurement
eval['day']<-gsub('(.*)(j)(\\d)(.*)', '\\3', eval$Method)

# Use the spectra set acquired at 'as25_3', also for the timeseries experiment. Use 'as25_3' because as25_1 and as25_2 have missing entries
eval['day']<-ifelse(grepl('as25_3.*', eval$day), '1', eval$day)
# set the 'day' variable to NA, if does not match 1-6
eval['day']<-ifelse(!eval$day %in% c('1', '2', '3','4', '5', '6'), NA, eval$day)

# Add a column 'samplepreparation' indicating which samplepreparation was used. (smear_1, smear_2 and smear_3 are summarised to 'smear'), etc
eval['sampleprep']<- ifelse(grepl('smear', eval$Method), 'smear', 
                                   ifelse(grepl('extraction', eval$Method), 'brukerextraktion', 
                                          ifelse(grepl('as25', eval$Method), 'as25', NA)))

# Add a column 'quantity' indicating at which dilution step a sample was measured
eval['quantity']<-ifelse(grepl('quantity', eval$Method), gsub('(.*)(1|4)(\\-\\d{1,5})((\\_|\\.).*)', '\\2\\3', eval$spectra), NA)

# Correct for typos in the dilution series experiment 
eval['quantity']<-gsub('4', '1', eval$quantity)
eval['quantity']<-gsub('4', '1', eval$quantity)
eval['quantity']<-gsub('1-1565', '1-15625', eval$quantity)

# Add identification by a marker based database. Each phylogenetic group has been analysed separately and the database outputs were summarised using the script 'summarise_marker_based_args.R'
marker.ID.Shimadzu<-read.csv('./08_Data/01_Spectra/04_outputs/spectra_test/marker_based_ID_shimadzu_spectra_test_clean.csv', stringsAsFactors = FALSE)
marker.ID.microflex<-read.csv('./08_Data/01_Spectra/04_outputs/spectra_test/marker_based_ID_micoflex_spectra_test_clean.csv', stringsAsFactors = FALSE)

# merge the database outputs for spectra acquired on two devices
marker.ID<-merge(marker.ID.microflex, marker.ID.Shimadzu, by=intersect(colnames(marker.ID.microflex), colnames(marker.ID.Shimadzu)), all = TRUE)

# extract the strainnumber
marker.ID['strainnumber']<-ifelse(nchar(marker.ID$strainnumber)== 1, paste0('0', marker.ID$strainnumber), marker.ID$strainnumber)

# Remove these columns, as these are already comvered in the 'eval' file
marker.ID$species_NGS <-NULL
marker.ID$genus_NGS <- NULL

# Add 'marker' tag to columnnames to make these unique and to be able to distinguish these output from the outputs from other databases
colnames(marker.ID)[colnames(marker.ID) %in% c("match_count","species_identified","Genus","n_species_identified","n_genera_identified","profiles_matched", "X.Subunit.Mass.","correct_genus_identified","correct_species_identified")] <- paste0('markerID.', colnames(marker.ID)[colnames(marker.ID) %in% c("match_count","species_identified","Genus","n_species_identified","n_genera_identified","profiles_matched", "X.Subunit.Mass.","correct_genus_identified","correct_species_identified")])

# Remove the '.1' tag from the filenames
marker.ID['sample']<-gsub('\\.1$', '', marker.ID$sample)

# Some spectra did not yield a marker based ID, because they were excluded by the marker based classifier (e.g. because less than the group specific threshold of marker masses were detected)
# Check how many did not yield a marker basd ID and whether this is within the expected range
check<-eval[eval$spectra %in% setdiff(eval$spectra, marker.ID$sample),]
table(check$MALDI, check$Method) # These are within the expected range of the identification using PAPMID not giving an ID. We also double checked, whether these were in the input

# Add '0' to all one digit strainnumbers
eval['strainnumber']<-ifelse(nchar(eval$strainnumber)== 1, paste0('0', eval$strainnumber), eval$strainnumber)

# Remove duplocates
marker.ID<-marker.ID[!duplicated(marker.ID$sample),]

# Merge the marker based ID to the eva file
eval<-merge(eval, marker.ID, by.x=c('spectra',intersect(colnames(eval), colnames(marker.ID))), by.y = c('sample',intersect(colnames(eval), colnames(marker.ID))), all.x = TRUE)

# The marker based database applied group specific threshold. If less tahn these were detected, the species identification is not reliable. 
# Set all species identifications below the threshold to 'no identification possible
# Use the following thresholds: 
# E. coli / Shigella: 15
# Enterobacter: 20
# Klebsiella: 15
# S.aureus 7: 
# all other: 10

# extract the 'true' species and genus, as determined by NGS
eval['species_NGS']<-gsub('([[:alpha:]]+)(\\_)([[:alpha:]]+)(\\_)(.*)', '\\3', eval$Strain)
eval['genus_NGS']<-gsub('([[:alpha:]]+)(\\_)([[:alpha:]]+)(\\_)(.*)', '\\1', eval$Strain)

# Add a 'theshold' column
eval['marker.theshold']<-ifelse(eval$genus_NGS %in% c('Escherichia', 'Shigella', 'Klebsiella'), 15,
                                ifelse(eval$genus_NGS == 'Enterobacter', 20,
                                       ifelse(eval$genus_NGS == 'Staphylococcus', 7, 10)))

# If a spectra did not yield a species identification by a marker based database, the 'match.count Score' is NA. Set these to 0 
eval['markerID.match_count']<-ifelse(is.na(eval$markerID.match_count), 0, as.numeric(as.character(eval$markerID.match_count)))                                  

# Set all identifications below the phylogenetic group specific threshold to 'no identification possible'
eval["markerID.species_identified"]<-ifelse(as.numeric(eval$markerID.match_count) < eval$marker.theshold, 'no identification possible',eval$markerID.species_identified)
eval["markerID.Genus"]<-ifelse( as.numeric(eval$markerID.match_count) < eval$marker.theshold, 'no identification possible',eval$markerID.species_identified)
eval["markerID.n_species_identified"]<-ifelse( as.numeric(eval$markerID.match_count) < eval$marker.theshold, 0,eval$markerID.n_species_identified)
eval["markerID.n_genera_identified"]<-ifelse( as.numeric(eval$markerID.match_count) < eval$marker.theshold, 0,eval$markerID.n_genera_identified)
eval["markerID.correct_species_identified"]<-ifelse( as.numeric(eval$markerID.match_count) < eval$marker.theshold, 'no ID',eval$markerID.correct_species_identified)
# Remove the threshold column, is not needed any more
eval$marker.theshold<-NULL

# Add Species identification by the VitekMS database (only available for spectra acquired on the Axima Confidence device). Spectra were received as .xlsx / .csv files and summarised using the script 'readout_VitekMS_ID.R'
# Import files
VitekMS_report_dir<-'./08_Data/01_Spectra/04_outputs/spectra_test/testspectra_axima/'
vitekMS_filnames <- list.files(VitekMS_report_dir)
vitekMS_list <- lapply(paste0(VitekMS_report_dir, vitekMS_filnames), read.csv, sep=';')
VitekMS_reports <- do.call(rbind, vitekMS_list)

# Remove '.txt' tag and harmonise spectra name
VitekMS_reports['spectra']<-gsub('\\.txt$', '', VitekMS_reports$samplename_new)
VitekMS_reports['spectra']<-gsub('\\.', '\\_', VitekMS_reports$spectra)

# Remove duplicated spectra. These have been analysed twice by the database, once as 'as25' and once as 'day 1'.
VitekMS_reports<-VitekMS_reports[!grepl('^j1\\-as25\\_', VitekMS_reports$spectra),]

# Harmonise spectra name nomenclature
eval$spectra<-gsub('\\-', '\\_', eval$spectra)
VitekMS_reports$spectra<-gsub('\\-', '\\_', VitekMS_reports$spectra)

# Check whether there are any differences
# setdiff(eval[eval$MALDI == 'Axima Confidence', 'spectra'], VitekMS_reports$spectra) # There are none

# Set these columns to zero as the are not needed for further analysis or already corered by the 'eval' file
VitekMS_reports$X<-NULL
VitekMS_reports$samplename_new<-NULL
VitekMS_reports$files<-NULL
VitekMS_reports$lab<-NULL
VitekMS_reports$position <- NULL

# Remove 'Referenz' spectra, these are the DHalpha E.coli calibration spectra
VitekMS_reports<-VitekMS_reports[!grepl('Referenz', VitekMS_reports$strainnumber),]

# Add '0' to all one digit strainnumber
VitekMS_reports$strainnumber<-ifelse(nchar(VitekMS_reports$strainnumber) < 2, paste0('0', VitekMS_reports$strainnumber), VitekMS_reports$strainnumber)

# Add 'noID' as option to whether the species was correctly identified or not
VitekMS_reports$correct_species_identified<-ifelse(VitekMS_reports$identifType %in% c('NoIdentification','NotEnoughPeaks','TooNoisy'), 'noID',
                                                   ifelse(VitekMS_reports$correct_species_identified == 0, FALSE,
                                                          ifelse(VitekMS_reports$correct_species_identified == 1, TRUE, NA)))

# Add 'VitekMS' all columnnames to make them unique and to be able to distinguish these from other database outputs
colnames(VitekMS_reports)[!colnames(VitekMS_reports) %in% colnames(eval)]<-paste0('VitekMS.DB', colnames(VitekMS_reports)[!colnames(VitekMS_reports) %in% colnames(eval)])

# Merge these species identifications to the 'eval' file
eval<-merge(eval, VitekMS_reports, by =c("strainnumber", "spectra", "Strain"), all.x = TRUE)

# Merge the columns which are present but not complete for any of the two datasets 
eval['species_NGS']<- ifelse(!is.na(eval$species_NGS.x), eval$species_NGS.x, eval$species_NGS.y)
eval['genus_NGS']<- ifelse(!is.na(eval$genus_NGS.x), eval$genus_NGS.x, eval$genus_NGS.x)
eval['include']<- ifelse(!is.na(eval$include.x), eval$include.x, eval$include.y)

# Remove redundant columns
eval$species_NGS.x<-NULL
eval$species_NGS.y<-NULL
eval$genus_NGS.x<-NULL
eval$genus_NGS.y<-NULL
eval$include.x<-NULL
eval$include.y<-NULL
eval$method <- NULL

# Set all species identifications by the MALDI Biotyper database, which yielded a log Score < 1.7 to noID
eval["brukerDB.correct_species_identified"]<-ifelse(eval$brukerDB.Score<1.7, 'no ID', eval$brukerDB.correct_species_identified)
eval["brukerDB.correct_genus_identified"]<-ifelse(eval$brukerDB.Score<1.7, 'no ID', eval$brukerDB.correct_genus_identified)


# Add which strains are assigned to which group
eval<-merge(groups, eval, by.x=c('Strain',"Numbering_Shipment_1"), by.y = c('Strain','strainnumber'), all.y = TRUE)
eval<-eval[!is.na(eval$Strain),]

# Add 'single species ID' column, indicating whether the cirect species was uniquely identified using a marker based approach
eval['markerID.correct_single_species']<-ifelse(eval$markerID.correct_species_identified ==TRUE & eval$markerID.n_species_identified =='1', TRUE, FALSE)
eval<-eval[!is.na(eval$spectra),]
eval<-eval[!duplicated(eval$spectra),]

# This spectra had a weird samplename and the measurement was repeated. Remove it
eval<-eval[!eval$spectra == 'j1-bruker.j1-aline-bruker-extraction2.aline - bruker extraction j1.20200424-1650.BRU-20',]

# Check how many strains are there and how many are missing per sample spectra set
table(eval$Method) 
table(eval$sampleprep) 

# Split up the 'Streptococcus' group into in 'viridans streptococci' and and 'non-viridans streptocci' as 'viridans streptococci' are of special interest and form a challenge for species identification by MALDI-TOF MS
eval['Group']<-as.factor(ifelse(eval$genus_NGS == 'Streptococcus' & grepl('infantis|pseudopneumoniae|pneumoniae|gordonii', as.character(eval$Strain)), 'Viridans streptococci', 
                                ifelse(eval$genus_NGS == 'Streptococcus' & grepl('equinus|gallolyticus|lutetiensis|dysgalactiae', as.character(eval$Strain)), 'Other Streptococci', as.character(eval$Group))))


# Add 'fraction of high peaks' as endpoint, indicatingwhat fraction of detected peaks had m/z values > 10'000. This could be a measure for how 'balanced' a spectrum is
eval['frac_high_peaks']<-eval$n_high_peaks / eval$n_peaks

# There has probably been as mixup betwwen strain Nr 34 and 35 in 'smear_1' and 'as_1' (indicated by the wrong species identifications). Exclude these
eval <- eval[!(eval$Numbering_Shipment_1 %in% c(34,35) & !is.na(eval$sampleprep) & eval$Method %in% c('smear_1', 'as25_1')),]

# First compare the spectra features and species identifications yielded for spectra acquired under different measurement conditions, eg. sample preparation protocol and days of incubation 
# the dilution series included only two strains per group and the species identification by the different databases are not analysed separately
# Evaluate the species identifications, per measurement condition, MALDI and database and summarise this to %
# Do this separately for all experiments as there us an overlap between as25 and the timeseries experiment
eval2_sampleprep <- eval[!is.na(eval$sampleprep),] %>%
  group_by(MALDI, sampleprep) %>%
  summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            
            prop.correct.marker.single = sum(grepl('TRUE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),             
            prop.wrong.marker.single = sum(grepl('FALSE', markerID.correct_single_species) & grepl('FALSE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.correct.marker.multi = sum(grepl('FALSE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.noID.marker = sum(grepl('FALSE', markerID.correct_single_species) & grepl('no ID', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            
            prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified))) 

eval2_timeseries <- eval[!is.na(eval$day),] %>%
  group_by(MALDI, day) %>%
  summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            
            prop.correct.marker.single = sum(grepl('TRUE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),             
            prop.wrong.marker.single = sum(grepl('FALSE', markerID.correct_single_species) & grepl('FALSE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.correct.marker.multi = sum(grepl('FALSE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.noID.marker = sum(grepl('FALSE', markerID.correct_single_species) & grepl('no ID', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            
            prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified))) 

# Merge the two sample_ID evaluation files
eval2<-merge(eval2_sampleprep, eval2_timeseries, by = intersect(colnames(eval2_sampleprep),colnames(eval2_timeseries)), all = T)

# Add them to the 'eval' file
eval<-merge(eval, eval2, by = intersect(colnames(eval), colnames(eval2)), all.x  = T)
# Gather to the long format 
eval_prop_long<-gather(eval2, "classification_eval_prop", "value", -c(MALDI, sampleprep, day))

# Gather all Endpoints to one column in order to plot all at once using facet_wrap
eval_long<-gather(eval[,c('max_mass_ribo', 'spectra',"n_peaks", "highest_peaks", "mass_90", "n_high_peaks", "frac_high_peaks","frac_peaks_repr", "n_ribos_detected", "n_predicted_su", "rel_amount_su_detected", "mean_dist", "mean_dist_ppm", "mean_intensity","median_intensity", "total_intensity_peaks","brukerDB.Score", "brukerDB.n_species", "brukerDB.n_genera", "brukerDB.n_species_over_2", "brukerDB.n_genera_over_2", "brukerDB.diff", "brukerDB.correct_genus_identified", "brukerDB.correct_species_identified",  "markerID.match_count",  "markerID.n_genera_identified",  "markerID.correct_genus_identified", "markerID.correct_species_identified", "markerID.correct_single_species", "VitekMS.DBproba_rank_1","VitekMS.DBscore_rank_1",  "VitekMS.DBdiff", "brukerDB.diff")], "Endpoint", "value", -spectra)
eval_long<-merge(eval[,c("spectra","Strain", "Numbering_Shipment_1", "Shape", "Gram", "aerobic", "Family", "Group", "spectra", "species_NGS", "genus_NGS", "MALDI", "Method", "run", "position", "day", "sampleprep", "quantity", "Run")], eval_long, by='spectra')

# Convert TRUE/FALSE to numerical binaries. This applies for the columns indicating wether the correct species/genus was identified
eval_long['value']<-ifelse(eval_long$value == 'TRUE', 1,
                           ifelse(eval_long$value %in% c('FALSE', 'no ID'), 0, eval_long$value))
# Convert all features to a numeric value
eval_long["value"]<-as.numeric(as.character(eval_long$value)) 

# Add 'DB' column indicating which databse was used
eval_prop_long['DB']<-ifelse(grepl('bruker', eval_prop_long$classification_eval_prop), 'bruker', 
                             ifelse(grepl('marker', eval_prop_long$classification_eval_prop), 'marker', 
                                    ifelse(grepl('VitekMS', eval_prop_long$classification_eval_prop), 'VitekMS', NA)))

# rename calssification eval prop column
eval_prop_long['eval']<-eval_prop_long$classification_eval_prop
eval_prop_long$classification_eval_prop<-NULL

# Convert sample preparations as factors
eval_prop_long$sampleprep <- factor(eval_prop_long$sampleprep, levels=c("smear", "as25", "brukerextraktion"))  

# Log transform the intensity data
# Store the untransformed data in a separate column
eval_long['value_no_log']<-eval_long$value

# In order to be able to plot the log transformed data, logtransform all intensities + 1. (log10 of 1 is 0, these values will be displayd as 0)
eval_long['value']<-ifelse(eval_long$Endpoint == 'total_intensity_peaks', log10(as.numeric(eval_long$value) + 1), as.numeric(eval_long$value))

# If no peaks were found at all, or if no ribsosmal markers were detected, many of the Ednpoints are 'NA'. These can be set to 0. (This is however not true for the measuremnet error, if no marker mass was detected, these are left as 'NA)
eval_long['value']<-ifelse(eval_long$Endpoint %in% c('total_intensity_peaks', 'mean_dist_ppm', 'frac_high_peaks', 'frac_peaks_repr', 'n_ribos_detected','median_intensity', 'total_intensity_peaks') & is.na(eval_long$value), 0, as.numeric(eval_long$value))

# subset data to plot only including the relevant Endpoints
plot_data<-eval_long[eval_long$Endpoint %in% c("n_ribos_detected","median_intensity","total_intensity_peaks", "frac_peaks_repr", "rel_amount_su_detected"),]

# Before merging the data to the evaluation of the species identification, rename the 'DB' to 'Endpoint'
names(eval_prop_long)[names(eval_prop_long) == 'DB'] <- 'Endpoint'
eval_prop_long$classification_eval_prop<-NULL

# Add the species evaluation to the data to plot
plot_data<-merge(plot_data, eval_prop_long, by = intersect(colnames(eval_prop_long), colnames(plot_data)), all = TRUE)

# write out
write.csv(plot_data, './08_Data/01_Spectra/04_outputs/plot_data_long.csv')

# Add an 'Experiment' column for easier subsetting into the respective datasets
plot_data['experiment']<-ifelse(!is.na(plot_data$sampleprep), 'Sample prepreration', 
                                ifelse(!is.na(plot_data$day), 'Timeseries', 
                                       ifelse(!is.na(plot_data$quantity), 'Dilution Series', NA)))

# In order to be able to include the spectra which have been acquired under the 'as25' method at day 1 in both, the 'Samplepreparation' dataset AND the 'timeseries dataset, duplicate these rows.
# As in the datasets as25_1 and # as25_2 there are spectra missing, use 'as25_3' for this task
# Subset the rows to duplicate 
add<-plot_data[plot_data$Method == 'as25_3' & !is.na(plot_data$sampleprep),]
# Change the 'Experiment'
add['experiment']<-'Timeseries'

# Add these back to the plot data
plot_data<-merge(plot_data, add, by = intersect(colnames(plot_data), colnames(add)), all = T)

# For comparison between datasets of the same strains acquired under different conditions, it is valid to compare the number of ribsosmal subunits detected
plot_data1<-plot_data[plot_data$Endpoint != 'rel_amount_su_detected',]

# For comparisona between the different phylogenetic groups, the relative amount of ribosomal subunit detected has to be considered, as the number of ribsomal marker masses with a predicted mass in the MALDI mass range may vary between the groups
plot_data2<-plot_data[plot_data$Endpoint != 'n_ribos_detected',]

# First compare between the different measuring conditions, and comparing the not normalise values
plot_data<-plot_data1

# Rename the sample prepeeration protocols from 'brukerextraction to 'simple protein extraction', 'as25' to '25% FA overlay' and 'smear' to 'direct smear'
plot_data$sampleprep<-gsub('as25','25% FA overlay',plot_data$sampleprep)
plot_data$sampleprep<-gsub('brukerextraktion','simple protein extraction', plot_data$sampleprep)
plot_data$sampleprep<-gsub('smear','direct smear', plot_data$sampleprep)

# Add the string 'Day' in from of the number indicating after how many das of incubation a colony was measured
plot_data$day<-ifelse(!is.na(plot_data$day), paste('Day', plot_data$day), NA)

# Rename the dilution factors
plot_data$quantity<-gsub('\\-','\\:', plot_data$quantity)

# Make two separate plots for spectra acquired on the Axima Confidence and on the microflex Biotyper device
# In these overall plot include exclusively marker based Identifications, as neither the VitekMS database nor the MALDI Biotyper database cover all species included in this study
plot_data<-plot_data[!plot_data$Endpoint %in% c('VitekMS','bruker'),]

# Subset data which has been acquired on the microflex Biotyper
plot_data_bruker<-plot_data[plot_data$MALDI == 'microflex Biotyper',]
# Sort the species identification evaluation according to 'correctness'. This will then set the order in which these are displayd in a stacked bar plot
plot_data_bruker$eval<- factor(plot_data_bruker$eval, levels = c( "prop.wrong.marker.single", 
                                                                               "prop.noID.marker", 
                                                                               "prop.correct.marker.multi", 
                                                                               "prop.correct.marker.single"))

# Define how to call the Endpoints in the plot
labels_endpoints<-c("Number of detected ribosomal marker peaks", "Median relative intensity of ribosomal marker peaks", "Sum of the intensity of all detected peaks", "Fraction of reproducibly detected peaks", "Evaluation of the PAPMID species identification")
# Introduce linebreaks for long names
labels_endpoints<-str_wrap(labels_endpoints, width = 15)

# Rename the Endpoints
plot_data_bruker$Endpoint<- factor(plot_data_bruker$Endpoint, levels = c("n_ribos_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr", "marker"), 
                                   labels = labels_endpoints)

# Strain 17 is missing from microflex brukerextraction 2, axima brukerextraction 2, strain 20 is missing from microflex brukerextraction_3 and for strain 46 there are only 3 spectra in axima smear 3. As strain 34 and 35 were mixed up in as_1 and smear_1, these were also removed
# In order to perform paired wilcoxon rank test, remove strain 17 from microflex smear2, as25_2, strain 20 from as25_3 and smear_3. Also remove strain 34 and 35 from bruker_1, as the are already missing from as_1 and smear_1
plot_data_bruker<-plot_data_bruker[!(plot_data_bruker$Numbering_Shipment_1 == '17' & plot_data_bruker$Method %in% c("smear_2","as25_2") |
                                      plot_data_bruker$Numbering_Shipment_1 == '20' & plot_data_bruker$Method %in% c("smear_3","as25_3") | 
                                       plot_data_bruker$Numbering_Shipment_1 %in%  c('34','35') & plot_data_bruker$Method  == "brukerextraction_1"),]

# Extract the median and Interquartile ranges (IQR) to cite them in the running text of the manuscript
# Per Endpoint and sampleprep
check_sampleprep<-plot_data_bruker[plot_data_bruker$experiment == "Sample prepreration",] %>% 
  group_by(sampleprep, Endpoint) %>%
  summarize(IQR = quantile(value_no_log, probs = c(0.25, 0.5, 0.75), na.rm = T))

# Per Endpoint and 'day' (incubation time)
check_day<-plot_data_bruker[plot_data_bruker$experiment == "Timeseries",]%>% 
  group_by(day, Endpoint) %>%
  summarize(IQR = quantile(value_no_log, probs = c(0.25, 0.5, 0.75), na.rm = T))

# Create fake points to have same ylim everywhere (these wont be seen, as alpha == 0)
ylim_points<-data.frame(sampleprep = '25% FA overlay', day = 'Day 1', quantity = '1:5', Endpoint = rep(c("n_ribos_detected", "frac_peaks_repr", "total_intensity_peaks", "median_intensity", "marker"),2), value = c(38,1,8.7,8,1,0,0,0,0,0), value1 = c(38,1,8.7,8,1,0,0,0,0,0))
ylim_points$Endpoint<-factor(levels = c("n_ribos_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr", "marker"), 
                             ylim_points$Endpoint, labels = labels_endpoints)

# Subset for spectra with barying sample preparations
plot_data_bruker_sample<-plot_data_bruker[plot_data_bruker$experiment == 'Sample prepreration',]
plot_data_bruker_sample<-plot_data_bruker_sample[!(is.na(plot_data_bruker_sample$sampleprep) | is.na(plot_data_bruker_sample$Endpoint)),]

# Define the order in which the different protocols are displayed
plot_data_bruker_sample['sampleprep']<-factor(plot_data_bruker_sample$sampleprep, levels = c('direct smear', '25% FA overlay', 'simple protein extraction'))


# For the species identification evaluation, use the following colorcode throughout all figures: 
#"prop.bruker.correct.single.over.two"    #009E73
#"prop.correct.marker.single"             #009E73
#"prop.correct.VitekMS.single"            #009E73

#"prop.bruker.correct.multi.over.two"     #56B4E9
#"prop.correct.marker.multi"              #56B4E9
#"prop.correct.VitekMS.lowdiscriminatory" #56B4E9

#"prop.bruker.correct.under.two"          #0072B2

#"prop.noID.bruker"                       #999999     
#"prop.noID.marker.single"                #999999
#"prop.noID.VitekMS"                      #999999

#"prop.wong.bruker.under.two"             #E69F00
#"prop.wrong.VitekMS.lowdiscriminatory"   #E69F00

#"prop.wong.bruker.over.two"              #D55E00
#"prop.wrong.marker.single"               #D55E00
#"prop.wrong.VitekMS.single"              #D55E00  



# Plot
f1 <- ggplot(plot_data_bruker_sample, aes(x=sampleprep, y=as.numeric(as.character(value)))) + facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE, method.args=list(p.adjust.method = "BH"),
                                           ref.group = 'direct smear',
                                           aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1)
f2 <- f1+geom_boxplot(data = plot_data_bruker_sample[grepl('Number.*marker.*',plot_data_bruker_sample$Endpoint),]) + geom_point(data=ylim_points,x='direct smear', alpha = 0)
f2b <- f2+geom_boxplot(data = plot_data_bruker_sample[grepl('Median.*intensity.*',plot_data_bruker_sample$Endpoint),])
f3 <- f2b+geom_boxplot(data = plot_data_bruker_sample[grepl('Sum.*intensity.*',plot_data_bruker_sample$Endpoint),])
f4 <- f3+geom_boxplot(data = plot_data_bruker_sample[grepl('Fraction.*reproducibly.*',plot_data_bruker_sample$Endpoint),])
f5 <- f4+geom_col(data = plot_data_bruker_sample[grepl('.*PAPMID.*',plot_data_bruker_sample$Endpoint),], aes(fill=eval), width = 0.7)
f7_samplprep_all <- f5+scale_fill_manual('Identification',
                          breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker",
                                     "prop.wrong.marker.single" ),
                          values = c("#009E73", "#0072B2", "#999999","#D55E00"), 
                          labels = c("Correct single species identification", "Correct multi species identification", "No identification possibe", "Wrong identification")) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 ylab('') +xlab('') 
f7_samplprep_all 

# Extract the legend. Returns a gtable
leg.marker <- get_legend(f7_samplprep_all+ theme(legend.position="bottom"))

pdf('./marker_ID_legend.pdf', width= 12, height = 1)
as_ggplot(leg.marker)
dev.off()

# Subset for spectra with varying incubation times
plot_data_bruker_timeseries<-plot_data_bruker[plot_data_bruker$experiment == 'Timeseries',]

# In order to perform paired wilcoxon rank test, remove strain 20 from day 2,3,4,5,6, #strain 32 is missing from day 6, remove it from all days
plot_data_bruker_timeseries<-as.data.frame(plot_data_bruker_timeseries[!(is.na(plot_data_bruker_timeseries$day) | is.na(plot_data_bruker_timeseries$Endpoint)),])
plot_data_bruker_timeseries <- plot_data_bruker_timeseries %>% filter((as.character(Numbering_Shipment_1) != '20')|is.na(Numbering_Shipment_1) ) 
plot_data_bruker_timeseries <- plot_data_bruker_timeseries %>% filter((as.character(Numbering_Shipment_1) != '32')|is.na(Numbering_Shipment_1) ) 

# Define the order in which the 'day' variable is displayes
plot_data_bruker_timeseries['day']<-factor(plot_data_bruker_timeseries$day, levels = c('Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6'))

f1 <- ggplot(plot_data_bruker_timeseries, aes(x=day, y=as.numeric(as.character(value))))+facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                           ref.group = 'Day 1',
                                           aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1)
f2 <- f1+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Number.*marker.*',plot_data_bruker_timeseries$Endpoint),]) + geom_point(data=ylim_points,x='Day 1', alpha = 0)
f2b <- f2+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Median.*intensity.*',plot_data_bruker_timeseries$Endpoint),])
f3 <- f2b+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Sum.*intensity.*',plot_data_bruker_timeseries$Endpoint),])
f4 <- f3+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Fraction.*reproducibly.*',plot_data_bruker_timeseries$Endpoint),])
f5 <- f4+geom_col(data = plot_data_bruker_timeseries[grepl('.*PAPMID.*',plot_data_bruker_timeseries$Endpoint),], aes(fill=eval), width = 0.7)
f7_timeseries_all <- f5+scale_fill_manual('Identification',
                                          breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker",
                                                     "prop.wrong.marker.single" ),
                                          values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 ylab('') +xlab('')
f7_timeseries_all

# combine the plots for varying samplepreparation and vay^ryig incubation time to one
combined_plots<-cowplot::plot_grid(f7_samplprep_all + theme(legend.position="none") + ggtitle('Sample Preparation') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14)),
                                   f7_timeseries_all + theme(legend.position="none") + ggtitle('Time Series') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14)),
                                   rel_widths=c(0.125, 0.2), ncol = 2, align = "h",axis = "b")


# Export the plot
pdf('./bruker_overall_nobrukerID_no_quantity.pdf', height = 16, width = 11.25)
plot_grid(combined_plots, leg.marker, ncol = 1, rel_heights = c(0.9, 0.1))
dev.off()

# Subset for spectra from the dilution experiment
plot_data_bruker_quantity<-plot_data_bruker[plot_data_bruker$experiment == 'Dilution Series' & !(plot_data_bruker$quantity %in% c('1:3125', '1:15625')),]
plot_data_bruker_quantity<-plot_data_bruker_quantity[!(is.na(plot_data_bruker_quantity$quantity) | is.na(plot_data_bruker_quantity$Endpoint)),]


# Add spectra acquired under the 'as25' method and of the same strains for direct comparison
plot_data_bruker_to_add<-plot_data_bruker_sample[plot_data_bruker_sample$Numbering_Shipment_1 %in% c('07', '08', '09', '10', '17','18', '21','23','26','27', '33', '34','39', '40', '45', '46') &
                                                   plot_data_bruker_sample$sampleprep == '25% FA overlay',]
# Set '25% FA overlay' as 'quantity' variable
plot_data_bruker_to_add$quantity<-ifelse(is.na(plot_data_bruker_to_add$quantity), '25% FA overlay', plot_data_bruker_to_add$quantity)
plot_data_bruker_quantity1<-merge(plot_data_bruker_quantity, plot_data_bruker_to_add, by=intersect(colnames(plot_data_bruker_to_add), colnames(plot_data_bruker_quantity)), all = T)

# Remove 17 in quantity 2 and 34 in quantity 2 in order to perform paired paired wilcoxon rank tests
plot_data_bruker_quantity1<-plot_data_bruker_quantity1[!(plot_data_bruker_quantity1$Numbering_Shipment_1 == '17' & plot_data_bruker_quantity1$Method == 'quantity_2' |
                                                         plot_data_bruker_quantity1$Numbering_Shipment_1 == '34' & plot_data_bruker_quantity1$Method == 'quantity_1'),]

plot_data_bruker_quantity1<-plot_data_bruker_quantity1[!is.na(plot_data_bruker_quantity1$Endpoint),]

# Define order in which the dilution factors should be displayed
plot_data_bruker_quantity1$quantity<-factor(plot_data_bruker_quantity1$quantity,levels = c('25% FA overlay','1:5','1:25','1:125','1:625'))

# Summarise data to median and IQR, per dilution factor and Endpoint
check_quantity<- plot_data_bruker_quantity1[!grepl('Evaluation', plot_data_bruker_quantity1$Endpoint),]%>% 
  group_by(quantity, Endpoint) %>%
  summarize(IQR = quantile(value_no_log, probs = c(0.25, 0.5, 0.75), na.rm = T))

# Plot
plot_data_bruker_quantity_all_groups <- ggplot(plot_data_bruker_quantity1[!grepl('Evaluation', plot_data_bruker_quantity1$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ ., scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                                                                          ref.group = '25% FA overlay',
                                                                                          aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill = quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14), legend.position = 'none') + 
  scale_fill_brewer(palette = "Blues", direction=-1)

# Export the plot
pdf('./bruker_quantity_all_groups.pdf', height = 12, width = 5)
plot_data_bruker_quantity_all_groups
dev.off()

# Duliting the sample might have taxon specific effects. Plot per phylogenetic group
# For this the Streptococco can be summarised, as no species identification is considered here
plot_data_bruker_quantity1$Group<-ifelse(grepl('treptoc', plot_data_bruker_quantity1$Group), "Streptococcus", as.character(plot_data_bruker_quantity1$Group))

# Define how to call groups, add how many strains are included
group_labels_strep_together <- c("Enterobacteriaceae (18 strains)", "Listeria (2 strains)", "Burkholderia (4 strains)", "Bordetella (3 strains)",
                 "Streptococcus (9 strains)", "Staphylococcus (3 strains)", "Actinobacteria (5 strains)", "Gram negative Anaerobes (3 strains)")
# Introduce line breaks if the group names are very long
group_labels_strep_together <- str_wrap(group_labels_strep_together, width = 15)

# for the quantity experiment, only 2 strains per group have been measured
group_labels_strep_together_qunatity <- c("Enterobacteriaceae\n(2 strains)", "Listeria\n(2 strains)", "Burkholderia\n(2 strains)", "Bordetella\n(2 strains)",
                                 "Streptococcus\n(2 strains)", "Staphylococcus\n(2 strains)", "Actinobacteria\n(2 strains)", "Gram negative Anaerobes\n(2 strains)")

group_labels_strep_together_qunatity <- str_wrap(group_labels_strep_together_qunatity, width = 15)

plot_data_bruker_quantity1$Group<-factor(plot_data_bruker_quantity1$Group, levels = c("Enterobacteriaceae", "Listeria", "Burkholderia", "Bordetella",
                                                                                      "Streptococcus", "Staphylococcus", "Actinobacteria", "Anaerob_Gram_negative"),
                                                                            labels = group_labels_strep_together_qunatity)
# Plot dilutuion series per group
plot_data_bruker_quantity_per_group <- ggplot(plot_data_bruker_quantity1[!grepl('Evaluation', plot_data_bruker_quantity1$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ Group, scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                     ref.group = '25% FA overlay',
                     label = "p.signif", vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill = quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14), legend.position = 'none') + 
  scale_fill_brewer(palette = "Blues", direction=-1)

 # Export the plot
pdf('./bruker_quantity_per_groups.pdf', height = 12, width = 16)
plot_data_bruker_quantity_per_group
dev.off()

# Plot dilutuion Enterobacteriaceae
plot_data_bruker_quantity1_entero<-plot_data_bruker_quantity1[plot_data_bruker_quantity1$Group == "Enterobacteriaceae\n(2 strains)" & !grepl('Evaluation', plot_data_bruker_quantity1$Endpoint),]
plot_data_bruker_quantity1_entero<-plot_data_bruker_quantity1_entero[!is.na(plot_data_bruker_quantity1_entero$quantity) & !is.na(plot_data_bruker_quantity1_entero$Endpoint),]
plot_data_bruker_quantity_enterobacteriaceae <- ggplot(plot_data_bruker_quantity1_entero[!grepl('Evaluation', plot_data_bruker_quantity1_entero$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ ., scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                     ref.group = '25% FA overlay',
                     label = "p.signif", vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill= quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=15), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = 'none') +
  scale_fill_brewer(palette = "Blues", direction=-1)
# Export the plot
pdf('./bruker_quantity_per_groups_entero.pdf', height = 9.6, width = 3.5)
plot_data_bruker_quantity_enterobacteriaceae
dev.off()

# plot dilutions burkholderia
plot_data_bruker_quantity1_burko<-plot_data_bruker_quantity1[plot_data_bruker_quantity1$Group == "Burkholderia (2\nstrains)" & !grepl('Evaluation', plot_data_bruker_quantity1$Endpoint),]
plot_data_bruker_quantity1_burko<-plot_data_bruker_quantity1_burko[!is.na(plot_data_bruker_quantity1_burko$quantity) & !is.na(plot_data_bruker_quantity1_burko$Endpoint),]
plot_data_bruker_quantity_burko <- ggplot(plot_data_bruker_quantity1_burko[!grepl('Evaluation', plot_data_bruker_quantity1_burko$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ ., scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                     ref.group = '25% FA overlay',
                     label = "p.signif", vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill= quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=15), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = 'none') +
  scale_fill_brewer(palette = "Blues", direction=-1)

# Export the plot
pdf('./bruker_quantity_per_groups_burko.pdf', height = 9.6, width = 3.5)
plot_data_bruker_quantity_burko
dev.off()

# plot dilutions viridans strep
plot_data_bruker_quantity1_strep<-plot_data_bruker_quantity1[plot_data_bruker_quantity1$Group == "Streptococcus\n(2 strains)" & !grepl('Evaluation', plot_data_bruker_quantity1$Endpoint),]
plot_data_bruker_quantity1_strep<-plot_data_bruker_quantity1_strep[!is.na(plot_data_bruker_quantity1_strep$quantity) & !is.na(plot_data_bruker_quantity1_strep$Endpoint),]
plot_data_bruker_quantity_strep <- ggplot(plot_data_bruker_quantity1_strep[!grepl('Evaluation', plot_data_bruker_quantity1_strep$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ ., scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                     ref.group = '25% FA overlay',
                     label = "p.signif", vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill= quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=15), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = 'none') +
  scale_fill_brewer(palette = "Blues", direction=-1)

# Export the plot
pdf('./bruker_quantity_per_groups_strep.pdf', height = 9.6, width = 3.5)
plot_data_bruker_quantity_strep
dev.off()


# Plot the same series of plots for spectra acquired on the Axima Confidence device
# Subset the data
plot_data_axima<-plot_data[plot_data == 'Axima Confidence',]

# Define the stacking order of the species evaluation
plot_data_axima$eval<- factor(plot_data_axima$eval, levels = c( "prop.wrong.marker.single", 
                                                                  "prop.noID.marker", 
                                                                  "prop.correct.marker.multi", 
                                                                  "prop.correct.marker.single"))

plot_data_axima$Endpoint<- factor(plot_data_axima$Endpoint, levels = c("n_ribos_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr", "marker"),
                                  labels = labels_endpoints)

# Create fake points to have same ylim everywhere (these wont be seen, as alpha == 0)
ylim_points<-data.frame(sampleprep = '25% FA overlay', day = 'Day 1', quantity = '1:5', Endpoint = rep(c("n_ribos_detected", "frac_peaks_repr", "total_intensity_peaks", "median_intensity", "marker"),2), value = c(38,1,8.75,60,1,0,0,0,0,0), value1 = c(38,1,8.75,8,1,0,0,0,0,0))
ylim_points$Endpoint<-factor(levels = c("n_ribos_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr", "marker"), 
                             ylim_points$Endpoint, labels = labels_endpoints)

# Strain 17 is missing from microflex brukerextraction 2, axima brukerextraction 2, strain 20 is missing from microflex brukerextraction_3 and for strain 46 there are only 3 spectra in axima smear 3. As there is strong indication that strains 34 and 35 were not of the correct species in smear_1 and as_1, these were excluded. 
# In order to perform paired wilcoxon rank test, remove strain 17 from microflex smear2, as25_2, strain 20 from as25_3 and smear_3 and strain 34 and 35 from brukerextraction_1
plot_data_axima<-plot_data_axima[!(plot_data_axima$Numbering_Shipment_1 == '17' & plot_data_axima$Method %in% c("smear_2","as25_2") |
                                     plot_data_axima$Numbering_Shipment_1 == '46' & plot_data_axima$Method %in% c("as25_2", "brukerextraction_2") & grepl('4$',plot_data_axima$spectra)|
                                     plot_data_axima$Numbering_Shipment_1 %in%  c('34','35') & plot_data_axima$Method  == "brukerextraction_1"),]

# Subset for the set of spectra which was acquired under varying sample preparation protocols
plot_data_axima_sample<-plot_data_axima[plot_data_axima$experiment == 'Sample prepreration',]
plot_data_axima_sample<-plot_data_axima_sample[!(is.na(plot_data_axima_sample$sampleprep) | is.na(plot_data_axima_sample$Endpoint)),]
# Define the order in which the sample protocols are displayed
plot_data_axima_sample['sampleprep']<-factor(plot_data_axima_sample$sampleprep, levels = c('direct smear', '25% FA overlay', 'simple protein extraction'))
f1 <- ggplot(plot_data_axima_sample, aes(x=sampleprep, y=as.numeric(as.character(value)))) + facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE, method.args=list(p.adjust.method = "BH"),
                                            ref.group = 'direct smear',
                                            aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1)
f2 <- f1+geom_boxplot(data = plot_data_axima_sample[grepl('Number.*marker.*',plot_data_axima_sample$Endpoint),]) + geom_point(data=ylim_points,x='direct smear', alpha = 0)
f2b <- f2+geom_boxplot(data = plot_data_axima_sample[grepl('Median.*intensity.*',plot_data_axima_sample$Endpoint),])
f3 <- f2b+geom_boxplot(data = plot_data_axima_sample[grepl('Sum.*intensity.*',plot_data_axima_sample$Endpoint),])
f4 <- f3+geom_boxplot(data = plot_data_axima_sample[grepl('Fraction.*reproducibly.*',plot_data_axima_sample$Endpoint),])
f5 <- f4+geom_col(data = plot_data_axima_sample[grepl('.*PAPMID.*',plot_data_axima_sample$Endpoint),], aes(fill=eval), width = 0.7)
f7_samplprep_all <- f5+scale_fill_manual('Identification',
                                         breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker",
                                                    "prop.wrong.marker.single" ),
                                         values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')
f7_samplprep_all

# Subset the spectra which were acquired after different incubation times
plot_data_axima_timeseries<-plot_data_axima[plot_data_axima$experiment == 'Timeseries',]
plot_data_axima_timeseries<-plot_data_axima_timeseries[!(is.na(plot_data_axima_timeseries$day) | is.na(plot_data_axima_timeseries$Endpoint)),]
# Define the order in which the days of incubation time should be displayed
plot_data_axima_timeseries['day']<-factor(plot_data_axima_timeseries$day, levels = c('Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6'))
f1 <- ggplot(plot_data_axima_timeseries, aes(x=day, y=as.numeric(as.character(value))))+facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                            ref.group = 'Day 1',
                                            aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1)
f2 <- f1+geom_boxplot(data = plot_data_axima_timeseries[grepl('Number.*marker.*',plot_data_axima_timeseries$Endpoint),]) + geom_point(data=ylim_points,x='Day 1', alpha = 0)
f2b <- f2+geom_boxplot(data = plot_data_axima_timeseries[grepl('Median.*intensity.*',plot_data_axima_timeseries$Endpoint),])
f3 <- f2b+geom_boxplot(data = plot_data_axima_timeseries[grepl('Sum.*intensity.*',plot_data_axima_timeseries$Endpoint),])
f4 <- f3+geom_boxplot(data = plot_data_axima_timeseries[grepl('Fraction.*reproducibly.*',plot_data_axima_timeseries$Endpoint),])
f5 <- f4+geom_col(data = plot_data_axima_timeseries[grepl('.*PAPMID.*',plot_data_axima_timeseries$Endpoint),], aes(fill=eval), width = 0.7)
f7_timeseries_all <- f5+scale_fill_manual('Identification',
                                          breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker",
                                                     "prop.wrong.marker.single" ),
                                          values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')
f7_timeseries_all


# Combine the plot comparing samplepreparations and incubation time for spectra acquired on the Axima Confidence device
combined_plots<-cowplot::plot_grid(f7_samplprep_all + theme(legend.position="none") + ggtitle('Sample Preparation') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14)),
                                   f7_timeseries_all + theme(legend.position="none") + ggtitle('Time Series') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14)),
                                   rel_widths=c(0.125, 0.2), ncol = 2, align = "h",axis = "b")

# Output the combine figure
pdf('./axima_overall_noVitekMS_no_quantity.pdf', height = 16, width = 11.25)
plot_grid(combined_plots, leg.marker, ncol = 1, rel_heights = c(0.9, 0.1))
dev.off()

# Subset for spectra which were acquired with different quantities of bacterial material in a dilution series experimen
plot_data_axima_quantity<-plot_data_axima[plot_data_axima$experiment == 'Dilution Series' & !(plot_data_axima$quantity %in% c('1:3125', '1:15625')),]
plot_data_axima_quantity<-plot_data_axima_quantity[!(is.na(plot_data_axima_quantity$quantity) | is.na(plot_data_axima_quantity$Endpoint)),]


# In order to compare the diluted samples to samples which have been acquired under the '25% FA overlay' method, add measurements of the same strains (2 per phylogenetic group), which were used for the dilution series
plot_data_axima_to_add<-plot_data_axima_sample[plot_data_axima_sample$Numbering_Shipment_1 %in% c('07', '08', '09', '10', '17','18', '21','23','26','27', '33', '34','39', '40', '45', '46') &
                                                 plot_data_axima_sample$sampleprep == '25% FA overlay',]
plot_data_axima_to_add$quantity<-ifelse(is.na(plot_data_axima_to_add$quantity), '25% FA overlay', plot_data_axima_to_add$quantity)
plot_data_axima_quantity1<-merge(plot_data_bruker_quantity, plot_data_axima_to_add, by=intersect(colnames(plot_data_axima_to_add), colnames(plot_data_bruker_quantity)), all = T)

# Strain 17 is missing from microflex brukerextraction 2, axima brukerextraction 2, strain 20 is missing from microflex brukerextraction_3 and for strain 46 there are only 3 spectra in axima smear 3. As there is strong indication that strains 34 and 35 were not of the correct species in smear_1 and as_1, these were excluded. 
# In order to perform paired wilkoxon rank test, remove strain 17 and strains 34 in quantity 2 and the 4th spectrum of strain 46 of quantity_1
plot_data_axima_quantity1<-plot_data_axima_quantity1[!(plot_data_axima_quantity1$Numbering_Shipment_1 == '17' & plot_data_axima_quantity1$Method == 'quantity_2' |
                                                         plot_data_axima_quantity1$Numbering_Shipment_1 == '34' & plot_data_axima_quantity1$Method == 'quantity_1' |
                                                         plot_data_axima_quantity1$Numbering_Shipment_1 == '46' & plot_data_axima_quantity1$Method == 'quantity_1' & grepl('12$|8$', plot_data_axima_quantity1$position)),]

plot_data_axima_quantity1<-plot_data_axima_quantity1[!is.na(plot_data_axima_quantity1$Endpoint),]
plot_data_axima_quantity1<-plot_data_axima_quantity1[!is.na(plot_data_axima_quantity1$Group),]

# Define the order in which the dilution steps are displayed
plot_data_axima_quantity1$quantity<-factor(plot_data_axima_quantity1$quantity,levels = c('25% FA overlay','1:5','1:25','1:125','1:625'))

# Plot
plot_data_axima_quantity_all_groups <- ggplot(plot_data_axima_quantity1[!grepl('Evaluation', plot_data_axima_quantity1$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ ., scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                     ref.group = '25% FA overlay',
                     aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill= quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=15), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = 'none') +
  scale_fill_brewer(palette = "Blues", direction=-1)

# Output the figure 
pdf('./axima_quantity_all_groups.pdf', height = 12, width = 5)
plot_data_axima_quantity_all_groups
dev.off()

# As we expect thet diluting the samples can have taxon specific impacts, plot the dilution series results per group
# As here, no species identification results are considered, viridans and non viridans streprococci can be analysed congruently, rename accordingly
plot_data_axima_quantity1$Group<-ifelse(grepl('Strep', plot_data_axima_quantity1$Group), "Streptococci", as.character(plot_data_axima_quantity1$Group))

plot_data_axima_quantity1<-plot_data_axima_quantity1[!is.na(plot_data_axima_quantity1$quantity),]
plot_data_axima_quantity1<-plot_data_axima_quantity1[!is.na(plot_data_axima_quantity1$Group),]

plot_data_axima_quantity1$Group<-factor(plot_data_axima_quantity1$Group, levels = c("Enterobacteriaceae", "Listeria", "Burkholderia", "Bordetella",
                                                                                      "Streptococci", "Staphylococcus", "Actinobacteria", "Anaerob_Gram_negative"),
                                         labels = group_labels_strep_together_qunatity)
# Plot per group wise
plot_data_axima_quantity_per_group <- ggplot(plot_data_axima_quantity1[!grepl('Evaluation', plot_data_axima_quantity1$Endpoint),], aes(x=quantity, y=as.numeric(as.character(value))))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   +  facet_grid(Endpoint ~ Group, scales = 'free_y') + 
  stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                     ref.group = '25% FA overlay',
                     label = "p.signif", vjust = 1) +
  geom_boxplot(position = position_dodge(width = 0.9), aes(fill= quantity)) +
  geom_vline(xintercept= 1.5, colour = "black", linetype = "dashed", size = 0.4) + ggtitle('Dilution Series') + theme(axis.text=element_text(size=15), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = 'none') +
  scale_fill_brewer(palette = "Blues", direction=-1)
# Output the figure
pdf('./axima_quantity_per_groups.pdf', height = 12, width = 16)
plot_data_axima_quantity_per_group
dev.off()

# In order to assess the impact of different sample preparation protocols and incubation time per phylogenetic group, evaluate the species identification results per group
# Evaluate the identification reults per group, maldi, database and sample preparation protocol
# Set 0, to be sure not to include previous results
eval2_sampleprep<-NULL
eval2_sampleprep <- eval[!is.na(eval$sampleprep),] %>%
  group_by(Group, MALDI, sampleprep) %>%
  summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            
            prop.correct.marker.single = sum(grepl('TRUE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),             
            prop.wrong.marker.single = sum(grepl('FALSE', markerID.correct_single_species) & grepl('FALSE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.correct.marker.multi = sum(grepl('FALSE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.noID.marker = sum(grepl('FALSE', markerID.correct_single_species) & grepl('no ID', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            
            prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified))) 

# Evaluate the identification reults per group, maldi, database and incubation time
# Set 0, to be sure not to include previous results
eval2_timeseries<-NULL
eval2_timeseries<- eval[!is.na(eval$day),] %>%
  group_by(Group, MALDI, day, quantity) %>%
  summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            
            prop.correct.marker.single = sum(grepl('TRUE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),             
            prop.wrong.marker.single = sum(grepl('FALSE', markerID.correct_single_species) & grepl('FALSE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.correct.marker.multi = sum(grepl('FALSE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.noID.marker = sum(grepl('FALSE', markerID.correct_single_species) & grepl('no ID', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            
            prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified))) 

# Delet file with the same name from the previous analysis
eval2<-NULL
eval2<-merge(eval2_sampleprep, eval2_timeseries, by = intersect(colnames(eval2_sampleprep),colnames(eval2_timeseries)), all = T)

# Remove the previouse computed species identification evaluation, es the following plots will excludively include the group wise comparisons
eval<-eval[,!grepl('^prop', colnames(eval))]

# Add the species identification evaluation to the 'eval' file
eval<-merge(eval, eval2, by = intersect(colnames(eval), colnames(eval2)), all.x  = T)

# Export this file. This will be ised for the analysis on which spectra features are good endpoints for spectra quality
write.csv(eval, './08_Data/01_Spectra/04_outputs/spectra_test/eval_sum_clean.csv')

# Convert the new species identification evaluation (group wise) to the log format
# Delete file with the same name from previous analysis
eval_prop_long<-NULL
eval_prop_long<-gather(eval2, "classification_eval_prop", "value", -c(Group, MALDI, sampleprep, quantity, day))

# Add 'DB' column indicating which database was used for the species identification
eval_prop_long['DB']<-ifelse(grepl('bruker', eval_prop_long$classification_eval_prop), 'bruker', 
                             ifelse(grepl('marker', eval_prop_long$classification_eval_prop), 'marker', 
                                    ifelse(grepl('VitekMS', eval_prop_long$classification_eval_prop), 'VitekMS', NA)))

# rename the species identification evaluation column to 'eval
eval_prop_long['eval']<-eval_prop_long$classification_eval_prop
eval_prop_long$classification_eval_prop<-NULL

# Subset for data which should be plotted
plot_data<-NULL
plot_data<-eval_long[eval_long$Endpoint %in% c("n_ribos_detected","median_intensity","total_intensity_peaks",  "frac_peaks_repr", "rel_amount_su_detected"),]
plot_data$value <- ifelse(is.na(plot_data$value), 0, plot_data$value)

# Rename the 'DB' column to 'Endpoint'
names(eval_prop_long)[names(eval_prop_long) == 'DB'] <- 'Endpoint'

# Add the evalation of the species identification evaluation results to the 'plot_data'
plot_data<-merge(plot_data, eval_prop_long, by = intersect(colnames(eval_prop_long), colnames(plot_data)), all = TRUE)

# Create two 'plot_files': one with the normalised fraction of ribsosomal subunits detected (plot_data2) to compare between the phylogenetic groups and one with the absolu numbers (plot_data) when plotting per group (same strains compared)
plot_data1<-plot_data[plot_data$Endpoint != "rel_amount_su_detected",]
plot_data2<-plot_data[plot_data$Endpoint != "n_ribos_detected",]
plot_data<-plot_data1

# Rename the measurement conditions
plot_data$sampleprep<-gsub('as25','25% FA overlay',plot_data$sampleprep)
plot_data$sampleprep<-gsub('brukerextraktion','simple protein extraction', plot_data$sampleprep)
plot_data$sampleprep<-gsub('smear','direct smear', plot_data$sampleprep)
plot_data$day<-ifelse(!is.na(plot_data$day), paste('Day', plot_data$day), NA)
plot_data$quantity<-gsub('\\-','\\:', plot_data$quantity)

plot_data2$sampleprep<-gsub('as25','25% FA overlay',plot_data2$sampleprep)
plot_data2$sampleprep<-gsub('brukerextraktion','simple protein extraction', plot_data2$sampleprep)
plot_data2$sampleprep<-gsub('smear','direct smear', plot_data2$sampleprep)
plot_data2$day<-ifelse(!is.na(plot_data$day), paste('Day', plot_data2$day), NA)
plot_data2$quantity<-gsub('\\-','\\:', plot_data2$quantity)

# First compare the spectra of the phylogenetic groups which were acquired under the 'as25' method on the two devices. Subset from the 'plot_data2' file, including the relative number of ribsosomal subunits detected
# In order to compare between the phylogenetic groups, do not constider Bruker and VitekMS databases, but exclusively focus one a marker based approach as not neither the VitekMS database, not the icroflex Biotyper database cover all species included in this study.
plot_data_bruker<-plot_data2[!plot_data2$Endpoint %in% c('VitekMS', 'bruker'),]

# Subset for spectra which were acquired on the microflex Biotyper
plot_data_bruker<-plot_data_bruker[plot_data_bruker$MALDI == 'microflex Biotyper',]
plot_data_bruker<-plot_data_bruker[!is.na(plot_data_bruker$Group),]

# For this overall comparison of the phylogenetic groups, subset for spectra whic were acquired under the 'as25' method
plot_data_bruker_as25<-plot_data_bruker[plot_data_bruker$sampleprep == '25% FA overlay',]
plot_data_bruker_as25<-plot_data_bruker_as25[!is.na(plot_data_bruker_as25$Group) & !is.na(plot_data_bruker_as25$Endpoint),]

# Define order of the stacking of the species identification evaluation
plot_data_bruker_as25$eval<- factor(plot_data_bruker_as25$eval, levels = c("prop.wrong.marker.single","prop.noID.marker", "prop.correct.marker.multi","prop.correct.marker.single"))

# Define how the endpoints whould be labelled in th eplot including the relative number of ribosomal marker masses detected
labels_endpoints_rel<-c("Fraction of detected ribosomal marker peaks", "Median relative intensity\nof ribosomal marker peaks", "Sum of the intensity of all detected peaks", "Fraction of reproducibly detected peaks", "Evaluation of the PAPMID species identification")
# Introduce line breaks for very long labels
labels_endpoints_rel<-str_wrap(labels_endpoints_rel, width = 15)
# Sort and rename the endpoints
plot_data_bruker_as25$Endpoint<- factor(plot_data_bruker_as25$Endpoint, levels = c("rel_amount_su_detected", "median_intensity","total_intensity_peaks", "frac_peaks_repr","marker"), 
                                        labels = labels_endpoints_rel)
# Define how the phylogenetic groups should be displayed, including the number of strains which were included
group_labels <- c("Enterobacteriaceae (18 strains)", "Listeria (2 strains)", "Burkholderia (4 strains)", "Bordetella (3 strains)",
                  "Viridans streptococci (4 strains)", "Other streptococci (5 strains)",  "Staphylococcus (3 strains)", "Actinobacteria (5 strains)", "Gram negative Anaerobes (3 strains)")
# Sort and rename the groups accordingly
plot_data_bruker_as25$Group<- factor(plot_data_bruker_as25$Group, levels = c("Enterobacteriaceae","Listeria","Burkholderia", "Bordetella", "Viridans streptococci","Other Streptococci", "Staphylococcus", "Actinobacteria","Anaerob_Gram_negative"), 
                                     group_labels)
# Create ylim points in order to have the same y axis range for all. these will not be displayed as alpha == 0
ylim_points.rel<-data.frame(sampleprep = '25% FA overlay', day = 'Day 1', quantity = '1:5', Endpoint = rep(c("rel_amount_su_detected", "frac_peaks_repr", "total_intensity_peaks", "median_intensity", "marker"),2), value = c(0.8,1,8.7,6,1,0,0,0,0,0), value1 = c(0.8,1,8.7,6,1,0,0,0,0,0))
ylim_points.rel$Endpoint<-factor(ylim_points.rel$Endpoint, levels = c("rel_amount_su_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr", "marker"), 
                             labels = labels_endpoints_rel)
# Harmonise Gram assignment
plot_data_bruker_as25$Gram<-gsub('Negative', 'negative', plot_data_bruker_as25$Gram)
plot_data_bruker_as25$Gram<-gsub('Positive|variable', 'positive', plot_data_bruker_as25$Gram)

# Test whether there is a significant sifference between the Gram positive and the Gram negative groups
wilcox.test(plot_data_bruker_as25[grepl('Actinobacteria|Listeria|streptococci|Staphylococci',plot_data_bruker_as25$Group),'value'], 
            plot_data_bruker_as25[grepl("negative|Bordetella|Burkholderia|Enterobacteriaceae", plot_data_bruker_as25$Group),'value'])


# Summarise median and IQR per gram and endpoint
check_gram<-plot_data_bruker_as25 %>% 
  group_by(Gram, Endpoint) %>%
  summarize(IQR = quantile(value_no_log, probs = c(0.25, 0.5, 0.75), na.rm = T))

# Summarise median and IQR per gram and endpoint
check_group<-plot_data_bruker_as25 %>% 
  group_by(Group, Endpoint) %>%
  summarize(IQR = quantile(value_no_log, probs = c(0.25, 0.5, 0.75), na.rm = T))

# Plot spectra which were accquired under the 'as25' method on the microflex Biotyer device
f1 <- ggplot(plot_data_bruker_as25, aes(x=Group, y=as.numeric(as.character(value)), ymin=0,ymax=value))+facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = FALSE,method.args=list(p.adjust.method = "BH"),
                                            ref.group = 'Enterobacteriaceae (18 strains)', 
                                            aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1)
f2 <- f1+geom_boxplot(data = plot_data_bruker_as25[grepl('Fraction.*marker.*',plot_data_bruker_as25$Endpoint),]) + geom_point(data=ylim_points.rel,x='Enterobacteriaceae (18 strains)', alpha = 0)
f2b <- f2+geom_boxplot(data = plot_data_bruker_as25[grepl('Median.*intensity.*',plot_data_bruker_as25$Endpoint),])
f3 <- f2b+geom_boxplot(data = plot_data_bruker_as25[grepl('Sum.*intensity.*',plot_data_bruker_as25$Endpoint),])
f4 <- f3+geom_boxplot(data = plot_data_bruker_as25[grepl('Fraction.*reproducibly.*',plot_data_bruker_as25$Endpoint),])
f5 <- f4+geom_col(data = plot_data_bruker_as25[grepl('.*PAPMID.*',plot_data_bruker_as25$Endpoint),], aes(fill=eval), width = 0.7)
f7 <- f5+scale_fill_manual('Identification',
                           breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker",
                                      "prop.wrong.marker.single" ),
                           values = c("#009E73", "#0072B2", "#999999","#D55E00"), 
                           labels = c("Correct single species identification", "Correct multi species identification", "No identification possibe", "Wrong identification"))  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = "right") +
  ylab('') +xlab('')   
f7_group_bruker<-f7

# Output the figure
pdf("./group.comparison.bruker.nobrukerID.rel.ribo.pdf", width = 12, height = 12) 
f7_group_bruker
dev.off() 

# In order to compare the different databases, plot the species identification evaluations per group and database. For spectra acquired on the micoflex Biotyper, include species identification using a marker based approach and by the Biotyper software
# Subset for species identification evaluation as only endpoints
plot_data_sampleID_bruker<-plot_data[plot_data$MALDI == 'microflex Biotyper' & plot_data$Endpoint %in% c('marker', 'bruker') & plot_data$sampleprep == '25% FA overlay' & !is.na(plot_data$sampleprep),]

# Define how the species identification should be labelled per database
labels_endpoints_bruker_eval<-c("Evaluation of the microflex Biotyper species identification", "Evaluation of the PAPMID species identification")
# Introduce linebreaks for too long names
labels_endpoints_bruker_eval<-str_wrap(labels_endpoints_bruker_eval, width = 15)

# Soft and label database outputs accordingly
plot_data_sampleID_bruker$Endpoint<- factor(plot_data_sampleID_bruker$Endpoint, levels = c("bruker", "marker"), labels = labels_endpoints_bruker_eval)

# Define in which order the groups should be displayed
plot_data_sampleID_bruker$Group<- factor(plot_data_sampleID_bruker$Group, levels = c("Enterobacteriaceae","Listeria","Burkholderia", "Bordetella", "Viridans streptococci","Other Streptococci", "Staphylococcus", "Actinobacteria","Anaerob_Gram_negative"), 
                                         group_labels)

# Define order of the stacking of the species identification evaluation per database
plot_data_sampleID_bruker$eval<- factor(plot_data_sampleID_bruker$eval, levels = c( "prop.wrong.marker.single", "prop.wong.bruker.over.two",
                                                                                    "prop.wong.bruker.under.two", 
                                                                                    "prop.noID.marker", "prop.noID.bruker",
                                                                                    "prop.bruker.correct.under.two", 
                                                                                    "prop.correct.marker.multi", "prop.bruker.correct.multi.over.two", 
                                                                                    "prop.correct.marker.single", "prop.bruker.correct.single.over.two"))
# Plot
f1 <- ggplot(plot_data_sampleID_bruker, aes(x=Group, y=as.numeric(as.character(value)), ymin=0,ymax=value)) + facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   
f5 <- f1+geom_col(data = plot_data_sampleID_bruker[grepl("PAPMID", plot_data_sampleID_bruker$Endpoint),], aes(fill=eval), width = 0.7) 
f6 <- f5+geom_col(data = plot_data_sampleID_bruker[grepl("microflex", plot_data_sampleID_bruker$Endpoint),], aes(fill=eval), width = 0.7) 
f7 <- f6+scale_fill_manual('Identification',
                           breaks = c("prop.wong.bruker.over.two","prop.wrong.marker.single","prop.wong.bruker.under.two","prop.noID.bruker","prop.noID.marker","prop.bruker.correct.under.two",      
                                      "prop.bruker.correct.multi.over.two","prop.correct.marker.multi","prop.bruker.correct.single.over.two","prop.correct.marker.single"),
                           values = c("#D55E00","#D55E00","#E69F00","#999999","#999999","#56B4E9", "#0072B2", "#0072B2", "#009E73","#009E73")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   

# Rename the plot
f7_group_bruker_speciesID<-f7


# In order to elucidate which measurement conditions are the best per group, plot the different preperation protocols / incubation time per group seperately. 
# First produce these group wise plots for spectra acquired on the microflex Biotyper device, then for spectra acquired on the Axima Confidence, then combine and export
# as here the strains compared are again the same, use the absolut number of peaks as endpoint
plot_data_bruker<-plot_data[!plot_data$Endpoint %in% c('VitekMS', 'bruker'),]
plot_data_bruker<-plot_data_bruker[plot_data_bruker$MALDI == 'microflex Biotyper',]
plot_data_bruker<-plot_data_bruker[!is.na(plot_data_bruker$Group),]

# Rename file to 'plot_data_bruker_all', before subsetting per group
plot_data_bruker_all<-plot_data_bruker
# Add 'experiment' column for easier subsetting
plot_data_bruker_all['experiment']<-ifelse(!is.na(plot_data_bruker_all$sampleprep), 'Sample prepreration', 
                                ifelse(!is.na(plot_data_bruker_all$day), 'Timeseries', 
                                       ifelse(!is.na(plot_data_bruker_all$quantity), 'Dilution Series', NA)))
plot_data_bruker_all$eval<-factor(plot_data_bruker_all$eval, levels = c("prop.wrong.marker.single", "prop.noID.marker", "prop.correct.marker.multi", "prop.correct.marker.single"))

# In order to be able to include the spectra which have been acquired under the 'as25' method at day 1 in both, the 'Samplepreparation' dataset AND the 'timeseries dataset, duplicate these rows.
# As in the datasets as25_1 and # as25_2 there are spectra missing, use 'as25_3' for this task
# Subset the rows to duplicate 
add<-plot_data_bruker_all[plot_data_bruker_all$Method == 'as25_3' & !is.na(plot_data_bruker_all$sampleprep),]
add['experiment']<-'Timeseries'

# add the duplicated rows to the 'plot_data' file
plot_data_bruker_all<-merge(plot_data_bruker_all, add, by = intersect(colnames(plot_data_bruker_all), colnames(add)), all = T)

#plot_data_bruker_all['which_variable']<-ifelse(plot_data_bruker_all$experiment == 'Sample prepreration', plot_data_bruker_all$sampleprep,
#                                    ifelse(plot_data_bruker_all$experiment == 'Timeseries', plot_data_bruker_all$day, 
#                                           ifelse(plot_data_bruker_all$experiment == 'Dilution Series', plot_data_bruker_all$quantity, NA)))

plot_data_bruker_all<-plot_data_bruker_all[!is.na(plot_data_bruker_all$Group),]

# Create ylim opints (with absolut number of ribsosomal marker masses detected)
ylim_points<-data.frame(sampleprep = '25% FA overlay', day = 'Day 1', quantity = '1:5', Endpoint = rep(c("n_ribos_detected", "frac_peaks_repr", "total_intensity_peaks", "median_intensity", "marker"),2), value = c(35,1,8.7,8,1,0,0,0,0,0), value1 = c(38,1,8.7,8,1,0,0,0,0,0))
ylim_points$Endpoint<-factor(levels = c("n_ribos_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr", "marker"), 
                             ylim_points$Endpoint, labels = labels_endpoints)

# In order to perform paired wilcoxon rank text, remove the spectra of strains whcih are missing in oher set
plot_data_bruker_all<-plot_data_bruker_all[!(plot_data_bruker_all$Numbering_Shipment_1 == '17' & plot_data_bruker_all$Method %in% c("smear_2","as25_2") |
                                               plot_data_bruker_all$Numbering_Shipment_1 == '20' & plot_data_bruker_all$Method %in% c("smear_3","as25_3")|
                                                plot_data_bruker_all$Numbering_Shipment_1 %in% c('34', '35') & plot_data_bruker_all$Method == "brukerextraction_1"),]


# Loop through each group, save and output the plots on varying sample preparation protocols and quantities
combines_plots_bruker<-list()
for (i in 1:length(unique(as.character(plot_data_bruker_all$Group)))){
  plot_data_bruker<-NULL
  plot_data_bruker<-plot_data_bruker_all[plot_data_bruker_all$Group == unique(as.character(plot_data_bruker_all$Group))[i],]
  plot_data_bruker$Endpoint<-factor(plot_data_bruker$Endpoint, levels = c("n_ribos_detected","median_intensity","total_intensity_peaks","frac_peaks_repr","marker"), labels = labels_endpoints)
  # subset sample prep
  plot_data_bruker_sample<-plot_data_bruker[plot_data_bruker$experiment == 'Sample prepreration',]
  plot_data_bruker_sample<-plot_data_bruker_sample[!(is.na(plot_data_bruker_sample$sampleprep) | is.na(plot_data_bruker_sample$Endpoint)),]
  plot_data_bruker_sample['sampleprep']<-factor(plot_data_bruker_sample$sampleprep, levels = c('direct smear', '25% FA overlay', 'simple protein extraction'))
  
  f1 <- ggplot(plot_data_bruker_sample, aes(x=sampleprep, y=as.numeric(as.character(value))))+facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                              ref.group = 'direct smear', 
                                              label = "p.signif", vjust = 1)
  f2 <- f1+geom_boxplot(data = plot_data_bruker_sample[grepl('Number.*marker.*',plot_data_bruker_sample$Endpoint),]) + geom_point(data=ylim_points,x='direct smear', alpha = 0)
  f2b <- f2+geom_boxplot(data = plot_data_bruker_sample[grepl('Median.*intensity.*',plot_data_bruker_sample$Endpoint),])
  f3 <- f2b+geom_boxplot(data = plot_data_bruker_sample[grepl('Sum.*intensity.*',plot_data_bruker_sample$Endpoint),])
  f4 <- f3+geom_boxplot(data = plot_data_bruker_sample[grepl('Fraction.*reproducibly.*',plot_data_bruker_sample$Endpoint),])
  f5 <- f4+geom_col(data = plot_data_bruker_sample[grepl('.*PAPMID.*',plot_data_bruker_sample$Endpoint),], aes(fill=eval), width = 0.7)
  f7_samplprep_all <- f5+scale_fill_manual('Identification',
                                           breaks = c("prop.wrong.marker.single", "prop.noID.marker","prop.correct.marker.multi", "prop.correct.marker.single"),
                                           values = c("#D55E00", "#999999", "#0072B2","#009E73")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   
  
  # subset timeseries
  plot_data_bruker_timeseries<-plot_data_bruker[plot_data_bruker$experiment == 'Timeseries',]
  plot_data_bruker_timeseries<-plot_data_bruker_timeseries[!(is.na(plot_data_bruker_timeseries$day) | is.na(plot_data_bruker_timeseries$Endpoint)),]
  plot_data_bruker_timeseries<-as.data.frame(plot_data_bruker_timeseries[!(is.na(plot_data_bruker_timeseries$day) | is.na(plot_data_bruker_timeseries$Endpoint)),])
  plot_data_bruker_timeseries <- plot_data_bruker_timeseries %>% filter((as.character(Numbering_Shipment_1) != '20')|is.na(Numbering_Shipment_1) ) 
  plot_data_bruker_timeseries <- plot_data_bruker_timeseries %>% filter((as.character(Numbering_Shipment_1) != '32')|is.na(Numbering_Shipment_1) ) 
  plot_data_bruker_timeseries['day']<-factor(plot_data_bruker_timeseries$day, levels = c('Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6'))
  
  f1 <- ggplot(plot_data_bruker_timeseries, aes(x=day, y=as.numeric(as.character(value))))+facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                              ref.group = 'Day 1', 
                                              label = "p.signif", vjust = 1)
  f2 <- f1+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Number.*marker.*',plot_data_bruker_timeseries$Endpoint),])  + geom_point(data=ylim_points,x='Day 1', alpha = 0)
  f2b <- f2+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Median.*intensity.*',plot_data_bruker_timeseries$Endpoint),])
  f3 <- f2b+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Sum.*intensity.*',plot_data_bruker_timeseries$Endpoint),])
  f4 <- f3+geom_boxplot(data = plot_data_bruker_timeseries[grepl('Fraction.*reproducibly.*',plot_data_bruker_timeseries$Endpoint),])
  f5 <- f4+geom_col(data = plot_data_bruker_timeseries[grepl('.*PAPMID.*',plot_data_bruker_timeseries$Endpoint),], aes(fill=eval), width = 0.7)
  f7_timeseries_all <- f5+scale_fill_manual('Identification',
                                            breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker", "prop.wrong.marker.single" ),
                                            values = c("#009E73", "#0072B2", "#999999","#D55E00"))  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   
  combined_plots<-cowplot::plot_grid(f7_samplprep_all + theme(legend.position="none") + ggtitle('microflex Biotyper\nSamplepreparation') + theme(axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=16), strip.text = element_text(size=13)),
                                     f7_timeseries_all + theme(legend.position="none") + ggtitle('microflex Biotyper\nTimeseries') + theme(axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=16), strip.text = element_text(size=13)),
                                     rel_widths=c(0.125, 0.2), ncol = 2, align = "h",axis = "b")
  combined_plots_legend<-cowplot::plot_grid(combined_plots, cowplot::get_legend(f7_timeseries_all), ncol = 2, rel_widths = c(0.8, 0.2))
  combines_plots_bruker[[unique(as.character(plot_data_bruker_all$Group))[i]]]<-combined_plots
  pdf(paste0('./bruker.noBrukerID.no.quantity', unique(as.character(plot_data_bruker_all$Group))[i], '.pdf'), height = 12, width = 15)
    print(plot_grid(combined_plots, leg.marker, ncol = 1, rel_heights = c(0.9, 0.1)))
  dev.off()
}


# Produce the same plot for apectra acquired on the Axima Confidence device. First overall group comparisons, acquired under the 'as25' method group wise, then species identification evaluation by the different databases, then comparing sample preparation and varying dilutions per group
# Subset for spectra which have been acquired on the Axima Confidence. 
# Do notinclude VitekMS / Biotyper species identification evaluations, as these do not include all species covered in this study
# Subset from the plot_data2 file including the relative number of ribsosmal subunits detected
plot_data_axima<-plot_data2[!plot_data2$Endpoint %in% c('bruker', 'VitekMS'),]
plot_data_axima<-plot_data_axima[plot_data_axima$MALDI == 'Axima Confidence',]

# Define in which order the evaluation should be stacked
plot_data_axima$eval<- factor(plot_data_axima$eval, levels = c("prop.wrong.marker.single", "prop.noID.marker", "prop.correct.marker.multi", "prop.correct.marker.single"))

plot_data_axima$Endpoint<- factor(plot_data_axima$Endpoint, levels = c("rel_amount_su_detected", "median_intensity","total_intensity_peaks", "frac_peaks_repr","marker"),
                                  labels = labels_endpoints_rel)

# For overall group comparison, subset for spectra acquired under the 'as25' method
plot_data_axima_as25<-plot_data_axima[grepl('25% FA overlay',plot_data_axima$sampleprep),]

# Define in which order the groups should be displayed
plot_data_axima_as25$Group<- factor(plot_data_axima_as25$Group, levels = c("Enterobacteriaceae","Listeria","Burkholderia", "Bordetella", "Viridans streptococci","Other Streptococci", "Staphylococcus", "Actinobacteria","Anaerob_Gram_negative"), 
                                    group_labels)
# Plot
f1 <- ggplot(plot_data_axima_as25, aes(x=Group, y=as.numeric(as.character(value)), ymin=0,ymax=value))+facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = FALSE,method.args=list(p.adjust.method = "BH"),
                                            ref.group = 'Enterobacteriaceae (18 strains)', 
                                            aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1)
f2 <- f1+geom_boxplot(data = plot_data_axima_as25[grepl('Fraction.*marker.*',plot_data_axima_as25$Endpoint),]) + geom_point(data=ylim_points.rel,x='Enterobacteriaceae (18 strains)', alpha = 0)
f2b <- f2+geom_boxplot(data = plot_data_axima_as25[grepl('Median.*intensity.*',plot_data_axima_as25$Endpoint),])
f3 <- f2b+geom_boxplot(data = plot_data_axima_as25[grepl('Sum.*intensity.*',plot_data_axima_as25$Endpoint),])
f4 <- f3+geom_boxplot(data = plot_data_axima_as25[grepl('Fraction.*reproducibly.*',plot_data_axima_as25$Endpoint),])
f5 <- f4+geom_col(data = plot_data_axima_as25[grepl('.*PAPMID.*',plot_data_axima_as25$Endpoint),], aes(fill=eval), width = 0.7)
f7 <- f5+scale_fill_manual('Identification',
                           breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker",
                                      "prop.wrong.marker.single" ),
                           values = c("#009E73", "#0072B2", "#999999","#D55E00"), 
                           labels = c("Correct single species identification", "Correct multi species identification", "No identification possibe", "Wrong identification"))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=15), strip.text = element_text(size=12), legend.position = "right") +
  ylab('') +xlab('')   
# rename the plot
f7_group_axima<-f7

# Output figure
pdf("./group.comparison.axima.ribo.rel.pdf",  width = 12, height = 12) 
f7_group_axima
dev.off() 

# Combine the plots of the overall group comparison for spectra acquired on the microflex Biotyper and the Axima Confidence device
combined_plots<-cowplot::plot_grid(f7_group_bruker +  ggtitle('microflex Biotyper') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14), legend.position = "none"),
                                   f7_group_axima +  ggtitle('Axima Confidence') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14), legend.position = "none"),
                                   rel_widths=c(0.5, 0.5), ncol = 2, align = "h",axis = "b")
# Output the combined figure
pdf("./group.comparison.both.only.markerID.rel.ribo.pdf", width = 16, height = 16) 
plot_grid(combined_plots, leg.marker, ncol = 1, rel_heights = c(0.9, 0.1))
dev.off() 

# Compare the species identification evaluation
# Subset for spectra acquired on the Axima Confidence under the 'as25' method and with the species identification evaluation as single endpoints
plot_data_sampleID_axima<-plot_data[plot_data$MALDI == 'Axima Confidence' & plot_data$Endpoint %in% c('marker', 'VitekMS') & plot_data$sampleprep == '25% FA overlay' & !is.na(plot_data$sampleprep),]

# Define in which order the groups should be dispplayed
plot_data_sampleID_axima$Group<- factor(plot_data_sampleID_axima$Group, levels =  c("Enterobacteriaceae","Listeria","Burkholderia", "Bordetella", "Viridans streptococci","Other Streptococci", "Staphylococcus", "Actinobacteria","Anaerob_Gram_negative"), 
                                        group_labels)

# Define how these should be labelles
labels_endpoints_axima_eval<-c("Evaluation of the VitekMS species identification", "Evaluation of the PAPMID species identification")
labels_endpoints_axima_eval<-str_wrap(labels_endpoints_axima_eval, width = 15)

# Rename the enpoints accordingly
plot_data_sampleID_axima$Endpoint<- factor(plot_data_sampleID_axima$Endpoint, levels = c("VitekMS", "marker"), labels = labels_endpoints_axima_eval)

# Define in which order the species identification evaluation sould be stacked
plot_data_sampleID_axima$eval<- factor(plot_data_sampleID_axima$eval, levels = c( "prop.wrong.marker.single", "prop.wrong.VitekMS.single",
                                                                              "prop.wrong.VitekMS.lowdiscriminatory", 
                                                                              "prop.noID.marker", "prop.noID.VitekMS",
                                                                              "prop.correct.VitekMS.lowdiscriminatory", 
                                                                              "prop.correct.marker.multi", 
                                                                              "prop.correct.marker.single", "prop.correct.VitekMS.single"))

# Plot
f1 <- ggplot(plot_data_sampleID_axima, aes(x=Group, y=as.numeric(as.character(value)), ymin=0,ymax=value)) + facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   
f5 <- f1+geom_col(data = plot_data_sampleID_axima[grepl("PAPMID", plot_data_sampleID_axima$Endpoint),], aes(fill=eval), width = 0.7) 
f6 <- f5+geom_col(data = plot_data_sampleID_axima[grepl("VitekMS", plot_data_sampleID_axima$Endpoint),], aes(fill=eval), width = 0.7) 
f7 <- f6+scale_fill_manual('Identification',
                           breaks = c("prop.wrong.VitekMS.single","prop.wrong.marker.single","prop.wrong.VitekMS.lowdiscriminatory","prop.noID.VitekMS","prop.noID.marker","prop.correct.VitekMS.lowdiscriminatory",      
                                      "prop.correct.marker.multi","prop.correct.VitekMS.single","prop.correct.marker.single"),
                           values = c("#D55E00","#D55E00","#E69F00","#999999","#999999","#0072B2","#0072B2", "#009E73","#009E73")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   
# Rename the plot
f7_group_axima_speciesID<-f7

# Combine the two species identification plots
combined_plots<-cowplot::plot_grid(f7_group_bruker_speciesID +  ggtitle('microflex Biotyper') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14), legend.position = "none"),
                                   f7_group_axima_speciesID +  ggtitle('Axima Confidence') + theme(axis.text=element_text(size=18), axis.title=element_text(size=24),title=element_text(size=18), strip.text = element_text(size=14), legend.position = "none"),
                                   rel_widths=c(0.5, 0.5), ncol = 2, align = "h",axis = "b")
# Output the combined figure
pdf("./group.comparison.species_ID_all_DB.pdf", width = 14, height = 7.2) 
combined_plots
dev.off() 


# Subset again for spectra acquired on the Axima Confidence device, this time from the plot_data file, including the absolut number of ribsosomal subunits detected
# Here exlcusivly include species identification results by a marker based database, as neither the microflex Biotyper, not the Axima Confidence database cover all databases included in this study
plot_data_axima<-plot_data[!plot_data$Endpoint %in% c('bruker', 'VitekMS'),]
plot_data_axima<-plot_data_axima[plot_data_axima$MALDI == 'Axima Confidence',]

# Define in which order the species identification evaluation should be stacked
plot_data_axima$eval<- factor(plot_data_axima$eval, levels = c("prop.wrong.marker.single", "prop.noID.marker", "prop.correct.marker.multi", "prop.correct.marker.single"))

# Define in which order and how the endpoints sould be labelled
plot_data_axima$Endpoint<- factor(plot_data_axima$Endpoint, levels = c("n_ribos_detected", "median_intensity","total_intensity_peaks", "frac_peaks_repr","marker"),
                                  labels = labels_endpoints)

# Strain 17 is missing from microflex brukerextraction 2, axima brukerextraction 2, strain 20 is missing from microflex brukerextraction_3 and for strain 46 there are only 3 spectra in axima smear 3. As there is strong indication that strains 34 and 35 were not of the correct species in smear_1 and as_1, these were excluded. 
# In order to perform paired wilkoxon rank test, remove strain 17 and strains 34 in quantity 2 and the 4th spectrum of strain 46 of quantity_1
plot_data_axima<-plot_data_axima[!(plot_data_axima$Numbering_Shipment_1 == '17' & plot_data_axima$Method %in% c("smear_2","as25_2") |
                                     plot_data_axima$Numbering_Shipment_1 == '46' & plot_data_axima$Method %in% c("as25_2", "brukerextraction_2") & grepl('4$',plot_data_axima$spectra)|
                                     plot_data_axima$Numbering_Shipment_1 %in%  c('34','35') & plot_data_axima$Method  == "brukerextraction_1"),]
# Before subsetting per group, rename the file to 'plot_data_axima_all
plot_data_axima_all<-plot_data_axima

# Add 'experiment column for easire subsetting
plot_data_axima_all['experiment']<-ifelse(!is.na(plot_data_axima_all$sampleprep), 'Sample prepreration', 
                                          ifelse(!is.na(plot_data_axima_all$day), 'Timeseries', 
                                                 ifelse(!is.na(plot_data_axima_all$quantity), 'Dilution Series', NA)))

# In order to be able to include the spectra which have been acquired under the 'as25' method at day 1 in both, the 'Samplepreparation' dataset AND the 'timeseries dataset, duplicate these rows.
# As in the datasets as25_1 and # as25_2 there are spectra missing, use 'as25_3' for this task
# Subset the rows to duplicate 
add<-plot_data_axima_all[plot_data_axima_all$Method == 'as25_3' & !is.na(plot_data_axima_all$sampleprep),]
add['experiment']<-'Timeseries'
add['day']<-'Day 1'

# Add the duplicated files to the dataset to plot
plot_data_axima_all<-merge(plot_data_axima_all, add, by = intersect(colnames(plot_data_axima_all), colnames(add)), all = T)
plot_data_axima_all<-plot_data_axima_all[!is.na(plot_data_axima_all$Group),]

# Loop through all groups and plot spectra features and species identifications per sample preperation protocol and with varying quantity
combines_plots_axima<-list()
for (i in 1:length(unique(as.character(plot_data_axima_all$Group)))){
  plot_data_axima<-NULL
  plot_data_axima<-plot_data_axima_all[plot_data_axima_all$Group == unique(as.character(plot_data_axima_all$Group))[i],]
  plot_data_axima_sample<-plot_data_axima[plot_data_axima$experiment == 'Sample prepreration',]
  plot_data_axima_sample<-plot_data_axima_sample[!(is.na(plot_data_axima_sample$sampleprep) | is.na(plot_data_axima_sample$Endpoint)),]
  plot_data_axima_sample['sampleprep']<-factor(plot_data_axima_sample$sampleprep, levels = c('direct smear', '25% FA overlay', 'simple protein extraction'))
  
  f1 <- ggplot(plot_data_axima_sample, aes(x=sampleprep, y=as.numeric(as.character(value))))+facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                              ref.group = 'direct smear', 
                                              label = "p.signif", vjust = 1)
  f2 <- f1+geom_boxplot(data = plot_data_axima_sample[grepl('Number.*marker.*',plot_data_axima_sample$Endpoint),]) + geom_point(data=ylim_points,x='direct smear', alpha = 0)
  f2b <- f2+geom_boxplot(data = plot_data_axima_sample[grepl('Median.*intensity.*',plot_data_axima_sample$Endpoint),])
  f3 <- f2b+geom_boxplot(data = plot_data_axima_sample[grepl('Sum.*intensity.*',plot_data_axima_sample$Endpoint),])
  f4 <- f3+geom_boxplot(data = plot_data_axima_sample[grepl('Fraction.*reproducibly.*',plot_data_axima_sample$Endpoint),])
  f5 <- f4+geom_col(data = plot_data_axima_sample[grepl('.*PAPMID.*',plot_data_axima_sample$Endpoint),], aes(fill=eval), width = 0.7)
  f7_samplprep_all <- f5+scale_fill_manual('Identification',
                                            breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker", "prop.wrong.marker.single" ),
                                           values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   
  
  # subset timeseries
  plot_data_axima_timeseries<-plot_data_axima[plot_data_axima$experiment == 'Timeseries',]
  plot_data_axima_timeseries<-plot_data_axima_timeseries[!(is.na(plot_data_axima_timeseries$day) | is.na(plot_data_axima_timeseries$Endpoint)),]
  plot_data_axima_timeseries<-as.data.frame(plot_data_axima_timeseries[!(is.na(plot_data_axima_timeseries$day) | is.na(plot_data_axima_timeseries$Endpoint)),])
  plot_data_axima_timeseries <- plot_data_axima_timeseries %>% filter((as.character(Numbering_Shipment_1) != '20')|is.na(Numbering_Shipment_1) ) 
  plot_data_axima_timeseries['day']<-factor(plot_data_axima_timeseries$day, levels = c('Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6'))
  
  f1 <- ggplot(plot_data_axima_timeseries, aes(x=day, y=as.numeric(as.character(value))))+facet_grid(Endpoint ~., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   + stat_compare_means(method = "wilcox.test", paired = TRUE,method.args=list(p.adjust.method = "BH"),
                                              ref.group = 'Day 1', 
                                              label = "p.signif", vjust = 1)
  f2 <- f1+geom_boxplot(data = plot_data_axima_timeseries[grepl('Number.*marker.*',plot_data_axima_timeseries$Endpoint),]) + geom_point(data=ylim_points,x='Day 1', alpha = 0)
  f2b <- f2+geom_boxplot(data = plot_data_axima_timeseries[grepl('Median.*intensity.*',plot_data_axima_timeseries$Endpoint),])
  f3 <- f2b+geom_boxplot(data = plot_data_axima_timeseries[grepl('Sum.*intensity.*',plot_data_axima_timeseries$Endpoint),])
  f4 <- f3+geom_boxplot(data = plot_data_axima_timeseries[grepl('Fraction.*reproducibly.*',plot_data_axima_timeseries$Endpoint),])
  f5 <- f4+geom_col(data = plot_data_axima_timeseries[grepl('.*PAPMID.*',plot_data_axima_timeseries$Endpoint),], aes(fill=eval), width = 0.7)
  f7_timeseries_all <- f5+scale_fill_manual('Identification',
                                            #                                            breaks = c("prop.wong.bruker.over.two","prop.wrong.marker.single","prop.wong.bruker.under.two","prop.noID.bruker","prop.noID.marker","prop.bruker.correct.under.two",      
                                            #                                                       "prop.bruker.correct.multi.over.two","prop.correct.marker.multi","prop.bruker.correct.single.over.two","prop.correct.marker.single"),
                                            #                                            values = c("#56B4E9", "#009E73", "#0072B2", "#56B4E9", "#009E73", "#999999","#999999","#D55E00","#E69F00","#D55E00")) + 
                                            breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker", "prop.wrong.marker.single" ),
                                            values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')   
  combined_plots<-cowplot::plot_grid(f7_samplprep_all + theme(legend.position="none") + ggtitle('Axima Confidence\nSamplepreparation') + theme(axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=16), strip.text = element_text(size=13)),
                                     f7_timeseries_all + theme(legend.position="none") + ggtitle('Axima Confidence\nTimeseries') + theme(axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=16), strip.text = element_text(size=13)),
                                     rel_widths=c(0.125, 0.2), ncol = 2, align = "h",axis = "b")
  combined_plots_legend<-cowplot::plot_grid(combined_plots, cowplot::get_legend(f7_timeseries_all), ncol = 2, rel_widths = c(0.8, 0.2))
  combines_plots_axima[[unique(as.character(plot_data_axima_all$Group))[i]]]<-combined_plots
  pdf(paste0('./axima.noVitekMSIDID.no.quantity', unique(as.character(plot_data_axima_all$Group))[i], '.pdf'), height = 12, width = 15)
  print(plot_grid(combined_plots, leg.marker, ncol = 1, rel_heights = c(0.9, 0.1)))
  dev.off()
}


# Loop through the group wise plot for spectra acquired on the microflex Biotyper and the Axima Confidence, combine and output these
for (names in names(combines_plots_bruker)){
  combined_plots <- cowplot::plot_grid(combines_plots_bruker[[names]], combines_plots_axima[[names]], 
                                       rel_widths=c(0.5, 0.5), ncol = 2, align = "h",axis = "b")
  pdf(paste0('./', names, 'bothDevices.pdf'), height = 12, width = 17)
  print(plot_grid(combined_plots, leg.marker, ncol = 1, rel_heights = c(0.9, 0.1)))
  dev.off()
}


# Include plot to show that with increased spectra quality, better identification accuracy and resolution is possible
# For this and in order to include species identification resulta from all databases, include exclusively strains of species which are covered by all three databases
# Exclusively include strains of the three subgroups viridans streptococci, enterobacter cloacae complex and Burkholderia cepacia complex as these pose challenges for bacterial species identification by MALDI-TOF MS
# For viridans streptococci, these are: S. pneumoniae (#24) and S. pseudopneumoniae (#29)
# For the Burkholderia cepacia complex these are: B. contaminans (#18), B. multivorans (#19) and B.cenocepacia (#20)
# For the Enterobacter cloace complex these are: E. hormaechei (#36), E. asburiae (#37)  and E. ludwigii (#38)
# Exclusivela include spectra which were acquired on Day 1 and with different sample preperation protocols (smear, as25, brukerextraction)
# Subset for the right strains
eval_better_ID<-eval[eval$Numbering_Shipment_1 %in% c('24', '29','18', '19', '20', '36','37', '38') &is.na(eval$quantity) & grepl('as25|smear|extraction',eval$Method),]

# log transform intensity endpoints
eval_better_ID$total_intensity_peaks<-log10(as.numeric(eval_better_ID$total_intensity_peaks))
eval_better_ID$median_intensity<-log10(as.numeric(eval_better_ID$total_intensity_peaks))

# replace NA with 0
eval_better_ID$total_intensity_peaks<-ifelse(is.na(eval_better_ID$total_intensity_peaks),0, as.numeric(as.character(eval_better_ID$total_intensity_peaks)))
eval_better_ID$median_intensity<-ifelse(is.na(eval_better_ID$median_intensity),0, as.numeric(as.character(eval_better_ID$median_intensity)))
eval_better_ID$n_ribos_detected<-ifelse(is.na(eval_better_ID$n_ribos_detected),0, as.numeric(as.character(eval_better_ID$n_ribos_detected)))

# Split the spectra of each grou in three intensity groups, according to the sum of intensity of all peaks in each spectrum. The total intensity value in each group is split into three equal parts and the spectra are assigned a 'ointensity rou' accordingly
eval_better_ID <- eval_better_ID %>% group_by(MALDI, Group) %>% mutate(intensity_group = cut(total_intensity_peaks, 3, labels = FALSE))
eval_better_ID$intensity_group<-paste0('Intensity group ', eval_better_ID$intensity_group)
ggplot(eval_better_ID, aes(x = total_intensity_peaks, color=intensity_group)) +
  geom_histogram(binwidth = 0.02) + facet_grid(.~MALDI)

# Evaluate species identification evaluation per group, maldi and ribogroup 
eval2_better_ID<- eval_better_ID %>%
  group_by(Group, MALDI, intensity_group) %>%
  summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            
            prop.correct.marker.single = sum(grepl('TRUE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),             
            prop.wrong.marker.single = sum(grepl('FALSE', markerID.correct_single_species) & grepl('FALSE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.correct.marker.multi = sum(grepl('FALSE', markerID.correct_single_species) & grepl('TRUE', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            prop.noID.marker = sum(grepl('FALSE', markerID.correct_single_species) & grepl('no ID', markerID.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerID.correct_species_identified)),
            
            prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
            prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified))) 

# Merge to the 'eval' file
eval_better_ID<-merge(eval_better_ID, eval2_better_ID, by = intersect(colnames(eval_better_ID), colnames(eval2_better_ID)), all.x  = T)

# Convert to long format
eval_better_ID_prop_long<-gather(eval2_better_ID, "classification_eval_prop", "value", -c(Group, MALDI, intensity_group))

# Include 'DB' column
eval_better_ID_prop_long['DB']<-ifelse(grepl('bruker', eval_better_ID_prop_long$classification_eval_prop), 'bruker', 
                             ifelse(grepl('marker', eval_better_ID_prop_long$classification_eval_prop), 'marker', 
                                    ifelse(grepl('VitekMS', eval_better_ID_prop_long$classification_eval_prop), 'VitekMS', NA)))
# Rename the column including the species identification evaluation to 'eval'
eval_better_ID_prop_long['eval']<-eval_better_ID_prop_long$classification_eval_prop
eval_better_ID_prop_long$classification_eval_prop <- NULL

# Define in which order the species identification evaluation should be stacked
eval_better_ID_prop_long$eval<- factor(eval_better_ID_prop_long$eval, levels = c("prop.wrong.VitekMS.single", "prop.wrong.marker.single", "prop.wong.bruker.over.two",
                                                             "prop.wrong.VitekMS.lowdiscriminatory", "prop.wong.bruker.under.two", 
                                                             "prop.noID.VitekMS", "prop.noID.marker", "prop.noID.bruker",
                                                             "prop.bruker.correct.under.two", "prop.correct.VitekMS.lowdiscriminatory",
                                                             "prop.correct.marker.multi", "prop.bruker.correct.multi.over.two", 
                                                             "prop.correct.VitekMS.single", "prop.correct.marker.single", "prop.bruker.correct.single.over.two"))


#Gather all Endpoints to one column in order to plot all at once using facet_wrap
eval_long_better_ID<-gather(eval_better_ID[,c("spectra","n_ribos_detected","median_intensity","total_intensity_peaks", "frac_peaks_repr","mean_dist_ppm","brukerDB.Score",  "markerID.match_count",  "VitekMS.DBscore_rank_1")], "Endpoint", "value", -spectra)
eval_long_better_ID<-merge(eval_better_ID[,c("spectra","Strain", "Numbering_Shipment_1", "intensity_group","Group", "species_NGS", "genus_NGS", "MALDI", "Method", "run", "position", "day", "sampleprep", "quantity", "Run")], eval_long_better_ID, by='spectra')

plot_data_betterID<-eval_long_better_ID

# Subset for the endpoints of interest, also including the scores of the respective databases

# Rename DB to 'Endpoint'
names(eval_better_ID_prop_long)[names(eval_better_ID_prop_long) == 'DB'] <- 'Endpoint'

# Merge the two datasets
plot_data_betterID<-merge(plot_data_betterID, eval_better_ID_prop_long, by = intersect(colnames(eval_better_ID_prop_long), colnames(plot_data_betterID)), all = TRUE)

# Because we use the sum of the intensity of all peaks to group the spectra, we can take out total int as row
plot_data_betterID<-plot_data_betterID[plot_data_betterID$Endpoint != 'total_intensity_peaks',]

#Subset for spectra acquired on the microflex device
plot_data_betterID_bruker<-plot_data_betterID[plot_data_betterID$MALDI == 'microflex Biotyper' & !grepl('VitekMS', plot_data_betterID$Endpoint),]
# Define in which order the groups should be displayed and how they should be labelled
plot_data_betterID_bruker$Group<- factor(plot_data_betterID_bruker$Group, levels =  c("Burkholderia", "Viridans streptococci", "Enterobacteriaceae"), 
                                         labels = c("Burkholderia\ncepacia complex\n(3 strains)", "viridans\nstreptococci\n(2 strains)", "Enterobacter\ncloace complex\n(3 strains)"))

# Define how the endpoints shoudl be displayed
labels_endpoints_bruker_betterID<-c("Number of detected ribosomal marker peaks","Median relative intensity of the ribosomal marker peaks","Fraction of reproducibly detected peaks","Measurement error [ppm]","microflex Biotyper log Score", "Evaluation of the microflex Biotype species identification", "PAPMID Match Count Score","Evaluation of the PAPMID species identification")

# Introduce line breaks for very long labels
labels_endpoints_bruker_betterID<-str_wrap(labels_endpoints_bruker_betterID, width = 15)


# Rename the endpoints accordingly
plot_data_betterID_bruker$Endpoint<- factor(plot_data_betterID_bruker$Endpoint, levels = c("n_ribos_detected","median_intensity", "frac_peaks_repr", "mean_dist_ppm","brukerDB.Score", "bruker", "markerID.match_count","marker"), labels = labels_endpoints_bruker_betterID)


# Define in which order the species identification evaluation should be stacked
plot_data_betterID_bruker$eval<- factor(plot_data_betterID_bruker$eval, levels = c( "prop.wrong.marker.single", "prop.wong.bruker.over.two",
                                                                                    "prop.wong.bruker.under.two", 
                                                                                    "prop.noID.marker", "prop.noID.bruker",
                                                                                    "prop.bruker.correct.under.two", 
                                                                                    "prop.correct.marker.multi", "prop.bruker.correct.multi.over.two", 
                                                                                    "prop.correct.marker.single", "prop.bruker.correct.single.over.two"))
# Plot
f1 <- ggplot(plot_data_betterID_bruker, aes(x=intensity_group, y=as.numeric(as.character(value)), ymin=0,ymax=value)) + facet_grid(Endpoint ~ Group, scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   
f2 <- f1+geom_boxplot(data = plot_data_betterID_bruker[grepl('Score.*',plot_data_betterID_bruker$Endpoint),])
f2b <- f2+geom_boxplot(data = plot_data_betterID_bruker[grepl('Median.*intensity.*',plot_data_betterID_bruker$Endpoint),])
f2c <- f2b+geom_boxplot(data = plot_data_betterID_bruker[grepl('Number.*ribo.*',plot_data_betterID_bruker$Endpoint),])
f2d <- f2c+geom_boxplot(data = plot_data_betterID_bruker[grepl('Measurement.*ppm.*',plot_data_betterID_bruker$Endpoint),])
#f3 <- f2d+geom_boxplot(data = plot_data_betterID_bruker[grepl('Sum.*intensity.*',plot_data_betterID_bruker$Endpoint),])
f4 <- f2d+geom_boxplot(data = plot_data_betterID_bruker[grepl('Fraction.*reproducibly.*',plot_data_betterID_bruker$Endpoint),])
f5 <- f4+geom_col(data = plot_data_betterID_bruker[grepl("Evaluation..*PAPMID", plot_data_betterID_bruker$Endpoint),], aes(fill=eval), width = 0.7) 
f6 <- f5+geom_col(data = plot_data_betterID_bruker[grepl("Evaluation.*microflex.*", plot_data_betterID_bruker$Endpoint),], aes(fill=eval), width = 0.7) 
f7 <- f6+scale_fill_manual('Identification',
                           breaks = c("prop.wong.bruker.over.two","prop.wrong.marker.single","prop.wong.bruker.under.two","prop.noID.bruker","prop.noID.marker","prop.bruker.correct.under.two",      
                                      "prop.bruker.correct.multi.over.two","prop.correct.marker.multi","prop.bruker.correct.single.over.two","prop.correct.marker.single"),
                           values = c("#D55E00","#D55E00","#E69F00","#999999","#999999","#56B4E9", "#0072B2","#0072B2", "#009E73","#009E73"), 
                           labels= c("Wrong identification", "Wrong identification","Wrong identification low score","No identification possibe","No identification possibe","Correct identification low score", 
                                     "Correct multi species identification","Correct multi species identification","Correct single species identification","Correct single species identification")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   
# rename the plot
bruker_better_ID<-f7+ theme(axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=16), strip.text = element_text(size=13), legend.position = "none")
bruker_better_ID

# Extract the legend. Returns a gtable
leg <- get_legend(bruker_better_ID+ theme(legend.position="bottom"))

pdf('./marker_ID_legend_brukerID.pdf', width= 12, height = 1)
as_ggplot(leg)
dev.off()

# Output the figure
pdf('./bruker_better_ID_quality_level_quality_level_all_features_intensity_groups.pdf', width= 8, height = 14)
bruker_better_ID
dev.off()

# Subset spectra which have been acquired on the Axima Confidence and with varying sample preperation protocols
plot_data_betterID_axima<-plot_data_betterID[plot_data_betterID$MALDI == 'Axima Confidence' & !grepl('bruker', plot_data_betterID$Endpoint),]
# Define in which order the groups and how the groups should be displayed
plot_data_betterID_axima$Group<- factor(plot_data_betterID_axima$Group, levels =  c("Burkholderia", "Viridans streptococci", "Enterobacteriaceae"), 
                                        labels = c("Burkholderia\ncepacia complex\n(3 strains)", "viridans\nstreptococci\n(2 strains)", "Enterobacter\ncloace complex\n(3 strains)"))
# Define how the endpoints should be displayed
labels_endpoints_axima_betterID<-c("Number of detected ribosomal marker peaks","Median relative intensity of the ribosomal marker peaks","Fraction of reproducibly detected peaks","Measurement error [ppm]","VitekMS Confidence Score", "Evaluation of the VitekMS species identification", "PAPMID Match Count Score","Evaluation of the PAPMID species identification")
# Introduce line breaks for very long labels
labels_endpoints_axima_betterID<-str_wrap(labels_endpoints_axima_betterID, width = 15)

# Te VitekMS gives no Confidence Score (NA), when the 'Identification Type' is 'No Identification' / 'Not enough peaks'. Set these to 0 (have been identified with 0 confidence)
plot_data_betterID_axima$value<-ifelse(is.na(plot_data_betterID_axima$value) & plot_data_betterID_axima$Endpoint == "VitekMS.DBscore_rank_1", 0, plot_data_betterID_axima$value)


# Rename the endpoints accordingly
plot_data_betterID_axima$Endpoint<- factor(plot_data_betterID_axima$Endpoint, levels = c("n_ribos_detected","median_intensity",  "frac_peaks_repr", "mean_dist_ppm","VitekMS.DBscore_rank_1", "VitekMS", "markerID.match_count", "marker"), labels = labels_endpoints_axima_betterID)



# Define in the order in which the species identification evaluation whould be stacked
plot_data_betterID_axima$eval<- factor(plot_data_betterID_axima$eval, levels = c( "prop.wrong.marker.single", "prop.wrong.VitekMS.single",
                                                                                  "prop.wrong.VitekMS.lowdiscriminatory", 
                                                                                  "prop.noID.marker", "prop.noID.VitekMS",
                                                                                  "prop.correct.VitekMS.lowdiscriminatory", 
                                                                                  "prop.correct.marker.multi", 
                                                                                  "prop.correct.marker.single", "prop.correct.VitekMS.single"))
plot_data_betterID_axima<-plot_data_betterID_axima[!is.na(plot_data_betterID_axima$Group),]

# Plot
f1 <- ggplot(plot_data_betterID_axima, aes(x=intensity_group, y=as.numeric(as.character(value)), ymin=0,ymax=value)) + facet_grid(Endpoint ~ Group, scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('')  
f2 <- f1+geom_boxplot(data = plot_data_betterID_axima[grepl('Score.*',plot_data_betterID_axima$Endpoint),])
f2b <- f2+geom_boxplot(data = plot_data_betterID_axima[grepl('Median.*intensity.*',plot_data_betterID_axima$Endpoint),])
f2c <- f2b+geom_boxplot(data = plot_data_betterID_axima[grepl('Number.*ribo.*',plot_data_betterID_axima$Endpoint),])
f2d <- f2c+geom_boxplot(data = plot_data_betterID_axima[grepl('Measurement.*ppm.*',plot_data_betterID_axima$Endpoint),])
f4 <- f2d+geom_boxplot(data = plot_data_betterID_axima[grepl('Fraction.*reproducibly.*',plot_data_betterID_axima$Endpoint),])
f5 <- f4+geom_col(data = plot_data_betterID_axima[grepl("Evaluation..*PAPMID.*", plot_data_betterID_axima$Endpoint),], aes(fill=eval), width = 0.7) 
f6 <- f5+geom_col(data = plot_data_betterID_axima[grepl("Evaluation.*VitekMS", plot_data_betterID_axima$Endpoint),], aes(fill=eval), width = 0.7) 
f7 <- f6+scale_fill_manual('Identification',
                           breaks = c("prop.wrong.VitekMS.single","prop.wrong.marker.single","prop.wrong.VitekMS.lowdiscriminatory","prop.noID.VitekMS","prop.noID.marker","prop.correct.VitekMS.lowdiscriminatory",      
                                      "prop.correct.marker.multi","prop.correct.VitekMS.single","prop.correct.marker.single"),
                           values = c("#D55E00","#D55E00","#E69F00","#999999","#999999","#0072B2","#0072B2", "#009E73","#009E73"), 
                           labels= c("Wrong identification", "Wrong identification","Wrong identification low score","No identification possibe","No identification possibe","Correct identification low score", 
                                     "Correct multi species identification","Correct single species identification","Correct single species identification")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')   
# Rename the plot
axima_better_ID<-f7+ theme(axis.text=element_text(size=16), axis.title=element_text(size=20),title=element_text(size=16), strip.text = element_text(size=13), legend.position = "right")

axima_better_ID
# Output the figure
pdf('./axima_better_ID_quality_level_quality_level_3_all_features_intensity_groups_intensity_groups.pdf', width= 11, height = 14)
axima_better_ID
dev.off()



# Calibration
# Evaluate how the time after calibration impacts the measurement error
# Read in all calibration eval files
cal_eval_path <- './08_Data/01_Spectra/04_outputs/spectra_test/calibration/'
files <- list.files(cal_eval_path, pattern ='bruker_flex.*sum.*\\.csv')

# Read in all calibration output files into list
cal_eval<-lapply(paste0(cal_eval_path,files), function(i){
  read.csv(i,stringsAsFactors = FALSE, sep=';')
})

# Loop through the list and add 'file column'. the filename contains information on the target, the device and the day after calibration
for (i in 1:length(cal_eval)){
  cal_eval[[i]]['file']<-files[i]
}

# merge all calibration outputs to one
cal_eval <- cal_eval %>% reduce(full_join, by = colnames(cal_eval[[1]]))
# Extract da, MALDI and target information from filename
cal_eval['day']<-gsub('((.*\\-)|(.*\\_))(d(1|2|3|4|5|6|7))((\\-.*)|(\\_.*))','Day \\5',cal_eval$file)
cal_eval['MALDI']<-gsub('((.*\\-)|(.*\\_))(m(1|2))((\\-.*)|(\\_.*))','MALDI \\5',cal_eval$file)
cal_eval['which_target']<-ifelse(grepl('extern-(1|2|3)', cal_eval$file), gsub('((.*\\-)|(.*\\_))(extern-(1|2|3))((\\-.*)|(\\_.*))','external \\5', cal_eval$file), 
                           ifelse(grepl('extern', cal_eval$file), 'external', 'internal'))

cal_eval['target']<-ifelse(grepl('external',cal_eval$which_target), 'external', 'internal')

# Add unique spectra name by concatenating spectra and MALDI
cal_eval['spectra_unique']<-paste(cal_eval$spectra, cal_eval$MALDI)

# Remove empty spectrum 
cal_eval<-cal_eval[cal_eval$spectra_unique != '12-cal-day2-1.0_C2 MALDI 1',]

# subset for the spectra which have been acquired after different times after calibration
jx_cal<- cal_eval[cal_eval$target == 'internal' & !grepl('target', cal_eval$file),]

# Gather all endpoints to one column in order to plot all at once using facet_wrap
jx_cal_long<-gather(jx_cal[,c("spectra_unique","Score1","mean_dist", "mean_dist_ppm")], "Endpoint", "value", -spectra_unique)
jx_cal_long<-merge(jx_cal[,c("spectra_unique","day", "MALDI", "target")], jx_cal_long, by='spectra_unique')

# Convert to numeric
jx_cal_long['value']<-as.numeric(as.character(gsub('\\,', '.', jx_cal_long$value)))

# Plot all days, (exclude Day 6, weird measurement)
ggplot(jx_cal_long[jx_cal_long$Endpoint =='mean_dist' & jx_cal_long$target =='internal' & jx_cal_long$day != 'Day 6',], aes(x = day, y = value, fill = MALDI)) +   
  geom_boxplot() + 
  xlab('') + ylab('mean distance from predicted markers [Da]') +
  labs(fill = "MALDI") + 
  geom_point(position=position_jitterdodge(),alpha=0.3)

# Summaise median and IQR for spectra acquired on different days after calibration on the same target plate which was used for calibration
check_cal<-jx_cal_long[jx_cal_long$target =='internal',] %>% 
  group_by(day, Endpoint) %>%
  summarize(IQR = quantile(value, probs = c(0.25, 0.5, 0.75), na.rm = T))

time.cal.ppm<-ggplot(jx_cal_long[jx_cal_long$Endpoint =='mean_dist_ppm' & jx_cal_long$target =='internal' &jx_cal_long$day %in% c('Day 1', 'Day 7'),], aes(x = day, y = value, fill = MALDI)) +   
  geom_boxplot() +
  xlab('') + ylab('Measurement error [ppm]') +
  labs(fill = "MALDI") +
  theme(legend.title=element_blank()) + scale_fill_brewer(palette="Dark2")
time.cal.ppm

pdf('./time.cal.ppm.pdf', height = 4, width = 5)
time.cal.ppm
dev.off()

# Test whether there is a significant difference in the measurement error in spectra acquired on day 1 and on day 7 after calibration
wilcox.test(jx_cal_long[jx_cal_long$Endpoint =='mean_dist_ppm' & jx_cal_long$target =='internal' & jx_cal_long$day == 'Day 1', 'value'], jx_cal_long[jx_cal_long$Endpoint =='mean_dist_ppm' & jx_cal_long$target =='internal' & jx_cal_long$day == 'Day 7','value'], paired = TRUE)

# Evaluate whether there is a difference in measurement error between the target plates which was used for calibration (internal) and other target plates (external)
mx_d1<-cal_eval[cal_eval$day== 'Day 1' & grepl('target', cal_eval$file),]

#Gather all endpoints to one column in order to plot all at once using facet_wrap
mx_d1_long<-gather(mx_d1[,c("spectra_unique","mean_dist","mean_dist_ppm")], "Endpoint", "value", -spectra_unique)
mx_d1_long<-merge(mx_d1[,c("spectra_unique","day", "MALDI", "target","which_target")], mx_d1_long, by='spectra_unique')

mx_d1_long['value']<-as.numeric(as.character(gsub('\\,', '.', mx_d1_long$value)))

# plot everything
ggplot(mx_d1_long[mx_d1_long$Endpoint =='mean_dist',], aes(x = target, y = value, fill = MALDI)) +   
  geom_boxplot() + 
  xlab('') + ylab('mean distance from predicted markers [Da]') +
  labs(fill = "MALDI") + 
  geom_point(position=position_jitterdodge(),alpha=0.3)

# plot everything
ggplot(mx_d1_long[mx_d1_long$Endpoint =='mean_dist_ppm',], aes(x = target, y = value, fill = MALDI)) +   
  geom_boxplot() + 
  xlab('') + ylab('mean distance from predicted markers [ppm]') +
  labs(fill = "MALDI") + 
  geom_point(position=position_jitterdodge(),alpha=0.3)

# plot everything
ggplot(mx_d1_long[mx_d1_long$Endpoint =='mean_dist',], aes(x = MALDI, y = value, fill = which_target)) +   
  geom_boxplot() + 
  xlab('') + ylab('mean distance from predicted markers [Da]') +
  labs(fill = "target") + 
  geom_point(position=position_jitterdodge(),alpha=0.3)
  
# plot everything
calibration_target <- ggplot(mx_d1_long[mx_d1_long$Endpoint =='mean_dist_ppm',], aes(x = MALDI, y = value, fill = which_target)) +   
  geom_boxplot() + 
  xlab('') + ylab('Mean distance from predicted marker masses [ppm]') +
  labs(fill = "target") +theme(legend.title = element_blank()) + scale_fill_brewer(palette="Accent")
calibration_target

wilcox.test(mx_d1_long[mx_d1_long$Endpoint =='mean_dist_ppm' & mx_d1_long$target == 'internal','value'], mx_d1_long[mx_d1_long$Endpoint =='mean_dist_ppm' & mx_d1_long$target == 'external','value'], paired = F)

# summarise median and interquartile range 
check_cal_target<-mx_d1_long %>% 
  group_by(target, Endpoint) %>%
  summarize(IQR = quantile(value, probs = c(0.25, 0.5, 0.75), na.rm = T))

check_cal_MALDI<-mx_d1_long %>% 
  group_by(MALDI, Endpoint) %>%
  summarize(IQR = quantile(value, probs = c(0.25, 0.5, 0.75), na.rm = T))

cal.figure<-cowplot::plot_grid(time.cal.ppm + theme(legend.position="bottom") + theme(legend.text=element_text(size=13)) + theme(axis.text=element_text(size=14), axis.title=element_text(size=14),title=element_text(size=18), strip.text = element_text(size=14)), 
          calibration_target + theme(legend.position="bottom")  + theme(legend.text=element_text(size=13)) + theme(axis.text=element_text(size=14), axis.title=element_text(size=14),title=element_text(size=18), strip.text = element_text(size=14)), 
          labels = c("A", "B", size=24),
          ncol = 2, nrow = 1)

pdf("./calibration_figure_1-and-7.pdf", width=12, height = 6)
cal.figure
dev.off()
