# Aline Cuénod, 2020
# This script summarises the marker based identifications, which have been performed per group
# The following arguments are required: 
# (i) path to the common dir of the different output files (e.g. './Species_identification/PAPMID/')
# (ii) path to the file 'Strains_numbering_first_Shihpment.csv'
# (iii) output path to the summary csv file

# Load packages
library('tidyverse')
library('dplyr')

# Define arguments
args = commandArgs(trailingOnly=TRUE)

# Import all PAPMID output files
files = list.files(path = args[1],        # directory to search within
  pattern = ".*export.*csv$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)

PAPMID_export = lapply(files, read.csv)  # read all the matching files

# Convert all columns to chracter
for (i in 1:length(PAPMID_export)){
  PAPMID_export[[i]][]<-lapply(PAPMID_export[[i]], as.character)
}

# Merge all PAPMID outputs
PAPMID_export <- PAPMID_export %>% reduce(full_join, by = colnames(PAPMID_export[[1]]))

# Keep only 1 row per spectra, filter for the one with the highest MatchCount Score
PAPMID_export_filtered<-as.data.frame(PAPMID_export %>% 
  group_by(Sample) %>% filter(DataCount == max(DataCount)) %>% ungroup())
  
# Rename Winkia neuii -> Actinomyces neii in order to compare
PAPMID_export_filtered['Genus']<-ifelse(PAPMID_export_filtered$Genus == 'Winkia' & grepl('neuii', PAPMID_export_filtered$Species), 'Actinomyces', PAPMID_export_filtered$Genus)

# In order to evaluate the species identification result, create 'species' column, not including subspecies info
PAPMID_export_filtered['species']<-gsub('\\-*subsp.*', '',PAPMID_export_filtered$Species)


# Count number of species and genus with highest match, add which species have been identified per sample and keep one entry per spectrum
PAPMID_export_filtered <- PAPMID_export_filtered %>%
  group_by(Sample) %>%
  mutate(n_species_identified = n_distinct(paste(Genus, species, sep=' '))) %>% 
  mutate(n_genera_identified = n_distinct(Genus)) %>%
  mutate(species_identified = paste(unique(paste(Genus, species, sep=' ')), collapse='-')) %>%
  mutate(genus_identified = paste(unique(Genus), collapse='-')) %>%
  filter(row_number() == 1)

# Remove empty or redundant columns columns
PAPMID_export_filtered <- PAPMID_export_filtered[!PAPMID_export_filtered$Sample=='',c("Sample","Run","species_identified","genus_identified","DataCount","X.Subunit.Mass.","n_species_identified","n_genera_identified")]

# Some 'qnt' files have been analysed with PAPMID and classifer, keep only classifier results 
PAPMID_export_filtered<- PAPMID_export_filtered[!grepl('BRU-(17|18|39|40)-1-125', PAPMID_export_filtered$Sample),]

# Import all classifier outputs. (Same approach, but subtyping module)
# List all files
files = list.files(path = args[1], 
                   pattern = ".*Identification Report.*csv$", # regex pattern, some explanation below
                   recursive = TRUE,          # search subdirectories
                   full.names = TRUE          # return the full path
)

# Define delimiter
delimiter <- ','

# Read all files
output = lapply(files, read.csv, sep = delimiter)  # read all the matching files

# Convert all columns to chracter and add genus
for (i in 1:length(output)){
  colnames(output[[i]])<-c("X","spectra","sample","match_count","species_identified","profiles_matched")
  output[[i]]['Genus']<-ifelse(strsplit(files[i], '\\/')[[1]][13] != "", strsplit(files[i], '\\/')[[1]][13], strsplit(files[i], '\\/')[[1]][14])
  output[[i]][]<-lapply(output[[i]], as.character)
}


# Merge all classifier outputs
output <- output %>% reduce(full_join, by = colnames(output[[1]]))
output$X<-NULL
output$spectra<-NULL

# some 'quantity' spectra have wrongly been analysed with the strep viridans classifier, although they belong to other genera, remove these
output<-output[!(output$Genus == 'S_viridans' & grepl('qty', output$sample)),]

# Add 'Genus' column
output['Genus']<-ifelse(output$Genus == 'S_aureus_complex', 'Staphylococcus', 
                        ifelse(output$Genus == 'S_viridans', 'Streptococcus', output$Genus))

# Count how many species and genera have been identified
output['species_identified']<-gsub('-like', '.like', output$species_identified)
output['species_identified']<-gsub('(like\\-)(\\d{1,2})', 'like.\\2', output$species_identified)
output['n_species_identified']<-str_count(output$species_identified, pattern = "-")+1
output['n_genera_identified']<-1

# Merge PAPMID and Classifier Outputs
output<-merge(output, PAPMID_export_filtered, by.x=c('sample','match_count', 'species_identified', 'Genus','n_species_identified','n_genera_identified'), by.y = c('Sample', 'DataCount', 'species_identified', 'genus_identified','n_species_identified','n_genera_identified'), all=TRUE)

# Add strainnumber
if (grepl('Spectratest', args[1])){
  output['sample2']<-gsub('(.*AXI\\-[[:digit:]]{1,2})(\\-wh)(\\-1\\-.*)', '\\1\\3', output$sample)# remove '-wh'
  output['sample2']<-gsub('(.*AXI\\-\\d{1,2}[[:alnum:]]*)(\\-\\d\\-\\d{1,5})(\\_mabr.*)', '\\1\\3', output$sample2) # remove indication of quantity
  pattern<- ifelse(any(grepl('BRU', output$sample)), '(.*BRU\\-(WP1\\-)*)(\\d{1,2})(.*)', '(.*)(\\-|\\_)(\\d{1,2}[[:alnum:]]*)(\\_mabr.*)') # two versions of pattern: one for buker and one for axima normal
  output['strainnumber']<-gsub(pattern, '\\3', output$sample2)
  output['strainnumber']<-gsub('[[:alpha:]]', '', output$strainnumber)
  
  output['strainnumber']<-ifelse(nchar(output$strainnumber)==1, paste0('0',output$strainnumber), output$strainnumber)
  output['sample2']<-NULL
  unique(output$strainnumber)
  } else{ 
    if (grepl('Participating_labs', args[1])) {
      output['sample']<-gsub('St\\.Gallen', 'St\\-Gallen', output$sample)
      output['strainnumber']<-gsub('([^\\.]+\\.)([[:digit:]]{1,2})(.*)', '\\2', output$sample)
    }
    }
    
output['strainnumber']<-ifelse(nchar(output$strainnumber)==1, as.character(paste0('0', output$strainnumber)), as.character(output$strainnumber))

# Remove entries with weird strainnumber
output<-output[!nchar(output$strainnumber)>2,]

# Remove 'contaminated'
output<-output[!grepl('conta', output$sample),]

# In order to evaluate whether the identification is correct, add specied by NGS
numbering<-read.csv2(args[2], sep=',')
numbering['species_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\3', numbering$Strain)
numbering['genus_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\1', numbering$Strain)
numbering['genus_NGS']<-gsub('Winkia', 'Actinomyces', numbering$genus_NGS)
numbering['Numbering_Shipment_1']<-ifelse(nchar(numbering$Numbering_Shipment_1)==1, paste0('0',numbering$Numbering_Shipment_1), numbering$Numbering_Shipment_1)

# Add strainnumber
output<-merge(output, numbering, by.x = 'strainnumber', by.y = 'Numbering_Shipment_1')


# Check if correct genus amongst identified genera
output['correct_genus_identified']<-mapply(grepl, output$genus_NGS, output$Genus)
output['correct_species_identified']<-mapply(grepl, output$species_NGS, output$species_identified)

write.csv(output, args[3], row.names = F)