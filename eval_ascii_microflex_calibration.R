# Aline Cuénod, 2020
# This script summarises the following information for the MALDI-TOF MS spectra, acquired of a microflex Biotyper system to assess the impact of the time after calibration on the measurement precision 
# Total number of peaks, number of phylogenetic marker masses detected (absolut and relativ) (every spectra is exclusively queried for the subunits of the relative genome), average distance of predicted to detected, mass with the highest m/z value, mass at the 90th percentile of m/z values, fraction of peaks > 10'000 Daltons.
# The following input arguments are required:
# (i) directory to the acquired spectra, including spectra as ASCII files and the species identification results by the MALDI_Biotyper database
# (ii) path to the file 'predicted_masses.csv' including all predicted masses and produced by the script 'import_GAP.R'
# (iii) path to the file 'Strains_numbering_first_Shihpment.csv'
# (iv) name and path of the summary output file (1 row per spectra)
# (v) name and file to the to the per subunits output file (one row per detected subunit)

# Load packages
library('MALDIquant')
library('MALDIquantForeign')
library('stringr')
library('dplyr')
library('tidyr')
library('pdftools')

# Define arguments
args = commandArgs(trailingOnly=TRUE)

# Import all brukerreports. These have been saved as PDF file in the same directory as the fid files
brukerreport_dir<-list.files(
  path = args[1],            # directory to search within
  pattern = "\\.pdf$", #
  recursive = TRUE,          #
  full.names = TRUE          #
)

# List all files
brukerreport<-lapply(brukerreport_dir, pdf_text)

# If more than one report was created (because of repetitive measurements), concatenate these
brukerreport[[1]]<-unlist(brukerreport, recursive = FALSE)

# Using regular expressions, extract pages with first two matches
tt<-ifelse(sum(grepl('Organism', brukerreport[[1]]))>1, paste(brukerreport[[1]][grepl('Organism', brukerreport[[1]])], collapse = ''), brukerreport[[1]][grepl('Organism', brukerreport[[1]])])

# Extract which targets these have been measured on
tt<-gsub('\\\n\\s+extern', 'extern', tt)
tt<-gsub('\\\n\\s+intern', 'intern', tt)

# Split text at every new line
tt2<-as.data.frame(strsplit(tt, "\\\n")[[1]][grepl('12-cal',strsplit(tt, "\\\n")[[1]])])
colnames(tt2)<-'text'
tt2['text']<-gsub('^\\s+','', tt2$text)

# Rename columns
tt3<-separate(tt2, text, c('position', 'col', 'Species1', 'Score1', 'Species2', 'Score2'), sep = "\\s{2,}")

# Add 'spectra column'
tt3['spectra']<-paste0(tt3$col, '.0_',tt3$position)

# Import spectra (picked peaks) as ASCIIS
asciis_dir<-list.files(
  path = args[1],            # directory to search within
  pattern = ".*[A-Z]{1}\\d{1,2}.txt$", #
  recursive = TRUE,          #
  full.names = TRUE          #
)

asciis = lapply(asciis_dir, read.delim)



spectra_name<-list()
peaks <- list()
peaks2<-list()
for (i in 1:length(asciis)){
  spectra_name<-as.character(asciis[[i]][4,])
  spectra_name<-gsub('\\\\+','\\.',unlist(gsub('(.*)(12\\-.*)(\\\\1$)','\\2',spectra_name)))
  peaks[[spectra_name]]<-separate(as.data.frame(asciis[[i]][-c(1:6),]), col="asciis[[i]][-c(1:6), ]", into = c("mass", "intensity_raw"), sep=" ") %>% unique() %>%
    mutate(mass =  as.numeric(as.character(mass))) %>%
    mutate(intensity_norm =  as.numeric(intensity_raw)/median(as.numeric(intensity_raw)))
  peaks2[[spectra_name]]<-createMassPeaks(peaks[[spectra_name]]$mass, peaks[[spectra_name]]$intensity_norm)
  }


# Query each spectrum for 
# (i) the total number of peaks, 
# (ii) the peak with the highes m/z value, 
# (iii) the peak with the m/z value at the 90th percentile, 
# (iv) the number of peaks > 10'000 Daltons and 

n_peaks<-list()
sample_name<-list()
highest_peaks<-list()
sample_with_position<-list()
mass_90<-list()
n_high_peaks<-list()
for (i in 1:length(peaks)){
  n_peaks[i]<-length(peaks[[i]]$mass)
  highest_peaks[i]<-ifelse(n_peaks[i] != 0, max(as.numeric(peaks[[i]]$mass)), NA)
  mass_90[i]<-ifelse(n_peaks[i] != 0, quantile(as.numeric(peaks[[i]]$mass),probs = .90), NA)
  n_high_peaks[i]<-ifelse(n_peaks[i] != 0, sum(as.numeric(peaks[[i]]$mass)>10000), NA)
}
# Summarize all to one dataframe
n_peaks<-data.frame(n_peaks = as.character(n_peaks), highest_peaks = as.character(highest_peaks), mass_90 = as.character(mass_90), spectra = names(peaks), n_high_peaks = as.character(n_high_peaks))
# All calibration measurements were of strain Nr. 12
n_peaks['strain_number']<-'12'

# Extract how many phylogenetic marker peaks have been detected, which have been detected and what was the average distance of the predicted to the detected mass
# Read in predicted masses
ribos<-read.csv2(args[2])

# Add strainnumber to file
numbering<-read.csv2(args[3], sep=',')
ribos<-merge(ribos, numbering, by = 'Strain', all.x = T)
ribos<-ribos[ribos$within_mass_range=='TRUE',]
ribos['Numbering_Shipment_1']<-ifelse(nchar(ribos$Numbering_Shipment_1)==1, paste0('0', ribos$Numbering_Shipment_1), ribos$Numbering_Shipment_1)


# set error
error <- 800
ppm <- error / 1000000

# Query each spectrum for the phylogenetic marker peaks predicted from its genome. 
# Create two dataframes: 
# (i)'per_su': here each line is one phylogenetic marker mass. Variables listed cover the spectra in which it was detected, the predicted and the detected mass and the intensity (realtiv)
# (ii) 'sum': here each line summarises one spectrum. Variables listed cover the number of 
#   - phylogenetic markers detected (absolut and normalised by dividing thriugh the number of predicted subunits for that strain), 
#   - their median relative intennsity, 
#   - the mean distance from the predicted to the detected mass (in Daltons and in ppm)
#   - the freaction of peaks which could reproducibly detected in more than half of the technical replicates. 

# Create empty lists
subunit_detected_all_in_one_list<-list()

subunit_detected_all_in_one_each_su<-data.frame()
mass_detected_all_in_one_each_su<-data.frame()
intensity_detected_all_in_one_each_su<-data.frame()
n_predicted_su_all <- list()
frac_peaks_repr<-list()

# Loop through each spectrum and strain
for (strain in 'Escherichia_coli_805237-12'){ # only strain included
  strain_number <- unique(ribos[ribos$Strain == strain, 'Numbering_Shipment_1'])
  peaks_of_interest<-peaks2
  subunit_detected_per_strain<-list()
  n_peak_of_interest <- n_peaks
  for (i in 1:length(peaks_of_interest)){
    subunit_detected<- list()
    mass_detected<- list()
    intensity_detected<- list()

    sample_with_position<-names(peaks_of_interest)[[i]]
    # add % subunits detected (not all strains have the same number of predicted masses)
    n_predicted_su <- length(ribos[ribos$Strain == strain, 'Subunit'])

    for (subunit in ribos[ribos$Strain == strain, 'Subunit']){ # check for the presence of a peak for each subunit (mass +/- error range)
      mass <- as.numeric(as.character(ribos[ribos$Strain == strain & ribos$Subunit == subunit, 'Mass']))
      subunit_temp<-if (any(peaks_of_interest[[i]]@mass<(mass+(mass*ppm)) &  peaks_of_interest[[i]]@mass>(mass-(mass*ppm)))){
        paste(subunit, mass) # If nothing  is detected include NA, so that if non is detected, still someting can be appended
      } else {
        NA
      }
      subunit_detected <- append(subunit_detected, subunit_temp)
      subunit_detected<- unlist(unique(subunit_detected[!is.na(subunit_detected)])) # remove where not detected

      mass_detected_temp <- if (!is.na(subunit_temp)){
        as.character(peaks_of_interest[[i]]@mass[peaks_of_interest[[i]]@mass<(mass+(mass*ppm)) &  peaks_of_interest[[i]]@mass>(mass-(mass*ppm))])
      } else {
        NA
      }
      mass_detected_temp <- paste(subunit, mass_detected_temp)
      mass_detected <- append(mass_detected, mass_detected_temp)
      mass_detected <- mass_detected[!grepl('NA$',mass_detected)] # remove where not detected

      intensity_detected_temp <- if (!is.na(subunit_temp)){
        as.character(peaks_of_interest[[i]]@intensity[peaks_of_interest[[i]]@mass<(mass+(mass*ppm)) &  peaks_of_interest[[i]]@mass>(mass-(mass*ppm))])
      } else {
        NA
      }
      intensity_detected_temp <- paste(subunit, intensity_detected_temp)
      intensity_detected <- append(intensity_detected, intensity_detected_temp)
      intensity_detected <- intensity_detected[!grepl('NA$',intensity_detected)] # remove where not detected

    }
    subunit_detected_all_in_one_list[[sample_with_position]] <- paste(subunit_detected, collapse = ',') # summarise
    #subset into different colonies measured
    if (all(grepl('\\d\\.', n_peak_of_interest$spectra))){
      for (col in unique(str_extract(n_peak_of_interest$spectra, '\\d\\.'))){
        col_peaks<-peaks_of_interest[grepl(col, names(peaks_of_interest))]
        n_peaks_col<-n_peak_of_interest[grepl(col, n_peak_of_interest$spectra),]
        if (sum(grepl(col, n_peak_of_interest$spectra) & as.numeric(as.character(n_peak_of_interest$n_peaks))>0)<2){ # at least two spectra have two have more than 0 peaks, otherwise no filtering of peaks is possible
          frac_peaks_repr[[sample_with_position]]<-0
        }
        else {
          for (j in 1:length(col_peaks)){
            frac_peaks_repr[[as.character(n_peaks_col$spectra[[j]])]]<-length(filterPeaks(binPeaks(col_peaks, method='relaxed', tolerance = ppm), minFrequency=0.51)[[j]]@mass) / length(col_peaks[[j]]@mass)
          }
        }
      }
    } else {
      frac_peaks_repr[[sample_with_position]]<-length(filterPeaks(binPeaks(peaks_of_interest, method='relaxed', tolerance = ppm), minFrequency=0.51)[[i]]@mass) / length(peaks_of_interest[[i]]@mass)
    }
    n_predicted_su_all[strain] <- n_predicted_su
    subunit_detected_all_in_one_each_su <- rbind(subunit_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(subunit_detected), " "))) %>% mutate(spectra = sample_with_position))
    mass_detected_all_in_one_each_su <- rbind(mass_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(mass_detected), " "))) %>% mutate(spectra = sample_with_position))
    intensity_detected_all_in_one_each_su <- rbind(intensity_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(intensity_detected), " "))) %>% mutate(spectra = sample_with_position))
  }
}

# Add columnnames
colnames(subunit_detected_all_in_one_each_su)<-c('Subunit', 'predicted_mass', 'spectra')
colnames(mass_detected_all_in_one_each_su)<-c('Subunit', 'detected_mass', 'spectra')
colnames(intensity_detected_all_in_one_each_su)<-c('Subunit', 'intensity', 'spectra')

# Merge to one 'each_su' dataframe
read_out_each_su<-merge(subunit_detected_all_in_one_each_su, mass_detected_all_in_one_each_su, by=c(intersect(colnames(subunit_detected_all_in_one_each_su), colnames(mass_detected_all_in_one_each_su))))
read_out_each_su<-merge(read_out_each_su, mass_detected_all_in_one_each_su, by=c(intersect(colnames(read_out_each_su), colnames(mass_detected_all_in_one_each_su))))
read_out_each_su<-merge(read_out_each_su, intensity_detected_all_in_one_each_su, by=c(intersect(colnames(read_out_each_su), colnames(intensity_detected_all_in_one_each_su))))

# add distance (measurement error) in Daltons and in ppm
read_out_each_su['dist']<-abs(as.numeric(as.character(read_out_each_su$predicted_mass)) - as.numeric(as.character(read_out_each_su$detected_mass)))
read_out_each_su['dist_ppm']<-(read_out_each_su$dist / as.numeric(as.character(read_out_each_su$predicted_mass))) *1000000

# If multiple peaks per subunit detected, only keep the one with higher intensity
read_out_each_su<-read_out_each_su %>% group_by(spectra,Subunit) %>%
  filter(intensity == max(intensity)) %>%
  filter(dist_ppm == min(dist_ppm))

# Summarise to 'sum' dataframe
read_out_each_su_sum<-read_out_each_su %>% group_by(spectra) %>%
  summarize(mean_dist = mean(as.numeric(as.character(dist))), mean_intensity = mean(as.numeric(as.character(intensity))), mean_dist_ppm = mean(as.numeric(as.character(dist_ppm))))

# Add fraction of reproducibly detected peaks to the 'n_peaks' dataframe
frac_peaks_repr<-do.call(rbind, frac_peaks_repr)
colnames(frac_peaks_repr)<-'frac_peaks_repr'
frac_peaks_repr<- cbind(spectra = rownames(frac_peaks_repr), data.frame(frac_peaks_repr, row.names=NULL))
n_peaks<-merge(n_peaks, frac_peaks_repr, by='spectra')

# Build dataframe from detected subuni data
ribos_detected<-do.call(rbind, subunit_detected_all_in_one_list)
colnames(ribos_detected)<-'subunits_detected'
ribos_detected<- cbind(spectra = rownames(ribos_detected), data.frame(ribos_detected, row.names=NULL))
ribos_detected['strain_number']<-'12'
numbering['Numbering_Shipment_1']<-gsub('(^\\d{1}$)','0\\1', numbering$Numbering_Shipment_1)
ribos_detected<-merge(ribos_detected, numbering, by.x = 'strain_number', by.y = 'Numbering_Shipment_1', all.x = T)

# Count subunits detected there is one more subunit detected that there are commas, so add 1 to count the subunits detected
ribos_detected['n_ribos_detected']<-ifelse(ribos_detected$subunits_detected == "",0,str_count(ribos_detected$subunits_detected, ',') + 1)
n_predicted_su_all<-t(as.data.frame(n_predicted_su_all))
colnames(n_predicted_su_all)<-'n_predicted_su'
rownames(n_predicted_su_all)<-gsub('\\.', '-', rownames(n_predicted_su_all))
rownames(n_predicted_su_all)<-gsub('\\-TS\\-', '\\(TS\\)', rownames(n_predicted_su_all))
ribos_detected<-merge(ribos_detected, n_predicted_su_all, by.x = 'Strain', by.y = "row.names", all=T)
ribos_detected['rel_amount_su_detected']<-ribos_detected$n_ribos_detected / ribos_detected$n_predicted_su

# Merge all
read_out<-merge(n_peaks, ribos_detected, by=intersect(colnames(n_peaks), colnames(ribos_detected)), all.x = T)
read_out<-merge(read_out, read_out_each_su_sum, by = intersect(colnames(read_out), colnames(read_out_each_su_sum)), all.x = T)
read_out<-merge(read_out, tt3, by = intersect(colnames(read_out), colnames(tt3)), all.x = T)

# Write output
write.csv2(read_out, args[4], row.names=FALSE)
write.csv2(read_out_each_su, args[5], row.names=FALSE)
