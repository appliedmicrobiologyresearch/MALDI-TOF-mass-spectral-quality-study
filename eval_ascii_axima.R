# Aline Cuénod, 2020
# This script summarises the following information for the MALDI-TOF MS spectra, acquired of a Axima Confidence system system: 
# Total number of peaks, number of phylogenetic marker masses detected (absolut and relativ) (every spectra is exclusively queried for the subunits of the relative genome), average distance of predicted to detected, mass with the highest m/z value, mass at the 90th percentile of m/z values, fraction of peaks > 10'000 Daltons.
# The following input arguments are required:
# (i) path to the file called 'positions.csv', a .csv file including the information which sample was measured on which run and position
# (ii) directory to the ascii files (unique sample name) eg.'./ASCII-Files-picked-peaks/Axima_Confidence/' # spectra were extracted as 'mzXml' files and converted to ASCII files. The '.mzXml files' exclusively cover m/z and 'intensity' values assigned by the Shimadzu Launchpad Software
# (iii) path to the file 'predicted_masses.csv' including all predicted masses and produced by the script 'import_GAP.R'
# (iv) path to the file 'Strains_numbering_first_Shihpment.csv'
# (v) name and path of the summary output file (1 row per spectra)
# (vi) name and file to the to the per subunits output file (one row per detected subunit)
# When running the scrip, a warning messege will be evoged through MALDIQuant which is triggered by emtpy spectra. This can be ignored. 

# Load packages
library('MALDIquant')
library('MALDIquantForeign')
library('stringr')
library('dplyr')
library('tidyr')

# Define arguments
args = commandArgs(trailingOnly=TRUE)

# Import 'positions' file indicating which spectrum belongs to which strain
positions<-read.csv(args[1], sep=',')
positions['spectra_name']<-gsub('-','_',tolower(positions$File))
positions['spectra_name']<-gsub('\\.txt*','',positions$spectra_name)

# Import all ascii files
asciis_dir<-list.files(
  path = args[2],        # directory to search within
  pattern = "*.*txt$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)

asciis = lapply(asciis_dir, read.delim)


# Normalise the intensity by deviding throuogh the median
# Create 'MassPeaks' object
spectra_name<-list()
peaks <- list()
peaks2<-list()
for (i in 1:length(asciis)){
  spectra_name<-as.character(asciis[[i]][1,])
  spectra_name<-gsub('.*Sample\\=','',spectra_name)
  spectra_name<-gsub('#.*', '', spectra_name)
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
# (v) the total intensity of all peaks.
n_peaks<-list()
sample_name<-list()
highest_peaks<-list()
sample_with_position<-list()
mass_90<-list()
n_high_peaks<-list()
total_intensity_peaks<-list()
for (i in 1:length(peaks)){
  n_peaks[i]<-length(peaks2[[i]]@mass)
  highest_peaks[i]<-ifelse(n_peaks[i] != 0, max(peaks2[[i]]@mass), NA)
  mass_90[i]<-ifelse(n_peaks[i] != 0, quantile(as.numeric(peaks2[[i]]@mass),probs = .90), NA)
  n_high_peaks[i]<-sum(peaks2[[i]]@mass>10000)
  total_intensity_peaks[i]<-sum(as.numeric(peaks[[i]]$intensity_raw))
}

# Summarize all to one dataframe
n_peaks<-data.frame(n_peaks = as.character(n_peaks), highest_peaks = as.character(highest_peaks), mass_90 = as.character(mass_90), spectra = names(peaks), n_high_peaks = as.character(n_high_peaks), total_intensity_peaks = as.character(total_intensity_peaks))

# Extract spectra name
n_peaks['spectra_name']<-gsub('-','_',tolower(n_peaks$spectra))
n_peaks['spectra_name']<-gsub('\\#.*','',n_peaks$spectra_name)
n_peaks['spectra_name']<-gsub('(.*)([[:alpha:]]{4,5}\\_[[:alnum:]]+\\_\\d{4}\\_[[:alnum:]]{2,4})(.*)', '\\2', n_peaks$spectra_name)
n_peaks['spectra_name']<-gsub('^smear\\_2\\_\\_', '', n_peaks$spectra_name)

# Add Position
n_peaks<-merge(n_peaks, positions, by='spectra_name')
n_peaks$File<-NULL

# Extract Strainnumber
n_peaks['strain_number']<-ifelse(grepl('(\\d+.*)(\\-\\d+)(\\-\\d+)', n_peaks$Sample), gsub('(\\d+.*)(\\-\\d+)(\\-\\d+)', '\\1',n_peaks$Sample), str_extract(n_peaks$Sample, '\\d{1,2}[[:alnum:]]*$'))
n_peaks['strain_number']<-str_extract(n_peaks$strain_number, '\\d{1,2}')
n_peaks['strain_number']<-ifelse(nchar(n_peaks$strain_number)==1, paste0('0', n_peaks$strain_number), n_peaks$strain_number)
n_peaks['run']<-gsub('(.*)(\\_\\d[[:alpha:]]\\d$)', '\\1', n_peaks$spectra_name)

# Exclude calibration spectra
n_peaks<-n_peaks[!grepl('REF|Ref', n_peaks$Sample),]
peaks2<-peaks2[!grepl('REF|Ref', n_peaks$Sample)]

# Extract how many phylogenetic marker peaks have been detected, which have been detected and what was the average distance of the predicted to the detected mass
# Read in predicted masses
ribos<-read.csv2(args[3], sep=';')

# Add strainnumber to file
numbering<-read.csv2(args[4], sep=',')
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
frac_peaks_repr_2<-list()

# Loop through each spectrum and strain

for (strain in unique(ribos$Strain)){ # Loop through all strains 
  strain_number <- unique(ribos[ribos$Strain == strain, 'Numbering_Shipment_1'])
  for (run in unique(n_peaks$run)){  # loop through all runs too. This is important for the 'fraction of reproducibly detected peaks' which should only be calculated for the technical replicates of the same experiment. 
    peaks_of_interest<-if (any(n_peaks$strain_number == strain_number & n_peaks$run == run)){ #only compare the spectra to the predicted ribos of that strain # add if statement, if strain is missing, no error
      peaks2[names(peaks2) %in% n_peaks[n_peaks$strain_number == strain_number & n_peaks$run == run,'spectra']] 
    } else {
      next
    }
    n_peak_of_interest<-if (any(n_peaks$strain_number == strain_number & n_peaks$run == run)){ #only compare the spectra to the predicted ribos of that strain # add if statement, if strain is missing, no error
      n_peaks[n_peaks$strain_number == strain_number & n_peaks$run == run,] 
    } else {
      next
    }
    subunit_detected_per_strain<-list()
    
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
      #subset into different dillutions (only for qnt experimnet)
      if (all(grepl('1\\-[0-9]+', n_peak_of_interest$spectra))){
        for (conc in unique(str_extract(n_peak_of_interest$spectra, '1\\-[0-9]+'))){
          conc_peaks<-peaks_of_interest[grepl(conc, names(peaks_of_interest))]
          n_peaks_conc<-n_peak_of_interest[grepl(conc, n_peak_of_interest$spectra),]
          if (sum(grepl(conc, n_peak_of_interest$spectra) & as.numeric(as.character(n_peak_of_interest$n_peaks))>0)<2){ # at least two spectra have two have more than 0 peaks, otherwise no filtering of peaks is possible
            frac_peaks_repr[[sample_with_position]]<-0
            frac_peaks_repr_2[[sample_with_position]]<-0
          }
          else {
            for (j in 1:length(conc_peaks)){
              if(length(conc_peaks[[j]]@mass)> 0){
                frac_peaks_repr[[as.character(n_peaks_conc$spectra[[j]])]]<-length(filterPeaks(binPeaks(conc_peaks, method='relaxed', tolerance = ppm), minFrequency=0.51)[[j]]@mass) / length(conc_peaks[[j]]@mass)
              } else {
                frac_peaks_repr[[as.character(n_peaks_conc$spectra[[j]])]]<-0
              }
              if(length(conc_peaks[[j]]@mass)> 0 & length(conc_peaks[[length(conc_peaks)+1-j]]@mass) > 0){# read out the reproducibility of two random spectra of the same technical replicate
                frac_peaks_repr_2[[as.character(n_peaks_conc$spectra[[j]])]]<-length(filterPeaks(binPeaks(conc_peaks[c(j, length(conc_peaks)+1-j)], method='relaxed', tolerance = ppm), minFrequency=0.51)[[1]]@mass) / length(conc_peaks[[j]]@mass)
              } else {
                frac_peaks_repr_2[[as.character(n_peaks_conc$spectra[[j]])]]<-0
              }
            }
          }
        }
      } else {
        if(length(peaks_of_interest[[i]]@mass)> 0){
          frac_peaks_repr[[sample_with_position]]<-length(filterPeaks(binPeaks(peaks_of_interest, method='relaxed', tolerance = ppm), minFrequency=0.51)[[i]]@mass) / length(peaks_of_interest[[i]]@mass)
        } else {
          frac_peaks_repr[[sample_with_position]]<-0
        }
        if(length(peaks_of_interest[[i]]@mass) > 0 & length(peaks_of_interest[[length(peaks_of_interest)+1 -i]])){
          frac_peaks_repr_2[[sample_with_position]]<-length(filterPeaks(binPeaks(peaks_of_interest[c(i, length(peaks_of_interest)+1-i)], method='relaxed', tolerance = ppm), minFrequency=0.51)[[1]]@mass) / length(peaks_of_interest[[i]]@mass)
        }
        else{
          frac_peaks_repr_2[[sample_with_position]]<-0
        }
      }
      n_predicted_su_all[strain] <- n_predicted_su
      subunit_detected_all_in_one_each_su <- rbind(subunit_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(subunit_detected), " "))) %>% mutate(spectra = sample_with_position))
      mass_detected_all_in_one_each_su <- rbind(mass_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(mass_detected), " "))) %>% mutate(spectra = sample_with_position))
      intensity_detected_all_in_one_each_su <- rbind(intensity_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(intensity_detected), " "))) %>% mutate(spectra = sample_with_position))
    }
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

# add distnace (measurement error) in Daltons and in ppm
read_out_each_su['dist']<-abs(as.numeric(as.character(read_out_each_su$predicted_mass)) - as.numeric(as.character(read_out_each_su$detected_mass)))
read_out_each_su['dist_ppm']<-(read_out_each_su$dist / as.numeric(as.character(read_out_each_su$predicted_mass))) *1000000

# If multiple peaks per subunit detected, only keep the one with higher intensity (relative to the mean)
read_out_each_su<-read_out_each_su %>% group_by(spectra,Subunit) %>%
  filter(intensity == max(intensity)) %>%
  filter(dist_ppm == min(dist_ppm))

# Summarise to 'sum' dataframe
read_out_each_su_sum<-read_out_each_su %>% group_by(spectra) %>%
  summarize(mean_dist = mean(as.numeric(as.character(dist))), mean_dist_ppm = mean(as.numeric(as.character(dist_ppm))), mean_intensity = mean(as.numeric(as.character(intensity))), median_intensity = median(as.numeric(as.character(intensity))), max_mass_ribo = max(as.numeric(as.character(predicted_mass))))

# Add fraction of reproducibly detected peaks to the 'n_peaks' dataframe
frac_peaks_repr_df<-do.call(rbind, frac_peaks_repr)
colnames(frac_peaks_repr_df)<-'frac_peaks_repr'
frac_peaks_repr_df<- cbind(spectra = rownames(frac_peaks_repr_df), data.frame(frac_peaks_repr_df, row.names=NULL))
frac_peaks_repr_df_2spectra<-do.call(rbind, frac_peaks_repr_2)
colnames(frac_peaks_repr_df_2spectra)<-'frac_peaks_repr_2spectra'
frac_peaks_repr_df_2spectra<- cbind(spectra = rownames(frac_peaks_repr_df_2spectra), data.frame(frac_peaks_repr_df_2spectra, row.names=NULL))
frac_peaks_repr_df<- merge(frac_peaks_repr_df, frac_peaks_repr_df_2spectra, by = 'spectra')

n_peaks<-merge(n_peaks, frac_peaks_repr_df, by='spectra')

# Build dataframe from detected subuni data
ribos_detected<-do.call(rbind, subunit_detected_all_in_one_list)
colnames(ribos_detected)<-'subunits_detected'
ribos_detected<- cbind(spectra = rownames(ribos_detected), data.frame(ribos_detected, row.names=NULL))
ribos_detected<-merge(ribos_detected, n_peaks[, c('spectra', 'strain_number')], by='spectra')

# Add strainnumber
numbering['Numbering_Shipment_1']<-gsub('(^\\d{1}$)','0\\1', numbering$Numbering_Shipment_1)
ribos_detected<-merge(ribos_detected, numbering, by.x = 'strain_number', by.y = 'Numbering_Shipment_1', all = T)

# Count subunits detected there is one more subunit detected that there are commas, so add 1 to count the subunits detected
ribos_detected['n_ribos_detected']<-ifelse(grepl('\\d',ribos_detected$subunits_detected), str_count(ribos_detected$subunits_detected, ',') + 1, 0)
n_predicted_su_all<-t(as.data.frame(n_predicted_su_all))
colnames(n_predicted_su_all)<-'n_predicted_su'
rownames(n_predicted_su_all)<-gsub('\\.', '-', rownames(n_predicted_su_all))
rownames(n_predicted_su_all)<-gsub('\\-TS\\-', '\\(TS\\)', rownames(n_predicted_su_all))
ribos_detected<-merge(ribos_detected, n_predicted_su_all, by.x = 'Strain', by.y = "row.names", all=T)
ribos_detected['rel_amount_su_detected']<-ribos_detected$n_ribos_detected / ribos_detected$n_predicted_su

# Merge all
read_out<-merge(n_peaks, ribos_detected, by=intersect(colnames(n_peaks), colnames(ribos_detected)), all.x = T)
read_out<-merge(read_out, read_out_each_su_sum, by = intersect(colnames(read_out), colnames(read_out_each_su_sum)), all.x = T)

# Write output
write.table(read_out, args[5], row.names=FALSE, dec='.', quote = F, sep = ';')
write.table(read_out_each_su, args[6], row.names=FALSE, dec='.', quote = F, sep = ';')
