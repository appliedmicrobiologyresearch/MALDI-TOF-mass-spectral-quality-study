# Aline Cuénod, 2020
# This script compares correctly to incorrectly identified mass spectra. 'Incorrectly identified spectra' include spectra for which the wrong species was assigned of for which no identification was possible). 
# We exclusively include spectra of species which are covered by the MALDI Biotyper, the Vitek MS and the PAPMID database. 

# Load packages
library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')
library('stringr')
library('ggpubr')

setwd('')

# Load data. The 'eval_sum_clean.csv' file includes all species identification of the spectra acquired for this study and was exported by the script 'spectra_comparison.R'
eval_testset<-read.csv('./08_Data/01_Spectra/04_outputs/spectra_test/eval_sum_clean.csv')
eval_testset['frac_high_peaks']<-eval_testset$n_high_peaks / eval_testset$n_peaks

# Set all NA values to 0. This often happens when no ribosomal marker peak / no peak at all was detected
# Do not set the mean distance to 0, if NA. If no marker peak was detected, the emasurement error cannot be assessed (is truly NA)
eval_testset$max_mass_ribo<-ifelse(is.na(eval_testset$max_mass_ribo), 0, as.numeric(as.character(eval_testset$max_mass_ribo)))
eval_testset$n_peaks<-ifelse(is.na(eval_testset$n_peaks), 0, as.numeric(as.character(eval_testset$n_peaks)))
eval_testset$highest_peaks<-ifelse(is.na(eval_testset$highest_peaks), 0, as.numeric(as.character(eval_testset$highest_peaks)))
eval_testset$mass_90<-ifelse(is.na(eval_testset$mass_90), 0, as.numeric(as.character(eval_testset$mass_90)))
eval_testset$n_high_peaks<-ifelse(is.na(eval_testset$n_high_peaks), 0, as.numeric(as.character(eval_testset$n_high_peaks)))
eval_testset$total_intensity_peaks<-ifelse(is.na(eval_testset$total_intensity_peaks), 0, as.numeric(as.character(eval_testset$total_intensity_peaks)))
eval_testset$frac_peaks_repr<-ifelse(is.na(eval_testset$frac_peaks_repr), 0, as.numeric(as.character(eval_testset$frac_peaks_repr)))
eval_testset$n_ribos_detected<-ifelse(is.na(eval_testset$n_ribos_detected), 0, as.numeric(as.character(eval_testset$n_ribos_detected)))
eval_testset$mean_intensity<-ifelse(is.na(eval_testset$mean_intensity), 0, as.numeric(as.character(eval_testset$mean_intensity)))
eval_testset$median_intensity<-ifelse(is.na(eval_testset$median_intensity), 0, as.numeric(as.character(eval_testset$median_intensity)))
eval_testset$frac_high_peaks<-ifelse(is.na(eval_testset$frac_high_peaks), 0, as.numeric(as.character(eval_testset$frac_high_peaks)))

# Gather to long format, with all identification evaluations in 'ID_correct' columns and the 'DB' column indicating the database used 
eval_long_ID<-gather(eval_testset[,c("spectra","markerID.correct_single_species","VitekMS.DBcorrect_species_identified", "brukerDB.correct_species_identified")], "DB", "ID_correct", -spectra)
eval_ID<-merge(eval_testset, eval_long_ID, by = 'spectra', all.x = T)

# Only include strains which are present in all databases (exclude all which are missing in one of the databases)
# When considering VitekMS / microfkex Biotyper DB, exclude the strains for which there is no species entry. The PAPMID DB covers all species included in this study
# Although Shigella are included in the VitekMS DB, exlude these as they are displayed as E. coli
# Strains 3,4,13,14,15,16,25,35,40,41 are missing in the VitekMS database and 3,4,13,14,15,16,18,35 are missingin the MALDI Biotyper database
# eval_ID['included_in_VitekMSDB']<-ifelse(eval_ID$Numbering_Shipment_1 %in% c(3,4,13,14,15,16,25,35,40,41), 0, 1)
# eval_ID['included_in_brukerDB']<-ifelse(eval_ID$Numbering_Shipment_1 %in% c(3,4,13,14,15,16,18,35), 0, 1)
exclude<-union(c(3,4,13,14,15,16,25,35,40,41),c(3,4,13,14,15,16,18,35))
eval_ID<-eval_ID[!eval_ID$Numbering_Shipment_1 %in% exclude,]

# Set all 'noID' to FALSE. 'Incorrectly identified spectra include thereafter not identified and wrongly identified spectra
eval_ID$ID_correct<-ifelse(eval_ID$ID_correct %in% c('noID','no ID'), FALSE, eval_ID$ID_correct)

# Rename 'ID correct'
eval_ID$ID_correct<-ifelse(eval_ID$ID_correct == FALSE, str_wrap('Incorrectly identified spectra', width = 15),eval_ID$ID_correct)
eval_ID$ID_correct<-ifelse(eval_ID$ID_correct == TRUE, str_wrap('Correctly identified spectra', width = 15), eval_ID$ID_correct)

# Exclude NA rows
eval_ID<-eval_ID[!is.na(eval_ID$spectra),]

# Harmonise group labelling
eval_ID$Group<-ifelse(grepl('treptococci', eval_ID$Group), 'Streptococcus', eval_ID$Group)
group_labels_strep_together <- c("Enterobacteriaceae\n(11 strains)", "Listeria\n(2 strains)", "Burkholderia\n(3 strains)", "Bordetella\n(3 strains)",
                                 "Streptococcus\n(8 strains)", "Staphylococcus\n(1 strain)", "Actinobacteria\n(5 strains)", "Gram negative\nAnaerobes\n(3 strains)")

eval_ID$Group<-factor(eval_ID$Group, levels = c("Enterobacteriaceae", "Listeria", "Burkholderia", "Bordetella",
                                                                                      "Streptococcus", "Staphylococcus", "Actinobacteria", "Anaerob_Gram_negative"),
                                         labels = group_labels_strep_together)

# rename DB
eval_ID$DB<-ifelse(grepl('marker',eval_ID$DB), 'PAPMID', 
                   ifelse(grepl('VitekMS', eval_ID$DB), 'VitekMS', 
                          ifelse(grepl('bruker', eval_ID$DB), 'Biotyper',NA)))

# Plot add n_reproducibly detected peaks (absolut)
eval_ID['n_peaks_repr']<-round(eval_ID$frac_peaks_repr* eval_ID$n_peaks)

# Exclude empty spectra
eval_ID_noempty<-eval_ID[eval_ID$n_peaks != 0 ,]
eval_ID_noempty$MALDI<-factor(eval_ID_noempty$MALDI, levels = c("microflex Biotyper","Axima Confidence"))

# Plot the number of reproducibly detected peaks (absolut) agains the total number of peaks
datacount_repr<-ggplot(eval_ID_noempty, aes(x=n_peaks, y=n_peaks_repr)) + facet_grid(Group~MALDI,scales = "free_x") +
  geom_point(alpha=0.1) +
  geom_smooth() +ylab('Number of peaks, which were detected in at least 3/4 of the technical replicates') +
  xlab('Total Number of Peaks')

pdf('./datacount_repr.pdf', width = 6, height = 10)
datacount_repr
dev.off()

# Plot the same and exclusively for spectra acquired on a micoflex Biotyper
datacount_repr_bruker<-ggplot(eval_ID_noempty[eval_ID_noempty$MALDI== 'microflex Biotyper',], aes(x=n_peaks, y=n_peaks_repr)) + facet_grid(.~Group) +
  geom_point(alpha=0.1) +
  geom_smooth() +ylab('Number of peaks, which were detected\nin at least 3/4 of the technical replicates') +
  xlab('Total Number of Peaks')

pdf('./datacount_repr_bruker.pdf', width = 12, height = 4)
datacount_repr_bruker
dev.off()

# Plot the same exlusively for Enterobacteriaceae and microflex spectra
datacount_repr_bruker_entero<-ggplot(eval_ID_noempty[eval_ID_noempty$MALDI== 'microflex Biotyper' & eval_ID_noempty$Group == 'Enterobacteriaceae\n(11 strains)',], aes(x=n_peaks, y=n_peaks_repr)) + facet_grid(.~Group) +
  geom_point(alpha=0.1) +
  geom_smooth() +ylab('Number of peaks, which were detected\nin at least 3/4 of the technical replicates') +
  xlab('Total Number of Peaks')

pdf('./datacount_repr_bruker_entero.pdf', width = 4, height = 4)
datacount_repr_bruker_entero
dev.off()

# Plot the number of reproducibly detected peaks against the sum of the intensity of all peaks
datacount_repr_bruker_intensity<-ggplot(eval_ID_noempty[eval_ID_noempty$MALDI== 'microflex Biotyper',], aes(x=total_intensity_peaks, y=n_peaks_repr)) + facet_grid(.~Group) +
  geom_point(alpha=0.1) +
  geom_smooth() +ylab('Number of peaks, which were detected\nin at least 3/4 of the technical replicates') +
  xlab('Sum of the intenisty of all detected peaks')

pdf('./datacount_repr_bruker_intensity.pdf', width = 12, height = 4)
datacount_repr_bruker_intensity
dev.off()

# Convert to the long format with all feature values in the 'value' column and the 'Endpoint' indicating which feature is displayed
eval_long<-gather(eval_ID_noempty[,c("spectra","MALDI", "DB", "ID_correct","n_peaks","mass_90", "n_ribos_detected", "mean_dist", "median_intensity", "total_intensity_peaks", "frac_high_peaks")], "Endpoint", "value", -c(spectra, MALDI, DB, ID_correct))
test<-compare_means(value~ID_correct, eval_long, method = "wilcox.test", paired = FALSE, group.by = c('MALDI','Endpoint','DB'))

# Compare the the median and the Interquartile Range od spectral features between correctly and incorrectly identified spectra
# Check specta acquired on a mircroflex Biotyper for numbers in the main text
eval_ID_bruker<-eval_ID_noempty[eval_ID_noempty$MALDI == 'microflex Biotyper' & eval_ID_noempty$n_peaks != 0 & !is.na(eval_ID_noempty$ID_correct),]
eval_ID_bruker$DB<-factor(eval_ID_bruker$DB, levels =c('Biotyper', 'PAPMID'))

# Number of ribosomal merker masses detected
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(n_ribos_detected, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra','n_ribos_detected'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra','n_ribos_detected'], paired = F)$p.value

# Median relative intensity of the ribosomal marker masses
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(median_intensity, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra','median_intensity'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra','median_intensity'], paired = F)$p.value

# Total Number of peaks
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(n_peaks, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra','n_peaks'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra','n_peaks'], paired = F)$p.value

# Sum of the intensity of all detected peaks
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(total_intensity_peaks, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra','total_intensity_peaks'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra','total_intensity_peaks'], paired = F)$p.value

# Measurement error
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(mean_dist_ppm, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra','mean_dist_ppm'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra','mean_dist_ppm'], paired = F)$p.value

# Fraction of high peaks
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(frac_high_peaks, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra' & eval_ID_bruker$DB == 'Biotyper','frac_high_peaks'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra' & eval_ID_bruker$DB == 'Biotyper','frac_high_peaks'], paired = F)$p.value

# Mass at the 90th percentile
eval_ID_bruker %>% group_by(ID_correct) %>% summarize(IQR = quantile(mass_90, probs = c(0.25, 0.5, 0.75), na.rm = T))
wilcox.test(eval_ID_bruker[eval_ID_bruker$ID_correct == 'Correctly\nidentified\nspectra' & eval_ID_bruker$DB == 'Biotyper','mass_90'], 
            + eval_ID_bruker[eval_ID_bruker$ID_correct == 'Incorrectly\nidentified\nspectra' & eval_ID_bruker$DB == 'Biotyper','mass_90'], paired = F)$p.value

# In order to plot the intensities, pseudolog transform them. Use log10 and replace tge '-Inf' values (resulting from values == 0) to 0
eval_ID_noempty$total_intensity_peaks<-log10(eval_ID_noempty$total_intensity_peaks)
eval_ID_noempty$total_intensity_peaks<-ifelse(eval_ID_noempty$total_intensity_peaks == '-Inf', 0, eval_ID_noempty$total_intensity_peaks)
eval_ID_noempty$median_intensity<-log10(eval_ID_noempty$median_intensity)
eval_ID_noempty$median_intensity<-ifelse(eval_ID_noempty$median_intensity == '-Inf', 0, eval_ID_noempty$median_intensity)


# Plot the mass spectral quality feature for each MALDI device and phylogenetic group separately, comparing correct to incorrectly identified spectra
# First plot spectra acquired on an Axima Confidence
eval_ID_axi<-eval_ID_noempty[eval_ID_noempty$MALDI == 'Axima Confidence' & eval_ID_noempty$n_peaks != 0 & !is.na(eval_ID_noempty$ID_correct),]
# Define in which order the databases should be plottet
eval_ID_axi$DB<-factor(eval_ID_axi$DB, levels =c('VitekMS', 'PAPMID'))

# Number of ribososmal marker peaks
axi_n_ribos<-ggplot(eval_ID_axi, aes(x=Group, y=n_ribos_detected, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Number of detected\nribosomal marker peaks') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Median relative intensity of the robsosmla marker peaks
axi_median_ribo_int<-ggplot(eval_ID_axi, aes(x=Group, y=median_intensity, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Median relative intensity of\nribosomal marker peaks (log10)') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Use an unpaired wilcoxon rank test to find significance level
wilcox.test(eval_ID_axi[eval_ID_axi$ID_correct == 'Correctly\nidentified\nspectra' & eval_ID_axi$DB =='VitekMS','n_ribos_detected'], 
            eval_ID_axi[eval_ID_axi$ID_correct == 'Incorrectly\nidentified\nspectra' & eval_ID_axi$DB =='VitekMS','n_ribos_detected'])

# Total number of peaks detected
axi_n_peaks<-ggplot(eval_ID_axi, aes(x=Group, y=n_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Total number of peaks') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442")) 

# Fraction of peaks > 10'000
axi_frac_high_peaks<-ggplot(eval_ID_axi, aes(x=Group, y=frac_high_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Fraction of peaks > 10,000 Daltons') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'bottom') + scale_fill_manual(values = c("#0072B2", "#F0E442"), name = '') 

# Sum of the Intensity of all peaks
axi_total_int<-ggplot(eval_ID_axi, aes(x=Group, y=total_intensity_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Sum of the intensity\nof all detected peaks (log10)') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Measurement error (here no phylogenetic group specific effect is expected, plot for all groups together)
axi_median_dist<-ggplot(eval_ID_axi, aes(x=ID_correct, y=mean_dist_ppm, fill = ID_correct)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Mean measurement error [ppm]') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# m/z value at the 90th percentile
axi_mass_90<-ggplot(eval_ID_axi, aes(x=Group, y=mass_90, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('m/z value at the 90th percentile') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Plot mass spectral quality features for spectra acquired on a microflex Biotyper
eval_ID_bruker<-eval_ID_noempty[eval_ID_noempty$MALDI == 'microflex Biotyper' & eval_ID_noempty$n_peaks != 0 & !is.na(eval_ID_noempty$ID_correct),]
# Define the order in which the databases whould be displayed
eval_ID_bruker$DB<-factor(eval_ID_bruker$DB, levels =c('Biotyper', 'PAPMID'))

# Number of ribosomal marker peaks detected
bruker_n_ribos<-ggplot(eval_ID_bruker, aes(x=Group, y=n_ribos_detected, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Number of detected\nribosomal marker peaks') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Median relative intensity of these
bruker_median_ribo_int<-ggplot(eval_ID_bruker, aes(x=Group, y=median_intensity, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Median relative intensity of\nribosomal marker peaks (log10)') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Total number of peaks
bruker_n_peaks<-ggplot(eval_ID_bruker, aes(x=Group, y=n_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Total number of peaks') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Fraction of high peaks
bruker_frac_high_peaks<-ggplot(eval_ID_bruker, aes(x=Group, y=frac_high_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Fraction of peaks > 10,000 Daltons') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'bottom') + scale_fill_manual(values = c("#0072B2", "#F0E442"), name = '')  

# Sum of the intensity of all detected peaks
bruker_total_int<-ggplot(eval_ID_bruker, aes(x=Group, y=total_intensity_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Sum of the intensity\nof all detected peaks (log10)') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Measurement error (here no phylogenetic group specific effect is expected, plot for all groups together)
bruker_median_dist<-ggplot(eval_ID_bruker, aes(x=ID_correct, y=mean_dist_ppm, fill = ID_correct)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Mean measurement error [ppm]') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# m/z value at the 90th percentile
bruker_mass_90<-ggplot(eval_ID_bruker, aes(x=Group, y=mass_90, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(DB~.) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('m/z value at the 90th percentile') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Plot the measurement error for the proker spectra
pdf('./which.measurement.error.bruker.pdf', width=3, height = 3)
bruker_median_dist
dev.off()

# Combine all mass specrtal features for spectra acquired on the microflex Biotyper 
combined_plots<-cowplot::plot_grid(bruker_n_ribos, bruker_median_ribo_int, 
                                   bruker_total_int, 
                                   bruker_n_peaks, bruker_mass_90, 
                                   bruker_frac_high_peaks, 
                                   ncol = 1, align = "v",axis = "b")
pdf('./which.Endpoints.all.bruker.per.group.pdf', width=9, height = 21)
combined_plots
dev.off()

# Combine all mass specrtal features for spectra acquired on the Axima Confidence
combined_plots<-cowplot::plot_grid(axi_n_ribos, axi_median_ribo_int, 
                                   axi_total_int, 
                                   axi_n_peaks, axi_mass_90, 
                                   axi_frac_high_peaks, 
                                   ncol = 1, align = "v",axis = "b")
pdf('./which.Endpoints.all.axi.per.group.pdf', width=9, height = 21)
combined_plots
dev.off()


# For the main text figure, exclusively plot the comparison between correctly and incorrectly identified spectra by the microflex Biotyper database
eval_ID_bruker_biotyper<-eval_ID_bruker[eval_ID_bruker$DB == 'Biotyper',]

bruker_n_ribos_biotyper<-ggplot(eval_ID_bruker_biotyper, aes(x=Group, y=n_ribos_detected, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Number of detected\nribosomal marker peaks') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

bruker_total_int_biotyper<-ggplot(eval_ID_bruker_biotyper, aes(x=Group, y=total_intensity_peaks, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  stat_compare_means(method = "wilcox.test", na.rm = T, paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Sum of the intensity\nof all detected peaks (log10)') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

bruker_median_ribo_int_biotyper<-ggplot(eval_ID_bruker_biotyper, aes(x=Group, y=median_intensity, fill= factor(ID_correct))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Median relative intensity\nof ribosomal\nmarker peaks (log10)') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

bruker_median_dist_biotyper<-ggplot(eval_ID_bruker_biotyper, aes(x=ID_correct, y=mean_dist_ppm, fill = ID_correct)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(label = ifelse(p < 0.0001,..p.signif.., ..p.format..)), vjust = 1) +
  xlab('') + ylab('Mean measurement\nerror [ppm]') +
  theme_bw() + theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), legend.position = 'none') + scale_fill_manual(values = c("#0072B2", "#F0E442"))

# Combine these Biotyper plota and export these
combined_plots<-cowplot::plot_grid(bruker_n_ribos_biotyper,
                                   bruker_median_ribo_int_biotyper,
                                   bruker_total_int_biotyper,
                                   ncol = 1, align = "v",axis = "b")
pdf('./which.Endpoints.biotyper.group.pdf', width=9, height = 9)
combined_plots
dev.off()

# Export measurement error for specrta identified by the biotyper
pdf('./which.measurement.error.bruker.biotyper.pdf', width=3, height = 3)
bruker_median_dist_biotyper
dev.off()

