# Aline Cuénod, 2020
## This file reads in a html file and exports a .csv file summarising the the bruker species identification (one line per spectrum). It requires the folloing two arguments 
## This file reads in a html file and exports a .csv file summarising the the bruker species identification (one line per spectrum). It requires the folloing three arguments 
# (i) input path to the html file
# (ii) .csv file saying whihc run belongs to which repetition
# (iii) Input file which shich strainnumber is which species per NGS.
# (iii) output file to the csv file

library(data.table)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

rawHTML <- paste(readLines(args[1]), collapse = '\n')
rawHTML <- unlist(strsplit(rawHTML,'\\\n\\s*<'))

# Extract, the identifies species, the score, the position and the run from the html, by grepping for the respective patterns
Species_all<-list()
Scores_all<-list()
species.names <- 1:10
scores.names<-1:10
col<-list()
position<-list()
run<-list()
for (i in 1:47){
  Species<-list()
  Score <-list()
  temp<-NULL
  selected_lines<-NULL
  selected_lines<-if(any(grepl(paste0("Mabritec.*",ifelse(nchar(i)==1, paste0('0', i), i),"(\\-|\\_|\\.)*[[:alnum:]]*\\\\0"), rawHTML))){
    which(grepl(paste0("Mabritec.*",ifelse(nchar(i)==1, paste0('0', i), i),"(\\-|\\_|\\.)*[[:alnum:]]*\\\\0"), rawHTML))
  }else{
    NA}
  selected_lines2<-if(any(grepl(paste0("Mabritec.*","BRU\\-", i,'[[:alnum:]]*$'), gsub('\\-1\\-.*','',rawHTML)))){
    which(grepl(paste0("Mabritec.*","BRU\\-", i,'[[:alnum:]]*$'), gsub('\\-1\\-.*','',rawHTML), rawHTML))
  }else{
    NA}
  selected_lines<-union(selected_lines, selected_lines2)
  selected_lines<-selected_lines[!is.na(selected_lines)]
 temp <- rawHTML[selected_lines]
 run<-append(run, gsub('(.*\\M.*\\\\)(2020.{4}\\-.{4})(\\\\.*)(\\\\)(\\d\\_[[:alnum:]]{2,3})(\\\\.*)', '\\2',temp))
 col <- append(col, gsub('(.*\\M.*\\\\2020.{4}\\-.{4}\\\\)(.*)(\\\\)(\\d\\_[[:alnum:]]{2,3})(\\\\.*)', '\\2',temp))
 position <- append(position, gsub('(.*\\M.*\\\\2020.{4}\\-.{4}\\\\)(.*)(\\\\)(\\d\\_[[:alnum:]]{2,3})(\\\\.*)', '\\4',temp))
 for (j in 1:10){
   temp2<-NULL
   temp2 <- rawHTML[selected_lines+23+((j-1)*5)]
   Species[[j]]<-list()
   Species[[j]]<-append(Species[[j]], gsub('(.*>)(.*)(<.*|\\+|\\s$)', '\\2',temp2))
   temp3<-NULL
   temp3 <- rawHTML[selected_lines+24+((j-1)*5)]
   Score[[j]]<-list()
   Score[[j]]<-append(Score[[j]], gsub('(.*>)(.*)(<.*)', '\\2',temp3))
   temp4<-NULL
   temp4 <- rawHTML[selected_lines+23+((j-1)*5)]
   }
 names(Species) <- species.names
 names(Score) <- scores.names
 Species_all[[i]]<-list()
 Species_all[[i]]<-append(Species_all[[i]], Species)
 Scores_all[[i]]<-list()
 Scores_all[[i]]<-append(Scores_all[[i]], Score)
}



# Build two dataframes, one including the identified species, rank and spectra, and the other one including the identification score, rank and spectra 
tt4_scores<-cbind(tt4, Scores_all_df)
tt4_scores_long<- tt4_scores %>% gather("Rank","Score", -c('run', 'col', 'position'))
tt4_species<-cbind(tt4, Species_all_df)
tt4_species_long<- tt4_species %>% gather("Rank","Species", -c('run', 'col', 'position'))

# Merger these
tt4_all<-merge(tt4_scores_long, tt4_species_long, by=intersect(colnames(tt4_species_long), colnames(tt4_scores_long)), all= F)
tt4_all['Score']<-as.numeric(ifelse(tt4_all$Species=='no peaks found', 0, as.character(tt4_all$Score)))
tt4_all<-tt4_all[!duplicated(tt4_all),]

# Add 'spectra' as separate column variable
tt4_all['spectra']<-paste(tt4_all$col, tt4_all$position, sep = '.')


# Add method argument fom 'run' file
run<-read.csv2(args[2], sep=',')
colnames(run)<-c('Method','run')

# Spectra which are labelled as 20200221-1731, should belong to run 20200221-1732, relabel
tt4_all['run']<-ifelse(tt4_all$run == '20200221-1731' & grepl('(0\\_(A|B|C|E|F)(5|6|7|8))',tt4_all$position), '20200221-1732', as.character(tt4_all$run))

tt5<-merge(tt4_all, run, by='run', all.x = T)

# the run 20200428-1650 includes measureents from different tests, specify and remove duplicates
tt5['Method']<-ifelse(tt5$run == '20200428-1650' & grepl('(0\\_C(9|10|11|12))|(0\\_D(7|8|11|12))',tt5$position), 'as25_3',
            ifelse(tt5$run == '20200428-1650' & grepl('(0\\_B(9|10|11|12))',tt5$position), 'smear_3',
                   ifelse(tt5$run == '20200428-1650' & grepl('(0\\_(C|E|F|G|H)(5|6|7|8))|(0\\_A(9|10|11|12))',tt5$position), 'quantity_3', as.character(tt5$Method))))

# Remove duplicates
tt5<-tt5[!duplicated(tt5),]

# Add strainnumber
tt5['strainnumber']<-gsub('\\-1\\-.*', '', tt5$col)
tt5['strainnumber']<-gsub('(.*\\-)(\\d{1,2}$)', '\\2', tt5$strainnumber)

# Add separate column for the identifies species and genus 
tt5['species_identified']<-sapply(strsplit(as.character(tt5$Species)," "), `[`, 2)
tt5['genus_identified']<-sapply(strsplit(as.character(tt5$Species)," "), `[`, 1)

# Count how many different species and genera have been identified with a score higher than 2, and what the difference between highest and second highest species is, if there are more than 1 
tt5_sum<-tt5 %>% 
  group_by(run, spectra) %>%
  mutate(n_species = n_distinct(species_identified)) %>%
  mutate(n_genera = n_distinct(genus_identified)) %>%
  mutate(n_species_over_2 = n_distinct(species_identified[Score > 2])) %>%
  mutate(n_genera_over_2 = n_distinct(genus_identified[Score > 2])) %>%
  group_by(run, spectra, species_identified) %>%
  slice(which.min(Rank)) %>%
  ungroup() %>%
  group_by(run, spectra) %>%
  mutate(diff = as.numeric(sort(as.numeric(Score), decreasing=T)[1] - sort(as.numeric(Score), decreasing=T)[2])) %>%
  filter(Rank == '1')

# Grep for the strainnumber in the spectra name
tt5_sum$strainnumber<-gsub('(BRU\\-)*(\\d{1,2})([[:alpha:]]*)', '\\2', tt5_sum$strainnumber)
tt5_sum$strainnumber<-ifelse(nchar(tt5_sum$strainnumber)==1, paste0('0', tt5_sum$strainnumber), tt5_sum$strainnumber)
 

# Add species as assigned by by NGS
numbering<-read.csv2(args[3], sep=',')
numbering['species_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\3', numbering$Strain)
numbering['genus_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\1', numbering$Strain)
numbering['genus_NGS']<-gsub('Winkia', 'Actinomyces', numbering$genus_NGS)
numbering['Numbering_Shipment_1']<-ifelse(nchar(numbering$Numbering_Shipment_1)==1, paste0('0', numbering$Numbering_Shipment_1), numbering$Numbering_Shipment_1)

tt5_sum<-merge(tt5_sum, numbering, by.x = 'strainnumber', by.y = 'Numbering_Shipment_1')


# Check, whether the microflex Biotyper assigned the correct (same as NGS) species and genus. 
tt5_sum['correct_genus_identified']<-mapply(grepl, tt5_sum$genus_NGS, tt5_sum$genus_identified)
tt5_sum['correct_species_identified']<-mapply(grepl, tt5_sum$species_NGS, tt5_sum$species_identified)

# Add 'include' column. Spectra which are identified as a wrong genus are most likely contaminations, and are excluded from further analyses.
# There are some exeptions to this rule: Shigella / E.coli; Klebsiella / Raoultella, Winkia / Actinomyces, Klebsiella aerogenes / Enterobacter aerogenes. This will applied in downstream scripts. 
tt5_sum['include']<-ifelse(tt5_sum$correct_genus_identified, TRUE, 
                           ifelse(tt5_sum$Score<1.7, TRUE, #include bad spectra, only inlcude if genus from with more than 1.7 confidence
                                  ifelse(tt5_sum$genus_identified == 'Escherichia' & tt5_sum$genus_NGS == 'Shigella', TRUE, FALSE)))

tt5_sum <- apply(tt5_sum,2,as.character)

# Export as CSV 
write.table(tt5_sum,args[4], row.names = F, dec='.', sep=';')
