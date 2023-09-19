library(data.table)
library(plyr)
library(scales)
### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt",header=TRUE)
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)

# Get important categories from my table into peter df
for(i in 1:nrow(petertab)){
  for(j in 1:nrow(sctab)){
    if(!is.element(petertab$`Isolate name`[i],sctab$ID)){
      petertab$heterozygosity[i] <- NA
      petertab$new_eco[i] <- NA
      petertab$basename[i] <- NA
    }
    else if(petertab$`Isolate name`[i] == sctab$ID[j]){
      petertab$heterozygosity[i] <- sctab$heterozygosity[j]
      petertab$new_eco[i] <- sctab$ecology_category[j]
      petertab$basename[i] <- sctab$basename[j]
    }
  }
}
# Remove rows with unassigned clades
petertab <- petertab[complete.cases(petertab$Clades),]
petertab <- petertab[complete.cases(petertab$basename),]

# Add a column with admixture information and remove admixed strains
petertab$admixture <- ifelse(petertab$Clades == "M. Mosaic" |
                            petertab$Clades == "8. Mixed origin" |
                            petertab$Clades == "7. Mosaic beer","Admixed","Non-admixed")
nonadpeter <- petertab[petertab$admixture != "Admixed",]
# Remove non-diploids
diploids <- nonadpeter[nonadpeter$Ploidy == 2,]
# Remove highly heterozygous strains
nonhet <- diploids[diploids$heterozygosity <= 0.003,]

filename <- c()
for(i in 1:length(nonhet$basename)){
  filename[i] <- paste(nonhet$basename[i],nonhet$Char_state[i],sep=" ")
}
filename
write(filename, "/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/char_state.txt",ncolumns = 1)
