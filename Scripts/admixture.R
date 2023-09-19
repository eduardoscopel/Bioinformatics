barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}


data <- read.table("admix_pilot_test.3.Q")
indTable <- read.table("pop3.ind", col.names=c("Sample","Sex","Pop"))
popGroups <- read.table("pop.group.txt", col.names=c("Spec","Broad","Pop"))
mergedAdmixtureTable = cbind(data, indTable)
mergedAdmWithPopGroups = merge(mergedAdmixtureTable, popGroups, by="Pop")
ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$Pop),]
ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$Pop),]
par(mar=c(10,4,4,4))
barplot(t(as.matrix(ordered[,2:5])), col=rainbow(3), border=NA,
        names.arg=barNaming(ordered$Pop), las=2)


data <- read.table("admix_pilot_test.8.Q")
indTable <- read.table("pop_summarized.ind", col.names=c("Sample","Sex","Pop"))
popGroups <- read.table("pop.group.txt", col.names=c("Spec","Pop","PopGroup"))
mergedAdmixtureTable = cbind(data, indTable)
mergedAdmWithPopGroups = merge(mergedAdmixtureTable, popGroups, by="Pop")
ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$Pop),]
ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$Pop),]
barplot(t(as.matrix(subset(ordered, select=V1:V8))), col=rainbow(8), border=NA)
