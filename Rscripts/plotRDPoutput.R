plotRDPmain <- function(filename)
{
	x <- read.table(filename, sep=',', header=TRUE);
	annt <- read.table('/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/EBB_Illumina/data.htseq.org/FHCRC_Hullar/FCD05FJ/IlluminaPE/sample.annotation.txt', sep='\t', header=TRUE);
	conf_cutoff <- 0.5;
	plotRDPphylum(x, conf_cutoff, annt, filename);
	plotRDPclass(x, conf_cutoff, annt, filename);
	plotRDPorder(x, conf_cutoff, annt, filename);
	plotRDPfamily(x, conf_cutoff, annt, filename);
	plotRDPgenus(x, conf_cutoff, annt, filename);

}

plotRDPphylum <- function(x, conf_cutoff, annt, filename)
{
	phylums <- as.vector(unique(x[x$phylumConf>=conf_cutoff,]$phylum));
	phylums <- phylums[!is.na(phylums)]; # we treat NA separately
	n <- length(phylums);
	CC <- rep(0, n+1);
	for(i in 1:n) { 
		CC[i] <- sum(x[!is.na(x$phylum)&x$phylum==phylums[i]&x$phylumConf>=conf_cutoff, ]$count); 
	}
	# add Unclassified (which actually includes NA)
	CC[n+1] <- sum(x$count) - sum(CC[1:n]);
	phylums <- c(phylums, "Unclassified");
	
	Y <- data.frame(names=phylums, counts=CC);
	o <- order(Y$counts, decreasing=TRUE);
	
	# write out genus data frame
	output_filename <- paste(filename, '.phylum', sep='');
	write.table(Y[o,], output_filename, sep='\t', row.names=FALSE);
	
	# make counts into percentages
	Y$counts <- Y$counts*100 / sum(Y$counts);
	
	# plot into png
	png(paste(output_filename,'.png',sep=''));
	par(oma=c(0,5,0,0), las=2);
	barplot(Y[o,]$counts, names.arg=Y[o,]$names, horiz=TRUE, xlim=c(0,100), xlab='Abundance (%)');
	id <- substr(filename, 3, 7);
	cur_annt <- annt[annt$Sample==id,];
	title(paste("Sample ", id, " Fat:", cur_annt$Fat, " (", cur_annt$Name, ") phyla", sep=""));
	dev.off();
}

plotRDPclass <- function(x, conf_cutoff, annt, filename)
{
	classs <- as.vector(unique(x[x$classConf>=conf_cutoff,]$class));
	classs <- classs[!is.na(classs)]; # we treat NA separately
	n <- length(classs);
	CC <- rep(0, n+1);
	for(i in 1:n) { 
		CC[i] <- sum(x[!is.na(x$class)&x$class==classs[i]&x$classConf>=conf_cutoff, ]$count); 
	}
	# add Unclassified (which actually includes NA)
	CC[n+1] <- sum(x$count) - sum(CC[1:n]);
	classs <- c(classs, "Unclassified");
	
	Y <- data.frame(names=classs, counts=CC);
	o <- order(Y$counts, decreasing=TRUE);
	
	# write out genus data frame
	output_filename <- paste(filename, '.class', sep='');
	write.table(Y[o,], output_filename, sep='\t', row.names=FALSE);
	
	# make counts into percentages
	Y$counts <- Y$counts*100 / sum(Y$counts);
	
	# plot into png
	png(paste(output_filename,'.png',sep=''));
	par(oma=c(0,5,0,0), las=2);
	barplot(Y[o,]$counts, names.arg=Y[o,]$names, horiz=TRUE, xlim=c(0,100), xlab='Abundance (%)');
	id <- substr(filename, 3, 7);
	cur_annt <- annt[annt$Sample==id,];
	title(paste("Sample ", id, " Fat:", cur_annt$Fat, " (", cur_annt$Name, ") class", sep=""));
	dev.off();
}

plotRDPorder <- function(x, conf_cutoff, annt, filename)
{
	orders <- as.vector(unique(x[x$orderConf>=conf_cutoff,]$order));
	orders <- orders[!is.na(orders)]; # we treat NA separately
	n <- length(orders);
	CC <- rep(0, n+1);
	for(i in 1:n) { 
		CC[i] <- sum(x[!is.na(x$order)&x$order==orders[i]&x$orderConf>=conf_cutoff, ]$count); 
	}
	# add Unclassified (which actually includes NA)
	CC[n+1] <- sum(x$count) - sum(CC[1:n]);
	orders <- c(orders, "Unclassified");
	
	Y <- data.frame(names=orders, counts=CC);
	o <- order(Y$counts, decreasing=TRUE);
	
	# write out genus data frame
	output_filename <- paste(filename, '.order', sep='');
	write.table(Y[o,], output_filename, sep='\t', row.names=FALSE);
	
	# make counts into percentages
	Y$counts <- Y$counts*100 / sum(Y$counts);
	
	# plot into png
	png(paste(output_filename,'.png',sep=''));
	par(oma=c(0,5,0,0), las=2);
	barplot(Y[o[1:10],]$counts, names.arg=Y[o[1:10],]$names, horiz=TRUE, xlim=c(0,100), xlab='Abundance (%)');
	id <- substr(filename, 3, 7);
	cur_annt <- annt[annt$Sample==id,];
	title(paste("Sample ", id, " Fat:", cur_annt$Fat, " (", cur_annt$Name, ") order top 10", sep=""));
	dev.off();
}

plotRDPfamily <- function(x, conf_cutoff, annt, filename)
{
	familys <- as.vector(unique(x[x$familyConf>=conf_cutoff,]$family));
	familys <- familys[!is.na(familys)]; # we treat NA separately
	n <- length(familys);
	CC <- rep(0, n+1);
	for(i in 1:n) { 
		CC[i] <- sum(x[!is.na(x$family)&x$family==familys[i]&x$familyConf>=conf_cutoff, ]$count); 
	}
	# add Unclassified (which actually includes NA)
	CC[n+1] <- sum(x$count) - sum(CC[1:n]);
	familys <- c(familys, "Unclassified");
	
	Y <- data.frame(names=familys, counts=CC);
	o <- order(Y$counts, decreasing=TRUE);
	
	# write out genus data frame
	output_filename <- paste(filename, '.family', sep='');
	write.table(Y[o,], output_filename, sep='\t', row.names=FALSE);
	
	# make counts into percentages
	Y$counts <- Y$counts*100 / sum(Y$counts);
	
	# plot into png
	png(paste(output_filename,'.png',sep=''));
	par(oma=c(0,5,0,0), las=2);
	barplot(Y[o[1:10],]$counts, names.arg=Y[o[1:10],]$names, horiz=TRUE, xlim=c(0,100), xlab='Abundance (%)');
	id <- substr(filename, 3, 7);
	cur_annt <- annt[annt$Sample==id,];
	title(paste("Sample ", id, " Fat:", cur_annt$Fat, " (", cur_annt$Name, ") family top 10", sep=""));
	dev.off();
}
	
plotRDPgenus <- function(x, conf_cutoff, annt, filename)
{
	genuss <- as.vector(unique(x[x$genusConf>=conf_cutoff,]$genus));
	genuss <- genuss[!is.na(genuss)]; # we treat NA separately
	n <- length(genuss);
	CC <- rep(0, n+1);
	for(i in 1:n) { 
		CC[i] <- sum(x[!is.na(x$genus)&x$genus==genuss[i]&x$genusConf>=conf_cutoff, ]$count); 
	}
	# add Unclassified (which actually includes NA)
	CC[n+1] <- sum(x$count) - sum(CC[1:n]);
	genuss <- c(genuss, "Unclassified");
	
	Y <- data.frame(names=genuss, counts=CC);
	o <- order(Y$counts, decreasing=TRUE);
	
	# write out genus data frame
	output_filename <- paste(filename, '.genus', sep='');
	write.table(Y[o,], output_filename, sep='\t', row.names=FALSE);
	
	# make counts into percentages
	Y$counts <- Y$counts*100 / sum(Y$counts);
	
	# plot into png
	png(paste(output_filename,'.png',sep=''));
	par(oma=c(0,5,0,0), las=2);
	# only listing top 10 for now!!!!
	barplot(Y[o[1:10],]$counts, names.arg=Y[o[1:10],]$names, horiz=TRUE, xlim=c(0,100), xlab='Abundance (%)');
	id <- substr(filename, 3, 7);
	cur_annt <- annt[annt$Sample==id,];
	title(paste("Sample ", id, " Fat:", cur_annt$Fat, " (", cur_annt$Name, ") genus top 10", sep=""));	
	dev.off();
}

filename <- commandArgs(TRUE);
print(paste("processing RDP table output", filename, "......", sep=''));
plotRDPmain(filename);

