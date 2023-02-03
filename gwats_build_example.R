####################################################################################################
# Brachy GWATS Builder
# Author: Jared Streich
# Created: 20-01-06
# Version 0.1.0
####################################################################################################


####################################################################################################
####################################### Load Packages ##############################################
####################################################################################################

# install.packages("seqGDS2VCF")
# library(seqGDS2VCF)
library(raster)

####################################################################################################
########################################### Start ##################################################
####################################################################################################

##### Set working directory
setwd("~/Desktop/projects/BrachyMTGenomes/")

SNP_Markers <- read.table("Brachy_mt_Markers_Bsylv_RefmtG_sampNames.ped")
SNP_Map <- read.table("Brachy_mt_Markers_Bsylv_RefmtG_sampNames.map")
rownames(SNP_Markers) <- SNP_Markers[,1]
SNP_Markers <- SNP_Markers[c(T,F)]
colnames(SNP_Markers) <- as.character(paste(SNP_Map[,2], SNP_Map[,4], sep = "_"))
SNP_Markers <- SNP_Markers[,7:ncol(SNP_Markers)]
head(SNP_Markers)

pca.b <- cmdscale(dist(SNP_Markers), k = 19)
rownames(pca.b) <- rownames(SNP_Markers)
plot(pca.b[,1], pca.b[,2], cex = 0.7, col = "aquamarine3", pch = 16)
text(pca.b[,1], pca.b[,2], label = rownames(pca.b), cex = 0.7, col = "grey30")


##### Set working directory
setwd("~/Desktop/projects/gwats/Brachypodium/")


##### Load Table
# md <- read.delim("SraRunTable.txt")

##### Convert to character data type
# md <- cbind(paste("BVZ_", as.character(md[,7]), ".fastq.gz.fasta", sep = ""), 
# 			paste("/home/ju0/Brachy/dist-mt-build/", as.character(md[,11]), 
# 			"_mtG.fasta", sep = ""), as.character(md[,12]), as.character(md[,14]))

##### Get list of commands to move fasta files
# write.table(cbind("mv", md[,1:2]), "",  quote = F, col.names = F, row.names = F)
# samps <- read.delim("")

coord <- read.delim("Brachy_META_NO_NAs.txt")
coord <- cbind(as.character(coord[,1]), coord[,5], coord[,6])
coord.r <- coord[,1]
coord <- cbind(as.numeric(coord[,2]), as.numeric(coord[,3]))
rownames(coord) <- coord.r
colnames(coord)[1:2] <- c("Latitude", "Longitude")
head(coord)
dim(coord)
match(rownames(coord), rownames(swil.1.plink))
input.points <- coord


SNP_Markers <- read.table("Brachy_479samps_36567vars_2020-01-13.ped")
SNP_Map <- read.table("Brachy_479samps_36567vars_2020-01-13.map")
rownames(SNP_Markers) <- SNP_Markers[,1]
# SNP_Markers[SNP_Markers == "A"] <- 0
# SNP_Markers[SNP_Markers == "T"] <- 1
SNP_Markers[1:10,1:8]
# head(SNP_Map)
# head(SNP_Markers)
dim(SNP_Markers)

SNP_Markers <- SNP_Markers[,7:ncol(SNP_Markers)]
SNP_Markers <- SNP_Markers[c(T,F)]
colnames(SNP_Markers) <- as.character(SNP_Map[,2])
Pheno_data <- coord
match(rownames(Pheno_data), rownames(SNP_Markers))

Pheno_data <- nm(Pheno_data, SNP_Markers)
SNP_Markers <- nm(SNP_Markers, Pheno_data)
SNP_Markers <- order_Rnames(SNP_Markers, Pheno_data)

Pheno_data <- Pheno_data[complete.cases(Pheno_data[,1]), ]
dim(Pheno_data)

library(rworldmap)
newmap <- getMap(resolution = "high")
plot(newmap, xlim = c(-10, 55), ylim = c(28, 48), asp = 1, 
	interior = T, col = "grey30", lwd = 0.1)
points(cbind(Pheno_data[,2], Pheno_data[,1]), cex = 1.5, pch = 16, col = "dodgerblue3")

tot <- Pheno_data[duplicated(Pheno_data[,2]), ]
nrow(tot)
####################################################################################################
######################################### Variable names ###########################################
####################################################################################################

# Increments between middle of months
inc <- c(31, 29, 30, 30, 31, 30, 31, 31, 30, 31, 30, 31)


# Set Directory to raster layers
setwd("/Volumes/Smithers/Rasters/")

# Give directory a variable name to be read in a loop
# !!!!!!!!!!!!! Check month order for normal or computer order !!!!!!!!!!!!!!!!!!!!!
txtfiles = list.files(pattern='*tif')
txtfiles

# Give directory a variable name to be read in a loop
txtfiles = list.files(pattern='swc_fr_*')
txtfiles

# Column name prefix
pheno.colheader <- "swc_fr_"
pheno.colheader

# Output file name
outfile <- paste(pheno.colheader,"GEMMA_day.txt", sep = "")
outfile

# Plink Name
plinkfile <- paste(pheno.colheader,"plink.txt", sep = "")
plinkfile

# Keep file for gemma
keepfile <- paste(pheno.colheader, "namesKeep.txt", sep = "")
keepfile



####################################################################################################
############# Create raw list of mlonthly phenotype data for day imputation ########################
####################################################################################################
input.points <- input.points[complete.cases(input.points[1,]), ]
ClimData <- input.points
dim(ClimData)

# Start Loop throug set and extract monthly values
for(i in 1:length(txtfiles)){
	# print name of current text file
	print(txtfiles[i])

	# Read in individual raster layers
	biolayer <- raster(txtfiles[i])
	
	# Extract BioClim Values per point	
	BioCol <- extract(biolayer, input.points[,2:1])

	# Paste New Column on Climate Data to ClimData Variable
	ClimData <- cbind(ClimData, BioCol)
	j <- i+2
	colnames(ClimData)[j] <- txtfiles[i]
	print(dim(ClimData))
	print(i)
	print(head(ClimData))
}# End i loop for all layers

head(ClimData)
min(ClimData[,1])
dim(ClimData)
####################################################################################################
####### Check to see the order of column extraction for human or computer numeric order ############
####################################################################################################

# Human numeric order of variables: 1,2,3,4,5,6,7,8,9,10,11,12
layer <- cbind(ClimData[,1:2], clim0 = ClimData[,14], ClimData[,3:14], clim13 = ClimData[,3])
layer <- layer[complete.cases(layer[,5]), ]
min(layer)
dim(layer)
head(layer)
layer <- layer[,3:ncol(layer)]
min(layer)

zer_adj <- function(x){
	if(min(x) <= 0){
		correction <- 0 - min(x)
		correction <- correction + 1
		x <- x + correction
		print("zero correction made")
		print(correction)
		}
	else{
		correction <- 0
		print("All Values above zero")
		}
	return(x)
	}
# Check adjustment
layer.ad <- zer_adj(layer)
min(layer.ad)
dim(layer.ad)
##### Plus Minus days for scan
plsmin <- 10

##### Create Interpolation Loop
for(l in 1:nrow(layer)){
	x <- as.matrix(layer[l,])
	for(i in 2:13){
		em <- x[i-1] # Earlier Month
		lm <- x[i]   # Later month
		if(i == 2){ # If Jan, get Dec Data
			em.0 <- x[13]
			lm.0 <- x[2]
			mms <- seq(em.0, lm.0, length.out = inc[i-1])
		}
		else{
			mmn <- seq(em, lm, length.out = inc[i-1])
			mms <- c(mms, mmn)
		}
	}
	mms <- c(mms, mms)
	for(j in 1:2){
		for(i in 1:(length(mms)/2)){
			ri <- mean(c(mms[i], mms[i+5], mms[i+10], mms[i+15]))
			if(i == 1){
				ri.p <- ri
			}
			else{
				ri.p <- c(ri.p, ri)
			}	
		}
		mms <- c(ri.p, ri.p)
	}	
	layer.n <- mms[16:(380+plsmin)]
	# plot(layer.n, cex = 0.5, pch = 16, col = "Green1")
	if(l == 1){
		layer.p <- layer.n
	}
	else{
		layer.p <- cbind(layer.p, layer.n)
	}
}

dim(layer.p)
head(layer.p)
layer.p <- layer.p - min(layer.ad)

# Observe image sorted by middle of year
layer.t <- t(layer.p)
dim(layer.t)
dim(coord)
dim(input.points)
rownames(layer.t) <- rownames(layer.ad)
layer.t <- layer.t[order(layer.t[,183]), ]
# cols <- palette(rainbow(200))
# cols <- cols[(length(cols)-40):1]

# Image of year sorted at middle date
# image(t(layer.t), col = cols)
# abline(v = 0.5, lwd = 2)
# rownames(layer.t) <- rownames(coord)
layer.t <- layer.t[,1:365]

# Create Phenotype data
setwd("~/Desktop/projects/gwats/Brachypodium/")
outfile
write.table(layer.t, file = outfile, sep = "\t", quote = F, row.names = F, col.names = F)

# write keep file
write.table(rownames(layer), file = keepfile, quote = F, col.names = F, row.names = F)

##### Create Plink based pheno files
plink.file <- cbind(rownames(layer.t), rownames(layer.t), rep(0, times = nrow(layer.t)), 
	rep(0, times = nrow(layer.t)), rep(0, times = nrow(layer.t)), layer.t)

plink.file <- plink.file[,1:370]
dim(plink.file)
colnames(plink.file) <- c("FID", "IID", "MATID", "PATID", "SEX", paste("Day_", c(1:365), sep = "") )
plink.file[1:10,1:10]
write.table(plink.file, plinkfile, quote = F, col.names = T, row.names = F)
