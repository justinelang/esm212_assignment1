#analyze ESM grassland community data
#species accumulation curves and similarity
####
#install and load package vegan
	install.packages("vegan")
	library(vegan)
	library(tidyr)
	library(dplyr)
####
#read presence-absence data
	grasspa <- read.csv("grassland_dat_pa.csv",header=TRUE)
	padat <- grasspa[,4:40]		#plot x species P/A matrix in format for vegan analysis
####
#read cover data
	grasscov <- read.csv("grassland_dat_cover.csv",header=TRUE)
	covdat <- grasscov[,4:40]	#plot x species cover matrix in format for vegan analysis
####
#create exact species accumulation curves for 4 plots using 10 subplots
	accum1 <- specaccum(padat[1:10,]) 	#plot 1, subplots 1-10
	accum2 <- specaccum(padat[13:22,])	#plot 2, subplots 1-10
	accum3 <- specaccum(padat[25:34,])	#plot 3, subplots 1-10
	accum4 <- specaccum(padat[37:46,])	#plot 4, subplots 1-10
##
#plot accumulation curves for each plot
	plot(accum1,ylim=c(0,20),xlab="subplots",ylab="species")
		lines(accum2,col="red")
		lines(accum3,col="blue")
		lines(accum4,col="green")
	 	text(10.1,9,"V")	#valley plot
	 	text(10.1,14,"N") #North plot
	 	text(10.1,16,"S") #South plot
	 	text(10.1,18,"N2") #second north plot
####
#plot single accumulation curve using 40 subplots
	accum_allsubplots <- specaccum(padat[c(1:10,13:22,25:34,37:46),])
	plot(accum_allsubplots,xlab="subplots",ylab="species")
#####
#####
#####
#calculate Bray-Curtis dissimilarity between plots
#First, aggregate subplot data
	p1 <- apply(padat[1:12,],2,max)
	p2 <- apply(padat[13:24,],2,max)
	p3 <- apply(padat[25:36,],2,max)
	p4 <- apply(padat[37:48,],2,max)
# create new matrix to calculate pairwise plot dissimilarity
	plots_spp <- rbind(p1,p2,p3,p4)
#produce dissimilarity matrix
	plots_bray_diss <- vegdist(plots_spp,method="bray")
#####
#####
#####
#Analyze cover data
#Calculate area-weighted cover, total area sampled = 10 *0.25 + 2*2.5 = 7.5 sq m
	p1_cov <- (2.5/7.5)*apply(covdat[1:10,],2,mean) + 
			(5/7.5)*apply(covdat[11:12,],2,mean)
	p2_cov <- (2.5/7.5)*apply(covdat[13:22,],2,mean) + 
			(5/7.5)*apply(covdat[23:24,],2,mean)
	p3_cov <- (2.5/7.5)*apply(covdat[25:34,],2,mean) + 
			(5/7.5)*apply(covdat[35:36,],2,mean)
	p4_cov <- (2.5/7.5)*apply(covdat[37:46,],2,mean) + 
			(5/7.5)*apply(covdat[47:48,],2,mean)
#
#####
#####
#####
#Analyze plot dissimilarity based on cover
	plots_cov <- rbind(p1_cov,p2_cov,p3_cov,p4_cov)
	plots_euclid_diss <- vegdist(plots_cov,method="euclidean")
	plots_hclust <- hclust(plots_euclid_diss,method="average")		#cluster analysis of plots (contrived)
		plot(plots_hclust)
	plots_hclust_pa <- hclust(plots_bray_diss,method="average")
		plot(plots_hclust_pa)
##
#####
#####
#####
#Analyze plot diversity
	plots_shannon <- diversity(plots_cov,index="shannon") 	#shannon H'
	plots_simpson <- diversity(plots_cov,index="simpson")	#simpson diversity
	plots_richness <- specnumber(plots_cov)					#species richness
	
	