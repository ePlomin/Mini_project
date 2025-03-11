################# Clean Orca Script ###############

# Prepartations ####

## Libs & Cleaning####
rm(list = ls())
library(hierfstat)
library(pegas)
library(adegenet)
library(vcfR)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

## Data import ####
vcf <- read.vcfR("C:/Users/fejar/Desktop/Uni/Göteborg/Conservation and the genetics of populations/Orca Project/R/Input/orca_data_3323.vcf")

# Create genind object 
data <- vcfR2genind(vcf)

## Add population & ecotype information:
#Read in Csv
pop_data <- read.csv("C:/Users/fejar/Desktop/Uni/Göteborg/Conservation and the genetics of populations/Orca Project/R/Input/orca_popinfo.csv", header = TRUE, stringsAsFactors = TRUE, sep = ";")

#Assign pops
pop(data) <- pop_data$pop_id

# Save Ecotype information
data@other$ecotype <- pop_data$ecotype

#Check if worked
pop(data) 
data@other$ecotype

# Data Handling ####

## Cleaning up dataset ####
#We have one wrongly assigned Individual that we want to change
data@other$ecotype[23] <- "O"

#And we wanna change the name of the Ecotypes for nicer visualization
data@other$ecotype <- as.character(data@other$ecotype)
data@other$ecotype[data@other$ecotype == "R"] <- "Resident"
data@other$ecotype[data@other$ecotype == "T"] <- "Transient"
data@other$ecotype[data@other$ecotype == "O"] <- "Offshore"
data@other$ecotype[data@other$ecotype == "A"] <- "Antarctic"
data@other$ecotype[is.na(data@other$ecotype)] <- "Icelandic"
data@other$ecotype <- factor(data@other$ecotype)
table(data@other$ecotype)


## Splitting ####
# We want to create 2 seperate datasets: one for FST-outlier loci assumed to be under selection and one without them were we assume all loci to behave neutrally

### Outlier detection ####
## Lets have a look at FST-outlier loci
#Calculate FST
stat <- basic.stats(genind2hierfstat(data))$perloc
# Plot results
ggplot(stat,aes(Fst))+ geom_histogram(fill="steelblue4") + 
  geom_vline(xintercept = quantile(stat$Fst,0.95),
             color="forestgreen",linetype ="dashed") +
  geom_vline(xintercept = quantile(stat$Fst,0.99),
             color="red",linetype ="dashed") + theme_bw()
# Everything beyond the redline (= beyond the 99% quantile) we consider FST-outliers likely under selection

### Outlier filtering & Dataset creation####
#Pool Outliers
outliers <- stat %>%
  rownames_to_column('locus') %>%  
  filter(Fst>quantile(Fst,0.99))
#Some name changes for downstream analysis
locus_key <- data.frame(loc1=locNames(data),loc2=colnames(genind2hierfstat(data))[-1])

#Dataset of only outliers (a = adaptive)
data_a <- data[loc = filter(locus_key,loc2 %in% outliers$locus)$loc1]

#Dataset without outliers (n = neutral)
`%notin%` <- Negate(`%in%`) 
data_n <- data[loc = filter(locus_key,loc2 %notin% outliers$locus)$loc1]



# Structure Analysis ####

## Neutral Structure ####
#First let's have a look at the dataset without the outlier loci. Since we filtered the loci we assumed to be under selection, all variation left in this dataset should only be due to Drift & Gene flow 

### PCA ####
x.orca_n <- tab(data_n, freq=TRUE, NA.method="mean")
pca.orca_n <- dudi.pca(x.orca_n, center=TRUE, scale=FALSE, scannf = FALSE, nf =4)
s.class(pca.orca_n$li, fac=pop(data_n), col=transp(funky(15),.6), axesel=FALSE, cstar=0, cpoint=3, clabel = 0.6, grid = FALSE)
  add.scatter.eig(pca.orca_n$eig[1:30],4,1,2, ratio=.3)
  title(xlab = "PC1: 20.9 %", ylab = "PC2: 6.6 %")

# As we can see most variation is explained by the first 4 PCs, to better visualize differences let's now turn to DAPC

### DAPC ####
# As the first 4 Pcs explain most variation (see PCA) we only retain those for downstream analysis
  
grp_n <- find.clusters(data_n, n.pca=4, max.n=9,scale=FALSE)
# Best-fit BIC values vary a lot here, somewhere between 5 & 9, however 5 Cluster gives the most consistent results, as at all other numbers of clusters individual clusters can change in each run (even with 9 clusters!), but at 5 the clusters always stay the same. Also 5 is the first value after the steepest drop-off in BIC-values (mostly) => We decided to go with 5 Clusters

dapc2 <- dapc(data_n, pop=grp_n$grp, n.pca=4)
# All 4 can be retained

#First DAPC
scatter(dapc2, xax=1, yax=2, grp=dapc2$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")),
        pch=16, bg="white",cstar = 1, cellipse = 1,clabel = 0.7)

#Now we can use the DAPC clusters on the PCA again:
s.class(pca.orca_n$li, fac=dapc2$grp, clab=1, col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")), csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
# We see that each cluster neatly encompasses an ecotype! For better visualization we can therefore now do this:

# We create a data frame where assigned cluster is listed against ecotype
cluster_assignments <- dapc2$grp
ecotype_info <- data@other$ecotype
cluster_to_ecotype_df <- data.frame(
  Cluster = cluster_assignments,
  Ecotype = ecotype_info
)

# We then can assign each cluster the name of the most frequent ecotype, which is the only ecotype that occurs in it as we've seen above
ecotype_summary <- table(cluster_to_ecotype_df$Cluster, cluster_to_ecotype_df$Ecotype)
most_frequent_ecotypes <- apply(ecotype_summary, 1, function(x) names(x)[which.max(x)])

# We can then replace the labels for the clusters and create a better understandable graph for the report
dapc2$grp <- factor(dapc2$grp, labels = most_frequent_ecotypes)
scatter(dapc2, xax=1, yax=2, grp=dapc2$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")),
        pch=16, bg="white",cstar = 0, cellipse = 2,clabel = 0,scree.da=TRUE,
        scree.pca=TRUE, leg=TRUE)

# And the same with the PCA with DAPC clusters
s.class(pca.orca_n$li, fac=dapc2$grp, clab=0.7, col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")), csta=1, cpoint=2, cellipse =1, xax=1, yax=2, grid = FALSE)
add.scatter.eig(pca.orca_n$eig[1:30],4,1,2, ratio=.3)
title(xlab = "PC1: 20.9 %", ylab = "PC2: 6.6 %")

## Adaptive Structure ####
# So here we only use the 34 loci with suspiciously high FST which we assume to be under selection, otherwise we proceed analogous to the neutral sructure

### PCA ####
x.orca_a <- tab(data_a, freq=TRUE, NA.method="mean")
pca.orca_a <- dudi.pca(x.orca_a, center=TRUE, scale=FALSE, scannf = FALSE, nf =4)
s.class(pca.orca_a$li, fac=pop(data_a), col=transp(funky(15),.6), axesel=FALSE, cstar=0, cpoint=2, clabel = 0.6, grid = FALSE)
add.scatter.eig(pca.orca_a$eig[1:10],3,1,2, ratio=.3)
title(xlab = "PC1: 69.8 %", ylab = "PC2: 10.4 %")
# Here most variation seems to be explained by the first 3 PCs! to better visualize differences let's now turn to DAPC

### DAPC ####
# As the first 3 Pcs explain most variation (see PCA) we only retain those for downstream analysis n.pca

grp_a <- find.clusters(data_a, n.pca=3, max.n=9,scale=FALSE)
# Even more difficult to say than for neutral structure. Because of consistency (see Neutral Structure < DAPC) we decided to go with 5 clusters again, where can als the slope starting to level off.

dapc3 <- dapc(data_a, pop=grp_a$grp, n.pca=3)
# All 3 can be retained

#First DAPC
scatter(dapc3, xax=1, yax=2, grp=dapc3$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")),
        pch=16, bg="white",cstar = 1, cellipse = 1,clabel = 0.7)

#Now we can use the DAPC clusters on the PCA again:
s.class(pca.orca_a$li, fac=dapc3$grp, clab=1, col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")), csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
# We see that each cluster neatly encompasses an ecotype # except 1 individual of one transient individual that becomes antarctic! For better visualization we can still do this?:

# We create a data frame where assigned cluster is listed against ecotype
cluster_assignments <- dapc3$grp
ecotype_info <- data@other$ecotype
cluster_to_ecotype_df <- data.frame(
  Cluster = cluster_assignments,
  Ecotype = ecotype_info
)

# We then can assign each cluster the name of the most frequent ecotype, which is the only ecotype that occurs in it as we've seen above
ecotype_summary <- table(cluster_to_ecotype_df$Cluster, cluster_to_ecotype_df$Ecotype)
most_frequent_ecotypes <- apply(ecotype_summary, 1, function(x) names(x)[which.max(x)])

# We can then replace the labels for the clusters and create a better understandable graph for the report
dapc3$grp <- factor(dapc3$grp, labels = most_frequent_ecotypes)
scatter(dapc3, xax=1, yax=2, grp=dapc3$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")),
        pch=16, bg="white",cstar = 0, cellipse = 2,clabel = 0,scree.da=TRUE,
        scree.pca=TRUE, leg=TRUE)

# And the same with the PCA with DAPC clusters
s.class(pca.orca_a$li, fac=dapc3$grp, clab=0.7, col=transp(c("forestgreen","dodgerblue4","deeppink","orange2","brown")), csta=1, cpoint=2, cellipse =1, xax=1, yax=2, grid = FALSE)
add.scatter.eig(pca.orca_a$eig[1:10],3,1,2, ratio=.3)
title(xlab = "PC1: 69.8 %", ylab = "PC2: 10.4 %")


## Zoom Neutral Structure: Resident Orcas ####
# Let's also have a closer look at whats going on in the neutral structure within the resident ecotype

# First we subset the data set again
resident_inds <- which(data_n@other$ecotype == "Resident")
data_zn <- data_n[resident_inds, ] # zn = zoom-in neutral

#then we run the rest as usual

### PCA ####
x.orca_zn <- tab(data_zn, freq=TRUE, NA.method="mean")
pca.orca_zn <- dudi.pca(x.orca_zn, center=TRUE, scale=FALSE, scannf = FALSE, nf =4)
s.class(pca.orca_zn$li, fac=pop(data_zn), col=transp(funky(15),.6), axesel=FALSE, cstar=0, cpoint=3, clabel = 0.7, grid =FALSE)
add.scatter.eig(pca.orca_zn$eig[1:35],2,1,2, ratio=.3)
title(xlab = "PC1: 11.5 %", ylab = "PC2: 7.0 %")
# The first 2 PCs should do. To better visualize differences let's now turn to DAPC

# As the first 2 Pcs explain most variation (see PCA) we only retain those for downstream analysis n.pca

### DAPC ####

grp_zn <- find.clusters(data_zn, n.pca=2, max.n=4,scale=FALSE)
# Again very difficult to say, 3 or 4???

dapc5 <- dapc(data_zn, pop=grp_zn$grp, n.pca=2)
# 2 can be retained

#First DAPC
scatter(dapc5, xax=1, yax=2, grp=dapc5$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2")),
        pch=16, bg="white",cstar = 1, cellipse = 1,clabel = 0.7)

#Now we can use the DAPC clusters on the PCA again:
s.class(pca.orca_zn$li, fac=dapc5$grp, clab=1, col=transp(c("forestgreen","dodgerblue4","deeppink","orange2")), csta=1, cpoint=2, cellipse =1, xax=1, yax=2, grid = FALSE)
add.scatter.eig(pca.orca_zn$eig[1:35],2,1,2, ratio=.3)
title(xlab = "PC1: 11.5 %", ylab = "PC2: 7.0 %")

#And tidying it up
cluster_assignments <- dapc5$grp
pop_info <- data_zn@pop
cluster_to_pop_df <- data.frame(
  Cluster = cluster_assignments,
  Population = pop_info
)

# We then can assign each cluster the name of the most frequent ecotype, which is the only ecotype that occurs in it as we've seen above
pop_summary <- table(cluster_to_pop_df$Cluster, cluster_to_pop_df$Population)
most_frequent_pop <- apply(pop_summary, 1, function(x) names(x)[which.max(x)])

# We can then replace the labels for the clusters and create a better understandable graph for the report
dapc5$grp <- factor(dapc5$grp, labels = most_frequent_pop)
s.class(pca.orca_zn$li, fac=dapc5$grp, clab=0.7, col=transp(c("forestgreen","dodgerblue4","deeppink","orange2")), csta=1, cpoint=2, cellipse =1, xax=1, yax=2, grid = FALSE)


# FST ####
# Now we wanna have a look at genetic distance between populations:

## Neutral FST ####
# This might tell us about the relationship between populations/ecotypes as neutral variation tells us about demographic happenings 

# getting pairwise FSTs (takes some time!)
data_fn <- genind2hierfstat(data_n)
matFst_n <- pairwise.WCfst(data_fn)
# Here we replace all the NA values with 0 for easier plotting
matFst_n <- matFst_n %>% 
  replace(., is.na(.), 0)

# Getting a Tree
cows.tree <- nj(matFst_n)
plot(cows.tree, type="phylogram", tip.col="black", font=2)
add.scale.bar(x=0, y=8.9, length=0.05, lwd=2, cex=1, font=1)

###get confidence intervals to check significance
# change population names to numbers
pop_num <- as.numeric(data_fn$pop)

#make a data frame with population numbers in the first column, followed by the data
data_fn_dat_num <- data.frame(pop_num,data_fn[,-1])

# Performs bootstrapping over loci of pairwise Fst
fst_boot <- boot.ppfst(data_fn_dat_num) # This might take a while
# lower CI
fst_boot$ll
# upper CI
fst_boot$ul

# None of the CIs overlap with 0 so they are all significantly above that


## Adaptive FST ####
# This might be a proxy for similar adaptations

data_fa <- genind2hierfstat(data_a)
matFst_a <- pairwise.WCfst(data_fa)
# Here we replace all the NA values with 0 for easier plotting
matFst_a <- matFst_a %>% 
  replace(., is.na(.), 0)

# Reorder the matrix so that same ecotypes are grouped together for better visualization
pop_to_ecotype <- c(
  "SR" = "Resident",
  "AR" = "Resident",
  "OS" = "Offshore",
  "AT" = "Transient",
  "CT" = "Transient",
  "ICE" = "Icelandic",
  "RUS" = "Resident",
  "BS" = "Resident",
  "MI" = "Antarctic"
)
ecotypes <- pop_to_ecotype[rownames(matFst_a)]
pop_order <- order(ecotypes)
matFst_a_reordered <- matFst_a[pop_order, pop_order]

# visualize the final heatmap 
pheatmap(matFst_a_reordered, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = TRUE)

###get confidence intervals to check significance
# change population names to numbers
pop_num <- as.numeric(data_fa$pop)

#make a data frame with population numbers in the first column, followed by the data
data_fa_dat_num <- data.frame(pop_num,data_fa[,-1])

# Performs bootstrapping over loci of pairwise Fst
fst_boot <- boot.ppfst(data_fa_dat_num) # This might take a while
# lower CI
fst_boot$ll
# upper CI
fst_boot$ul

# lower CI overlaps 0 between 7/8  2/7  2/8: 2 = AR, 7= RUS, 8 =BS
# => No significant FSTs between Rus,BS,AR,  
# BUT significant FST (even tho ~0.001) between AT & CT


# GenDiv (Ho, He, FIS) ####

## Felix' approach (table) ####
# Let's calculate the genetic diversity per pop statistics in the neutral dataset (cause outlier loci might give wrong signal)

# Get the unique populations
populations <- unique(data_n@pop)  

# Initialize empty vectors for Ho, Hs, and FIS per population
ho_he_fis_per_pop <- sapply(populations, function(pop) {
  # Subset the data for that population
  pop_data <- data_n[data_n@pop == pop, ]
  
  # Calculate Ho, Hs, and FIS per locus for that population (use the subset of loci)
  pop_stats <- basic.stats(genind2hierfstat(pop_data))  # Convert to hierfstat format
  
  ho_for_pop <- pop_stats$Ho  # Get Ho per locus
  hs_for_pop <- pop_stats$Hs  # Get Hs per locus (heterozygosity within population)
  
  # Calculate the mean Ho and Hs across loci for this population
  mean_ho <- mean(ho_for_pop, na.rm = TRUE)
  mean_hs <- mean(hs_for_pop, na.rm = TRUE)
  
  # Calculate FIS as (Hs - Ho) / Hs
  mean_fis <- (mean_hs - mean_ho) / mean_hs
  
  return(c(mean_ho, mean_hs, mean_fis))  # Return Ho, Hs, and FIS
})

# Convert the result into a matrix for easier viewing
ho_he_fis_matrix <- t(ho_he_fis_per_pop)  # Transpose so populations are rows
colnames(ho_he_fis_matrix) <- c("Ho", "Hs", "FIS")

# Add population names as row names
rownames(ho_he_fis_matrix) <- populations

# Create a final data frame with Ho, Hs, and FIS for each population
final_table <- data.frame(ho_he_fis_matrix)

# Print the final table
print(final_table)

#Safe
write.csv(final_table, "stats_table.csv", row.names = FALSE, sep = "/")

## Emil's approach (graph) ####

### HETEROZYGOSITY ####



Ho_per_locus <- apply(data_n@tab, 2, function(x) {
  het_count <- sum(x == 1, na.rm = TRUE)  # Count heterozygous calls
  total <- sum(!is.na(x))  # Total non-missing genotypes
  het_count / total  # Observed heterozygosity
})


Ho_per_individual <- apply(data_n@tab, 1, function(x) {
  het_count <- sum(x == 1, na.rm = TRUE)  # Count heterozygous loci
  total <- sum(!is.na(x))  # Total loci
  het_count / total  # Ho per individual
})




Ho_by_pop <- data.frame(Population = pop(data_n), Ho = Ho_per_individual) %>%
  group_by(Population) %>%
  summarise(Mean_Ho = mean(Ho, na.rm = TRUE))

print(Ho_by_pop)

#He_summary$Mean_Ho <- data.frame(Population = pop(data_n), Ho = Ho_per_individual) %>%
group_by(Population) %>%
  summarise(Mean_Ho = mean(Ho, na.rm = TRUE))



### Visualizing Heterozygosity ###

#quick overview:
summary(Ho_by_pop)

Ho_by_pop

#Variance within populations:
Ho_by_pop %>% summarise(Mean = mean(Mean_Ho), SD = sd(Mean_Ho))

#compare distributions of observed heterozygosity across populations
ggplot(Ho_by_pop, aes(x = Population, y = Mean_Ho, fill = Population)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  theme_minimal() +
  labs(title = "Observed Heterozygosity Across Populations",
       x = "Population",
       y = "Mean Observed Heterozygosity (Ho)")






### Expected Heterozygosity ###

genind_hier <- genind2hierfstat(data_n)
het_values <- basic.stats(genind_hier)

# Extract Expected Heterozygosity (He)
He_by_pop <- data.frame(Population = unique(pop(data_n)),
                        Mean_He = colMeans(het_values$Hs, na.rm = TRUE))

# Compute mean He per population
He_summary <- aggregate(Mean_He ~ Population, data = He_by_pop, mean, na.rm = TRUE)
colnames(He_summary) <- c("Population", "Mean_He")

# Compute mean observed heterozygosity per population
Ho_summary <- aggregate(Mean_Ho ~ Population, data = Ho_by_pop, mean, na.rm = TRUE)

# Merge the two summaries by Population
He_summary <- merge(He_summary, Ho_summary, by = "Population")


# Print results
print(He_summary)







### COMPARING He & Ho ###

#PLot of Ho
ggplot(Ho_by_pop, aes(x = Population, y = Mean_Ho, fill = Population)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Observed Heterozygosity by Population", y = "Mean Ho", x = "Population")




#Scatterplot Ho vs He

het_summary <- merge(Ho_by_pop, He_by_pop, by = "Population")


ggplot(het_summary, aes(x = Mean_He, y = Mean_Ho, color = Population)) +
  geom_point(size = 4) +  
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid") +  # 1:1 line
  theme_minimal() +
  labs(title = "Expected vs Observed Heterozygosity", 
       x = "Expected Heterozygosity", 
       y = "Observed Heterozygosity")

het_summary




#Check if normal distribution:
ggplot(He_summary, aes(x = Mean_He)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "blue", alpha = 0.5) +
  geom_density(color = "red", size = 1) +
  theme_minimal() +
  labs(title = "Histogram & Density Plot of Expected Heterozygosity (He)")

#OR SHAPIRO TEST:
shapiro.test(He_summary$Mean_He)


##Statistical test for Hexp

#Statistical comparison:
kruskal.test(Mean_He ~ Population, data = He_summary)

# Pairwise comparisons
pairwise.wilcox.test(He_summary$Mean_He, He_summary$Population, p.adjust.method = "bonferroni")


##Statistical test for Hobs

#Statistical comparison:
kruskal.test(Mean_Ho ~ Population, data = Ho_by_pop)

# Pairwise comparisons
pairwise.wilcox.test(Ho_by_pop$Mean_Ho, Ho_by_pop$Population, p.adjust.method = "bonferroni")

# Admixture (sNMF plots) ####

#### ANCESTRY #####

# Ancestry for all pops ###
#setup
#first, create a data frame containing both the population (harbour) and genetic group assignment for each individual
pop.data <- data.frame(Population=pop_data$pop_id, Group=pop_data$ecotype)

#then, summarise data2 and count the number of each genetic group in each population
count <- pop.data %>% count(Population, Group) %>% data.frame()
head(count)


#genofile <- vcf2geno("C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/orca_data_3323.vcf", "genofile.geno")


#first, create a hierfstat object
data_dat <- genind2hierfstat(data_n)

#then use this object to write a ".str" file (originally used in the STRUCTURE software, which can also calculate ancestry coefficients)
write.struct(data_dat, ilab=indNames(data_n), pop = data_n$pop, fname = "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data.str")

#lastly, convert this file into a ".geno" file
struct2geno(input.file = "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data.str", ploidy = 2, FORMAT = 2, extra.column = 2, extra.row = 0) 

#create a string with the filename, then read the file
genofile <- "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data.str.geno"
cic.geno <- read.geno(genofile)

project <- snmf(genofile, K = 2:9, entropy = TRUE, repetitions = 2, project = "new")

plot(project, col = transp("steelblue4"), pch = 19)

#calculate the ancestry of each individual for K=5
cic.snmf <- snmf(genofile, K=5, project="new")

#extract the probability matrix for K=5
qmatrix <- Q(cic.snmf, K=5)

barplot(t(qmatrix), 
        col=c("forestgreen","dodgerblue4","deeppink","orange2", "red", "blue", "black"), 
        border=NA, space=0, 
        xlab="Individuals", 
        ylab = "Ancestry")


#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix_new <- data.frame(ID=rownames(pop.data),pop.data,qmatrix)

#then, we restructure the data a little, so that we get two columns: 1. Prob: the ancestry probability; and 2. Variable: variable name, name of the cluster

#create empty list
qmatrix_new_list <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix_new)){ 
  qmatrix_new_list[[i-3]] <- qmatrix_new[,c(1:3,i)] %>% mutate("Var"=rep(colnames(qmatrix_new)[i],nrow(qmatrix_new)))
  colnames(qmatrix_new_list[[i-3]]) <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix_new <- bind_rows(qmatrix_new_list)

#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
for(i in c("AR","SR","RU","BS","AT","CT","OS","IC","MI")){
  max_anc <- filter(qmatrix_new,ID==i) %>% filter(Prob==max(Prob))
  qmatrix_new$Variable <- replace(qmatrix_new$Variable,qmatrix_new$Variable==max_anc$Variable,max_anc$Kmeans_cluster)
}


#plot the ancestry coefficients in a nicer-looking barplot
anc_plot <- ggplot(qmatrix_new, aes(factor(ID), Prob, fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(x = "", title = , y = "Probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  ggtitle("Admixture Plot for all ecotypes")+
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5)) +
  scale_fill_manual(name="Cluster",values=c("forestgreen","dodgerblue4","purple","orange","red", "blue", "black","pink","yellow"))
anc_plot




### Ancestry only residents ###
#setup
#first, create a data frame containing both the population (harbour) and genetic group assignment for each individual
pop.data <- data.frame(Population=pop_data$pop_id, Group=pop_data$ecotype)

#then, summarise data2 and count the number of each genetic group in each population
count <- pop.data %>% count(Population, Group) %>% data.frame()
head(count)

pop.data2 <- pop.data %>% filter(Group == "R")
#then, summarise data2 and count the number of each genetic group in each population
count2 <- pop.data2 %>% count(Population, Group) %>% data.frame()
head(count2)


#first, create a hierfstat object
data_dat <- genind2hierfstat(data_n2)

#then use this object to write a ".str" file (originally used in the STRUCTURE software, which can also calculate ancestry coefficients)
write.struct(data_dat, ilab=indNames(data_n2), pop = data_n2$pop, fname = "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data_n2.str")

#lastly, convert this file into a ".geno" file
struct2geno(input.file = "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data_n2.str", ploidy = 2, FORMAT = 2, extra.column = 2, extra.row = 0) 

#create a string with the filename, then read the file
genofile2 <- "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data_n2.str.geno"
cic.geno2 <- read.geno(genofile2)

project2 <- snmf(genofile2, K = 2:4, entropy = TRUE, repetitions = 2, project = "new")

plot(project2, col = transp("steelblue4"), pch = 19)

#calculate the ancestry of each individual for K=4
cic.snmf2 <- snmf(genofile2, K=3, project="new")

#extract the probability matrix for K=4
qmatrix2 <- Q(cic.snmf2, K=3)


barplot(t(qmatrix2), 
        col=c("forestgreen","dodgerblue4","deeppink","orange2"), 
        border=NA, space=0, 
        xlab="Individuals", 
        ylab = "Ancestry")


#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix_new2 <- data.frame(ID=rownames(pop.data2),pop.data2,qmatrix2)



#then, we restructure the data a little, so that we get two columns: 1. Prob: the ancestry probability; and 2. Variable: variable name, name of the cluster

#create empty list
qmatrix_new_list2 <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix_new2)){ 
  qmatrix_new_list2[[i-3]] <- qmatrix_new2[,c(1:3,i)] %>% mutate("Var"=rep(colnames(qmatrix_new2)[i],nrow(qmatrix_new2)))
  colnames(qmatrix_new_list2[[i-3]]) <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix_new2 <- bind_rows(qmatrix_new_list2)

#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
existing_ids <- qmatrix_new2$ID
for(i in c("AR", "SR", "RUS", "BS")) {
  if(i %in% existing_ids) {
    max_anc2 <- filter(qmatrix_new2, ID == i) %>% filter(Prob == max(Prob, na.rm = TRUE))
    if(nrow(max_anc2) > 0) {
      qmatrix_new2$Variable <- replace(qmatrix_new2$Variable, qmatrix_new2$Variable == max_anc2$Variable, max_anc2$Kmeans_cluster)
    }
  }
}



#plot the ancestry coefficients in a nicer-looking barplot
anc_plot2 <- ggplot(qmatrix_new2, aes(factor(ID), Prob, fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(x = "", title = , y = "Probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  ggtitle("Admixture Plot for resident ecotype")+
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5)) +
  scale_fill_manual(name="Cluster",values=c("forestgreen","dodgerblue4","purple","orange"))
anc_plot2









### Ancestry non-residents ###
#setup
#first, create a data frame containing both the population (harbour) and genetic group assignment for each individual
pop.data <- data.frame(Population=pop_data$pop_id, Group=pop_data$ecotype)

#then, summarise data2 and count the number of each genetic group in each population
count <- pop.data %>% count(Population, Group) %>% data.frame()
head(count)

pop.data3 <- pop.data %>% filter(Group != "R" | is.na(Group))
#then, summarise data3 and count the number of each genetic group in each population
count3 <- pop.data3 %>% count(Population, Group) %>% data.frame()
head(count3)

#first, create a hierfstat object
data_dat3 <- genind2hierfstat(data3)

#then use this object to write a ".str" file (originally used in the STRUCTURE software, which can also calculate ancestry coefficients)
write.struct(data_dat3, ilab=indNames(data3), pop = data3$pop, fname = "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data3.str")

#lastly, convert this file into a ".geno" file
struct2geno(input.file = "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data3.str", ploidy = 2, FORMAT = 2, extra.column = 2, extra.row = 0) 

#create a string with the filename, then read the file
genofile3 <- "C:/Users/emilt/OneDrive/Documents/BIO418/Mini_Project/Raw_Data/data3.str.geno"
cic.geno3 <- read.geno(genofile3)

project3 <- snmf(genofile3, K = 2:18, entropy = TRUE, repetitions = 2, project = "new")

plot(project3, col = transp("steelblue4"), pch = 19)

#calculate the ancestry of each individual for K=5
cic.snmf3 <- snmf(genofile3, K=5, project="new")

#extract the probability matrix for K=5
qmatrix3 <- Q(cic.snmf3, K=5)

barplot(t(qmatrix3), 
        col=c("forestgreen","dodgerblue4","purple","orange","red"), 
        border=NA, space=0, 
        xlab="Individuals", 
        ylab = "Ancestry")



#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix_new3 <- data.frame(ID=rownames(pop.data3),pop.data3,qmatrix3)

#create empty list
qmatrix_new_list3 <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix_new3)){ 
  qmatrix_new_list3[[i-3]] <- qmatrix_new3[,c(1:3,i)] %>% mutate("Var"=rep(colnames(qmatrix_new3)[i],nrow(qmatrix_new3)))
  colnames(qmatrix_new_list3[[i-3]]) <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix_new3 <- bind_rows(qmatrix_new_list3)

#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
existing_ids <- qmatrix_new3$ID
for(i in c("AR", "SR", "RUS", "BS","ICE")) {
  if(i %in% existing_ids) {
    max_anc3 <- filter(qmatrix_new3, ID == i) %>% filter(Prob == max(Prob, na.rm = TRUE))
    if(nrow(max_anc3) > 0) {
      qmatrix_new3$Variable <- replace(qmatrix_new3$Variable, qmatrix_new3$Variable == max_anc3$Variable, max_anc3$Kmeans_cluster)
    }
  }
}

#plot the ancestry coefficients in a nicer-looking barplot
anc_plot3 <- ggplot(qmatrix_new3, aes(factor(ID), Prob, fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(x = "", y = "Probability") +
  ggtitle("Admixture Plot for NON-resident ecotype")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5)) +
  scale_fill_manual(name="Cluster",values=c("forestgreen","dodgerblue4","purple","orange2","red"))
anc_plot3




anc_plot
anc_plot2
anc_plot3










