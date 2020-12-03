## MATH 459 - Statistical Bioinformatics
## Final Project Pipeline
## Hayden Brown & Bryan Edmonson

library('snpStats')
library(hash)
library(ggplot2)
library(tidyverse)
library(mice)
library('VIM')
library(lattice)

################### STEP 1: DATA IMPORTATION #######################
data(families)

# CHARTS PART 1: NUMBER OF MEMBERS PER FAMILY
fam_uniq <- unique(pedData$familyid)
length(fam_uniq) # 756 unique families

f_hash <- hash()
for (family in fam_uniq) {
  famlist <- subset(pedData, familyid == family)
  num_members <- nrow(famlist)
  f_hash[[family]] <- num_members
}
members <- values(f_hash, USE.NAMES = FALSE)
mem_data <- as.data.frame(members)

number_ticks <- function(n) {function(limits) pretty(limits, n)}
ggplot(mem_data, aes(members)) + 
  geom_histogram(binwidth = 1, fill=rainbow(11)) + 
  scale_x_continuous(breaks=number_ticks(11)) + 
  scale_y_continuous(breaks=number_ticks(11)) + 
  ggtitle("Number of Members per Family") +
  geom_text(stat="count", aes(label=..count..), vjust=-1)


# CHARTS PART 2: PROPORTION OF NA VALUES

data_complete <- na.omit(genotypes)
row.names(data_complete) <- 1:nrow(data_complete)
comp_df <- as.data.frame(data_complete)
num_comp <- nrow(comp_df)
num_NA <- nrow(genotypes) - num_comp

pie_df <- data.frame(Category = c("Complete", "Contains NA Values"), "freq" = c(num_comp, num_NA))
require(grid)
require(gridExtra)
bar <- ggplot(pie_df, aes (x=Category, y = freq)) +
  geom_bar(stat = 'identity', aes(fill = Category)) +
  theme_classic() +
  ggtitle("Number of Genomes by NA Value Status") +
  geom_text(aes(label = freq), vjust=1)

pie <- ggplot(pie_df, aes (x="", y = freq, fill = factor(Category))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste(round(freq / sum(freq) * 100, 1), "%")),
            position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of Genomes with NA Values") + 
  coord_polar("y")
grid.arrange(bar,pie,ncol=1)


# CHARTS PART 3: PARENT VS CHILD PROPORTION VS THAT W NA VALUES
genmat <- as(genotypes, 'numeric')
data_na <- genmat[!complete.cases(genmat),]
data_na <- as.data.frame(data_na)
na_ids <- row.names(data_na)

na_members <- numeric()
for (id in na_ids) {
  member <- pedData[id,]
  
  if(is.na(member$father) && is.na(member$mother)) {
    status <- "parent"
  } else {
    status <- "child"
  }
  na_members <- c(na_members, status)
}
na_members <- as.data.frame(na_members)

ggplot(na_members, aes(na_members)) + 
  geom_bar(color='black', aes(fill=na_members)) + 
  ggtitle("Incomplete Genomes by Family Member Status") +
  theme(legend.position = "none") +
  geom_text(stat="count", aes(label=..count..), vjust=1)

# CHARTS PART 4: ELIMINATING FAMILIES WITH MISSING DATA

# removing families with missing parent data
row.names(na_members) <- na_ids

rm.na_parent <- numeric()
rm.na_all_children <- numeric()
rm.na_some_children <- numeric()
for (id in na_ids) {
  temp <- pedData[id,]
  family <- temp$familyid

  # all parents with NA values has their family removed
  if (na_members[id,] == "parent"){
    rm.na_parent <- c(rm.na_parent, family)
  }
  # children with NA values are a little more complicated
  else {
    # find parent ids
    famtemp <- pedData[pedData$familyid == family,]
    parent_ids <- numeric()
    for(j in rownames(famtemp)) {
      member <- pedData[id,]
      if(is.na(member$father) && is.na(member$mother)){
        parent_ids <- c(parent_ids, j)
      }
    }
    # if parent in na_ids, skip (they will be filtered when parent is checked)
    if(parent_ids[1] %in% na_ids || parent_ids[2] %in% na_ids){
      next
    } else {
      # find all siblings 
      sibling_ids <- numeric()
      for(k in rownames(famtemp)) {
        member <- pedData[id,]
        if(!is.na(member$father) && k != id) {
          sibling_ids <- c(sibling_ids, k)
        }
      }
      # if any sibling not in na_ids, add to na_some_children
      for (l in sibling_ids) {
        if (!(l %in% na_ids)) {
          rm.na_some_children <- c(rm.na_some_children, family)
          break
        }
      }
      # else add to na_all_children
      if(!(family %in% rm.na_some_children)) {
      rm.na_all_children <- c(rm.na_all_children, family)
      }
    }
  }
}
# remove repeats from lists if any
rm.na_parent <- unique(rm.na_parent)
rm.na_all_children <- unique(rm.na_all_children)
rm.na_some_children <- unique(rm.na_some_children)

for (fam in rm.na_parent) {
  if (fam %in% rm.na_all_children) {
    rm.na_all_children <- rm.na_all_children[rm.na_all_children != fam]
  }
  if (fam %in% rm.na_some_children) {
    rm.na_some_children <- rm.na_some_children[rm.na_some_children != fam]
  }
}
num.parent <- length(rm.na_parent)
num.all <- length(rm.na_all_children)
num.some <- length(rm.na_some_children)
num.complete <- length(fam_uniq) - num.parent - num.all - num.some

# plot time
fam_df <- data.frame(Category = c("Complete", "Some Children Removed", "All Children Removed", "Parent Removed"), "freq" = c(num.complete, num.some, num.all, num.parent))
fam_df$Category <- factor(fam_df$Category, levels = c("Complete", "Some Children Removed", "All Children Removed", "Parent Removed"))

fam_df %>%
  mutate(Legend = ifelse(Category %in% c("Complete","Some Children Removed"), "Family Kept", "Family Removed")) %>%
  ggplot(aes(x=Category, y=freq)) + 
  ggtitle("Genome Families by Members With NA Values") +
  geom_bar(stat='identity', aes(fill = Legend)) +
  scale_fill_manual(name = "Legend", values = c("seagreen1", "red2")) +
  geom_text(aes(label = freq), vjust=-0.5)


# CHARTS DONE: CREATE FINAL DATA SET

fam_cut_list <- c(rm.na_all_children, rm.na_parent)
fam_list <- rownames(pedData[!pedData$familyid %in% fam_cut_list,])
length(fam_list) # 673 genomes

final_data <- genmat[fam_list,]
final_data <- na.omit(final_data) # 508 remaining genomes
final_data <- as.data.frame(final_data)

b_and_a <- data.frame(cat= c("Before", "After"), num_families = c(756, 172),num_genomes = c(3017, 508))
b_and_a$cat <- factor(b_and_a$cat, levels = c("Before", "After"))

require(grid)
require(gridExtra)
p1 <- ggplot(b_and_a, aes(x = cat, y = num_families)) + 
      geom_bar(stat = 'identity', aes(fill = cat)) + 
      ggtitle("Number of Families") +
      theme(legend.position = "none") +
      geom_text(aes(label = num_families), vjust=-0.5)

p2 <- ggplot(b_and_a, aes(x = cat, y = num_genomes)) + 
      geom_bar(stat = 'identity', aes(fill = cat)) +
      ggtitle("Number of Genomes") +
      theme(legend.position = "none") +
      geom_text(aes(label = num_genomes), vjust=-0.5)

grid.arrange(p1,p2,ncol=2, top=textGrob("Before vs. After Data Cleaning",gp=gpar(fontsize=20,font=3)))


################# STEP 2: RANDOM MISSING VALUES ####################


set.seed(1)

# takes random sample 15% of data set and converts to NA values
miss <- final_data
miss <- as.data.frame(miss)
while(sum(is.na(miss) == TRUE) < (nrow(miss) * ncol(miss) * 0.15)){
  miss[sample(nrow(miss),1), sample(ncol(miss),1)] <- NA
}

final_data == miss # identical except for the new NA values

# visualize pattern of missing data

require(mice)
miss_sample <- head(miss, 40)
md.pattern(miss_sample) # find a way to make this usable?

require('VIM')
aggr_plot <- aggr(miss_sample, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
# ^ make this usable too!

tempData <- mice(miss,meth='pmm',seed=101)

pmm <- tempData
summary(pmm)

pmm$imp
pmm$meth

completedData <- complete(pmm,1)
densityplot(x=pmm, data= ~rs91126+rs62927+rs79960+rs19348+rs99786+
              rs52628+rs6699+rs12373+rs35215+rs41229+rs86267+
              rs23261+rs69208+rs16483+rs8558+rs55762+rs8124+
              rs72056+rs82369+rs97686+rs77065+rs53106+rs37378+
              rs83832+rs35431+rs61158+rs32410+rs85906+rs83977+
              rs24527+rs73721+rs36088+rs32998+rs5566+rs98256+
              rs29479+rs42938+rs32018+rs39483+rs42367+rs87640+rs98918)

num_NAs <- sum(sapply(miss, function(x) sum(is.na(x))))
numtrues <- colSums(completedData==final_data)
num_false <- sum(sapply(numtrues, function(x) 508-x ), na.rm=T)

success_rate = (1- (num_false / num_NAs)) * 100 # 93.31706 %

success_matrix = as.data.frame(completedData == final_data)

image(as.matrix(success_matrix), col = c("red", "blue"), main = "Correctness of Imputed Data Set vs. Original")

library(reshape)
plotDat <- reshape::melt(success_matrix)
ggplot(plotDat, aes(X1, X2, fill = c("red2", "seagreen1")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ggtitle("ggplot geom_tile")

