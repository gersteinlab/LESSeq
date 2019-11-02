################################################################################
###   R code to perform Fisher exact test, log-linear model+ likelihood ratio test, Wilcoxon rank sum test to detect differential alternative splicing
################################################################################



################################################################################
###
### Fisher exact test (when only one replicate is available for each condition)
###



###
### One Tab delimited data file (DATA_FILE_COUNT) is required, containing raw counts of reads compatible with each of the two forms of local events 
###
### Data file format
### Row - 1st row:header, from 2nd row: raw read count values (the two forms of each local event should be in consecutive rows)
### Column - 1st column:ID of each form of local event, 2nd and 3rd column: the two samples 



### Specify data file paths
file_count <- DATA_FILE_COUNT
file_output <- OUTPUT_FILE_NAME


### Read data files
data_val <- read.delim(file_count,sep="\t",header=T,row.names=1,as.is=T) 
data_val <- as.matrix(data_val)

### Conduct test

fisher.p <- bonp <- bhp <- matrix(0,nrow(data_val)/2,1)

for(j in seq(2,nrow(data_val),2)){
	
	print(j/2)
		
	fisher_mat <- matrix(c(round(as.numeric(unlist(data_val[j-1,1]))),round(as.numeric(unlist(data_val[j-1,2]))),round(as.numeric(unlist(data_val[j,1]))),round(as.numeric(unlist(data_val[j,2])))),nrow=2)
	fisher.p[j/2]<-fisher.test(fisher_mat,alternative="two.sided")$p.value		
}

### Write output
output_mat <- cbind(rownames(data_val[seq(2,nrow(data_val),2),]),fisher.p,p.adjust(fisher.p,method="bonferroni"),p.adjust(fisher.p,method="BH"))
colnames(output_mat) <- c("ID","rawP","bonP","bhP")		
write.table(output_mat, file_output,quote=F,row.names=F,sep="\t")
		


################################################################################
###
### Log linear model + likelihood ratio test (when two or more replicates are available for each condition)
###
library(multtest)
library(epicalc)
library(lmtest)
library(xtable)
library(MASS)



###
### Two Tab delimited data files are required, one with raw counts of reads compatible with each of the two forms of local events (DATA_FILE_COUNT_ONE), the other with raw counts of all reads mapped to the local events (DATA_FILE_COUNT_ALL)
### The two files should have the same dimensions, as well as same row and column names
###
### Data file format
### Row - 1st row:header, from 2nd row: raw read count values
### Column - 1st column:ID of each form of local event,from 2nd column: samples (all replicates from condition1 should come together first, followed by all replicates from condition2)


### Specify data file paths and number of replicates in each condition
file_count <- DATA_FILE_COUNT_ONE
file_count_norm <- DATA_FILE_COUNT_ALL
condition1_num <- CONDITION1_REPLICATES_NUMBER
condition2_num <- CONDITION2_REPLICATES_NUMBER
file_output <- OUTPUT_FILE_NAME


### Read data files
mat <- read.delim(file_count,sep="\t",header=T,row.names=1,as.is=T) 
offset_val <- read.delim(file_count_norm,sep="\t",header=T,row.names=1,as.is=T) 
mat<-as.matrix(mat)
offset_val<-as.matrix(offset_val)


### Load test function

lrt_test <- function(mat,offset_val,condition1_num,condition2_num){
	
    tissue <- c(rep(1,condition1_num),rep(2,condition2_num))
    rep <- c(seq(1:condition1_num),seq(1:condition2_num))		
    
    pval <- stat <-  matrix(0,nrow(mat),1)
	
	for(d in 1:nrow(mat)){

		
		if(is.na(mat[d,1])==F){
			cat(d,"\n")
			
			m2 <- glm(formula = (as.numeric(unlist(mat[d, ]))) ~ tissue + rep, family="poisson",offset=offset_val[d,]) 
			m1 <- glm(formula = (as.numeric(unlist(mat[d, ]))) ~ rep, family="poisson",offset=offset_val[d,]) 
			
			
			if(is.na(lrtest(m1,m2)[2,5])==F){
				stat[d]<-lrtest(m1,m2)[2,4]
				pval[d] <- lrtest(m1,m2)[2,5]
			}
			
			
			if(is.na(lrtest(m1,m2)[2,5])==T){
				pval[d] <- 1
			}
			
		}
		if(is.na(mat[d,1])==T){
			stat[d] <- NA
			pval[d] <- NA
		}
	}
	
	return(cbind(stat,pval))
}


### Conduct test and write output
output_lrt <- lrt_test(round(mat)+1,log(round(offset_val)+1),condition1_num,condition2_num)
output_mat <- cbind(rownames(mat),output_lrt,p.adjust(output_lrt[,2],method="bonferroni"),p.adjust(output_lrt[,2],method="BH"))
colnames(output_mat) <- c("ID","LRT_statistics","rawP","bonP","bhP")
write.table(output_mat,file_output,quote=F,row.names=F,sep="\t")


################################################################################
###
### Wilcoxon rank sum test (when many replicates are available for each condition)
###



###
### One Tab delimited data file (DATA_FILE_REL) is required, containing estimated relative expression levels of each of the two forms of local events 
###
### Data file format
### Row - 1st row:header, from 2nd row: estimated relative expression values
### Column - 1st column:ID of each form of local event,from 2nd column: samples (all replicates from condition1 should come together first, followed by all replicates from condition2)



### Specify data file paths and number of replicates in each condition
file_rel <- DATA_FILE_REL
condition1_num <- CONDITION1_REPLICATES_NUMBER
condition2_num <- CONDITION2_REPLICATES_NUMBER
file_output <- OUTPUT_FILE_NAME


### Read data files
mat <- read.delim(file_rel,sep="\t",header=T,row.names=1,as.is=T) 
mat<-as.matrix(mat)

### Load Wilcoxon rank sum test function
my.wilcox.test.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

### Conduct test and write output
test<-matrix(nrow=nrow(mat),ncol=2)
for (i in 1:nrow(test)) {
	test[i,1]<-mean(mat[i,1:condition1_num])-mean(mat[i,condition1_num+(1:condition2_num)])
	test[i,2]<-my.wilcox.test.p.value(mat[i,1:condition1_num],mat[i,condition1_num+(1:condition2_num)])
}
test<-cbind(rownames(mat),test,p.adjust(test[,2],method="bonferroni"),p.adjust(test[,2],method="BH"))
colnames(test)<-c("ID","Diff","rawP","bonP","bhP")
write.table(test,file_output,quote=F,sep="\t",row.names=F)















