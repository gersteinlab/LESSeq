################################################################################
###   R code to identify local events from splicing graphs
################################################################################



### A splicing graph for each gene is required, which is generated by "classify" program and named "LOCAL_EVENT_GROUP_ID.matrix"
### Specify the directory containing splicing graphs (need to include "/" at the end)
prefix <- SPLICING_GRAPH_DIRECTORY


### Specify the path to gene grouping file that was used to generate the gene splicing graphs (UCSC_GENE2ISOFORM format)
group <- LOCAL_EVENT_GROUP_FILE


### Specify the prefix of output local events annotation files
out_prefix <- OUT_PREFIX


### Name output annotation files 
ES.1 <- paste(out_prefix,"ES.interval",sep="")
ES.2 <- paste(out_prefix,"ES.map",sep="")
RI.1 <- paste(out_prefix,"RI.interval",sep="")
RI.2 <- paste(out_prefix,"RI.map",sep="")
A5SS.1 <- paste(out_prefix,"A5SS.interval",sep="")
A5SS.2 <- paste(out_prefix,"A5SS.map",sep="")
A3SS.1 <- paste(out_prefix,"A3SS.interval",sep="")
A3SS.2 <- paste(out_prefix,"A3SS.map",sep="")
MXE.1 <- paste(out_prefix,"MXE.interval",sep="")
MXE.2 <- paste(out_prefix,"MXE.map",sep="")
AFE.1 <- paste(out_prefix,"AFE.interval",sep="")
AFE.2 <- paste(out_prefix,"AFE.map",sep="")
ALE.1 <- paste(out_prefix,"ALE.interval",sep="")
ALE.2 <- paste(out_prefix,"ALE.map",sep="")
T3.1 <- paste(out_prefix,"T3.interval",sep="")
T3.2 <- paste(out_prefix,"T3.map",sep="")


### Read data files and generate local events
select<-read.table(group,as.is=T)
count<-table(select[,1])
select<-names(count[count>1])
ES.count<-RI.count<-A5SS.count<-A3SS.count<-MXE.count<-AFE.count<-ALE.count<-T3.count<-1

for (j in 1:length(select)) {
	print (paste("processing gene: ",select[j],sep=""))

	file<-paste(prefix,select[j],".matrix",sep="")
	cor<-readLines(file,n=1)
	cor.split<-unlist(strsplit(cor,"\t",fixed=T))
	chr<-cor.split[1]
	strand<-cor.split[2]
	pos<-as.integer(unlist(strsplit(gsub("[^0-9]"," ",cor.split[3])," {1,}"))[-1])
	matrix<-read.table(file,skip=1,as.is=T)
	
	if (ncol(matrix)<3) {next}
	else {

		usage<-apply(matrix,2,sum)
		k<-length(pos)

		for (i in 2: (length(usage)-1)) {
    		if (nrow(matrix)==usage[i-1] & nrow(matrix)==usage[i+1] & nrow(matrix)>usage[i]) {
       			if ((pos[i*2-1]-pos[i*2-2])>0 & (pos[i*2+1]-pos[i*2])>0) {
       				write(c(paste(select[j],i,"1",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],3,paste(pos[i*2-3],pos[i*2-1],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2],pos[i*2+2],sep=",")),ES.1,append=T,sep="\t",ncolumns=8) 
       				write(c(paste(select[j],i,"2",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],2,paste(pos[i*2-3],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2+2],sep=",")),ES.1,append=T,sep="\t",ncolumns=8) 
					write(c(ES.count,paste(select[j],i,"1",sep="|")),ES.2,append=T,sep="\t",ncolumns=2)
					write(c(ES.count,paste(select[j],i,"2",sep="|")),ES.2,append=T,sep="\t",ncolumns=2)
					ES.count<-ES.count+1
				}
				if (pos[i*2-1]-pos[i*2-2]==0 & pos[i*2+1]-pos[i*2]==0) {
       				write(c(paste(select[j],i,"1",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],3,paste(pos[i*2-3],pos[i*2-1],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2],pos[i*2+2],sep=",")),RI.1,append=T,sep="\t",ncolumns=8) 
       				write(c(paste(select[j],i,"2",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],2,paste(pos[i*2-3],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2+2],sep=",")),RI.1,append=T,sep="\t",ncolumns=8) 
					write(c(RI.count,paste(select[j],i,"1",sep="|")),RI.2,append=T,sep="\t",ncolumns=2)
					write(c(RI.count,paste(select[j],i,"2",sep="|")),RI.2,append=T,sep="\t",ncolumns=2)
					RI.count<-RI.count+1
				}
				if ((strand=="+" & pos[i*2-1]-pos[i*2-2]==0 & pos[i*2+1]-pos[i*2]>0) | (strand=="-" & pos[i*2-1]-pos[i*2-2]>0 & pos[i*2+1]-pos[i*2]==0)) {
       				write(c(paste(select[j],i,"1",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],3,paste(pos[i*2-3],pos[i*2-1],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2],pos[i*2+2],sep=",")),A5SS.1,append=T,sep="\t",ncolumns=8) 
       				write(c(paste(select[j],i,"2",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],2,paste(pos[i*2-3],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2+2],sep=",")),A5SS.1,append=T,sep="\t",ncolumns=8) 
					write(c(A5SS.count,paste(select[j],i,"1",sep="|")),A5SS.2,append=T,sep="\t",ncolumns=2)
					write(c(A5SS.count,paste(select[j],i,"2",sep="|")),A5SS.2,append=T,sep="\t",ncolumns=2)
					A5SS.count<-A5SS.count+1					
				}
				if ((strand=="+" & pos[i*2-1]-pos[i*2-2]>0 & pos[i*2+1]-pos[i*2]==0) | (strand=="-" & pos[i*2-1]-pos[i*2-2]==0 & pos[i*2+1]-pos[i*2]>0)) {
       				write(c(paste(select[j],i,"1",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],3,paste(pos[i*2-3],pos[i*2-1],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2],pos[i*2+2],sep=",")),A3SS.1,append=T,sep="\t",ncolumns=8) 
       				write(c(paste(select[j],i,"2",sep="|"),chr,strand,pos[i*2-3],pos[i*2+2],2,paste(pos[i*2-3],pos[i*2+1],sep=","),paste(pos[i*2-2],pos[i*2+2],sep=",")),A3SS.1,append=T,sep="\t",ncolumns=8) 
					write(c(A3SS.count,paste(select[j],i,"1",sep="|")),A3SS.2,append=T,sep="\t",ncolumns=2)
					write(c(A3SS.count,paste(select[j],i,"2",sep="|")),A3SS.2,append=T,sep="\t",ncolumns=2)
					A3SS.count<-A3SS.count+1					
				}				
       		}
		}
		
		if (ncol(matrix)>=4) {		
			for (i in 4: length(usage)) {
   	 			if (nrow(matrix)==usage[i-3] & nrow(matrix)==usage[i] & sum(matrix[,i-2]==matrix[,i-1])==0) {
					if ((pos[i*2-1]-pos[i*2-2])>0 & (pos[i*2-3]-pos[i*2-4])>0 & (pos[i*2-5]-pos[i*2-6])>0) {
       				write(c(paste(select[j],i,"1",sep="|"),chr,strand,pos[i*2-7],pos[i*2],3,paste(pos[i*2-7],pos[i*2-5],pos[i*2-1],sep=","),paste(pos[i*2-6],pos[i*2-4],pos[i*2],sep=",")),MXE.1,append=T,sep="\t",ncolumns=8) 
       				write(c(paste(select[j],i,"2",sep="|"),chr,strand,pos[i*2-7],pos[i*2],3,paste(pos[i*2-7],pos[i*2-3],pos[i*2-1],sep=","),paste(pos[i*2-6],pos[i*2-2],pos[i*2],sep=",")),MXE.1,append=T,sep="\t",ncolumns=8) 
					write(c(MXE.count,paste(select[j],i,"1",sep="|")),MXE.2,append=T,sep="\t",ncolumns=2)
					write(c(MXE.count,paste(select[j],i,"2",sep="|")),MXE.2,append=T,sep="\t",ncolumns=2)
					MXE.count<-MXE.count+1					
       				}
       			}
			}		
		}
		
		if (nrow(matrix)==usage[3] & sum(matrix[,1]==matrix[,2])==0 & pos[3]-pos[2]>0 & pos[5]-pos[4]>0) {
			if (strand=="+") {
       			write(c(paste(select[j],"AFE","1",sep="|"),chr,strand,pos[3],pos[6],2,paste(pos[3],pos[5],sep=","),paste(pos[4],pos[6],sep=",")),AFE.1,append=T,sep="\t",ncolumns=8) 
       			write(c(paste(select[j],"AFE","2",sep="|"),chr,strand,pos[1],pos[6],2,paste(pos[1],pos[5],sep=","),paste(pos[2],pos[6],sep=",")),AFE.1,append=T,sep="\t",ncolumns=8) 
				write(c(AFE.count,paste(select[j],"AFE","1",sep="|")),AFE.2,append=T,sep="\t",ncolumns=2)
				write(c(AFE.count,paste(select[j],"AFE","2",sep="|")),AFE.2,append=T,sep="\t",ncolumns=2)
				AFE.count<-AFE.count+1	
			}
			if (strand=="-") {
       			write(c(paste(select[j],"ALE","1",sep="|"),chr,strand,pos[3],pos[6],2,paste(pos[3],pos[5],sep=","),paste(pos[4],pos[6],sep=",")),ALE.1,append=T,sep="\t",ncolumns=8) 
       			write(c(paste(select[j],"ALE","2",sep="|"),chr,strand,pos[1],pos[6],2,paste(pos[1],pos[5],sep=","),paste(pos[2],pos[6],sep=",")),ALE.1,append=T,sep="\t",ncolumns=8) 
				write(c(ALE.count,paste(select[j],"ALE","1",sep="|")),ALE.2,append=T,sep="\t",ncolumns=2)
				write(c(ALE.count,paste(select[j],"ALE","2",sep="|")),ALE.2,append=T,sep="\t",ncolumns=2)
				ALE.count<-ALE.count+1				
			}
		}
		
		if(nrow(matrix)==usage[length(usage)-2] & sum(matrix[,length(usage)-1]==matrix[,length(usage)])==0 & pos[k-1]-pos[k-2]>0 & pos[k-3]-pos[k-4]>0) {
			if (strand=="+") {
       			write(c(paste(select[j],"ALE","1",sep="|"),chr,strand,pos[k-5],pos[k-2],2,paste(pos[k-5],pos[k-3],sep=","),paste(pos[k-4],pos[k-2],sep=",")),ALE.1,append=T,sep="\t",ncolumns=8) 
       			write(c(paste(select[j],"ALE","2",sep="|"),chr,strand,pos[k-5],pos[k],2,paste(pos[k-5],pos[k-1],sep=","),paste(pos[k-4],pos[k],sep=",")),ALE.1,append=T,sep="\t",ncolumns=8) 
				write(c(ALE.count,paste(select[j],"ALE","1",sep="|")),ALE.2,append=T,sep="\t",ncolumns=2)
				write(c(ALE.count,paste(select[j],"ALE","2",sep="|")),ALE.2,append=T,sep="\t",ncolumns=2)
				ALE.count<-ALE.count+1				
			}
			if (strand=="-") {
       			write(c(paste(select[j],"AFE","1",sep="|"),chr,strand,pos[k-5],pos[k-2],2,paste(pos[k-5],pos[k-3],sep=","),paste(pos[k-4],pos[k-2],sep=",")),AFE.1,append=T,sep="\t",ncolumns=8) 
       			write(c(paste(select[j],"AFE","2",sep="|"),chr,strand,pos[k-5],pos[k],2,paste(pos[k-5],pos[k-1],sep=","),paste(pos[k-4],pos[k],sep=",")),AFE.1,append=T,sep="\t",ncolumns=8) 
				write(c(AFE.count,paste(select[j],"AFE","1",sep="|")),AFE.2,append=T,sep="\t",ncolumns=2)
				write(c(AFE.count,paste(select[j],"AFE","2",sep="|")),AFE.2,append=T,sep="\t",ncolumns=2)
				AFE.count<-AFE.count+1	
			}
		}
		
		if(strand=="+" & nrow(matrix)==usage[length(usage)-1] & nrow(matrix)>usage[length(usage)] & pos[k-1]-pos[k-2]==0) {
       			write(c(paste(select[j],"T3","1",sep="|"),chr,strand,pos[k-3],pos[k],1,pos[k-3],pos[k]),T3.1,append=T,sep="\t",ncolumns=8) 
       			write(c(paste(select[j],"T3","2",sep="|"),chr,strand,pos[k-3],pos[k-2],1,pos[k-3],pos[k-2]),T3.1,append=T,sep="\t",ncolumns=8) 
				write(c(T3.count,paste(select[j],"T3","1",sep="|")),T3.2,append=T,sep="\t",ncolumns=2)
				write(c(T3.count,paste(select[j],"T3","2",sep="|")),T3.2,append=T,sep="\t",ncolumns=2)
				T3.count<-T3.count+1				   								
		}
		if (strand=="-" & nrow(matrix)==usage[2] & nrow(matrix)>usage[1] & pos[3]-pos[2]==0) {
       			write(c(paste(select[j],"T3","1",sep="|"),chr,strand,pos[1],pos[4],1,pos[1],pos[4]),T3.1,append=T,sep="\t",ncolumns=8) 
       			write(c(paste(select[j],"T3","2",sep="|"),chr,strand,pos[3],pos[4],1,pos[3],pos[4]),T3.1,append=T,sep="\t",ncolumns=8) 
				write(c(T3.count,paste(select[j],"T3","1",sep="|")),T3.2,append=T,sep="\t",ncolumns=2)
				write(c(T3.count,paste(select[j],"T3","2",sep="|")),T3.2,append=T,sep="\t",ncolumns=2)
				T3.count<-T3.count+1				   								
		}

	}

}










