## This script will load raw microarray data from Illumina's Bead Studio and perform pre-processing for GSEA analysis
## Output files will include five sets
## 1. TB and Media for TST+ (TBAM_POS)
## 2. TB and Media for TST- (TBAM_NEG)
## 3. TB-Media for TST+ and TST- (TBMM)
## 4. Media only for TST+ and TST- (MEDIA)
## 5. TB only for TST+ and TST- (TB)

library(lumi)
library(gdata)
library(genefilter)

## Users should download data raw data from GSE76873 and extract results into a data frame containing probe names as rows and sample names as columns.
## Otherwise, the raw data output file from Illumina BeadStudio is the inputFile
## Load Data and generate 'lumiResults' data.frame containing raw intensity values
inputFile<-"2010.07.01.chethan.FinalReport_noRNAday5.txt"
decodeFile<-"targetsForSignalIntensityFiltering2.txt"  ## Available in SeshadriLab Github folder
d<-lumiR(inputFile)
dnorm<-lumiExpresso(d, bg.correct = FALSE, variance.stabilize=FALSE, normalize.param = list(method='quantile'))
lumiResults<-data.frame(dnorm@featureData$ProbeID, dnorm@featureData$TargetID, exprs(dnorm), check.names=F)
names(lumiResults)[1]<-c("ProbeID")
names(lumiResults)[2]<-c("TargetID")
names(lumiResults)<-trim(names(lumiResults))
names(lumiResults)

## Subset control data and perform signal intensity filtering
threshold<-0.75                                                # signal intensity threshold is 75%-ile of control probes 
dcontrol<-d@controlData
sampleName<-sampleNames(d)
table(dcontrol$controlType)
dcontrol<-dcontrol[dcontrol$controlType=="NEGATIVE",]          # Subset only negative control prbes
dcontrol<-dcontrol[,sampleName]          	                   # Drop ProbeID and TargetID
quantile.ctrl<-apply(dcontrol,2,quantile,probs=threshold)      # Obtain 75%-ile of control probes for each sample
quantile.ctrl<-log2(quantile.ctrl)
sfData<-lumiResults
sfData<-data.frame(sfData[3:ncol(sfData)],row.names=sfData$ProbeID,check.names=FALSE)
sfData<-log2(sfData)                                           # Rename rows and log2 transform expression values
cohort<-read.delim(decodeFile,stringsAsFactors=FALSE)          # Load file decoding sample and experimental cohort
cohort$Group<-as.factor(cohort$Group)
experiment<-data.frame(levels(cohort$Group),check.names=FALSE)
names(experiment)<-"Group"
experiment$Group<-as.character(experiment$Group)
experiment                                                     # There are four groups of samples (TSTPOS.media, TSTPOS.tb, TSTNEG.media, TSTNEG.tb)
sfD<-sfData                                                    # Backup data prior to signal intensity filter
for(i in 1:ncol(sfData)){                                      # Apply signal intensity filter
	cutoff<-quantile.ctrl[i]
	sfData[i][sfData[i] < cutoff,]<-NaN
	}

## Generate four data frames by grouping experimental samples using filtered data
for(i in 1:nrow(experiment)){
	subgroup<-cohort[cohort$Group==experiment[i,1],]
	subgroup.exprs<-sfData[,names(sfData)%in%subgroup$Sample]
	assign(as.character(experiment[i,1]),subgroup.exprs)
	}

percentPassToRetainProbe<-0.75                                 # Minimum % samples per group that must pass signal intensity filter to keep probe
filter.stats<-matrix(data=NA,nrow=nrow(d),ncol=nrow(experiment))
minSamples<-vector(mode="numeric",length=nrow(experiment))
signalFilteredNamesList<-as.list(1:nrow(experiment))
for(i in 1:nrow(experiment)){
	x<-get(experiment[i,1])                                    # Get expression values for group[i]
	minSamples[i]<-ceiling(ncol(x)*percentPassToRetainProbe)   # calculate number of samples needed to pass
	x[["Filter"]]<-with(x,ifelse(x[,c(1:ncol(x))] > 0, x$Filter<-1, x$Filter<-0))
	x$Sum<-apply(x$Filter,1,sum,na.rm=TRUE)                    # 'Filter' is a matrix within x which states whether an expression value > 0
	filter.stats[,i]<-x$Sum
	x<-subset(x,x$Sum >= minSamples[i])                        # Subset data according to filter
	x<-data.frame(x[1:(ncol(x)-2)], check.names=FALSE)         # Drop 'Filter and 'Sum' after filtering
	signalFilteredName<-paste(experiment[i,],"signalFiltered",sep="")
	assign(signalFilteredName,x)                               # Generates four filtered data frames
	signalFilteredNamesList[i]<-signalFilteredName
	}
	
## Plot distribution applied to signal intensity filter per group.  Data shows cutoff is conservative (more probes included than needed)
filter.stats2<-data.frame(filter.stats)
names(filter.stats2)<-experiment$Group
par(mfrow=c(2,2))
for(i in 1:ncol(filter.stats2)){
	hist(filter.stats[,i],main=experiment$Group[i],xlab="Cummulative Samples",ylab="Unique Probes")
	abline(v=minSamples[i],col="red")
	}
dev.off()

## Reformat signal intensity filtered dataframes to include ProbeIDs in first column and merge into a single dataframe
negtb<-data.frame(row.names(TST_NEG.tbsignalFiltered), TST_NEG.tbsignalFiltered)
names(negtb)[1]<-"ProbeID"
negmedia<-data.frame(row.names(TST_NEG.mediasignalFiltered), TST_NEG.mediasignalFiltered)
names(negmedia)[1]<-"ProbeID"
postb<-data.frame(row.names(TST_POS.tbsignalFiltered), TST_POS.tbsignalFiltered)
names(postb)[1]<-"ProbeID"
posmedia<-data.frame(row.names(TST_POS.mediasignalFiltered), TST_POS.mediasignalFiltered)
names(posmedia)[1]<-"ProbeID"
neg.m<-merge(negtb,negmedia,by.x="ProbeID",by.y="ProbeID",all.x=TRUE,all.y=TRUE)          # Union is larger because not all probes overlap
pos.m<-merge(postb,posmedia,by.x="ProbeID",by.y="ProbeID",all.x=TRUE,all.y=TRUE)
m<-merge(neg.m,pos.m,by.x="ProbeID",by.y="ProbeID",all.x=TRUE,all.y=TRUE)                 # Merged data contains all probes that pass signal intensity filter
m2<-sfD[row.names(sfD)%in%m$ProbeID,]                                                     # Subset data from backup object prior to signal intensity filter
m2<-2^m2																				  # Undo log2 conversion

## NO VARIANCE FILTERING!

## Collapse probes to genes (N=13970)
m3<-data.frame(row.names(m2),m2,check.names=FALSE)
pg<-lumiResults[1:2]
m4<-merge(pg,m3,by.x="ProbeID",by.y="ProbeID")
m4$TargetID<-as.character(m4$TargetID)
m4<-m4[order(m4$TargetID),]

abs.max<-function(x){                                                                     # Function abs.max will return the greatest magnitude regardless of sign
	x2<-abs(x)
	index<-grep(max(x2),x2)
	return(x[index])
	}

collapsed<-unlist(with(m4, tapply(m4[,3], full.lim3$TargetID, abs.max),check.rows=FALSE))
i <- 4
while(i<=ncol(m4)){
	temp<-unlist(with(m4, tapply(m4[,i], m4$TargetID, abs.max)))
	collapsed<-cbind(collapsed,temp)
	i <- i+1
	}
m5<-data.frame(row.names(collapsed),collapsed,check.names=FALSE,row.names=NULL)
colnames<-names(m4)
names(m5)<-colnames[2:length(colnames)]
names(m5)[1]<-"Gene"

m5[1:20,1:5]
lumiResults[lumiResults$TargetID=="ABCA1",c("ProbeID","TargetID","I-","I-TB")]             # Manual Double Check

## Group replicate samples together using filtered data.
for(i in 1:nrow(experiment)){
	subgroup<-cohort[cohort$Group==experiment[i,1],]
	subgroup.exprs<-m5[,names(m5)%in%subgroup$Sample]
	assign(paste(as.character(experiment[i,1]),"2",sep=""),subgroup.exprs)
	}

TBAM.POS<-cbind(m5$Gene,TST_POS.media2,TST_POS.tb2)
names(TBAM.POS)[1]<-"Gene"
TBAM.NEG<-cbind(m5$Gene,TST_NEG.media2,TST_NEG.tb2)
names(TBAM.NEG)[1]<-"Gene"
MEDIA<-cbind(m5$Gene,TST_NEG.media2,TST_POS.media2)
names(MEDIA)[1]<-"Gene"
TB<-cbind(m5$Gene,TST_NEG.tb2,TST_POS.tb2)
names(TB)[1]<-"Gene"

## Perform TB-media calculation for TST_POS and TST_NEG and generate TBMM dataframe.
tbmm.neg<-TST_NEG.tb2-TST_NEG.media2
tbmm.neg<-as.data.frame(tbmm.neg)
names(tbmm.neg)<-sub("TB","TBMM",names(tbmm.neg))
tbmm.pos<-TST_POS.tb2-TST_POS.media2
tbmm.pos<-as.data.frame(tbmm.pos)
names(tbmm.pos)<-sub("TB","TBMM",names(tbmm.pos))
names(tbmm.pos)

TBMM<-cbind(m5$Gene,tbmm.neg,tbmm.pos)
names(TBMM)[1]<-"Gene"

# Write data to file
write.table(TBMM,file="20140805_TBMM.gct",row.names=FALSE,quote=FALSE,sep="\t")
write.table(TBAM.POS,file="20140805_TBAMPOS.gct",row.names=FALSE,quote=FALSE,sep="\t")
write.table(TBAM.NEG,file="20140805_TBAMNEG.gct",row.names=FALSE,quote=FALSE,sep="\t")
write.table(MEDIA,file="20140805_MEDIA.gct",row.names=FALSE,quote=FALSE,sep="\t")
write.table(TB,file="20140805_TB.gct",row.names=FALSE,quote=FALSE,sep="\t")

## Output files need to be further modified to be consistent with GSEA input

## End Code	
