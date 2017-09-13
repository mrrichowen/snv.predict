#' snv.predict
#'
#' This function predicts snvs in single cell RNA-seq data by creating a model based on true positives provided by the user and false positives from synthetic spike-ins (currently only ERCCs are supported). It will create a directory in the working  directory called snv.predict.results containing a table of snv predictions.
#' @param nuseq_classifiers List of named data frames as output from nuseq-sc.
#' @param truths List of named data.frames with snv information. Easiest to use the convenience function 'convert_to_pos' and add a fifth column with snvs annotated as true_positive (1), false_positive (0) and unknown (NA).
#' @keywords
#' @export
#' @examples
#' snv.predict(classifiers,truths)

snv.predict <- function(nuseq_classifiers,truths){

  require(caret)
  require(ROCR)
  require(randomForest)

data <- nuseq_classifiers

fun <- function(x){
  x <- na.omit(x)
}

data <- lapply(data,fun)

print("loading ground truths and checking samples match")

x <- truths

a <- NULL
z <- NULL
for(i in 1:length(x)){
  a <- x[[i]]
  z <- rbind(z,a)
}

gt <- z

print("generating new classifiers")

fun <- function(x){

  #Normalise counts by total depth
  x$alt_counts_norm <- (x$alt_counts+0.01)/(x$total_counts+0.01)
  x$alt_forward_counts_norm <- (x$alt_forward_counts+0.01)/(x$total_forward_counts+0.01)
  x$alt_reverse_counts_norm <- (x$alt_reverse_counts+0.01)/(x$total_reverse_counts+0.01)
  x$ref_counts_norm <- (x$ref_counts+0.01)/(x$total_counts+0.01)
  x$ref_forward_counts_norm <- (x$ref_forward_counts+0.01)/(x$total_forward_counts+0.01)
  x$ref_reverse_counts_norm <- (x$ref_reverse_counts+0.01)/(x$total_reverse_counts+0.01)

  #Ratios
  x$alt_bqual_ratio <- (x$alt_bqual+0.01)/(sum(x$alt_bqual)+0.01)
  x$alt_mqual_ratio <- (x$alt_mqual+0.01)/(sum(x$alt_mqual)+0.01)
  x$alt_distance_ratio <- (x$alt_distance+0.01)/(sum(x$alt_distance)+0.01)

  x$ref_bqual_ratio <- (x$ref_bqual+0.01)/(sum(x$ref_bqual)+0.01)
  x$ref_mqual_ratio <- (x$ref_mqual+0.01)/(sum(as.numeric(x$ref_mqual))+0.01)
  x$ref_distance_ratio <- (x$ref_distance+0.01)/(sum(x$ref_distance)+0.01)

  x$total_bqual_ratio <- (x$total_bqual+0.01)/(sum(x$total_bqual)+0.01)
  x$total_mqual_ratio <- (x$total_mqual+0.01)/(sum(as.numeric(x$total_mqual))+0.01)
  x$total_distance_ratio <- (x$total_distance+0.01)/(sum(x$total_distance)+0.01)

  x <- as.data.frame(x,stringsAsFactors=F)

  return(x)
}

x <- lapply(data,fun)

#Make training matrix
a <- NULL
z <- NULL
for(i in 1:length(x)){
  a <- x[[i]]
  a$sample_id <- names(x)[i]
  z <- rbind(z,a)
}
cm <- z

print("deduplicating and adding ground truths")

x <- paste(cm$chrom,cm$coord,cm$alt,sep=":")
gt <- na.omit(gt)
gt <- gt[gt$label==1,]
y <- paste(gt$chrom,gt$coord,gt$alt,sep=":")
calls <- x%in%y
cm[calls==TRUE,"snp_class"] <- "true_positive"
cm[grep("ERCC",cm$chrom),"snp_class"] <- "false_positive"
dups <- !duplicated(paste(cm$chrom,cm$coord,cm$alt,sep=":"))
data <- cm
cm <- cm[dups,]

print("filtering classifiers")

cor.m <- cor(cm[,c(5:(ncol(cm)-2))])
hc <- findCorrelation(cor.m, cutoff=0.75)

if(length(hc)<20){
  print("too many high correlating classifiers, using all")
  hc <- colnames(cm[,c(5:(ncol(x)-2))])
}

if(length(hc)>19){
  print("highly correlated classifiers removed, printing remaining classifiers used for modelling")
  print(colnames(cm[,hc+4]))
}

#Train the model
cm <- na.omit(cm)
y <- as.factor(cm$snp_class)
x <- cm[,hc+4]

print("tuning model")

bestmtry <- NULL
mtry <- NULL
for(i in seq(500,2500,500)){
  bestmtry <- as.data.frame(tuneRF(x, y, stepFactor=1.5,trace=F,plot=FALSE, improve=0.05, ntree=i))
  bestmtry$ntree <- i
  mtry <- rbind(mtry,bestmtry)
}
mtry <- mtry[which.min(mtry$OOBError),]
model <- randomForest(x,y,ntree=mtry$ntree,mtry=mtry$mtry)

print("predicting snvs")
x <- as.matrix(data[,5:(ncol(data)-2)])

model.prob <- predict(model, type="prob", newdata=x, probability=TRUE)

x <- as.data.frame(cbind(data[,c(1:4,ncol(data))],model.prob))

print("saving outputs")

write.table(x,file=paste0(getwd(),"/predicted_snps"),row.names=F,quote=F,sep="\t")
write.table(plot,file=paste0(getwd(),"/plot_ROC"),row.names=F,quote=F,sep="\t")

print("output data saved")

}