#' snv.predict
#'
#' This function predicts snvs in single cell RNA-seq data by creating a model based on true positives provided by the user and false positives from synthetic spike-ins (currently only ERCCs are supported). It will create a table of snv predictions.
#' @param nuseq_classifiers List of named data frames as output from nuseq-sc.
#' @param truths List of named data.frames with snv information. Easiest to use the convenience function 'convert_to_pos' and add a fifth column with snvs annotated as true_positive (1), false_positive (0) and unknown (NA).
#' @keywords
#' @export
#' @examples
#' snv.predict(classifiers,truths)

snv.predict <- function(snv.summary,classifiers=NULL){

  require(caret)
  require(ROCR)
  require(randomForest)

x <- na.omit(snv.summary) 
y <- as.factor(x$snp_class)
  
if(!is.null(classifiers)){
    x <- x[,classifiers]
  }
  
if(is.null(classifiers)){
    x <- x[,]
  }
  
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
