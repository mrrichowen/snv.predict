#' snv.cross.predict
#'
#' This function makes a snv prediction model on one set of data and then applies it to another, typically one without the ground truth information required to build a cell type specific model.
#' @param train_data Output from snv.bench. A data frame containing all derived snv classifiers with ground truths.
#' @param test_data Output from snv.bench. A data frame containing all derived snv classifiers (where training_set=FALSE, though not essential).
#' @param classifiers Output from snv.bench. A vector to keep low correlating classifiers only, all used if NULL.If applied to the training data, it will need to be used here.
#' @param threshold Defaults to 0.5. Values 0 to 1 can be input to tune predictions for specificity or sensitivity.
#' @keywords
#' @export
#' @examples
#' snv.cross.predict(train_data,test_data,classifiers,threshold)

snv.cross.predict <- function(train_data,test_data,classifiers=NULL,threshold=NULL){

  require(caret)
  require(ROCR)
  require(randomForest)

  x <- na.omit(train_data)
  x <- x[!duplicated(x$mutation),]
  names <- x[,c(1:4,ncol(x)-2,ncol(x)-1)]
  x <- x[,5:(ncol(x)-3)]

  if(is.null(classifiers)){
    classifiers <- as.data.frame(colnames(x),stringsAsFactors = F)
  }

  x <- x[,classifiers[,1]]
  y <- as.factor(names$snp_class)
  x <- as.matrix(x)

    bestmtry <- NULL
    mtry <- NULL
    for(i in seq(500,2500,500)){
      bestmtry <- as.data.frame(tuneRF(x, y, stepFactor=1.5,trace=F, plot=FALSE, improve=0.05, ntree=i))
      bestmtry$ntree <- i
      mtry <- rbind(mtry,bestmtry)
    }
  
  mtry <- mtry[which.min(mtry$OOBError),]
  model <- randomForest(x,y,ntree=mtry$ntree,mtry=mtry$mtry)

  x <- na.omit(test_data)
  x <- x[!duplicated(x$mutation),]
  names <- x[,c(1:4,ncol(x)-2,ncol(x)-1)]
  x <- x[,5:(ncol(x)-3)]
  x <- x[,classifiers[,1]]
  x <- as.matrix(x)

  model.prob <- predict(model, type="prob", newdata=x, probability=TRUE)
  z <- as.data.frame(model.prob)
  z <- cbind(names,model.prob)

  write.table(z,file="cross_predictions",row.names=F,quote=F,sep="\t")
}