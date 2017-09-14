#' snv.performance
#'
#' This function takes the output files from snv.benchmark, makes a snv prediction model and tests it on hold out data. The output is a confusion matrix, a set of snv predictions and a table for plotting a ROC curve (all applied to test data only, AUC is output in STDOUT).
#' @param snv_data Output from snv.benchmark. A data frame containing all derived snv classifiers with ground truth data appended.
#' @param test_index A vector to determine which samples to train the model on (1) and then test (2). This can be user input to allow consistency between modelling methods if required.
#' @param classifiers Defaults to use all classifiers. We recommend providing a vector of the low correlating classifiers to use in building the prediction model. This is already output from snv.bench.
#' @param threshold Defaults to 0.5. Values 0 to 1 can be input to tune predictions for specificity or sensitivity.
#' @param model.method Defaults to standard random forest, other options are "tuned_randomforest", "svm", and "tuned_svm".
#' @keywords
#' @export
#' @examples
#' snv.performance(snv_data,classifiers,test_index=NULL,threshold=NULL,model.method="randomforest")

snv.performance <- function(snv_data,classifiers=NULL,test_index=NULL,threshold=NULL,model.method="randomforest"){

  require(ROCR)
  require(ggplot2)
  require(cowplot)
  require(caret)

  if(is.null(threshold)){
    threshold <- 0.5
  }

  data <- na.omit(snv_data)
  dups <- !duplicated(data$mutation)

  data <- data[dups,]

  if(is.null(test_index)){
    set.seed(2)
    ind <- sample(2, nrow(data), replace = TRUE, prob=c(0.5, 0.5))
    write.table(ind,file="test_index",quote=F,row.names=F,sep="\t")
  }

  if(!is.null(test_index)){
    ind <- test_index
    ind <- ind[dups,]
  }

  x <- data[ind==1,]
  names <- x[,c(1:4,ncol(x)-2,ncol(x))]
  x <- x[,5:(ncol(x)-3)]
  x <- x[,classifiers[,1]]
  y <- as.factor(names$snp_class)
  x <- as.matrix(x)

  #run model
  if(model.method=="randomforest"){
    require(randomForest)
    model <- randomForest(x,y)
  }

  if(model.method=="tuned_randomforest"){
    require(randomForest)
    bestmtry <- NULL
    mtry <- NULL
    for(i in seq(500,2500,500)){
      bestmtry <- as.data.frame(tuneRF(x, y, stepFactor=1.5,trace=F, improve=0.05, ntree=i))
      bestmtry$ntree <- i
      mtry <- rbind(mtry,bestmtry)
    }
    mtry <- mtry[which.min(mtry$OOBError),]
    model <- randomForest(x,y,ntree=mtry$ntree,mtry=mtry$mtry)
  }

  if(model.method=="svm"){
    require(e1071)
    model <- svm(x,y,probability=T)
  }

  if(model.method=="tuned_svm"){
    require(e1071)
    svm.tune <- tune(svm,x,y,ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),tunecontrol = tune.control(sampling = "fix"))
    model <- svm(x,y, cost=svm.tune[[1]][[2]], gamma=svm.tune[[1]][[1]],probability=TRUE)
  }

  #Predict with test data
  x <- data[ind==2,]
  names <- x[,c(1:4,ncol(x)-2,ncol(x))]
  x <- x[,5:(ncol(x)-3)]

  if(is.null(classifiers)){
    classifiers <- colnames(x)
  }

  x <- x[,classifiers[,1]]
  y <- as.factor(names$snp_class)
  x <- as.matrix(x)

  model.prob <- predict(model, type="prob", newdata=x, probability=TRUE)

  if(grepl("randomforest",model.method)==TRUE){
    model.pred <- prediction(model.prob[,2], names$snp_class)
  }

  if(grepl("svm",model.method)==TRUE){
    model.pred <- prediction(attr(model.prob, "probabilities")[,1], y)
  }

  model.cm <- confusionMatrix(as.factor(model.pred@predictions[[1]]>threshold),y=="true_positive")

  z <- as.data.frame(model.prob)
  z$snp_class <- names$snp_class
  z$label <- "test"

  model.perf <- performance(model.pred,"tpr","fpr")
  auc <- performance(model.pred, measure = "auc")
  print(auc)

  plot <- as.data.frame(cbind(model.perf@x.values[[1]],model.perf@y.values[[1]]))
  colnames(plot) <- c("FPR","TPR")

  write.table(z,file="predictions",row.names=F,quote=F,sep="\t")
  write.table(plot,file="plot_ROC",row.names=F,quote=F,sep="\t")
  write.table(model.cm[[2]],file="confusion_matrix",quote=F,sep="\t")
}
