#' snv.bench
#'
#' This is a convenience function that takes the output files from nuseq-sc and combines them with a list of data frames with annotated snvs for input into other snv.predict functions. It will output a single snv.summary file concatenating all results and a vector of the low correlating classifiers recommended for subsequent analysis.
#' @param nuseq_classifiers List of named data frames as output from nuseq-sc.
#' @param truths List of data.frames with snv information. Easiest to use the convenience function 'convert_to_pos' and add a fifth column ('labels') with snvs annotated as true_positive (1), false_positive (0) and unknown (NA).
#' @param training_set If TRUE, only snvs with ground truth data will be output. This is useful if this data is to be used for training a model to be applied to a different dataset with snv.cross.predict.
#' @keywords
#' @export
#' @examples
#' snv.bench(nuseq_classifiers,truths,training_set=FALSE)

snv.bench <- function(nuseq_classifiers,truths,training_set=FALSE){

  require(caret)

  data <- truths

  fun <- function(x){
    x <- na.omit(x)
    x$mut <- paste(x$chrom,x$coord,x$alt,sep=":")
  }

  data <- lapply(data,fun)

  a <- NULL
  z <- NULL
  for(i in 1:length(data)){
    a <- data[[i]]
    z <- c(z,a)
  }

  mut_filter <- z

  data <- nuseq_classifiers

  fun <- function(x){
    x <- na.omit(x)
    x <- as.data.frame(x,stringsAsFactors=F)
  }

  data <- lapply(data,fun)

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

    return(x)
  }

  x <- lapply(data,fun)
  data <- x

  #Make one matrix for all data and filter mutations
  a <- NULL
  z <- NULL
  for(i in 1:length(x)){
    a <- x[[i]]
    a$sample_id <- names(x)[i]
    a$mutation <- paste(a$chrom,a$coord,a$alt,sep=":")
    a[,!apply(is.na(a[,1:ncol(a)-1]), 2, any)]
    a[a$mutation%in%mut_filter,"snp_class"] <- "true_positive"
    a[grep("ERCC",a$chrom),"snp_class"] <- "false_positive"
    z <- rbind(z,a)
  }

  x <- z

  data <- x

  if(training_set==TRUE){
    data <- na.omit(data)
    write.table(data,file="snv.summary",quote=F,row.names=F,sep="\t")
  }

  if(training_set==FALSE){
    write.table(data,file="snv.summary",quote=F,row.names=F,sep="\t")
  }

  x <- data
  x <- x[,c(5:(ncol(x)-3))]
  print(colnames(x))
  cm <- cor(x)
  hc <- findCorrelation(cm, cutoff=0.75)

  if(length(hc)<20){
    print("too many high correlating classifiers, using all")
  }

  if(length(hc)>19){
    print("highly correlated classifiers removed, remaining classifiers recommended for model")
    x <- x[,hc]
    classifiers <- colnames(x)
    write.table(classifiers,file="low_correlation_classifiers",quote=F,row.names=F,sep="\t")
  }
}