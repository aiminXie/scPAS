
#' scPAS : A tool for identifying Phenotype-Associated cell Subpopulations from single-cell sequencing data by integrating bulk data
#'
#' @param bulk_dataset Matrix. Bulk expression matrix of related disease. Each row represents a gene and each column represents a sample. The input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.
#' @param sc_dataset Matrix or seurat object. Single-cell RNA-seq expression matrix of related disease. Each row represents a gene and each column represents a sample. A Seurat object that contains the preprocessed data and constructed network is preferred. Otherwise, a cell-cell similarity network is constructed based on the input matrix.Otherwise, the raw count expression matrix will be processed by using Seurat's default parameters. See run_Seurat for details.
#' @param phenotype Phenotype annotation of each bulk sample. It can be a continuous dependent variable,
#' binary group indicator vector, or clinical survival data:
#'   \itemize{
#'   \item Continuous dependent variable. Should be a quantitative vector for \code{family = gaussian}.
#'   \item Binary group indicator vector. Should be either a 0-1 encoded vector or a factor with two levels for \code{family = binomial}.
#'   \item Clinical survival data. Should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable,
#'   with '1' indicating event (e.g.recurrence of cancer or death), and '0' indicating right censored.
#'   The function \code{Surv()} in package survival produces such a matrix.
#'   }
#' @param assay Name of Assay to get.
#' @param tag Names for each phenotypic group. Used for logistic regressions only.
#' @param nfeature Numeric. The Number of features to select as top variable features in sc_dataset. Top variable features will be used to intersect with the features of bulk_dataset. Default is NULL.All features will be used.
#' @param imputation Logical. imputation or not.
#' @param imputation_method Character. Name of alternative method for imputation.
#' @param alpha Numeric. Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector.
#' If \code{alpha = NULL}, a default searching vector is used. The range of alpha is in \code{[0,1]}. A larger alpha lays more emphasis on the l1 norm.
#' @param network_class  The source of feature-feature similarity network. By default this is set to \code{sc} and the other one is \code{bulk}.
#' @param cutoff Numeric. Cutoff for the percentage of the scPAS selected features in total features. This parameter is used to restrict the number of the
#' scPAS selected features. A cutoff less than \code{50\%} (default \code{20\%}) is recommended depending on the input data.
#' @param family Character. Response type for the regression model.  It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification, or \code{family = cox} for Cox regression.
#' @param FDR.threshold Numeric. FDR value threshold for identifying phenotype-associated cells.
#' The default is 0.05.
#' @param independent Logical. The background distribution of risk scores is constructed independently of each cell.
#'
#' @return This function returns a Seurat object with the following components added to :
#'   \item{scPAS_para}{ A list contains the final model parameters added to misc.}
#'   \item{PAS result}{ A data frame containing risk scores (scPAS_RS), normalized risk scores (scPAS_NRS), p-value (scPAS_Pvalue) , adjusted p-value (scPAS_FDR) cell classification labels (scPAS) added to metaData.}
#'
#' @import Seurat Matrix preprocessCore
#'
#' @export
scPAS <- function(bulk_dataset, sc_dataset, phenotype,assay = 'RNA', tag = NULL,nfeature = NULL,imputation=T,imputation_method=c('KNN','ALRA'),
                    alpha = NULL,network_class=c('SC','bulk'),independent=T, family = c("gaussian","binomial","cox"),permutation_times=2000,
                    FDR.threshold = 0.05){
  library(Seurat)
  library(Matrix)
  library(preprocessCore)
 # library(APML0)
  network_class=match.arg(network_class)
  family=match.arg(family)
  imputation_method =match.arg(imputation_method)
  DefaultAssay(sc_dataset) <- assay

  if(any(class(sc_dataset)=='Seurat')){

    if(is.null(nfeature)){
      common <- intersect(rownames(bulk_dataset),rownames(sc_dataset))
    }
    else if(is.numeric(nfeature) & length(nfeature)==1){
      sc_dataset <- FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = F,nfeatures = nfeature)

      common <- intersect(rownames(bulk_dataset),VariableFeatures(sc_dataset))
      common <- common[!grepl(pattern ='^RP[LS]',common)]
      common <- common[!grepl(pattern ='^MT-',common)]
    }else if(is.character(nfeature) & length(nfeature)>1){
      common <- intersect(rownames(bulk_dataset),rownames(sc_dataset))
      common <- intersect(common,nfeature)
      common <- common[!grepl(pattern ='^RP[LS]',common)]
      common <- common[!grepl(pattern ='^MT-',common)]
    }
  }else{
    print("Step 0:The single-cell data is not a Seurat object, and a default Seurat pipeline will be run.")
    sc_dataset <- run_Seurat(sc_dataset)
    if(is.null(nfeature)){
      common <- intersect(rownames(bulk_dataset),rownames(sc_dataset))
      common <- common[!grepl(pattern ='^RP[LS]',common)]
      common <- common[!grepl(pattern ='^MT-',common)]
    }
    else if(is.numeric(nfeature)){
      sc_dataset <- FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = F,nfeatures = nfeature)
      common <- intersect(rownames(bulk_dataset),VariableFeatures(sc_dataset))
    }else if(is.character(nfeature) & length(nfeature)>1){
      common <- intersect(rownames(bulk_dataset),rownames(sc_dataset))
      common <- intersect(common,nfeature)
      common <- common[!grepl(pattern ='^RP[LS]',common)]
      common <- common[!grepl(pattern ='^MT-',common)]
    }
  }
  if (length(common) == 0) {
    stop("There is no common genes between the given single-cell and bulk samples.")
  }

  print("Step 1:Quantile normalization of bulk data.")
  Expression_bulk <- preprocessCore::normalize.quantiles(as.matrix(bulk_dataset))
  rownames(Expression_bulk) <- rownames(bulk_dataset)
  colnames(Expression_bulk) <- colnames(bulk_dataset)
  Expression_bulk <- Expression_bulk[common,]


  if(imputation){
    sc_dataset <- imputation(sc_dataset,assay = assay,method = imputation_method)
    assay <- DefaultAssay(sc_dataset)
  }

  print("Step 2: Extracting single-cell expression profiles....")
  sc_exprs <- GetAssayData(object = sc_dataset, assay = assay,slot = 'data')

  Expression_cell <- sc_exprs
  rownames(Expression_cell) <- rownames(sc_exprs)
  colnames(Expression_cell) <- colnames(sc_exprs)
  Expression_cell <- Expression_cell[common,]
  sc_exprs <- NULL
  bulk_dataset <- NULL

  x <- t(Expression_bulk)


  # Construct a gene network
  if(network_class=='bulk'){
    print("Step 3: Constructing a gene-gene similarity by bulk data....")
    cor.m <- cor(x)

  }else{
    print("Step 3: Constructing a gene-gene similarity by single cell data....")
    cor.m <- sparse.cor(t(Expression_cell))
  }
  cor.m[which(cor.m<0)] <- 0
  SNN <- FindNeighbors(1-cor.m,distance.matrix=T)
  Network <- as.matrix(SNN$snn)
  diag(Network) <- 0
  Network[which(Network >0.2)] <- 1
  Network[which(Network <= 0.2)] <- 0

  print("Step 4: Optimizing the network-regularized sparse regression model....")
  if (family == "binomial"){
    y <- as.numeric(phenotype)
    z <- table(y)
    print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
    print("Perform logistic regression on the given phenotypes:")

  }
  if (family == "gaussian"){
    y <- as.numeric(phenotype)
    print("Perform linear regression on the given phenotypes:")
  }
  if (family == "cox"){
    y <- as.matrix(phenotype)
    if (ncol(y) != 2){
      stop("The size of survival data is wrong. Please check inputs and selected regression type.")
    }else{
      print("Perform cox regression on the given clinical outcomes:")
    }
  }

  if (is.null(alpha)){
    alpha <- c(0.001,0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  }
  lambda <- c()
  for (i in 1:length(alpha)){
    set.seed(123)
    fit0 <- APML0(x = x,y = y,family = family,penalty = 'Net',Omega = Network,alpha = alpha[i],nlambda = 100,nfolds = min(10,nrow(x)))
    fit1 <- APML0(x = x,y = y,family = family,penalty = 'Net',Omega = Network,alpha = alpha[i],lambda = fit0$lambda.min)
    lambda <- c(lambda,fit0$lambda.min)
    if (family == "binomial"){
      Coefs <- as.numeric(fit1$Beta[2:(ncol(x)+1)])
    }else{
      Coefs <- as.numeric(fit1$Beta)
    }
    Feature1 <- colnames(x)[which(Coefs > 0)]
    Feature2 <- colnames(x)[which(Coefs < 0)]
    percentage <- (length(Feature1) + length(Feature2)) / ncol(x)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("lambda = %s", fit0$lambda.min))
    print(sprintf("scPAS identified %d rick+ features and %d rick- features.", length(Feature1), length(Feature2)))
    print(sprintf("The percentage of selected feature is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))

    if (percentage < cutoff){
      break
    }
    cat("\n")
  }
  print("|**************************************************|")

  print("Step 5: calculating quantified risk scores....")
  names(Coefs) <- colnames(x)
  scaled_exp <- Seurat:::FastSparseRowScale(Expression_cell,display_progress = F)

  colnames(scaled_exp) <- colnames(Expression_cell)
  rownames(scaled_exp) <- rownames(Expression_cell)
  scaled_exp <- as(scaled_exp, "sparseMatrix")
  risk_score <- crossprod(scaled_exp,Coefs)

  print(paste0("Step 6: qualitative identification by permutation test program with ", as.character(permutation_times), " times random perturbations"))

  set.seed(12345)

  randomPermutation <- sapply(1:permutation_times,FUN = function(x){
    set.seed(1234+x)
    sample(Coefs,length(Coefs),replace = F)
  })
  randomPermutation <- as(randomPermutation, "sparseMatrix")
  risk_score.background <- crossprod(scaled_exp,randomPermutation)
  if(independent){
    mean.background <- rowMeans(risk_score.background)
    sd.background <- apply(risk_score.background,1,sd)
  }else{
    mean.background <- mean(as.matrix(risk_score.background))
    sd.background <- sd(as.matrix(risk_score.background))
  }



  Z <- (risk_score[,1]-mean.background)/sd.background

  p.value <- pnorm(q = abs(Z),mean = 0,sd = 1,lower.tail = F)
  q.value <- p.adjust(p = p.value,method = 'BH')
  risk_score_data.frame <- data.frame(cell=colnames(Expression_cell),
                                      raw_score=risk_score[,1],
                                      Z.statistics=Z,
                                      p.value=p.value,
                                      FDR=q.value
  )



  risk_score_data.frame$cell_lable <- ifelse(Z>0 & q.value<=FDR.threshold,'scPAS+',ifelse(Z<0 & q.value<=FDR.threshold,'scPAS-','0'))

  sc_dataset@misc$scPAS_para <- list(alpha = alpha[1:i], lambda = lambda, family = family,Coefs=Coefs,bulk=x,phenotype=y,Network=Network)

  sc_dataset <- AddMetaData(sc_dataset, metadata = risk_score_data.frame$raw_score, col.name = "scPAS_RS")

  sc_dataset <- AddMetaData(sc_dataset, metadata = risk_score_data.frame$Z.statistics, col.name = "scPAS_NRS")

  sc_dataset <- AddMetaData(sc_dataset, metadata = risk_score_data.frame$p.value, col.name = "scPAS_Pvalue")
  sc_dataset <- AddMetaData(sc_dataset, metadata = risk_score_data.frame$FDR, col.name = "scPAS_FDR")
  sc_dataset <- AddMetaData(sc_dataset, metadata = risk_score_data.frame$cell_lable, col.name = "scPAS")


  print("Finished.")
  return(sc_dataset)

}

#' Preprocess the single-cell raw data using functions in the \code{Seurat} package
#'
#' This function provide a simplified-version of Seurat analysis pipeline for single-cell RNA-seq data. It contains the following steps in the pipeline:
#' \itemize{
#'    \item Create a \code{Seurat} object from raw data.
#'    \item Normalize the count data present in a given assay.
#'    \item Identify the variable features.
#'    \item Scales and centers features in the dataset.
#'    \item Run a PCA dimensionality reduction.
#'    \item Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
#'    \item Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.
#'    \item Run t-distributed Stochastic Neighbor Embedding (t-SNE) dimensionality reduction on selected features.
#'    \item Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique.
#' }
#'
#' @param counts A \code{matrix}-like object with unnormalized data with cells as columns and features as rows.
#' @param project Project name for the \code{Seurat} object.
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#' @param normalization.method Method for normalization.
#'   \itemize{
#'   \item LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   This is then natural-log transformed using log1p.
#'   \item CLR: Applies a centered log ratio transformation.
#'   \item RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   No log-transformation is applied. For counts per million (CPM) set \code{scale.factor = 1e6}.
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Choose one of :
#'   \itemize{
#'   \item vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
#'   Then standardizes the feature values using the observed mean and expected variance (given by the fitted line).
#'   Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
#'   \item mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function)
#'   for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates
#'   z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong
#'   relationship between variability and average expression.
#'   \item dispersion (disp): selects the genes with the highest dispersion values
#'   }
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param dims_Neighbors Dimensions of reduction to use as input.
#' @param dims_TSNE Which dimensions to use as input features for t-SNE.
#' @param dims_UMAP Which dimensions to use as input features for UMAP.
#' @param meta.data meta data of single cell data.
#' @param verbose Print output.
#'
#' @return A \code{Seurat} object containing cell-cell similarity network, t-SNE and UMAP representations.
#' @import Seurat
#' @export
run_Seurat <- function(counts, project = "Single_Cell", min.cells = 400, min.features = 200,meta.data =NULL,
                       normalization.method = "LogNormalize", scale.factor = 10000,
                       selection.method = "vst", resolution = 0.6,
                       dims_Neighbors = 1:10, dims_TSNE = 1:10, dims_UMAP = 1:10,
                       verbose = TRUE){
  library(Seurat)
  if('matrix'  %in% class(counts) | 'dgCMatrix' %in% class(counts)){
    data <- CreateSeuratObject(counts = counts, project = project, min.cells = min.cells, min.features = min.features,meta.data =meta.data)
  }else if('Seurat' %in% class(counts)){
    data <- counts
  }else{
    stop("The class of scRNA-seq data is wrong. Please input a count matrix or Seurat object")
  }
  data <- NormalizeData(object = data, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
  data <- FindVariableFeatures(object = data, selection.method = selection.method, verbose = verbose)
  data <- ScaleData(object = data, verbose = verbose)
  data <- RunPCA(object = data, features = VariableFeatures(data), verbose = verbose)
  data <- FindNeighbors(object = data, dims = dims_Neighbors, verbose = verbose)
  data <- FindClusters( object = data, resolution = resolution, verbose = verbose)
  data <- RunTSNE(object = data, dims = dims_TSNE)
  data <- RunUMAP(object = data, dims = dims_UMAP, verbose = verbose)

  return(data)
}

#' The function of imputaion.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#' @param method The method for imputation. The default is 'RNA'.
#'
#' @return  A seurat object after imputaion.
#'
#' @import Matrix
#'
#'
#'
imputation <- function(obj,assay='RNA',method=c('KNN','ALRA')){
  method=match.arg(method)
  if(method=='KNN'){
    print("Step2: Imputation of missing values in single cell RNA-sequencing data with KNN....")
    obj <- imputation_KNN(obj = obj,assay = assay)
  }else if(method=='ALRA'){
    print("Step2: Imputation of missing values in single cell RNA-sequencing data with ALRA")
    obj <- imputation_ALRA(obj = obj,assay = assay)
  }else{
    warnings(paste0('The ',method, 'method does not exist, so imputaion is invalid!'))
  }
  return(obj)
}


#' A method for imputation of missing values in single cell RNA-sequencing data based on ALRA.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#'
#' @return  A seurat object after imputaion.
#'
#' @import Matrix
#'
#'
#'
imputation_ALRA <- function(obj,assay='RNA'){
  library(ALRA)
  library(Matrix)
  library(Seurat)
  data <- GetAssayData(object = obj, assay = assay,slot = 'data')
  data_alra <- t(alra(t(as.matrix(data)))[[3]])
  colnames(data_alra) <- colnames(data)
  data_alra <- Matrix(data_alra, sparse = T)

  obj[["imputation"]] <- CreateAssayObject(data = data_alra )
  DefaultAssay(obj) <- "imputation"
  return(obj)
}

#' A method for imputation of missing values in single cell RNA-sequencing data based on the average expression value of nearest neighbor cells.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#' @param LogNormalized Whether the data is LogNormalized.
#'
#' @return A seurat object after imputaion.
#'
#' @import Matrix
#'
#'
#'
imputation_KNN <- function (obj,assay='RNA', LogNormalized = T)
{
  library(Matrix)
  exp_sc <- GetAssayData(object = obj, assay = assay,slot = 'data')
  nn_network <- obj@graphs[[paste0(assay, "_nn")]]
  #nn_network <- obj@graphs$RNA_nn
  if (!is(object = exp_sc, class2 = "sparseMatrix")) {
    exp_sc <- as(exp_sc, "sparseMatrix")
  }
  if (!is(object = nn_network, class2 = "sparseMatrix")) {
    nn_network <- as(nn_network, "sparseMatrix")
  }
  if (LogNormalized) {
    exp_sc <- as(exp(exp_sc) - 1, "sparseMatrix")
  }
  network_count <- as(Diagonal(x = 1/rowSums(nn_network)),
                      "sparseMatrix")
  exp_sc_mean <- tcrossprod(x = tcrossprod(x = exp_sc, y = nn_network),
                            y = network_count)
  if (LogNormalized) {
    exp_sc_mean <- log1p(exp_sc_mean)
  }
  colnames(exp_sc_mean) <- colnames(exp_sc)
  obj[["imputation"]] <- CreateAssayObject(data = exp_sc_mean )
  DefaultAssay(obj) <- "imputation"
  return(obj)
}


#' A function compute the correlation of a sparse matrix.
#'
#' @param x Matrix. Normalized single cell expression profile extracted from Seurat object.
#'
#' @return A correlation matrix.
#'
#'
#'
#'
sparse.cor <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element

  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans

  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
                crossprod(t(Ex))*(n-length(ii))
  )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
}



#' scPAS.prediction: A function that uses the scPAS model to make predictions on independent data
#'
#' @param model seurat object. A seurat object containing the scPAS model
#' @param test.data Matrix or seurat object. Single-cell RNA-seq expression matrix of related disease. Each row represents a gene and each column represents a sample. A Seurat object that contains the preprocessed data and constructed network is preferred.
#' @param assay Name of Assay to get.
#' @param FDR.threshold Numeric. FDR value threshold for identifying phenotype-associated cells.
#' The default is 0.05.
#'
#'@return A seurat object or data frame containing the forecast results.
#'
#' @export
scPAS.prediction <- function(model, test.data,assay='RNA', FDR.threshold=0.05,imputation=F,imputation_method= 'KNN',independent=T){


  model <- model@misc$scPAS_para

  if(any(class(test.data)=='Seurat')){
    if(imputation){
      test.data <- imputation(test.data,assay = assay,method = imputation_method)
      assay <- DefaultAssay(test.data)
    }
    test.exp <- GetAssayData(object = test.data, assay = assay,slot = 'data')
    Expression_cell <- test.exp
    rownames(Expression_cell) <- rownames(test.exp)
    colnames(Expression_cell) <- colnames(test.exp)


  }else{
    test.exp <- as.matrix(test.data)
    Expression_cell <-  as(preprocessCore::normalize.quantiles(as.matrix(test.exp)),'dgCMatrix')
    rownames(Expression_cell) <- rownames(test.exp)
    colnames(Expression_cell) <- colnames(test.exp)
  }

  Coefs <- model$Coefs
  common <- intersect(names(Coefs),rownames(Expression_cell))

  if(sum(Coefs!=0)<20){
    stop("There are too few valid features and the test data may not be suitable for the model!")
  }


  Coefs <- Coefs[common]
  Expression_cell <- Expression_cell[common,]
  scaled_exp <- Seurat:::FastSparseRowScale(Expression_cell,display_progress = F)
  colnames(scaled_exp) <- colnames(Expression_cell)
  rownames(scaled_exp) <- rownames(Expression_cell)
  #scaled_exp <- as(scaled_exp, "sparseMatrix")

  scaled_exp[which(is.na(scaled_exp))] <- 0
  risk_score <- crossprod(scaled_exp,Coefs)

  set.seed(12345)

  randomPermutation <- sapply(1:2000,FUN = function(x){
    set.seed(1234+x)
    sample(Coefs,length(Coefs),replace = F)
  })
  randomPermutation <- as(randomPermutation, "sparseMatrix")
  risk_score.background <- crossprod(scaled_exp,randomPermutation)

  if(independent){
    mean.background <- rowMeans(risk_score.background)
    sd.background <- apply(risk_score.background,1,sd)
  }else{
    mean.background <- mean(as.matrix(risk_score.background))
    sd.background <- sd(as.matrix(risk_score.background))
  }

  Z <- (risk_score[,1]-mean.background)/sd.background


  p.value <- pnorm(q = abs(Z),mean = 0,sd = 1,lower.tail = F)
  q.value <- p.adjust(p = p.value,method = 'BH')

  risk_score_data.frame <- data.frame(sample=colnames(Expression_cell),
                                      scPAS_RS=risk_score[,1],
                                      scPAS_NRS=Z,
                                      scPAS_Pvalue=p.value,
                                      scPAS_FDR=q.value
  )
  risk_score_data.frame$scPAS <- ifelse(Z>0 & q.value<=FDR.threshold,'scPAS+',ifelse(Z<0 & q.value<=FDR.threshold,'scPAS-','0'))



  if(any(class(test.data)=='Seurat')){

    test.data <- AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_RS, col.name = "scPAS_RS")

    test.data <- AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_NRS, col.name = "scPAS_NRS")

    test.data <- AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_Pvalue, col.name = "scPAS_Pvalue")
    test.data <- AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_FDR, col.name = "scPAS_FDR")
    test.data <- AddMetaData(test.data, metadata = risk_score_data.frame$scPAS, col.name = "scPAS")
    return(test.data)
  }else{
    return(risk_score_data.frame)
  }



}
