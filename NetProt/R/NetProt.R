# Copyright (C) 2017  Wilson Wen Bin GOH
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#library(testthat)
#library(devtools)
#library(roxygen2)
#to build a vignette directory
#use_vignette("NetProt")
#cmd + shift + k to execute the vignette generation
#run devtools::document() to use roxygen2 for namespace generation

library("e1071")
library("genefilter")
library("pheatmap")
library("limma")
library("pvclust")
library("vioplot")

#' @title This is the RCC renal cancer control data to be included in my package
#'
#' @name RCC
#' @docType data
#' @author Wilson Goh \email{goh.informatics@gmail.com}
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#' @keywords data
NULL

#' @title This is the RC renal cancer data to be included in my package
#'
#' @name RC
#' @docType data
#' @author Wilson Goh \email{goh.informatics@gmail.com}
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#' @keywords data
NULL

#' @title This is the CR colorectal cancer data to be included in my package
#'
#' @name CR
#' @docType data
#' @author Wilson Goh \email{goh.informatics@gmail.com}
#' @references \url{http://www.nature.com/nm/journal/v21/n4/full/nm.3807.html}
#' @keywords data
NULL

#generic functions --- used by several methods
#' @title A test function for the t-test
#' @name my.t.test.p.value
my.t.test.p.value <- function(...)
{
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#' @title Gene Fuzzy Scoring (GFS)
#' @name gfs
#' @param vec a vector of values corresponding to one sample
#' @return new_vec A vector of GFS transformed values
#' @export
gfs <- function(vec)
{
  #alpha_1 <- x * length(vec)
  #alpha_2 <- y * length(vec)
  alpha_1 <- 0.1 * length(vec)
  alpha_2 <- 0.2 * length(vec)

  ranks <- rank(-vec)

  new_vec <- vec

  new_vec[ranks[1:round(alpha_1)]] <- 1
  new_vec[ranks[round(alpha_2)+1:length(vec)]] <- 0
  new_vec[ranks[round(alpha_1)+1:(round((alpha_2 - alpha_1)/5))]] <- 0.8 #bin 1
  new_vec[ranks[(round(alpha_1 + (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.6 #bin 2
  new_vec[ranks[(round(alpha_1 + 2 * (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.4 #bin 3
  new_vec[ranks[(round(alpha_1 + 3 * (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.2 #bin 4
  new_vec[ranks[(round(alpha_1 + 4 * (alpha_2 - alpha_1)/5))+1:(round((alpha_2 - alpha_1)/5))]] <- 0.1 #bin 5

  return(new_vec)
}


#' @title Gene No Fuzzy Scoring (GNFS)
#' @name gnfs
#' @param vec a vector of values corresponding to one sample
#' @return new_vec A vector of GNFS transformed values
#' @export
gnfs <- function(vec)
{
  #alpha_1 <- x * length(vec)
  #alpha_2 <- y * length(vec)
  alpha_1 <- 0.1 * length(vec)
  #alpha_2 <- 0.2 * length(vec)

  ranks <- rank(-vec)

  new_vec <- vec

  new_vec[ranks[1:round(alpha_1)]] <- 1
  new_vec[ranks[round(alpha_1)+1:length(vec)]] <- 0

  return(new_vec)
}

#' @title Matrix Binarization
#' @name binarize_mat
#' @param vec a data vector
#' @return new_mat A binarized data vector
#' @export
binarize_mat <- function(vec)
{
  vec[!is.na(vec)] <- 1
  vec[is.na(vec)] <- 0
  return(vec)
}

#' @title This generates a vector of p-values for a data matrix using the two-sample t-test
#' @name standard_t_test
#' @param x data matrix
#' @param y significance level
#' @return id_output A list of 2 p-value vectors, x vs y and y vs x
#' @export
standard_t_test <- function(x,y)
{

  ttest_pvals_1 <- c()
  ttest_pvals_2 <- c()

  #lets play with the first matrix
  classes_1 <- as.factor(matrix(unlist(strsplit(colnames(x[[1]]), "_")), ncol=ncol(x[[1]]))[1,])      #factors of the rank weight matrix
  cancer_class_1 <- x[[1]][, grep(levels(classes_1)[1], colnames(x[[1]]))]
  normal_class_1 <- x[[1]][, grep(levels(classes_1)[2], colnames(x[[1]]))]

  for (i in 1:nrow(x[[1]]))
  {
    test_1 <- my.t.test.p.value(normal_class_1[i,], cancer_class_1[i,], alternative="greater",var.equal = FALSE,paired = FALSE)
    ttest_pvals_1 <- append(ttest_pvals_1, test_1)
  }
  ttest_pvals_1[is.na(ttest_pvals_1)] <- 1
  names(ttest_pvals_1) <- rownames(x[[1]])

  #out1 <- x[[1]][names(ttest_pvals_1)[ttest_pvals_1 <= y],]

  #############
  classes_2 <- as.factor(matrix(unlist(strsplit(colnames(x[[2]]), "_")), ncol=ncol(x[[2]]))[1,])      #factors of the rank weight matrix
  cancer_class_2 <- x[[2]][, grep(levels(classes_2)[1], colnames(x[[2]]))]
  normal_class_2 <- x[[2]][, grep(levels(classes_2)[2], colnames(x[[2]]))]

  for (i in 1:nrow(x[[2]]))
  {
    test_2 <- my.t.test.p.value(cancer_class_2[i,], normal_class_2[i,], alternative="greater",var.equal = FALSE,paired = FALSE)
    ttest_pvals_2 <- append(ttest_pvals_2, test_2)
  }
  ttest_pvals_2[is.na(ttest_pvals_2)] <- 1
  names(ttest_pvals_2) <- rownames(x[[2]])

  #out2 <- x[[2]][names(ttest_pvals_2)[ttest_pvals_2 <= y],]

  id_output <- as.list(c())
  #id_output <- union(names(ttest_pvals_1)[ttest_pvals_1 <= y], names(ttest_pvals_2)[ttest_pvals_2 <= y])
  id_output[[1]] <- names(ttest_pvals_1)[ttest_pvals_1 <= y]
  id_output[[2]] <- names(ttest_pvals_2)[ttest_pvals_2 <= y]
  return(id_output)
}

#' @title This generates a vector of p-values for a data matrix using the two-sample t-test
#' @name qpsp_standard_t_test
#' @param x data matrix
#' @param y significance level
#' @return id_output A p-value vectors
#' @export
qpsp_standard_t_test <- function(x, y) #x is the data matrix, and y is the significance level
{
  classes <- as.factor(matrix(unlist(strsplit(colnames(x), "_")), ncol=ncol(x))[1,])      #factors of the rank weight matrix
  normal_vals_mean <- rowMeans(x[, grep(levels(classes)[1], colnames(x))])
  cancer_vals_mean <- rowMeans(x[, grep(levels(classes)[1], colnames(x))])

  x <- t(scale(t(x)))

  ttest <- rowttests(as.matrix(x), classes,  tstatOnly = FALSE)
  ttest_pvals <- ttest$p.value
  ttest_pvals[is.na(ttest_pvals)] <- 1

  names(ttest_pvals) <- rownames(x)
  out <- x[names(ttest_pvals)[ttest_pvals <= y],]
  return(out)
}

#' @title Generates a weight matrix by applying the generate weight vector function to all columns
#' @name generate_weight_matrix
#' @param x data matrix
#' @return weight_matrix A weight matrix
#' @export
generate_weight_matrix <- function(x)
{
  weight_matrix <- apply(x, 2, generate_weight_vector)
  rownames(weight_matrix) <- rownames(x)
  return(weight_matrix)
}


#' @title Generates a weight vector
#' @name generate_weight_vector
#' @param x data matrix
#' @return output A weight matrix
#' @export
generate_weight_vector <- function(x)
{
  output <- c()
  max <- 1
  first_percentile <- (0.9 * max)
  second_percentile_1 <- (0.875 * max)
  second_percentile_2 <- (0.85 * max)
  second_percentile_3 <- (0.825 * max)
  second_percentile_4 <- (0.80 * max)

  #we are going from 1 to 0.8, 0.6,0.4,0.2, and 0s for all else
  for (i in 1:length(x))
  {
    if (x[i] >= first_percentile)
    {
      output <- append(output, "1")
    }
    else if (x[i] >= second_percentile_1)
    {
      output <- append(output, "0.8")
    }
    else if (x[i] >= second_percentile_2)
    {
      output <- append(output, "0.6")
    }
    else if (x[i]  >= second_percentile_3)
    {
      output <- append(output, "0.4")
    }
    else if (x[i]  >= second_percentile_4)
    {
      output <- append(output, "0.2")
    }
    else
    {
      output <- append(output, "0")
    }
  }

  output <- as.numeric(output)
  return(output)
}

#' @title Generates a rank matrix
#' @name generate_rank_matrix
#' @param x data matrix
#' @return output A rank matrix
#' @export
generate_rank_matrix <- function(x)
{
  rank_matrix <- apply(x, 2, rank)
  #rank_matrix <- (nrow(x) - rank_matrix)/nrow(x)
  rank_matrix <- (rank_matrix)/nrow(x)

  rownames(rank_matrix) <- rownames(x)
  colnames(rank_matrix) <- colnames(x)

  return(rank_matrix)  #dont forget to assign the column and row names later
}


#' @title Generates a modulated weight matrix
#' @name subset_weights
#' @param x data matrix
#' @return output A modulated weight matrix list
#' @export
subset_weights <- function(x) #this provides the weights for each class by calculating the average weight for that particular gene for that class and is then used for multiplying against the original matrix
{
  output <- c()
  output <- as.list(output)
  classes <- as.factor(matrix(unlist(strsplit(colnames(x), "_")), ncol=ncol(x))[1,])    	#factors of the rank weight matrix
  normal_vals_mean <- rowMeans(x[, grep(levels(classes)[2], colnames(x))]) #this provides the vote for the particular gene for class A
  cancer_vals_mean <- rowMeans(x[, grep(levels(classes)[1], colnames(x))]) #this provides the vote for the particular gene for class B

  #colnames(modulator_matrix_A) <- NULL
  #colnames(modulator_matrix_B) <- NULL
  new_mat_A <- (normal_vals_mean * x)
  new_mat_B <- (cancer_vals_mean * x)

  colnames(new_mat_A) <- colnames(x)
  colnames(new_mat_B) <- colnames(x)
  rownames(new_mat_A) <- rownames(x)
  rownames(new_mat_B) <- rownames(x)

  output[[1]] <- new_mat_A
  output[[2]] <- new_mat_B
  #output <- append(output, new_mat_A)
  #output <- append(output, new_mat_B)
  return (output)
}


#ESSNET methods
#' @title This generates all paired differences, deltas, for each gene between samples in class A and B
#' @name internal_substraction
#' @param group_A samples in group A
#' @param group_A samples in group B
#' @return all_diffs A matrix of deltas between class A and B
#' @export
internal_substraction <- function(group_A, group_B) #returns an array of subtraction vectors for all genes/proteins
{
  all_diffs <- as.data.frame(c())

  size = length(group_A)
  sel<-cbind(rep(1:size,each=size),rep(1:size,times=size))

  for (x in 1:nrow(group_A))
  {
    A <- group_A[x,]
    B <- group_B[x,]
    delta=as.numeric(A[sel[,1]]-B[sel[,2]])
    all_diffs <- rbind(all_diffs, delta)
    #append the diff to the all diffs here
  }
  row.names(all_diffs) <- row.names(group_A)
  return(all_diffs)
}

#' @title Based on the output of internal_substraction, this calculates a p-value for each complex
#' @name essnet
#' @param diff_matrix A matrix of deltas derived from internal_substraction
#' @param complete_matrix The original expression matrix
#' @param overlap_threshold The required minimal overlap between the proteomics data and the network
#' @return p_val_vec A vector of p-values, one for each complex
#' @export
essnet <- function(diff_matrix, complex_matrix, overlap_threshold)
{
  p_val_vec <- c()
  complex_names <- names(complex_matrix) #contains the names of the complexes

  for (i in 1:length(complex_matrix))
  {
    overlap_size <- sum(as.numeric(rownames(diff_matrix) %in% complex_matrix[[i]])) #sums the number of trues
    complex_components <- complex_matrix[[i]]
    if (overlap_size <= overlap_threshold)
    {
      p_val_vec <- append(p_val_vec, 0)
    }
    else
    {
      list_vector <- c()
      subset <- diff_matrix[rownames(diff_matrix) %in% complex_components,]
      for (j in 1:nrow(subset))
      {
        list_vector <- append(list_vector, as.numeric(subset[j,]))
      }
      #p_val_vec <- append(p_val_vec, t.test(list_vector)$p.value)

      #size <- length(list_vector)
      #size <- sqrt(length(list_vector))
      #size <- overlap_size * sqrt(length(list_vector))
      size <- sqrt(overlap_size) * sqrt(length(list_vector)) #here we reduce the effect of k by taking square root of k
      sterr<-(sd(list_vector)/sqrt(size-1))
      stat<-mean(list_vector)/sterr
      t.pval <- sum(abs(rt(1000,size-1))>abs(stat))/1000
      p_val_vec <- append(p_val_vec, t.pval)
    }
  }
  names(p_val_vec) <- complex_names
  return(p_val_vec)
}


#SNET methods

#' @title Generates a vector of snet weights for a given sample
#' @name snet_weight_vector
#' @param x An expression vector
#' @return output A vector of weights
#' @export
snet_weight_vector <- function(x)
{
  output <- c()
  max <- 1
  first_percentile <- (0.9 * max)

  #we are going from 1 to 0.8, 0.6,0.4,0.2, and 0s for all else
  for (i in 1:length(x))
  {
    if (x[i] >= first_percentile)
    {
      output <- append(output, "1")
    }
    else
    {
      output <- append(output, "0")
    }
  }
  output <- as.numeric(output)
  return(output)
}

#' @title Generates a matrix of snet weights for a given data matrix
#' @name snet_generate_weight_matrix
#' @param x An expression matrix
#' @return output A matrix of weights
#' @export
snet_generate_weight_matrix <- function(x)
{
  weight_matrix <- apply(x, 2, snet_weight_vector)
  rownames(weight_matrix) <- rownames(x)
  return(weight_matrix)
}

#FSNET methods
#' @title Based on the rank weight matrix and the list of complexes, give a matrix of scores per complex
#' @name fsnet
#' @param rank_weight_matrix A matrix of weighted ranks
#' @param complex_matrix The original expression matrix
#' @param overlap_threshold The required minimal overlap between the proteomics data and the network
#' @return output A list of two matrices of fsnet scores, x vs y, and y vs x, for complexes for all samples
#' @export
fsnet <- function(rank_weight_matrix, complex_matrix, overlap_threshold) #this one returns a matrix of the fsnet scores of samples vs subnets which can then be evaluated by a fsnet_pval calculation approach
{
  fsnet_matrix_A <- c() #this will vectors of fsnet scores that can be assimilated into a matrix
  fsnet_matrix_B <- c() #this will vectors of fsnet scores that can be assimilated into a matrix

  for (i in 1:length(complex_matrix))
  {
    overlap_size <- sum(as.numeric(rownames(rank_weight_matrix[[1]]) %in% complex_matrix[[i]])) #sums the number of trues
    complex_components <- complex_matrix[[i]]
    if (overlap_size <= overlap_threshold)
    {
      fsnet_matrix_A <- rbind(fsnet_matrix_A,rep(0, ncol(rank_weight_matrix[[1]])))
      fsnet_matrix_B <- rbind(fsnet_matrix_B,rep(0, ncol(rank_weight_matrix[[2]])))
    }
    else
    {
      subset_A <- rank_weight_matrix[[1]][rownames(rank_weight_matrix[[1]]) %in% complex_components,]
      subset_B <- rank_weight_matrix[[2]][rownames(rank_weight_matrix[[2]]) %in% complex_components,]
      #here if we want to calculate the class modulator we can do so by splitting the classes by names and calculating the average normalized rank
      #subset <- subset_weights(subset)
      complex_sums_A <- colSums(subset_A) #otherwise, we simply add the rank_weight values up.
      complex_sums_B <- colSums(subset_B) #otherwise, we simply add the rank_weight values up.
      fsnet_matrix_A <- rbind(fsnet_matrix_A, complex_sums_A)
      fsnet_matrix_B <- rbind(fsnet_matrix_B, complex_sums_B)
    }
  }
  rownames(fsnet_matrix_A) <- names(complex_matrix)
  rownames(fsnet_matrix_B) <- names(complex_matrix)
  colnames(fsnet_matrix_A) <- colnames(rank_weight_matrix[[1]])
  colnames(fsnet_matrix_B) <- colnames(rank_weight_matrix[[2]])

  output <- as.list(c())
  output[[1]] <- fsnet_matrix_A
  output[[2]] <- fsnet_matrix_B
  return(output)
}

#PFSNET methods
#' @title Uses the original FSNET matrix but applies PFSNET's calculation method
#' @name pfsnet_theoretical_t_test
#' @param x FSNET matrix
#' @param y test significance level
#' @return out p-values based on pfsnet
#' @export
pfsnet_theoretical_t_test <- function(x, y) #where x is the fsnet matrix and y is the significance level of test
{
  #x <- t(scale(t(x)))
  #lets play with the first matrix
  #classes <- as.factor(matrix(unlist(strsplit(colnames(x[[1]]), "_")), ncol=ncol(x[[1]]))[1,])      #factors of the rank weight matrix
  #normal_class <- x[[1]][, grep(levels(classes)[1], colnames(x[[1]]))]
  #cancer_class <- x[[1]][, grep(levels(classes)[2], colnames(x[[1]]))]

  class_A <- x[[1]] #matrix 1
  class_B <- x[[2]] #matrix 2
  pvals <- c()

  for (i in 1:nrow(class_A)) #for each complex
  {
    all_pairs <- expand.grid(class_A[i,],class_B[i,]) #generate all pairs
    differences <- all_pairs[,1] - all_pairs[,2]

    size <- sqrt(length(differences))
    sterr <- sd(differences)/(size - 1)
    stat <- mean(differences)/sterr
    t.pval <- sum(abs(rt(1000,size-1))>abs(stat))/1000
    pvals <- append(pvals, t.pval)
  }
  pvals[is.na(pvals)] <- 1
  names(pvals) <- rownames(x[[1]])
  out <- pvals[pvals <= y]
  return(out)
}

#PPFSNET methods
#' @title Uses the original FSNET matrix but applies PPFSNET's calculation method
#' @name ppfsnet_theoretical_t_test
#' @param x FSNET matrix
#' @param y test significance level
#' @return out p-values based on pfsnet
#' @export
ppfsnet_theoretical_t_test <- function(x, y) #where x is the fsnet matrix and y is the significance level of test
{
  #first we generate two matrices of all paired differences

  #lets play with the first matrix
  classes <- as.factor(matrix(unlist(strsplit(colnames(x[[1]]), "_")), ncol=ncol(x[[1]]))[1,])      #factors of the rank weight matrix

  pvals <- c()
  for (i in 1:nrow(x[[1]]))
  {
    #cancer_class_A <- x[[1]][i, grep(levels(classes)[1], colnames(x[[1]]))]
    normal_class_A <- x[[1]][i, grep(levels(classes)[2], colnames(x[[1]]))]

    cancer_class_B <- x[[2]][i, grep(levels(classes)[1], colnames(x[[2]]))]
    #normal_class_B <- x[[2]][i, grep(levels(classes)[2], colnames(x[[2]]))]

    #all_pairs_A <- expand.grid(cancer_class_A, cancer_class_B)
    #all_pairs_B <- expand.grid(normal_class_A, normal_class_B)

    all_pairs_B_vs_A <- expand.grid(cancer_class_B, normal_class_A)
    all_pairs_B_vs_A <- abs(all_pairs_B_vs_A[,1] - all_pairs_B_vs_A[,2])
    #differences_A <-all_pairs_A[,1] - all_pairs_A[,2]
    #differences_B <-all_pairs_B[,1] - all_pairs_B[,2]

    #all_pairs <- expand.grid(differences_A, differences_B)

    #differences <- all_pairs[,1] - all_pairs[,2]
    differences <- all_pairs_B_vs_A

    size <- sqrt(length(differences)/2)
    sterr <- sd(differences)/(size - 1)
    stat <- mean(differences)/sterr
    t.pval <- sum(abs(rt(1000,size-1))>=abs(stat))/1000
    pvals <- append(pvals, t.pval)
  }
  pvals[is.na(pvals)] <- 1
  names(pvals) <- rownames(x[[1]])
  out <- pvals[pvals <= y]
  return(out)
}



#HE
#' @title The hypergeometric enrichment pipeline
#' @name he
#' @param x protein expression matrix
#' @param y complex vector
#' @param z class factors
#' @return hyp_pval A vector of p-values, one for each complex
#' @export
he <- function(x, y, z) # where x is the data matrix, y is the complex vector, and z is the factor
{
  hyp_pval <- c()

  for (n in 1:length(y))
  {
    test_ttest <- rowttests(as.matrix(x), z)
    pvals <- test_ttest$p.value
    pvals <- p.adjust(pvals, method="bonferroni")
    sig_prot_ttest <- rownames(x)[which(pvals <= 0.05)]

    q <- length(intersect(y[[n]], sig_prot_ttest)) #x is the intersection
    k <- length(y[[n]]) #size of the complex
    N <- nrow(x) #all the proteins
    m <- length(sig_prot_ttest)

    pval <- 1- phyper(q-1, m, N-m, k)
    hyp_pval <- append(hyp_pval, pval)
  }
  return(hyp_pval)
}


#GSEA
#' @title The original GSEA algorithm based on the Kolmogorov-Smirnov
#' @name gsea
#' @param data expression matrix
#' @param complex_vector A list of complexes
#' @param factors levels indicating the class of samples in the expression matrix
#' @return gsea_pval A vector of p-values; one for each complex
#' @export
gsea <- function(data, complex_vector, factors)
{
  gsea_pval <- c()
  pvals <- rowttests(as.matrix(data), factors)$p.value
  names(pvals) <- rownames(data)
  ranks <- rank(pvals)
  #first lets calculate all the t-scores then convert to ranks

  for (x in 1:length(complex_vector))
  {
    if (length(ranks[which(names(ranks) %in% complex_vector[[x]])]) >= 1)
    {
      ks_pval <- ks.test(ranks[which(names(ranks) %in% complex_vector[[x]])], ranks[which(!names(ranks) %in% complex_vector[[x]])])$p.value
      gsea_pval <- append(gsea_pval, ks_pval)
    }
    else
    {
      gsea_pval <- append(gsea_pval, 1)
    }
  }
  #gsea_pval <- p.adjust(gsea_pval, method="fdr") #turn this off if do not want to do multiple test correction
  return(gsea_pval)
}

#QPSP methods
#' @title QPSP rank matrix generation function
#' @name qpsp_generate_rank_matrix
#' @param x expression matrix
#' @return rank_matrix A matrix of weighted ranks
#' @export
qpsp_generate_rank_matrix <- function(x)
{
  rank_matrix <- apply(x, 2, rank)
  #rank_matrix <- rank_matrix/nrow(x)
  rank_matrix <- (rank_matrix)/nrow(x)

  rownames(rank_matrix) <- rownames(x)
  colnames(rank_matrix) <- colnames(x)

  return(rank_matrix)  #dont forget to assign the column and row names later
}

#' @title QPSP rank matrix generation function
#' @name qpsp
#' @param rank_weight_matrix A matrix of weights based on ranks
#' @param complex_list A list of complexes
#' @return output A vector of p-values
#' @export
qpsp <- function(rank_weight_matrix, complex_list)
{
  qpsp_matrix <- c()

  for (j in 1:length(complex_list))
  {
    #print(j)
    if (length(rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% complex_list[[j]])]) > 1)
    {
      qpsp_matrix <- rbind(qpsp_matrix, colSums(rank_weight_matrix[rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% complex_list[[j]])],]))
    }
    else if (length(rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% complex_list[[j]])]) == 1)
    {
      qpsp_matrix <- rbind(qpsp_matrix, rank_weight_matrix[rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% complex_list[[j]])],])
    }
    else
    {
      qpsp_matrix <- rbind(qpsp_matrix, c(rep(0, ncol(rank_weight_matrix))))
    }

  }

  rownames(qpsp_matrix) <- names(complex_list)

  output <- qpsp_standard_t_test(qpsp_matrix, 0.05)
  return(output)
}

#Simulate datasets
#' @title Proteomics data simulation function
#' @name generate_proteomics_sim
#' @param data A one class data matrix
#' @param fact Factors
#' @param prop Proportion of variables to insert effect size
#' @param rounds Number ofsimulations
#' @param effect_size A vector of intended effect sizes
#' @return prop A set of simulated datasets
#' @export
generate_proteomics_sim <- function(data, fact, prop, rounds, effect_size) #data is one class data, factor assigns class labels to data, prop is the proportion of variables to assign effect size, rounds is number of simulated data to geenrate, effect_size is a vector of effect sizes to sample from
{
  for (i in 1:rounds)
  {
    print(i)
    norm <- which(fact == levels(fact)[1])
    cancer <- which(fact == levels(fact)[2])

    null_mat <- matrix(rep(1,nrow(data)*ncol(data)), nrow= nrow(data), ncol=ncol(data)) #the multiplication factor
    null_null_mat <- matrix(rep(0,nrow(data)*ncol(data)), nrow= nrow(data), ncol=ncol(data)) #determines the adjustment factor

    diff <- sample(1:nrow(data), round(prop * nrow(data))) #choose 20%
    null_null_mat[diff, cancer] <- replicate(length(diff), sample(effect_size, 1))

    null_mat <- null_mat + null_null_mat

    data_adjust <- data * null_mat

    factor_vec <- rep("ns", nrow(data))
    factor_vec[diff] <- 'sig'

    data_adjust <- cbind(rep(rounds, nrow(data)), data_adjust, factor_vec)
    new_colnames <- c("GeneLength", as.vector(fact), "sig")
    colnames(data_adjust) <- new_colnames

    write.table(data_adjust, file=paste("sim.",i, ".txt",sep=""), sep="\t", row.names=T, col.names=T, quote=F)
  }
}

#Simulate complex
#' @title Pseudo-complex simulation function
#' @name generate_pseudo_cplx
#' @param tdata A one class data matrix
#' @param fact Factors
#' @param purity Proportion of differential proteins in complex
#' @return prop A set of simulated pseudo-complexes
#' @export
generate_pseudo_cplx <- function(tdata, fact, purity) #pseudo complex generator
{
  complex_vector <- c()
  data <- tdata
  data <- data[,-1]
  data <- data[,-9]
  #colnames(data) <- c(paste('N_', 1:4), paste('C_', 1:4))

  sig_prot <- rownames(tdata[which(tdata[,10]=='sig'),])
  ns_prot <- sample(rownames(tdata[which(tdata[,10]=='ns'),]), length(rownames(tdata[which(tdata[,10]=='sig'),])))

  cor_sig_clust <- hclust(dist(tdata[sig_prot,2:9]), method="ward.D")
  cor_ns_clust <- hclust(dist(tdata[ns_prot,2:9]), method="ward.D")

  sig_prot <- sig_prot[cor_sig_clust$order]
  ns_prot <- ns_prot[cor_ns_clust$order]

  sig_prot_complex <- split(sig_prot, sort(c(rep(rep(1:50), 10))))
  names(sig_prot_complex) <- paste("sig_",names(sig_prot_complex), sep='')

  ns_prot_complex <- split(ns_prot, sort(c(rep(rep(1:50), 10))))
  names(ns_prot_complex) <- paste("ns_",names(ns_prot_complex), sep='')

  for (j in 1:length(sig_prot_complex)) ######purity tests ---- turn this off if not wanted~! This substitutes the n% with proteins from NS
  {
    sig_prot_complex[[j]][1:round(((1 - purity)*length(sig_prot_complex[[j]])))] <- ns_prot_complex[[j]][1:round(((1 - purity)*length(sig_prot_complex[[j]])))]
  }

  complex_vector <- append(sig_prot_complex, ns_prot_complex)

  return(complex_vector)
}
