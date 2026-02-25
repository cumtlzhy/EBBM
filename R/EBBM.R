#' @title A novel biclustering algorithm for mining m6A co-methylation patterns based on beta-binomial distribution and data screening strategy
#'
#' @description  This package enables the discovery of co-methylation patterns for the epigenetic transcription modification m6A methylation. By inputting raw m6A sequencing reads from IP samples and input samples, it automatically uncovers co-methylation patterns hidden within the data.
#'
#'
#' @param a Ip matrix
#' @param b input matrix
#' @param iteration Number of initial iterations
#' @param burn_in Number of burnings
#' @param bic_row_labels Labels for sites in the patterns found, all with initial values of 0
#' @param number_bicluster Number of preset patterns
#'
#' @return the patterns hidden in the data
#' @export
#'
#' @examples
#' rm(list = ls(all = TRUE))
#' library(EBBM)
#' data("ip")
#' data("input")
#' number_bicluster <- 6
#' bic_row_labels <- c(rep(0,dim(ip)[1]))
#' patterns_found_by_EBBM <- multiple_GSB(a = ip,b = input,iteration = 1000,burn_in = 500,bic_row_labels = bic_row_labels,number_bicluster = number_bicluster)
#'


###################----------mulitiple bicluster model-----------######################

multiple_GSB <- function(a,b,iteration,burn_in,bic_row_labels,number_bicluster){
  out = list()
  for (c in 1:number_bicluster) {
    GSB_single_object <- GSB_SINGLE(a = a,b =b,iteration = iteration,burn_in = burn_in, bic_row_labels = bic_row_labels)

    assign(paste( c,"bicluster", sep="_"),GSB_single_object)
    out[[paste( c,"bicluster", sep="_")]] <- get(paste( c,"bicluster", sep="_"))
    save(out,file =  "./out")
    if( sum(GSB_single_object$row_labels)==0|sum(GSB_single_object$col_labels)==0 ){
      break;
    }else{
      bic_row_labels <- bic_row_labels + GSB_single_object$row_labels # Update bicluster labels
    }



  }
  return(out)

}
####################--------Single-bicluster discovery model-----------######################

GSB_SINGLE <- function(a,b,iteration,burn_in,bic_row_labels){
  a_b <- a + b
  n <- dim(a)[1]
  m<- dim(a)[2]

  matrix_pro_row <- matrix(0,n,1)# construct a matrix to save the value of problity of each row of iteration
  matrix_pro_col <- matrix(0,m,1)
  matrix_labels_row <- matrix(0,n,1)# construct a matrix to save the value of labels of each row of iteration
  matrix_labels_col <- matrix(0,m,1)

  row_theta <- c(2,3)
  col_theta <- c(4,10)
  pro_row <- rbeta(n,row_theta[1],row_theta[2])
  pro_col <- rbeta(m,col_theta[1],col_theta[2])
  row_labels <- c(rbinom(n = dim(a)[1], size = 1, prob = pro_row))
  col_labels <- c(rbinom(n = dim(a)[2], size = 1, prob = pro_col))

  matrix_pro_row[,1] <- pro_row
  matrix_pro_col[,1] <- pro_col
  matrix_labels_row[,1] <- row_labels
  matrix_labels_col[,1] <- col_labels

  bic_pro_theta <- c(9,4)
  bag_pro_theta <- c(3,10)
  bic_pro <- rbeta(1,bic_pro_theta[1],bic_pro_theta[2])
  bic_pro_chain <- bic_pro

  bag_pro <- rbeta(1,bag_pro_theta[1],bag_pro_theta[2])
  bag_pro_chain <- bag_pro
  k <- 1
  update_pro_rowi <-0
  update_pro_colj<-0
  likelihood <-0
  weighted_matrix_bic <- log((a+0.01/(b+0.01)))
  weighted_matrix_bic <- weighted_matrix_bic - min(weighted_matrix_bic+1)
  weighted_matrix_bag <- log((b+0.01/(a+0.01)))
  weighted_matrix_bag <- weighted_matrix_bag - min(weighted_matrix_bag+1)
  t <- 1
  repeat{

    for (i in 1:n) {
      # label of i equal 1
      bic_row1 <- i # row index of bicluster
      bic_col1 <- c(which(col_labels==1))    # column index of bicluster
            denstiy_bic <- dbinom(a[i,bic_col1],a_b[i,bic_col1],bic_pro)
      #####------Weighted correction, modifying only elements belonging to the bicluster------#######

      denstiy_bic[which(denstiy_bic==0)] <- .Machine$double.xmin
      denstiy_bag <- dbinom(a[i,bic_col1],a_b[i,bic_col1],bag_pro)
      denstiy_bag[which(denstiy_bag==0)] <- .Machine$double.xmin

      f.1 <- sum(log(denstiy_bic/denstiy_bag))
      v_ibar <- length(which(row_labels[-i]==1)) # the number of 1 in the row lables exclude i
      f.2 <-log((v_ibar + row_theta[1])/(n-v_ibar+row_theta[2]-1))
      odd_i <- f.1 + f.2
      update_pro_rowi[i] <- 1/(1+exp(-odd_i))

      if(bic_row_labels[i]== 0){
        row_labels[i] <- rbinom(1,1,update_pro_rowi[i])

      }else{
        row_labels[i] <- 0
      }


    }



    matrix_labels_row <- cbind(matrix_labels_row,row_labels)
    matrix_pro_row <- cbind(matrix_pro_row,update_pro_rowi)  # save the problity of label of i-th row

    # #  draw bic_pro 和bag_pro
    decision <- ((a[c(which(row_labels==1)),c(which(col_labels==1))]+0.01)/(b[c(which(row_labels==1)),c(which(col_labels==1))]+0.01))>1.5
    aa <- a[c(which(row_labels==1)),c(which(col_labels==1))]

    a_sum_update_bic_row <- sum(aa[which(decision==TRUE)])

    if(a_sum_update_bic_row==0){
      break
    }


    bb <- b[c(which(row_labels==1)),c(which(col_labels==1))]
    b_sum_update_bic_row <- sum(bb[which(decision==TRUE)])


    bic_pro <- rbeta(1,(bic_pro_theta[1] +  a_sum_update_bic_row),(bic_pro_theta[2]+b_sum_update_bic_row))

    if(bic_pro < 0.6){
      break
    }

    a_sum_update_bag_row <- sum(a) - a_sum_update_bic_row
    b_sum_update_bag_row <- sum(b) - b_sum_update_bic_row
    bag_pro <- rbeta(1,(bag_pro_theta[1]+a_sum_update_bag_row),(bag_pro_theta[2]+b_sum_update_bag_row))

    for (j in 1:m) {
      bic_col2 <- j
      bic_row2 <- c(which(row_labels==1))

      density_j_bic <- dbinom(a[bic_row2,j],a_b[bic_row2,j],bic_pro)
      ##############------Weighted Adjustment------##########

      density_j_bic[which(density_j_bic==0)] <- .Machine$double.xmin

      density_j_bag <- dbinom(a[bic_row2,j],a_b[bic_row2,j],bag_pro)

      density_j_bag[which(density_j_bag==0)] <- .Machine$double.xmin


      g.1 <- sum(log(density_j_bic/density_j_bag))

      w_jbar <- length(which(col_labels[-j]==1))
      g.2 <- log((w_jbar+col_theta[1])/(m-w_jbar+col_theta[2]-1))
      odd_j <- g.1 + g.2
      update_pro_colj[j] <- 1/(1+exp(-odd_j))

      j_label <- rbinom(1,1,update_pro_colj[j])
      col_labels[j] <- j_label

    }
    matrix_labels_col <- cbind(matrix_labels_col,col_labels)
    matrix_pro_col <- cbind(matrix_pro_col,update_pro_colj)
    # draw bic_pro 和bag_pro

    a_sum_update_bic_col <- sum(a[c(which(row_labels==1)),c(which(col_labels==1))])
    b_sum_update_bic_col <- sum(b[c(which(row_labels==1)),c(which(col_labels==1))])
    bic_pro <- rbeta(1,(bic_pro_theta[1]+a_sum_update_bic_col),(bic_pro_theta[2]+b_sum_update_bic_col))

    if(bic_pro < 0.6){
      break
    }

    bic_pro_chain[k+1] <- bic_pro

    a_sum_update_bag_col <- sum(a) - a_sum_update_bic_col
    b_sum_update_bag_col <- sum(b) - b_sum_update_bic_col
    bag_pro <- rbeta(1,(bag_pro_theta[1]+a_sum_update_bag_col),(bag_pro_theta[2]+b_sum_update_bag_col))

    bag_pro_chain[k+1] <- bag_pro
    # Compute the log-likelihood for each iteration
    l_bi <- mapply(function(x,y) log(dbinom(x,y,bic_pro)) , a[c(which(row_labels==1)),c(which(col_labels==1))],
                   a_b[c(which(row_labels==1)),c(which(col_labels==1))])
    l_bi[which(l_bi == -Inf)] <- -.Machine$double.xmin

    if(length(l_bi)==0){
      l_bi <- 0
    }
    likelihood.1 <- sum(l_bi)

    l_ba <- mapply(function(x,y) log(dbinom(x,y,bag_pro)) , a[c(which(row_labels==0)),c(which(col_labels==0))],
                   a_b[c(which(row_labels==0)),c(which(col_labels==0))])

    l_ba[which(l_ba == -Inf)] <- -.Machine$double.xmin
    if(length(l_ba)==0){
      l_ba <- 0
    }
    likelihood.2 <- sum(l_ba)
    likelihood[k] <- sum(likelihood.1,likelihood.2)
    #Convergence Test
    if((k == iteration & MCE(samp = bic_pro_chain[(iteration-500+1):iteration],v = 50)< 0.01) | k == (iteration + 500*t)){
      break


    }

    k <- k+1


  }

  if(k > iteration){
    iteration <- k
  }else{
    iteration <- iteration
  }

  # Calculate the mean


  if(a_sum_update_bic_row==0 | k < iteration){
    update_row_labels <-  rep(0,n)
    update_col_labels <-  rep(0,m)

  }
  else{
    update_row_labels <-  apply(matrix_labels_row[,(burn_in+2):iteration], 1, sum)/(iteration-burn_in)
    update_col_labels <-  apply(matrix_labels_col[,(burn_in+2):iteration], 1, sum)/(iteration-burn_in)


  }

  # Calculate the expected value

  # Quantile Threshold Determination
  update_row_labels[which(update_row_labels >= 0.75)] <- 1
  update_row_labels[which(update_row_labels < 0.75)] <- 0
  update_col_labels[which(update_col_labels >= 0.75)] <- 1
  update_col_labels[which(update_col_labels < 0.75)] <- 0
  out <- list()
  out[['row_labels']] <- update_row_labels
  out[['col_labels']] <- update_col_labels
  out[['matrix_row_labels']] <- matrix_labels_row
  out[['matrix_col_labels']] <- matrix_labels_col
  out[['matrix_pro_row']] <- matrix_pro_row
  out[['matrix_pro_col']] <- matrix_pro_col
  out[["bic_pro"]] <- bic_pro_chain
  out[["bag_pro"]] <- bag_pro_chain
  out[["log_likelihood"]] <- likelihood
  out[["iteration"]] <- iteration

  return(out)

}

####MCE Monte Carlo error#####
MCE <- function(samp,v){
  b <- 0
  mean_b <- 0
  k <- length(samp)/v
  for (i in 1:k) {

    mean_b[i] <- mean(samp[((i-1)*v+1):(i*v)])

  }
  mce <- sqrt(var(mean_b)/(k-1))
  return(mce)

}




