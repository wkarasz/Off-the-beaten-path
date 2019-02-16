# Parallel Matrix Multiplication
# William Karasz
# 2018-08-08
# https://www.cs.utexas.edu/users/rvdg/papers/SSUMMA.ps

library("rdriver")
library("rlogger")

strassen <- function(A,B){
  # Pad with zeros if necessary
  n0 <- log2(dim(A)[1])
  if (n0%%1==0) {
    # We have a valid matrix in the form 2^n
    additional_rows_cols <- 0
    total_rows_cols <- 2^n0;
    
  } else {
    # Determine number of rows and columns for padding
    total_rows_cols <- 2^ceil(n0);
    additional_rows_cols <- total_rows_cols-dim(A)[1]
  }
  
  A <- cbind(A,matrix(integer(1),nrow=dim(A)[1],ncol=additional_rows_cols))
  A <- rbind(A,matrix(integer(1),nrow=additional_rows_cols,ncol=total_rows_cols))
  
  B <- cbind(B,matrix(integer(1),nrow=dim(B)[1],ncol=additional_rows_cols))
  B <- rbind(B,matrix(integer(1),nrow=additional_rows_cols,ncol=total_rows_cols))
  
  
  n1 <- log2(dim(A)[1]);
  
  
  # Conformally divide matrix
  A11 <- A[1:2^(n1-1),1:2^(n1-1)]
  A12 <- A[1:2^(n1-1),(2^(n1-1)+1):dim(A)[1]]
  A21 <- A[(2^(n1-1)+1):dim(A)[1],1:2^(n1-1)]
  A22 <- A[(2^(n1-1)+1):dim(A)[1],(2^(n1-1)+1):dim(A)[1]]
  
  B11 <- B[1:2^(n1-1),1:2^(n1-1)]
  B12 <- B[1:2^(n1-1),(2^(n1-1)+1):dim(B)[1]]
  B21 <- B[(2^(n1-1)+1):dim(B)[1],1:2^(n1-1)]
  B22 <- B[(2^(n1-1)+1):dim(B)[1],(2^(n1-1)+1):dim(B)[1]]
  
  if (n1 == 1) {
    cs <- ServiceFactory.createService("Rmma");
    rs <- ResultSet(cs)
    cs$submit(rs,"mma", list(A11+A22,B11+B22))
    cs$submit(rs,"mma", list(A21+A22,B11))
    cs$submit(rs,"mma", list(A11,B12-B22))
    cs$submit(rs,"mma", list(A22,B21-B11))
    cs$submit(rs,"mma", list(A11+A12,B22))
    cs$submit(rs,"mma", list(A21-A11,B11+B12))
    cs$submit(rs,"mma", list(A12-A22,B21+B22))
    LL <- rs$waitForAll(as.integer(100))
    print(paste(c("LL: ", LL[['value']]), collapse=""))
    LS <- LL[which(sapply(LL, `[[`, "status") == 0)]	# completed only
    #avg <- mean(sapply(LS, `[[`, "value"))              # average
    cs$waitUntilInactive()
    cs$destroy()
    delete(rs)
    delete(cs)
    
    M1 <- LS[[1]]$value
    M2 <- LS[[2]]$value
    M3 <- LS[[3]]$value
    M4 <- LS[[4]]$value
    M5 <- LS[[5]]$value
    M6 <- LS[[6]]$value
    M7 <- LS[[7]]$value
  } else {
    M1 <- strassen(A11+A22,B11+B22);
    M2 <- strassen(A21+A22,B11);
    M3 <- strassen(A11,B12-B22);
    M4 <- strassen(A22,B21-B11);
    M5 <- strassen(A11+A12,B22);
    M6 <- strassen(A21-A11,B11+B12);
    M7 <- strassen(A12-A22,B21+B22);
  }
  
  C <- cbind(M1+M4-M5+M7, M3+M5)
  C <- rbind(C,cbind(M2+M4, M1-M2+M3+M6))
  C <- C[1:dim(A)[1],1:dim(A)[1]]
  return(C)
  
}





# Submit a random square matrix for evaluation

n<-4
NCols<-2^n
NRows<-NCols 

myMat<-matrix(runif(NCols*NRows), ncol=NCols) 

A<-myMat
B<-myMat

C<-strassen(A,B)

C
