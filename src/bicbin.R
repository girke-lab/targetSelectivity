# Purpose: find top binary clusters in code
# Code from Uitert, M. V., Meuleman, W. & Wessels, L. Biclustering Sparse Binary Genomic Data. Journal of Computational Biology 15, 1329?1345 (2008).

# modifications from fastbreak (https://code.google.com/p/fastbreak/):
#   this file has been modified to make it usable. BicBin now returns a list with the score, x and y
#   values instead of just the score

BicBin <- function(G, alpha, beta, p, proc_genes=TRUE) {
    M <- nrow(G) # TFs
    N <- ncol(G) # Genes
    
    max_scores <- c();
    
    # Selected TFs
    x <- rep.int(0, M);
    x[rnorm(M) < 0] <- 1
    xmatrix = c()
    
    # Selected Genes
    y <- rep.int(0, N); 
    y[rnorm(N) < 0] <- 1
    ymatrix = c()
    
    while (length(max_scores) < 2 || diff(max_scores)[length(max_scores)-1] != 0) {
        if (proc_genes) {
            n <- sum(y)
            
            s <- G %*% y # Rowsums over selected Genes (i.e., y == 1)
            ss <- sort(s, decr=TRUE, index.return=TRUE);
            
            C <- score(k <- cumsum(ss$x), m <- 1:M, n, p, alpha, beta)
            
            hits <- 1:which.max(C)
            x <- rep.int(0, M); x[ss$ix[hits]] <- 1 
            
            proc_genes <- FALSE;
        } else {
            m <- sum(x)
            
            s <- as.matrix(t(x) %*% G) # Colsums over selected TFs (i.e., x == 1)
            ss <- sort(s, decr=TRUE, index.return=TRUE);
            
            C <- score(k <- cumsum(ss$x), m, n <- 1:N, p, alpha, beta)
            
            hits <- 1:which.max(C)
            y <- rep.int(0, N); y[ss$ix[hits]] <- 1 
            
            proc_genes <- TRUE;
        }
        
        max_scores <- c(max_scores, max(C))
        xmatrix <-rbind(xmatrix,x)
        ymatrix <-rbind(ymatrix,y)
    }
    
    
    maxscore<-max(max_scores)
    return = list(
        score = maxscore,
        x = xmatrix[max_scores == maxscore,][1,],
        y = ymatrix[max_scores == maxscore,][1,]
    )
}

score <- function(k, m, n, p, alpha, beta) {
    mnp <- m*n*p
    res1 <- (k - mnp)^2 / (m^alpha * n^beta * (k+mnp))
    res2 <- (k - mnp)^2 / (3 * m^(1+alpha) * n^(1+beta) * p)
    
    hits <- k-(2*mnp) > 0
    res <- c(res1[hits], res2[!hits])
    
    # Same as: if (k - mnp > mnp) { res1 } else { res2 }
}

###################### Example #########################

#library(SparseM)
#mat <- matrix(sample(c(rep.int(0, 100), 1), 800*25000, replace=TRUE), 800, 25000)
#mat <- as.matrix.csr(mat)
#mat[10:15, 10:15] <- 1
#M <- nrow(mat) # TFs
#N <- ncol(mat) # Genes
#alpha <- 0.5
#beta <- 0.5
#p <- sum(attr(mat, "ra")) / (M*N)
#
#max <- 0;
#for (i in 1:300) {
#  print(i)
#  res <- BicBin(mat, alpha, beta, p, proc_genes=TRUE)
#  if (res > max) { max <- res }
#  res <- BicBin(mat, alpha, beta, p, proc_genes=FALSE)
#  if (res > max) { max <- res }
#}