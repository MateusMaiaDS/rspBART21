n_ <- 10
ones <- matrix(1,nrow = n_)
diag_ <- diag(nrow = n_)

determinant(diag_+tcrossprod(ones),logarithm = FALSE)$modulus^(-1/2)
(n_+1)^(-1/2)
