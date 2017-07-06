# tajima d

#' @title Tajima's pi estimator of theta
#' @export
#' @keywords internal
theta_w <- function(sfs_1d) {
    n <- length(sfs_1d) - 1
    i <- 1:(n-1)
    a <- sum(1/i)
    s <- sum(sfs_1d[i])
    return(s/a)
}

#' @title Watterson's estimator of theta
#' @export
#' @keywords internal
theta_pi <- function(sfs_1d) {
    n <- length(sfs_1d) - 1
    i <- 1:(n-1)
    inv <- n - i
    a <- 1/choose(n, 2)
    e <- sum(a * i * inv * sfs_1d[i])
    return(e)
}

#' @title Fay and Wu's estimator of theta
#' @export
#' @keywords internal
theta_h <- function(sfs_1d) {
    n <- length(sfs_1d) - 1
    i <- 1:(n-1)
    a <- 1/choose(n, 2)
    e <- sum(a * i**2 * sfs_1d[i])
    return(e)
}

#' @title Convert multidimensional joint SSF to multiple unidimensional SFSs
#' @export
convert_to_1dSFS <- function(x, relativize = TRUE) {
    if (is.vector(x)) {
        nama <- names(x)
        dim(x) <- c(1, length(x))
        colnames(x) <- nama
    }

    X <- sapply(colnames(x), function(z) (strsplit(gsub("ac.", "", z), "\\.")))
    X <- apply(t(simplify2array(X)), 2, as.integer)

    Z <- apply(x, 1, function(x) {
        sapply(1:ncol(X), function(j) {
            y = as.vector(by(x, X[,j], sum))
			y[1] = 0
			if (relativize) {
				sm = sum(y)
				if (sm!=0) y = y/sm
			}
			return(y)
        })
    })
    Z <- t(Z)

    colnames(Z) <- c(sapply(1:ncol(X), function(j) paste0(j,"_",unique(X[,j]))))
	
    return(Z)
}

#' @title Tajima's D test
#' @export
#' @keywords internal
tajima <- function(sfs_nD_matrix, internal = TRUE, do_H = FALSE) {

    if (is.vector(sfs_nD_matrix)) {
        nama <- names(sfs_nD_matrix)
        dim(sfs_nD_matrix) <- c(1, length(sfs_nD_matrix))
        colnames(sfs_nD_matrix) <- nama
    }

    X <- sapply(colnames(sfs_nD_matrix), function(x) (strsplit(gsub("ac.", "", x), "\\.")))
    X <- apply(t(simplify2array(X)), 2, as.integer)

    if (internal) {
        Z <- apply(sfs_nD_matrix, 1, function(x) {
            sapply(1:ncol(X), function(j) {
                y <- as.vector(by(x, X[,j], sum))
                d <- theta_pi(y) - theta_w(y)
				h <- theta_pi(y) - theta_h(y)
                n <- length(y) - 1
                S <- sum(y)
                i <- 1:(n-1)
                a1 <- sum(1/i)
                a2 <- sum(1/(i**2))
                b1 <- (n+1)/(3*(n-1))
                b2 <- (2*(n**2+n+3))/(9*n*(n-1))
                c1 <- b1 - 1/a1
                c2 <- b2 - (n+2)/(a1*n) + a2/(a1**2)
                e1 <- c1/a1
                e2 <- c2/(a1**2 + a2)
                Var <- e1*S + e2*S*(S-1)
                c(d/sqrt(Var), h)
            })
        })
    } else {
        Z <- apply(sfs_nD_matrix, 1, function(x) {
            sapply(1:ncol(X), function(j) {
                y <- as.vector(by(x, X[,j], sum))
                d <- theta_pi(y) - theta_w(y)
				h <- theta_pi(y) - theta_h(y)
				c(d, h)
            })
        })
    }
    return(t(Z))
}



