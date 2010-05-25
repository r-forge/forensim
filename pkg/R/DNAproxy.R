DNAproxy <- function (tab, x) {
    if (length(tab[[1]]) == 6) {
        mix <- TRUE
    }
    else {
        mix <- FALSE
    }
    DNAproxy.loc <- function(tabloc, x) {
        I <- intersect(tabloc$Expected[tabloc$D==0], tabloc$Observed[tabloc$D==0]) 
        if (length(I) != 0) {
            n <- rep(0, length(I))
            names(n) <- I
            H <- tabloc[, "Heights"]
            names(H) <- tabloc[, "Expected"]
            h <- n
            s <- n
            if (mix) {
                x1 <- tabloc[tabloc$Expected %in% I, x]
                x2 <- tabloc[tabloc$Expected %in% I, names(tabloc)[2:3][which(!(names(tabloc)[2:3]) %in% 
                  x)]]
                for (i in 1:length(I)) {
                  h[I[i]] <- H[I[i]]
                  if (x1[i] == 0) {
                    n[I[i]] <- 0
                    s[I[i]] <- 0
                  }
                  if (x1[i] != 0) {
                    if (x1[i] == 1 & x2[i] == 0) {
                      n[I[i]] <- 1
                      s[I[i]] <- 1
                    }
                    if (x1[i] == 1 & x2[i] > 1) {
                      n[I[i]] <- 1
                      s[I[i]] <- 0
                    }
                    if (x1[i] == 2 & x2[i] == 0) {
                      n[I[i]] <- 2
                      s[I[i]] <- 1
                    }
                    if (x1[i] == 1 & x2[i] == 1) {
                      n[I[i]] <- 1
                      s[I[i]] <- 0
                    }
                    if (x1[i] == 2 & x2[i] > 1) {
                      n[I[i]] <- 2
                      s[I[i]] <- 0
                    }
                  }
                }
                cbind.data.frame(deno = sum(n * s), num = sum(h * s))
            }
            else {
                for (i in 1:length(I)) {
                  if (length(tabloc$Expected) == 2) {
                    n[I[i]] <- 1
                  }
                  if (length(tabloc$Expected) == 1) {
                    n[I[i]] <- 2
                  }
                  h[I[i]] <- H[I[i]]
                }
                cbind.data.frame(deno = sum(n), num = sum(h))
            }
        }
        else {
            cbind.data.frame(deno = 0, num = 0)
        }
    }
    ind <- sum(unlist(sapply(tab, DNAproxy.loc, x)[2, ]))/sum(unlist(sapply(tab, DNAproxy.loc, x)[1, ]))
    if (is.na(ind)) {
        ind <- 0
    }
    ind
}

tabDNAproxy <- function (x, y = NULL, geno, tabcsv) 
{
    tabo <- recordDrop(x, y, geno, tabcsv)
    if (length(tabo[[1]]) == 6) {
        mix <- TRUE
    }
    else {
        mix <- FALSE
    }
    m <- names(geno)
    perloc <- vector("list", length(m))
    names(perloc) <- m
    if (mix) {
        ind1 <- DNAproxy(tabo, paste("c", x, sep = ""))
        ind2 <- DNAproxy(tabo, paste("c", y, sep = ""))
        for (M in names(tabo)) {
            Dloc <- NULL
            Hestim <- NULL
            for (h in 1:nrow(tabo[[M]])) {
                Dloc <- c(Dloc, tabo[[M]][h, "D"])
                Hestim <- c(Hestim, sum(tabo[[M]][h, 2:3] * c(ind1, 
                  ind2)))
            }
            perloc[[M]] <- cbind.data.frame(Dloc, Hestim)
        }
    }
    else {
        ind <- DNAproxy(tabo)
        for (M in names(tabo)) {
            Dloc <- NULL
            Hestim <- NULL
            for (h in 1:nrow(tabo[[M]])) {
                Dloc <- c(Dloc, tabo[[M]][h, "D"])
                if (length(tabo[[M]]$Expected) == 1) {
                  state <- 2
                }
                if (length(tabo[[M]]$Expected) == 2) {
                  state <- 1
                }
                Hestim <- c(Hestim, sum(state * ind))
            }
            perloc[[M]] <- cbind.data.frame(Dloc, Hestim)
        }
    }
    perloc
}
