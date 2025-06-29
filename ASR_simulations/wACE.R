'''
library(ape)
ace_weighted <- function (x, phy, type = "continuous", method = if (type == "continuous") "REML" else "ML",
    CI = TRUE, model = if (type == "continuous") "BM" else "ER",
    scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1, use.expm = FALSE,
    use.eigen = TRUE, marginal = FALSE, weights = 1.)
{
    if (!inherits(phy, "phylo"))
        stop("object phy is not of class phylo")
    if (is.null(phy$edge.length))
        stop("tree has no branch lengths")
    type <- match.arg(type, c("continuous", "discrete"))
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1)
        stop("phy is not rooted AND fully dichotomous.")
    if (length(x) != nb.tip)
        stop("length of phenotypic and of phylogenetic data do not match.")
        
    # Ensure weights are properly set up
    if (length(weights) == 1) {
        weights <- rep(weights, nb.tip)
    } else if (length(weights) != nb.tip) {
        stop("length of weights vector must match number of tips.")
    }
    
    if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label)) {
            x <- x[phy$tip.label]
            # Ensure weights are properly ordered if they have names
            if (!is.null(names(weights)) && all(names(weights) %in% phy$tip.label)) {
                weights <- weights[phy$tip.label]
            }
        } else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
    }
    
    obj <- list()
    if (kappa != 1)
        phy$edge.length <- phy$edge.length^kappa
    if (type == "continuous") {
        switch(method, REML = {
            minusLogLik <- function(sig2) {
                if (sig2 < 0) return(1e+100)
                V <- sig2 * vcv(phy)
                distval <- mahalanobis(x, center = mu, cov = V)
                logdet <- sum(log(eigen(V, symmetric = TRUE,
                  only.values = TRUE)$values))
                # Apply weights to the likelihood calculation
                (sum(weights * (log(2 * pi) + logdet/nb.tip + distval/nb.tip)))/2
            }
            mu <- rep(ace(x, phy, method = "pic")$ace[1], nb.tip)
            out <- nlm(minusLogLik, 1, hessian = TRUE)
            sigma2 <- out$estimate
            se_sgi2 <- sqrt(1/out$hessian)
            tip <- phy$edge[, 2] <= nb.tip
            minus.REML.BM <- function(p) {
                x1 <- p[phy$edge[, 1] - nb.tip]
                x2 <- numeric(length(x1))
                x2[tip] <- x[phy$edge[tip, 2]]
                x2[!tip] <- p[phy$edge[!tip, 2] - nb.tip]
                w <- numeric(length(x1))
                w[tip] <- weights[phy$edge[tip, 2]]
                w[!tip] <- 1
                # Apply weights to the tips
                -(-sum(w * (x1 - x2)^2/phy$edge.length)/(2 * sigma2) -
                  nb.node * log(sigma2))
            }
            out <- nlm(function(p) minus.REML.BM(p), p = rep(mu[1],
                nb.node), hessian = TRUE)
            obj$resloglik <- -out$minimum
            obj$ace <- out$estimate
            names(obj$ace) <- nb.tip + 1:nb.node
            obj$sigma2 <- c(sigma2, se_sgi2)
            if (CI) {
                se <- .getSEs(out)
                tmp <- se * qt(0.025, nb.node)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
        }, pic = {
            if (model != "BM") stop("the pic method can be used only with model = BM.")
            phy <- reorder(phy, "postorder")
            phenotype <- numeric(nb.tip + nb.node)
            phenotype[1:nb.tip] <- if (is.null(names(x))) x else x[phy$tip.label]
            contr <- var.con <- numeric(nb.node)

            ans <- .C(C_pic, as.integer(nb.tip), as.integer(phy$edge[,
                1]), as.integer(phy$edge[, 2]), as.double(phy$edge.length),
                as.double(phenotype), as.double(contr), as.double(var.con),
                as.integer(CI), as.integer(scaled))
            obj$ace <- ans[[5]][-(1:nb.tip)]
            names(obj$ace) <- nb.tip + 1:nb.node
            if (CI) {
                se <- sqrt(ans[[7]])
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
        }, ML = {
            if (model == "BM") {
                tip <- phy$edge[, 2] <= nb.tip
                dev.BM <- function(p) {
                  if (p[1] < 0) return(1e+100)
                  x1 <- p[-1][phy$edge[, 1] - nb.tip]
                  x2 <- numeric(length(x1))
                  x2[tip] <- x[phy$edge[tip, 2]]
                  x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
                  w <- numeric(length(x1))
                  w[tip] <- weights[phy$edge[tip, 2]]
                  w[!tip] <- 1
                  # Apply weights to the likelihood calculation
                  -2 * (-sum(w * (x1 - x2)^2/phy$edge.length)/(2 *
                    p[1]) - nb.node * log(p[1]))
                }
                out <- nlm(function(p) dev.BM(p), p = c(1, rep(mean(x),
                  nb.node)), hessian = TRUE)
                obj$loglik <- -out$minimum/2
                obj$ace <- out$estimate[-1]
                names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
                se <- .getSEs(out)
                obj$sigma2 <- c(out$estimate[1], se[1])
                if (CI) {
                  tmp <- se[-1] * qt(0.025, nb.node)
                  obj$CI95 <- cbind(obj$ace + tmp, obj$ace -
                    tmp)
                }
            }
        }, GLS = {
            if (is.null(corStruct)) stop("you must give a correlation structure if method = 'GLS'.")
            if (class(corStruct)[1] == "corMartins") M <- corStruct[1] *
                dist.nodes(phy)
            if (class(corStruct)[1] == "corGrafen") phy <- compute.brlen(attr(corStruct,
                "tree"), method = "Grafen", power = exp(corStruct[1]))
            if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
                dis <- dist.nodes(attr(corStruct, "tree"))
                MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
                M <- dis[as.character(nb.tip + 1), MRCA]
                dim(M) <- rep(sqrt(length(M)), 2)
            }
            one2n <- 1:nb.tip
            varAY <- M[-one2n, one2n]
            varA <- M[-one2n, -one2n]
            DF <- data.frame(x)
            V <- corMatrix(Initialize(corStruct, DF), corr = FALSE)
            # Incorporate weights into the inverse variance-covariance matrix
            W <- diag(weights)
            invV <- solve(V)
            invVW <- W %*% invV
            o <- gls(x ~ 1, DF, correlation = corStruct, weights = weights)
            GM <- o$coefficients
            obj$ace <- drop(varAY %*% invVW %*% (x - GM) + GM)
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                se <- sqrt((varA - varAY %*% invVW %*% t(varAY))[cbind(1:nb.node,
                  1:nb.node)])
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
        })
    }
    else {
        if (method != "ML")
            stop("only ML estimation is possible for discrete characters.")
        if (any(phy$edge.length < 0))
            stop("some branches have negative length")
        if (!is.factor(x))
            x <- factor(x)
        nl <- nlevels(x)
        lvls <- levels(x)
        x <- as.integer(x)
        if (is.character(model)) {
            rate <- matrix(NA, nl, nl)
            switch(model, ER = np <- rate[] <- 1, ARD = {
                np <- nl * (nl - 1)
                rate[col(rate) != row(rate)] <- 1:np
            }, SYM = {
                np <- nl * (nl - 1)/2
                sel <- col(rate) < row(rate)
                rate[sel] <- 1:np
                rate <- t(rate)
                rate[sel] <- 1:np
            })
        }
        else {
            if (ncol(model) != nrow(model))
                stop("the matrix given as 'model' is not square")
            if (ncol(model) != nl)
                stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
            rate <- model
            np <- max(rate)
        }
        index.matrix <- rate
        tmp <- cbind(1:nl, 1:nl)
        index.matrix[tmp] <- NA
        rate[tmp] <- 0
        rate[rate == 0] <- np + 1
        liks <- matrix(0, nb.tip + nb.node, nl)
        TIPS <- 1:nb.tip
        liks[cbind(TIPS, x)] <- 1
        if (anyNA(x))
            liks[which(is.na(x)), ] <- 1
        phy <- reorder(phy, "postorder")
        Q <- matrix(0, nl, nl)
        e1 <- phy$edge[, 1]
        e2 <- phy$edge[, 2]
        EL <- phy$edge.length
        if (use.eigen) {
            dev <- function(p, output.liks = FALSE) {
                if (any(is.nan(p)) || any(is.infinite(p)))
                  return(1e+50)
                comp <- numeric(nb.tip + nb.node)
                Q[] <- c(p, 0)[rate]
                diag(Q) <- -rowSums(Q)
                decompo <- eigen(Q)
                lambda <- decompo$values
                GAMMA <- decompo$vectors
                invGAMMA <- solve(GAMMA)
                for (i in seq(from = 1, by = 2, length.out = nb.node)) {
                  j <- i + 1L
                  anc <- e1[i]
                  des1 <- e2[i]
                  des2 <- e2[j]
                  v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*%
                    invGAMMA %*% liks[des1, ]
                  v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*%
                    invGAMMA %*% liks[des2, ]
                  # Apply weights to the likelihood calculation for tips
                  if (des1 <= nb.tip) {
                    v.l <- v.l^weights[des1]
                  }
                  if (des2 <= nb.tip) {
                    v.r <- v.r^weights[des2]
                  }
                  v <- v.l * v.r
                  comp[anc] <- sum(v)
                  liks[anc, ] <- v/comp[anc]
                }
                if (output.liks)
                  return(liks[-TIPS, , drop = FALSE])
                dev <- -2 * sum(log(comp[-TIPS]))
                if (is.na(dev))
                  Inf
                else dev
            }
        }
        else {
            if (!requireNamespace("expm", quietly = TRUE) &&
                use.expm) {
                warning("package 'expm' not available; using function 'matexpo' from 'ape'")
                use.expm <- FALSE
            }
            E <- if (use.expm)
                expm::expm
            else matexpo
            dev <- function(p, output.liks = FALSE) {
                if (any(is.nan(p)) || any(is.infinite(p)))
                  return(1e+50)
                comp <- numeric(nb.tip + nb.node)
                Q[] <- c(p, 0)[rate]
                diag(Q) <- -rowSums(Q)
                for (i in seq(from = 1, by = 2, length.out = nb.node)) {
                  j <- i + 1L
                  anc <- e1[i]
                  des1 <- e2[i]
                  des2 <- e2[j]
                  v.l <- E(Q * EL[i]) %*% liks[des1, ]
                  v.r <- E(Q * EL[j]) %*% liks[des2, ]
                  # Apply weights to the likelihood calculation for tips
                  if (des1 <= nb.tip) {
                    v.l <- v.l^weights[des1]
                  }
                  if (des2 <= nb.tip) {
                    v.r <- v.r^weights[des2]
                  }
                  v <- v.l * v.r
                  comp[anc] <- sum(v)
                  liks[anc, ] <- v/comp[anc]
                }
                if (output.liks)
                  return(liks[-TIPS, , drop = FALSE])
                dev <- -2 * sum(log(comp[-TIPS]))
                if (is.na(dev))
                  Inf
                else dev
            }
        }
        out <- nlminb(rep(ip, length.out = np), function(p) dev(p),
            lower = rep(0, np), upper = rep(1e+50, np))
        obj$loglik <- -out$objective/2
        obj$rates <- out$par
        oldwarn <- options("warn")
        options(warn = -1)
        out.nlm <- try(nlm(function(p) dev(p), p = obj$rates,
            iterlim = 1, stepmax = 0, hessian = TRUE), silent = TRUE)
        options(oldwarn)
        
        obj$index.matrix <- index.matrix
        if (CI) {
            lik.anc <- dev(obj$rates, TRUE)
            if (!marginal) {
                Q[] <- c(obj$rates, 0)[rate]
                diag(Q) <- -rowSums(Q)
                for (i in seq(to = 1, by = -2, length.out = nb.node)) {
                  anc <- e1[i] - nb.tip
                  des1 <- e2[i] - nb.tip
                  if (des1 > 0) {
                    P <- matexpo(Q * EL[i])
                    tmp <- lik.anc[anc, ]/(lik.anc[des1, ] %*%
                      P)
                    lik.anc[des1, ] <- (tmp %*% P) * lik.anc[des1,
                      ]
                  }
                  j <- i + 1L
                  des2 <- e2[j] - nb.tip
                  if (des2 > 0) {
                    P <- matexpo(Q * EL[j])
                    tmp <- lik.anc[anc, ]/(lik.anc[des2, ] %*%
                      P)
                    lik.anc[des2, ] <- (tmp %*% P) * lik.anc[des2,
                      ]
                  }
                  lik.anc <- lik.anc/rowSums(lik.anc)
                }
            }
            rownames(lik.anc) <- nb.tip + 1:nb.node
            colnames(lik.anc) <- lvls
            obj$lik.anc <- lik.anc
        }
    }
    if (!is.null(phy$node.label)) {
        if (!is.null(obj$ace))
            names(obj$ace) <- phy$node.label
        if (!is.null(obj$CI95))
            rownames(obj$CI95) <- phy$node.label
        if (!is.null(obj$lik.anc))
            rownames(obj$lik.anc) <- phy$node.label
    }
    obj$call <- match.call()
    # Store the weights in the output
    obj$weights <- weights
    class(obj) <- "ace"
    obj
}

