#
#Spatio-temporal model with shared and specific effects
#

library(INLA)
library(spdep)
library(spacetime)

#Load data
load("joint_st_dismap_sim_data.RData")

#Data together
d <- data.frame(CPROV = rep(stmapa.OyE.sim@data$CPROV, 3))
d$ANYO <- rep(stmapa.OyE.sim@data$ANYO, 3)
d$Expected <- with(stmapa.OyE.sim@data, c(Esp.OC, Esp.Eso, Esp.Sto))
d$Observed <- with(stmapa.OyE.sim@data, c(Obs.OC, Obs.Eso, Obs.Sto))


#Number of space-time units per disease
n.st <- nrow(stmapa.OyE.sim@data)
d$DISEASE <- rep(c("Bucal_Cavity", "Esophagus", "Estomach"), each = n.st)

#Index for spatial provinces CPROV ranges from 1 to 50
d$sp.idx <- rep(1:47, 19)

#Index for time
d$tm.idx <- d$ANYO - 1995

#Add indices for effects
d$disease <- rep(1:3, each = n.st)
d$intercept <- as.factor(d$disease)

#Dummy indices for space and time
d$s.dummy <- NA
d$t.dummy <- NA

#Create spacial indices for specific effects
d$s.1 <- NA
d$s.1 [d$disease == 1] <- d$sp.idx [d$disease == 1]
d$s.2 <- NA
d$s.2 [d$disease == 2] <- d$sp.idx [d$disease == 2]
d$s.3 <- NA
d$s.3 [d$disease == 3] <- d$sp.idx [d$disease == 3]

#Create temporal indices for specific effects
d$t.1 <- NA
d$t.1 [d$disease == 1] <- d$tm.idx [d$disease == 1] 
d$t.2 <- NA
d$t.2 [d$disease == 2] <- d$tm.idx [d$disease == 2]
d$t.3 <- NA
d$t.3 [d$disease == 3] <- d$tm.idx [d$disease == 3]

#Spatial adjacency matrix
W.sp <- as(nb2mat(AdySP, style = "B"), "Matrix")


#Spatio-temporal adjancency
W.tm <- Diagonal(19, x = 0)

W.tm[1, 1 + 1] <- 1
W.tm[19, 19 -1] <- 1
for(i in 2:18) {
    W.tm[i, i - 1] <- 1
    W.tm[i, i + 1] <- 1
}

#Prior for coefficients of copied effects
#  IMPORTANT: In paper, beta ~logNormal (0, prec = 1/5.9))
#   In paper, E[beta] = exp(0 + 5.9/2); Var[Beta] = [exp(5.9) - 1]*exp(2*0+5.9)
#   In THIS R CODE, use Normal prior with the same mean and variance
#     mean = 19.10595; prec = 1/132887.3
#prior.beta.s = list(prior = "normal", param = c(19.10595, 1/132887.3),
#  fixed = FALSE)
#prior.beta.t = list(prior = "normal", param = c(19.10595, 1/132887.3),
#  fixed = FALSE)

#Centered at the mode: exp(0 - 5.9)
#prior.beta.s = list(prior = "normal", param = c(0.002739445, 1/132887.3),
#  fixed = FALSE)
#prior.beta.t = list(prior = "normal", param = c(0.002739445, 1/132887.3),
#  fixed = FALSE)

#Prior to match the results in the paper
prior.beta.s = list(prior = "normal", param = c(0, 1 / 5.9),
  fixed = FALSE, initial = 0.01)
#prior.beta.s = list(prior = 
#  "expression:
#   logdens = -theta - 0.5 * log(5.9) + ((-(theta - 0)^2) / (2 * 5.9));
#   return(logdens)", fixed = FALSE)

prior.beta.t = list(prior = "normal", param = c(0, 1 / 5.9),
  fixed = FALSE, initial = 0.01)

#
# Prior for precisions
#
#prior.prec = list(initial = 0, fixed = TRUE)
# Gamma prior on precision
#file.saved <- "st_model_GAMMA.RData"
#prior.prec = list(prior = "loggamma", param = c(0.01, 0.01), initial = 0)

#Flat prior on sigma: Ugarte et al. (2018):
file.saved <- "st_model_UNIFORM.RData"
prior.prec = list(prior = 
  "expression:
   logdens = -log_precision / 2;
   return(logdens)",
  initial = 0)


# Half-Cuachy prior on sigma
#file.saved <- "st_model_HCAUCHY.RData"
#log(2) -log(pi) - log(25) = -3.670459
#prior.prec = list(prior = 
#  "expression:
#   scale_param = 25;
#   logdens = -3.670459 -log(1 + (exp(-theta) / (scale_param^2))) - (theta/2);
#   return(logdens)",
#  initial = 0)

#Indices for spatial  disease-specific effects
d$sp.idx1 <- d$s.1
d$sp.idx2 <- d$s.2
d$sp.idx3 <- d$s.3

#Indices for temporal  disease-specific effects
d$tm.idx1 <- d$t.1
d$tm.idx2 <- d$t.2
d$tm.idx3 <- d$t.3




# Scale models?
inla.scale <- FALSE

formula = Observed ~ -1 +
    intercept +
    f(sp.idx1, model = "besag",  scale.model = inla.scale,
      graph = W.sp, hyper = list(prec = prior.prec)) +
    f(sp.idx2, model = "besag",  scale.model = inla.scale,
      graph = W.sp, hyper = list(prec = prior.prec)) +
    f(sp.idx3, model = "besag",  scale.model = inla.scale,
      graph = W.sp, hyper = list(prec = prior.prec)) +
    f(s.dummy, model = "besag",  scale.model = inla.scale, graph = W.sp,
      hyper = list(prec = prior.prec))+
    f(s.1, copy = "s.dummy",  range = c(0, Inf),
      hyper = list(beta = prior.beta.s)) +
    f(s.2, copy = "s.dummy",  range = c(0, Inf),
      hyper = list(beta = prior.beta.s)) +
    f(s.3, copy = "s.dummy",  range = c(0, Inf),
      hyper = list(beta = prior.beta.s)) +
    f(tm.idx1, model = "besag", scale.model = inla.scale,
      graph = W.tm, hyper = list(prec = prior.prec)) + 
    f(tm.idx2, model = "besag", scale.model = inla.scale,
      graph = W.tm, hyper = list(prec = prior.prec)) + 
    f(tm.idx3, model = "besag", scale.model = inla.scale,
      graph = W.tm, hyper = list(prec = prior.prec)) + 
    f(t.dummy, model = "besag", scale.model = inla.scale,
      graph = W.tm, hyper = list(prec = prior.prec)) +
    f(t.1, copy = "t.dummy",  range = c(0, Inf),
      hyper = list(beta = prior.beta.t)) +
    f(t.2, copy = "t.dummy",  range = c(0, Inf),
      hyper = list(beta = prior.beta.t)) +
    f(t.3, copy = "t.dummy",  range = c(0, Inf),
      hyper = list(beta = prior.beta.t))


r = inla(formula,  data = d,  family = "poisson", control.predictor=list(compute=TRUE),
  E = Expected, verbose = TRUE,
  control.compute = list(config = TRUE),
  control.inla(strategy = "laplace", npoints = 21))
  #, diff.logdens = 15, dz = 0.25))

#  r <- inla.hyperpar(r, verbose = TRUE, diff.logdens = 5)
r <- inla.rerun(r)

summary(r)

#Display spatial effects
spain <- SpatialPolygonsDataFrame(stmapa.OyE.sim@sp,
  data.frame(SHARED = r$summary.random$s.dummy[, "mean"]),
  match.ID = FALSE )
#Specific
spain$SP.BUCAL <- r$summary.random$sp.idx1[, "mean"]
spain$SP.ESO <- r$summary.random$sp.idx2[, "mean"]
spain$SP.ESTO <- r$summary.random$sp.idx3[, "mean"]

spplot(spain, c("SHARED", "SP.BUCAL", "SP.ESO", "SP.ESTO"),
  at = seq(-0.4, 0.5, by = 0.1), col.regions = gray(10:1/11))

# Posterior means of total spatial effect
spain$TOTSP.BUCAL <- r$summary.random$sp.idx1[, "mean"] + 
  r$summary.random$s.1[, "mean"]
spain$TOTSP.ESO <- r$summary.random$sp.idx2[, "mean"] + 
  r$summary.random$s.2[, "mean"]
spain$TOTSP.ESTO <- r$summary.random$sp.idx3[, "mean"] + 
  r$summary.random$s.3[, "mean"]

dev.new()
spplot(spain, c("TOTSP.BUCAL", "TOTSP.ESO", "TOTSP.ESTO"),
  at = seq(-0.4, 0.5, by = 0.1), col.regions = gray(10:1/11))

#Temporal trends
dev.new()
par(mfrow = c(2, 2))
plot(1996:2014, r$summary.random$t.dummy[, "mean"], ylim = c(-0.4, 0.4))
plot(1996:2014, r$summary.random$tm.idx1[, "mean"], ylim = c(-0.4, 0.4))
plot(1996:2014, r$summary.random$tm.idx2[, "mean"], ylim = c(-0.4, 0.4))
plot(1996:2014, r$summary.random$tm.idx3[, "mean"], ylim = c(-0.4, 0.4))

# Temporal total effect (post. means)
dev.new()
plot(1996:2014, r$summary.random$tm.idx1[, "mean"] +
  r$summary.random$t.1[, "mean"], ylim = c(-0.3, 0.3), type = "l")
lines(1996:2014, r$summary.random$tm.idx2[, "mean"] +
  r$summary.random$t.2[, "mean"], lty = 2)
lines(1996:2014, r$summary.random$tm.idx3[, "mean"] +
  r$summary.random$t.3[, "mean"], lty = 3)


# Use inla.posterior.samples() to obtain samples from the joint post.
# dsitribution of hyperparameters
samples <- inla.posterior.sample(1000, r)
betas <- sapply(1:1000, function(x) (samples[[x]]$hyperpar)[9:14])
betas <- as.data.frame(t(betas))
rm(samples)
plot(betas)
cor(betas)


#Use the posterior density of the hyperparameters
xx <- (r$joint.hyper)
probs <- xx[, 15]
probs <- exp(probs - max(probs))
probs <- probs / sum(probs)

#save(file = "file.saved.Rdata", list = ls())

