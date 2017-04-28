
## 28 Jan 2015
## Functions to solve for initial conditions for counts of individuals in each category
library(deSolve)

## Function to calculate proportion of individuals in each compartment
getProportions <- function(params){
  A <- matrix(
    c( params["epsilonI0"] - (params["muI0"] + params["diagI0"] + params["gammaI0I1"]), params["epsilonI1"], params["epsilonI2"], params["epsilonJ0"], params["epsilonJ1"], params["epsilonJ2"],
       params["gammaI0I1"], -(params["muI1"] + params["diagI1"] + params["gammaI1I2"]), 0, 0, 0, 0,
       0, params["gammaI1I2"], -(params["muI2"] + params["diagI2"]), 0, 0, 0,
       params["diagI0"], 0, 0, -(params["muJ0"] + params["gammaJ0J1"]), 0, 0,
       0, params["diagI1"], 0, params["gammaJ0J1"], -(params["muJ1"] + params["gammaJ1J2"]), 0,
       0, 0, params["diagI2"], 0, params["gammaJ1J2"], -params["muJ2"]
    ),
    nrow = 6,
    byrow = TRUE
  )
  vecs <- eigen(A)$vectors
  vals <- eigen(A)$values
  index <- which.min(abs(vals))
  v <- vecs[,index]
  props <- v/sum(v)
  names(props) <- c("pI0", "pI1", "pI2", "pJ0", "pJ1", "pJ2")
  props
}

## Function to calculate absolute counts of individuals in each compartment
## given the total number of diagnosed individuals 
getCounts <- function(props, sumJ){
  propsJ <- props[c("pJ0", "pJ1", "pJ2")]
  propsJ <- propsJ/sum(propsJ) 
  numJ <- propsJ * sumJ
  propsI <- props[c("pI0", "pI1", "pI2")]
  pI <- sum(propsI)
  sumI <- sumJ/(1 - pI) - sumJ
  numI <- propsI * (sumI + sumJ)
  counts <- c(numI, numJ)
  names(counts) <- c("numI0", "numI1", "numI2", "numJ0", "numJ1", "numJ2")
  round(counts)
}

## An alternative function for initializing counts based on a constant flow into the I0 class
## and based on all rates into and out of compartments known
## This approach does not attempt to be self consistent with the epsilons

## Specify the model
model <- function(elapsed_time, y, params){
  with(as.list(c(y, params)), {
    dy1 <- h - (muI0 + diagI0 + gammaI0I1)*y[1]                  # I0
    dy2 <- gammaI0I1*y[1] - (muI1 + diagI2 + gammaI1I2)*y[2]           # I1
    dy3 <- gammaI1I2*y[2] - (muI2 + diagI3)*y[3]                   # I2 
    dy4 <- diagI1*y[1] - (muJ0 + gammaJ0J1)*y[4]                   # J0
    dy5 <- gammaJ0J1*y[4] + diagI2*y[2] - (muJ1 + gammaJ1J2)*y[5]     # J1
    dy6 <- gammaJ1J2*y[5] + diagI3*y[3] - muJ2*y[6]                # J2
    list(c(dy1, dy2, dy3, dy4, dy5, dy6))
  })
}

## Objective function to find the growth rate
objective.fn <- function(growth_rate, parms, n_diagnosed, n_yrs, model){
  parms <- c(parms, h = growth_rate)
  yini <- c(y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0, y6 = 0)
  times <- seq(0, n_yrs, 0.01)
  out <- ode(y = yini, times, model, parms)
  observed <- sum(out[nrow(out), c("y4","y5","y6")])
  abs(n_diagnosed - observed)
}

get_initial_counts_h_constant <- function(params, elapsed_time, n_diagnosed, initial_counts){
  ## Specify the model
  m <- model(elapsed_time, initial_counts, params)
  ## Find h such that n_diagnosed is correct
  res <- optimize(f = objective.fn, interval = c(0, 100), parms = parms, n_diagnosed = n_diagnosed, n_yrs = n_yrs, model = m)
  ## Write h to the parameter vector
  params <- c(params, h = res$minimum)
  ## Solve the differential equations numerically to find counts of each class at t0
  out <- ode(y = initial_counts, times, model, parms)
  ## Return counts
  counts <- round(out[nrow(out), -1])
  names(counts) <- c('I0','I1','I2','J0','J1','J2')
  counts
}













