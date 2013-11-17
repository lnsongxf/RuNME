# External Libraries
library('ipoptr')

# Set Problem Parameters

nAgents <- 25
nGoods  <- 25

# a
a <- rep(1,nAgents*nGoods)
#a <- runif(nAgents*nGoods)

# nu
nu <- rep(-2,nAgents*nGoods)
nu <- -1 * runif(nAgents*nGoods)

# lambda (Social Weights)
lambda <- rep(1,nAgents)

# endowment
#e <- rep(1,nAgents*nGoods)
e <- runif(nAgents*nGoods)

# calculate total amounts of each good (for constraint)
ej <- colSums(matrix(e, nAgents, nGoods, byrow=TRUE))

# Objective Function
eval_f <- function(x) {
    # Utility of each agent from each good
    util <- a * (x ^ (1 + nu)) / (1 + nu)
    
    # Sum weighted utility
    welfare <- sum(lambda * rowSums(matrix(x,nAgents,nGoods,byrow=TRUE)))
    
    return(-welfare)
}

# Gradient
eval_grad_f <- function(x) {

    l = rep(lambda, each=nGoods)

    return(- l * a * (x ^ nu))
}

# Constraint
eval_g <- function(x) {

    diff <- matrix(x,nAgents,nGoods,byrow=TRUE)
    diffs <- colSums(diff)

    return(diffs)
}

eval_jac_g_structure <- rep(list(c(seq(1,nAgents*nGoods))),nGoods)
eval_jac_g <- function(x) {

    jac = NULL
    for (i in 1:nGoods){
    	x = rep(0,nGoods)
	x[i] = 1

    	jac = append(jac, rep(x, nAgents))
    }

    return(jac)
}

print(length(eval_jac_g(e)))
print(length(unlist(eval_jac_g_structure)))

# Initial Guess = Endowment
x0 <- e

# Bounds
lb <- rep(0,nGoods * nAgents)
ub <- rep(Inf,nGoods * nAgents)

# Constraint Bounds
constraint_lb <- rep(0,nGoods)
constraint_ub <- ej

# Options
opts <- list("max_iter"=1000)

res <- ipoptr(x0=x0,
              eval_f=eval_f,
              eval_grad_f=eval_grad_f,
              eval_g=eval_g,
              eval_jac_g_structure=eval_jac_g_structure,
              eval_jac_g=eval_jac_g,
              lb = lb,
              ub = ub,
              constraint_lb = constraint_lb,
              constraint_ub = constraint_ub,
	      opts = opts)

print(res)
