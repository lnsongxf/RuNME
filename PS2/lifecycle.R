# External Libraries
library('ipoptr')

ipoptr_solve <- function(gamma, eta, beta, B, wmax, amin, tauL, tauK, r, T, period=1){

  # Adjust for period other than 1 year
  if (period != 1){
    r <- (1 + r) ^ period - 1
    beta <- (1 + beta) ^ period - 1
  }


  # Construct Wage Function
  wage <- function(t) {
    return(1 - (4/T) * (1-wmax) * t + (4/(T*T)) * (1-wmax) * t^2)
  }


  # Objective Function
  eval_f <- function(x) {

    # Utility at each time t
    c <- (x[1:T] ^ (1 - gamma)) / (1 - gamma) # consumption
    l <- B * (x[(T+1):(2*T)] ^ (1 + eta)) / (1 + eta) # labor
    
    u <- c + l

    # Create betas
    betas <- beta ^ c(1:T)

    discounted.u <- sum(betas * u)

    return(-discounted.u)

  }

  # Gradient of Objective Function
  eval_grad_f <- function(x) {
  
    grad = rep(0,3*T)

    for (i in 1:T) {
      # c_t  
      grad[i] = beta^i * (x[i]^(- gamma))

      # l_t
      grad[T+i] = beta^i * B * (x[T+i]^eta)
    }

    return(grad)
  }

  # Constraints
  eval_g <- function(x) {

    g = rep(0,T)
    
    g[1] = x[T+1] * wage(0) * (1 - tauL) - x[1] - x[2*T+1]

    g[T] = x[T+T] * wage(T-1) * (1 - tauL) + x[2*T + T] * (1 + r * (1 - tauK)) - x[T]

    for (i in 2:T-1) {
      g[i] = x[T+i] * wage(i-1) * (1 - tauL) + x[2*T + i] * (1 + r * (1 - tauK)) - x[i] - x[2*T + i] 
    } 

    return(g)

  }

  # Jacobian of Constraint
  eval_jac_g_structure <- rep(list(c(seq(1, 3*T))), T)
  print(eval_jac_g_structure)
  eval_jac_g <- function(x) {
    jac = NULL
    for (t in 1:T){
      tmp = rep(0,3*T)

      # time t constraint depends on c_t, l_t, A_t, A_t-1

      # wage
      tmp[T+t] = wage(t) * (1 - tauL)

      # consumption
      tmp[t] = -1

      # Savings
      if (t != T){
        tmp[2*T+t] = -1
      }

      # Savings from Previous Pd
      if (t != 1){
        tmp[2*T+t-1] = (1 + r * (1 - tauK))
      }

      jac = append(jac, tmp)
    }

    return(jac)
  }

  # Initial Guess
  x0 <- append(rep(1,2*T), rep(0,T))

  # Bounds

  lb <- append(rep(0,2*T), append(rep(amin,T-1), 0))
  ub <- append(rep(Inf, 3*T-1), 0)
  #print(lb)
  #print(ub)
  # Constraint Bounds
  constraint_lb <- rep(0, T)
  constraint_ub <- rep(0, T)

  # Options
  opts <- list("max_iter"=1000)
  #return(0)
  res <- ipoptr(x0=x0,
		eval_f = eval_f,
		eval_grad_f = eval_grad_f,
		eval_g = eval_g,
		eval_jac_g_structure = eval_jac_g_structure,
		eval_jac_g = eval_jac_g,
		lb = lb,
		ub = ub,
		constraint_lb = constraint_lb,
		constraint_ub = constraint_ub,
		opts = opts)

		       
  print(res)
  return(res)
}	

test1 <- function(){
  gamma = 1.5
  eta = 3
  beta = 0.9
  B = -1.

  wmax = 10.
  amin = 0.
  r = 0

  tauL = 0.
  tauK = 0.

  T = 3
  period = 1

  ipoptr_solve(gamma=gamma,
	       eta=eta,
	       beta=beta,
	       B=B,
	       wmax=wmax,
	       amin=amin,
	       r=r,
	       tauL=tauL,
	       tauK=tauK,
	       T=T,
	       period=period)  
}

test0 <- function(){
  T = 3
  gamma = 1.5
  eta = 3
  beta = 0.9
  B = -1

  x = c(5,1,1,1,1,1,0,0,0)
  
  u <- (x[1:T] ^ (1 - gamma)) / (1 - gamma) + B * (x[T+1:2*T] ^ (1 + eta)) / (1 + eta)


  betas <- beta ^ c(1:T)

  discounted.u <- sum(betas * u)

  return(discounted.u)
}

main <- function(){

  # Run tests
  #test0()
  test1()

}

main()