#!/bin/python

############
# PS2
# ECON 41800
# 
# R YANG
# 2013
############

## External Libraries

# numerical/matrix library
import numpy

# optimization library
import scipy.optimize

def testSuite():
    #test0()
    test1()
    #test2()
    #test3()
    return 0

# Utility Function test
"""
def test0():
    
    gamma = 3
    nu = 5.
    beta = 1
    B = -1
    T = 3

    u = utility(gamma, nu, beta, B, T)
    mu = marg_utility(gamma, nu, beta, B, T)

    c = [1,1,1]
    l = [2,1,1]
    a = [0,0,0]

    x = []
    x = x + c
    x = x + l
    x = x + a
    x = numpy.array(x)

    print u(x)
    print mu(x)
"""
# Single period test

def test1():

    # Define parameters
    gamma = 3
    nu = 3
    beta = 0.1
    B = -1.

    wmax = 1.
    amin = 0.
    r = 0.

    tauL = 0.
    tauK = 0.

    T = 4
    period = 1

    lifeCycle(gamma, nu, B, beta, wmax, amin, tauL, tauK, r, T, period=1)
    return 0


## Utility Parameters
# gamma
# nu
# beta - annual discount rate
# B

## Wage Parameters
# wmax - max wage (occurs at t=T/2)

## Tax Parameters
# tauL - proportional tax rate on labor
# tauK - proportional tax rate on interest 

## Other Parameters
# amin
# r - annual interest rate
# T - number of time periods
# period - number of years per period (can be decimal)

## Builds Utility function with given parameters
## NOTE: Returns negative utility for convenient use with minimization algorithms

def utility(gamma, nu, beta, B, T):

    print (gamma, nu, beta, B, T)

    # Construct array of discounts
    betas = beta ** numpy.array(xrange(T))


    #consumption = lambda x: (x[0:T]**(1 - gamma)) / (1 - gamma)
    #labor = lambda x: (x[T:2*T]**(1 + nu)) / (1 + nu)

    #u = lambda x: -numpy.dot(betas, consumption(x) + B * labor(x))

    u = lambda x: -numpy.dot(betas, (x[0:T]**(1 - gamma)) / (1 - gamma) + B * (x[T:2*T]**(1 + nu)) / (1 + nu))

    return u

def marg_utility(gamma, nu, beta, B, T):
    
    betas = beta ** numpy.array(xrange(T))

    consumption = lambda x: - betas * x[0:T]**(-gamma)
    labor = lambda x: -B * betas * x[T:(2*T)]**(nu)
    savings = lambda x: numpy.zeros(T)

    mu = lambda x: numpy.concatenate((consumption(x), labor(x), savings(x)))
    return mu

def foo(x):
    print x
    return x

## Finds path of consumption, labor, and savings which maximizes utility

def lifeCycle(gamma, nu, B, beta, wmax, amin, tauL, tauK, r, T, period=1):

    if period != 1:
        ## If period not yearly, adjust interest and discount rates
        # period = 1 -> 1 year per period

        # adjust interest (r)
        r = (1 + r) ** (period) - 1

        # adjust discount factor (beta) 
        # ** double check this **
        beta = (1 + beta) ** (period) - 1
    
    ## Construct wage as function of time
    # w(t) = 1 - (4 / T)(1 - wmax)t + (4 / T^2)(1 - wmax)t^2
    #wage = lambda t: 1 - (4 / T) * (1 - wmax) * t + (4 / (T * T)) * (1 - wmax) * t * t
    wage = numpy.array(range(T))
    wage = 1 - (4/T) * (1 - wmax) * wage + (4 / (T**2)) * (1 - wmax) * (wage ** 2)
    print "WAGE:", wage

    ## Construct utility function (c,l) -> u
    ## utility = lambda cl: (cl[0]**(1-gamma))/(1-gamma) + (B * cl[1]**(1+nu))/(1-gamma)

    ## Note: x variable = [c_1 ... c_T, l_1 ... l_T, A_1 ... A_T] = 3T by 1 matrix
    ## c_t = x[t]
    ## l_t = x[T + t]
    ## A_T = x[2*T + t]

    ## Set up objective function
    # max \sum \beta^i u(c,l)
    obj = utility(gamma, nu, beta, B, T)

    ## Set up jacobian of objective function
    jac = marg_utility(gamma, nu, beta, B, T)

    ## Set up constraints

    # w_t * l_t * (1 - tauL) + A_{t-1} * (1 + r * (1 - tauK)) = c_t + A_t

    constraints = []
    cons = []

    for t in xrange(T):
        if T == 1:
            # No savings in 1 pd case - this is kind of hacky
            f = lambda x: wage[t] * x[T] * (1 - tauL) - x[0] + 0 * (foo(x))
            constraints.append({'type':'eq', 'fun':f})
            cons.append(f)
            break

        if t == 0:
            # In first period, no initial assets
            f = lambda x: wage[t] * x[T] * (1 - tauL) - x[0] - x[2*T]
            constraints.append({'type':'eq', 'fun':f})
            cons.append(f)
            continue

        if t == T-1:
            # In last period, no final assets
            f = lambda x: wage[t] * x[T+t] * (1 - tauL) + x[2*T + t - 1] * (1 + r * (1 - tauK)) - x[t]
            constraints.append({'type':'eq', 'fun':f})
            cons.append(f)
            continue

        if t > 0 and t < T-1:
            f = lambda x: wage[t] * x[T+t] * (1 - tauL) + x[2*T + t - 1] * (1 + r * (1 - tauK)) - x[t] - x[2*T+t]
            constraints.append({'type':'eq', 'fun':f})
            cons.append(f)
            continue

    ## Set up initial guess (all 0)
    x0 = numpy.array([0] * 3 * T)

    ## Set up bounds
    bounds = []

    # minimal consumption, labor
    # c_t, l_t >= 0 for all t

    bounds = bounds + (2*T)*[(0.0,1e10)]

    # minimal capital (max borrowing)
    # A_t >= amin for all t

    bounds = bounds + T*[(amin,1e10)]

    #bounds = []
    #bounds = bounds + (3*T)*[(0,10)]
    print bounds

    ## Select optimization algorithm
    solver = 'SLSQP'

    ## Solve for optimal consumption, labor, and savings
    res = scipy.optimize.fmin_slsqp(obj, x0, fprime = jac, eqcons = cons, f_eqcons = None, bounds=bounds, iter=1000, iprint=2)

    print "RES"
    print res

    return 0

if __name__ == "__main__":
    print "Running Tests"
    testSuite()

