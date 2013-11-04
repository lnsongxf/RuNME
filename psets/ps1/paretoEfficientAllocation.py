# External Libraries

import numpy
import scipy.optimize
import time

def foo(x):
    print x
    time.sleep(0.5)
    return x

def efficientAlloc(n, m, nu, a, e, l):

    """
    ** Computes Solution to social planner's problem
    ** with m goods and n agents

    ** Problem: Max Sum lambda_i u_i(x)
    ** Nonlinear Optimization Problem

    ** Parameters
    * nu - n x m matrix, < 0
    * a  - n x m matrix, > 0
    * e  - n x m matrix, > 0 (initial endowment)
    * l  - n x 1 matrix, > 0 (social weights)

    ** Output
    * allocation - n x m matrix
    *            - allocation[i,j] = agent i's allocation of good j

    """

    # 1. Set Objective Function
    #obj = lambda x : print(x)
    obj = lambda x : -numpy.sum(a * (x ** (nu + 1)) / (nu + 1))
    #obj = lambda x : -numpy.sum(a * (numpy.reshape(foo(x),(n,m)) ** (1 + nu)) / (1 + nu))
    #obj = lambda x : -numpy.sum(a * ((numpy.reshape(x,(n,m)) ** (1 + nu)) / (1 + nu)))
    #obj = lambda x: -numpy.sum(foo(x) ** numpy.random.rand(25))

    # 2. Set Initial Guess = Initial Endowment
    x0 = numpy.zeros(n*m)
    numpy.copyto(x0, e)

    print x0
    print e
    print nu

    # 3. Choose Solver
    method = 'SLSQP'

    # 4. Set Equality Constraints and Bounds
    consType = "eq"
    consFun = lambda x : numpy.sum(numpy.reshape(x,(n,m)) - e, 0) #== numpy.zeros([m,1])
    consFun = lambda x : numpy.sum(numpy.reshape(x - e, (n,m)))

    cons = ({'type': consType, 'fun': consFun})

    bounds = [(0,None)] * (n*m)

    # 5. Perform Optimization
    res = scipy.optimize.minimize(obj, x0, method=method, bounds=bounds)#constraints=cons)#, bounds=bounds)

    print res

    return 0




def main():

    # Number of Agents
    n = 5

    # Number of Goods
    m = 1
    
    ## Generate Random Parameters

    # param(i*m + j) = param[i,j]

    # Generate Preferences
    nu = -numpy.random.rand(n*m)
    a = numpy.random.rand(n*m)

    # Generate Endowment
    e = numpy.random.rand(n*m)

    # Generate Social Weights
    l = numpy.random.rand(n)

    efficientAlloc(n, m, nu, a, e, l)

    return 0

if __name__ == "__main__":
    main()
