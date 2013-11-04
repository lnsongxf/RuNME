# External Libraries
import numpy

def cournotEqmSelection(cost, demand, iters=100, verbose=False):
    """
    ** Compute equilibrium to game with n players,
    ** stable under DAgMirN criterion (Gauss-Jacobi iteration)


    ** Parameters

    * cost - nx1 array of slope of each player's marginal cost
    * demand - 1x2 array of demand function parameters
    *        - qty(price) = demand[0] - demand[1] * price
    * iters  - maximum number of iterations

    ** Output

    * qty - nx1 array of quantity produced by each firm

    """
    
    # 0. Number of Players
    N = len(cost)
    print N, "players"

    print "Marginal Cost"
    print cost
    
    print "Demand Parameters"
    print demand

    # 1. Set up initial guess (qty=0 for all players)
    qty = numpy.zeros([N,1])
    print "Initial Guess"
    print qty

    # 2. Define best-response function
    # q_i = (d_0 - sum q_j) / (2 + c_i * d_1)
    # Can be written as Ax = b, where
    # A_ii = 2 + c_i * d_1,
    # A_ij = 1
    # b_i = d_0

    # 3. Run Gauss-Jacobi

    print "\nBeginning Gauss-Jacobi Iterations"
    print "================================="

    converge = numpy.zeros([N,1])
    converged = False

    for i in xrange(100):

        for n in xrange(N):
            # Best Response Function of Firm N
            # given that every other firm m produces qty[m]
            qtyN = (demand[0] - numpy.sum(qty) + qty[n]) / (2 + cost[n] * demand[1])

            if verbose:
                print qtyN

            converge[n] = qtyN == qty[n]
            qty[n] = qtyN
            
        # Check if converged
        if sum(converge) == N:
            print "Converged after", i, "iterations"
            print "=================================\n"

            converged = True
            break
    
    if converged == False:
        print "Stopped after", iters, "iterations"
        print "=================================\n"

    # 4. Calculate Price and Display Results
    price = demand[0] - demand[1] * sum(qty)

    print "Quantity"
    print qty

    print "Price"
    print price

    return (qty, price)

def main():

    #rcost = numpy.array([[1,3]]).T

    #rdemand = numpy.array([0.5,0.5])

    # Set number of players
    n = 2
    
    # Set max number of iterations
    iters = 100

    rcost = numpy.random.rand(n,1)
    rdemand = numpy.random.rand(n,1)

    cournotEqmSelection(rcost, rdemand, iters=iters)

    return 0

if __name__ == "__main__":
    main()
