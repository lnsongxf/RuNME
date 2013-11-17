# External Libraries
import numpy
import time
import math

def cournotEqmSelection(cost, demand, iters=100, tol=1e-10, algo="gj"):
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
    #qty = numpy.zeros(N)
    qty = numpy.repeat(-10e9,N)
    print "Initial Guess"
    print qty

    # 2. Define best-response function
    # q_i = (d_0 - sum q_j) / (2 + c_i * d_1)
    # Can be written as Ax = b, where
    # A_ii = 2 + c_i * d_1,
    # A_ij = 1
    # b_i = d_0

    # 3. Run Gauss-Jacobi

    print "\nBeginning Iterations"
    print "================================="

    converge = numpy.zeros([N,1])
    converged = False

    for i in xrange(100):

        Q = numpy.sum(math.e ** qty)

        if algo == "gj":
            # Gauss-Jacobi Step
            
            qtyN = map(math.log,(demand[0] - Q + math.e ** qty) / (2 + cost * demand[1]))
            #qtyN = (qtyN > numpy.zeros(N)) * qtyN

            converge = abs(qtyN - qty)

            qty = numpy.copy(qtyN)

        if algo == "gs":
            # Gauss-Seidel Step
            for n in xrange(N):
                # Best Response Function of Firm N
                # given that every other firm m produces qty[m]

                qtyN = math.log((demand[0] - Q + math.e ** qty[n]) / (2 + cost[n] * demand[1]))

                converge[n] = abs(qtyN - qty[n])

                qty[n] = qtyN
        
        # Check if converged
        if max(converge) < tol:
            print "CONVERGED: Converged after", i, "iterations"
            print converge
            print "=================================\n"

            converged = True
            break

    if converged == False:
        print "REACHED MAX ITERATIONS: Stopped after", iters, "iterations"
        print "=================================\n"

    # 4. Calculate Price and Display Results
    price = (demand[0] - sum(math.e ** qty)) / demand[1]

    print "Quantity"
    print math.e ** qty

    print "Price"
    print price

    return (qty, price)

def main():

    # Set number of players
    n = 5
    # Set max number of iterations
    iters = 100

    rcost = 10 * numpy.random.rand(n)

    rdemand = numpy.array([100 * numpy.random.rand(), 10 * numpy.random.rand()])


    cournotEqmSelection(rcost, rdemand, iters=iters, algo='gs')

    return 0

if __name__ == "__main__":
    main()
