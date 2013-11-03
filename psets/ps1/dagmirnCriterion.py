# External Libraries
import numpy

def cournotEqmSelection(cost, demand):
    """
    ** Compute equilibrium to game with n players,
    ** stable under DAgMirN criterion (Gauss-Jacobi iteration)


    ** Parameters

    * cost - nx1 array of slope of each player's marginal cost
    * demand - 1x2 array of demand function parameters
    *        - qty(price) = demand[0] - demand[1] * price

    ** Output

    * qty - nx1 array of quantity produced by each firm

    """
    
    # 0. Number of Players
    N = len(cost)
    print N, "players"

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

    # 3. Run Gauss-Jacobi (Just iterations now)
    converge = numpy.zeros([N,1])
    iters = 100
    for i in xrange(100):

        for n in xrange(N):
            qtyN = (demand[0] - numpy.sum(qty) + qty[n]) / (2 + cost[n] * demand[1])
            print qtyN
            converge[n] = qtyN == qty[n]
            qty[n] = qtyN
            
        # Check if converged
        if sum(converge) == N:
            print "DONE"
            iter = i
            break
    
    price = demand[0] - demand[1] * sum(qty)

    print "Result after", i, "iterations"
    print qty
    print price

    return qty


def main():

    cost = numpy.array([[1,3]]).T

    demand = numpy.array([5,1])

    cournotEqmSelection(cost, demand)

    return 0

if __name__ == "__main__":
    main()
