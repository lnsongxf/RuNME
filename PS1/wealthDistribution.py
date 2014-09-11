# External Libraries
import numpy
import matplotlib.pyplot as ppt

import sys

def iterate(wage, wealth, p, p0, p1, alpha, gamma):

    nextpd = []

    if wage == 0:
        nextpd.append((0, (1 - gamma) * wealth, p0))
        nextpd.append((1, (1 - gamma) * wealth + (1 - alpha), 1 - p0))
    if wage == 1:
        nextpd.append((0, (1 - gamma) * wealth, 1 - p1))
        nextpd.append((1, (1 - gamma) * wealth + (1 - alpha), p1))

    return nextpd

def wealthMovement(grid, prob, alpha, gamma, pi0, pi1):

    # Get N
    N = len(grid)

    newProb = numpy.zeros([N,2])

    # For w = 0
    # For each initial grid entry
    for i in xrange(N):
        nxt = iterate(0, grid[i], prob[i,0], pi0, pi1, alpha, gamma)

        n1 = numpy.argmax(grid >= nxt[0][1]) - 1 # grid slot for w=0 entry
        n2 = numpy.argmax(grid >= nxt[1][1]) - 1 # grid slot for w=1 entry

        newProb[n1,0] += prob[i,0] * nxt[0][2]
        newProb[n2,1] += prob[i,0] * nxt[1][2]

    # For w = 1
    for i in xrange(N):
        nxt = iterate(1, grid[i], prob[i,1], pi0, pi1, alpha, gamma)

        n1 = numpy.argmax(grid >= nxt[0][1]) # grid slot for w=0 entry
        n2 = numpy.argmax(grid >= nxt[1][1]) # grid slot for w=1 entry

        newProb[n1,0] += prob[i,1] * nxt[0][2]
        newProb[n2,1] += prob[i,1] * nxt[1][2]

    return newProb

def wealthDistribution(alpha, gamma, pi0, pi1, N, iters=300):
    """
    ** Computes ergodic distribution

    * Every period, A_t = (1-alpha) * w + (1-gamma) * A_{t-1}

    ** Parameters

    alpha - proportion of wage consumed, in [0,1]
    gamma - proportion of assets consumed, in [0,1]

    p0    - P(w_{t+1} = 0 | w_t = 0)
    p1    - P(w_{t+1} = 1 | w_t = 1)

    N     - Number of points on grid

    ** Output

    prob[:,0] - vector of probabilities for w = 0
    prob[:,1] - vector of probabilities for w = 1

    """

    # 1. Compute maximum wealth, Amax
    amax = (1 - alpha) / gamma
    #print "Amax = ", amax

    # 2. Create grid on [0, amax]
    resolution = float(amax) / N
    grid = numpy.linspace(0,amax,num=N)
    prob = numpy.zeros([N,2])

    # 3. Construct Markov Chain

    # Initial distribution: wage = 1, assets at 0
    prob[0,0] = 1.0

    # Iterate through process
    for i in xrange(iters):
        prob = wealthMovement(grid, prob, alpha, gamma, pi0, pi1)

    # Plot W = 0
    ppt.plot(grid,prob[:,0])

    # Plot W = 1
    ppt.plot(grid,prob[:,1])

    # Show Graph
    #ppt.show()

    # Save Graph
    ppt.savefig(str(alpha) + "_" + str(gamma) + "_" + str(pi0) + "_" + str(pi1) + "_" + str(N) + "_" + str(iters) + "_plot.png")

    # Clear Graph
    ppt.close()

    return prob

def main():

    a_ = [0.1,0.5,0.9]
    g_ = [0.1,0.5,0.9]
    p0_ = [0.1,0.5,0.9]
    p1_ = [0.1,0.5,0.9]
    N_ = [100,500,1000]
    iters = [10,100,1000]

    # Test for each of parameter combination

    wealthDistribution(0.9, 0.1, 0.1, 0.5, 10000, 100)
    sys.exit()

    for a in a_:
        for g in g_:
            for p0 in p0_:
                for p1 in p1_:
                    for N in N_:

                        print [a,g,p0,p1,N]
                        try:
                            for i in iters:
                                wealthDistribution(a,g,p0,p1,N,iters=i)
                        except:
                            continue

    return 0

if __name__ == "__main__":
    main()
