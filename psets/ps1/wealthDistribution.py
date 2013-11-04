# External Libraries
import numpy
import matplotlib.pyplot as ppt

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

        #if prob[i,0] != 0:
        #    print i
        #    print nxt


        n1 = numpy.argmax(grid >= nxt[0][1]) - 1 # grid slot for w=0 entry
        n2 = numpy.argmax(grid >= nxt[1][1]) - 1 # grid slot for w=1 entry

        #if prob[i,0] != 0:
        #    print "P",prob[i,0]
        #    print nxt[0][2]
        #    print nxt[1][2]
            
        newProb[n1,0] += prob[i,0] * nxt[0][2]
        newProb[n2,1] += prob[i,0] * nxt[1][2]

    # For w = 1
    for i in xrange(N):
        nxt = iterate(1, grid[i], prob[i,1], pi0, pi1, alpha, gamma)

        #if prob[i,1] != 0:
        #    print i
        #    print nxt


        n1 = numpy.argmax(grid >= nxt[0][1]) # grid slot for w=0 entry
        n2 = numpy.argmax(grid >= nxt[1][1]) # grid slot for w=1 entry
        """
        if prob[i,0] != 0:
            print "P",prob[i,0]
            print nxt[0][2]
            print nxt[1][2]
        """ 
        newProb[n1,0] += prob[i,1] * nxt[0][2]
        newProb[n2,1] += prob[i,1] * nxt[1][2]


    #print "Newprob:\n", newProb

    return newProb

def wealthDistribution(alpha, gamma, pi0, pi1, N):
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

    dist0 - vector of probabilities for w = 0
    dist1 - vector of probabilities for w = 1

    """

    # 1. Compute maximum wealth, Amax
    amax = (1 - alpha) / gamma

    print "Amax = ", amax

    # 2. Create grid on [0, amax]
    resolution = float(amax) / N
    grid = numpy.linspace(0,amax,num=N)
    prob = numpy.zeros([N,2])




    # 3. Construct Markov Chain

    # Initial distribution: wage = 1, assets at 0
    prob[0,0] = 1.0

    #print grid
    #print prob
    
    for i in xrange(15):
        prob = wealthMovement(grid, prob, alpha, gamma, pi0, pi1)

    #print prob

    # Plot W = 0
    #ppt.plot(grid,prob[:,0])

    # Plot W = 1
    #ppt.plot(grid,prob[:,1])
    #ppt.show()

    print "SUM: ", sum(prob)
    
    return 0

def main():

    a = 0.1
    g = 0.1
    p0 = 0.5
    p1 = 0.5
    N = 10000

    wealthDistribution(a,g,p0,p1,N)

    return 0

if __name__ == "__main__":
    main()
