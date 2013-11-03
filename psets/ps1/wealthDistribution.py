# External Libraries
import numpy



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

    # 3. Construct Markov Chain
    
    return 0

def main():

    a = 0.5
    g = 0.5
    p0 = 0.5
    p1 = 0.5
    N = 10

    wealthDistribution(a,g,p0,p1,N)

    return 0

if __name__ == "__main__":
    main()
