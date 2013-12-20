# teststochdp.py

from stochdp import *
import numpy as np
import matplotlib.pyplot as ppt

## Sample runs of functions provided in
## module stochdp.py

##############################################
# Case 1.
# One period with bequest, no labor (wages=0),
# where u(c_1, b) = c_1^0.5 + 0.9 * b^0.5
##############################################

#####
# Note: the utility function supplied should
# take two arguments, x and k0
# where x is a k-length array of the controls
# and k0 is the single state variable
#
# for example, if we were solving a problem
# where an agent chooses savings and labor,
# we might let x[0] = k_{t+1} and x[1] = l
#####

# x - controls
# y - state variables
def utilityDefault(x, y):
	k1, l = x
	k0, w = y

	if l >= 1 or k0**0.33 + (0.9 * k0) - k1< 0:
		return -1

	return (k0**0.33 + (0.9 * k0) - k1 + (w * l))**0.5

def bequestValueDefault(x):
	return x**0.5

def productionDefault(k0, l):
	return k0 ** 0.5 + (1 - 0.1) * k0

def statewageDefault(state, wage):
	return wage

opts = {}

## Model Parameters

# Utility and Production Functions
opts['utility'] = utilityDefault
opts['bequest'] = bequestValueDefault
opts['production'] = productionDefault
opts['statewage'] = statewageDefault
opts['beta'] = 0.9

# Mean wages over time
# opts['wages'] = None

# Number of periods
opts['T'] = 1

# Number states describing
# change in wages
# opts['states'] = [0,1,2]
# opts['trans'] = mTrans

## Optimization Parameters

# Degree of Chebyshev Polynomial
# and number of nodes
opts['deg'] = 5.
opts['pts'] = 40.

# Range of points
opts['init'] = 0.1
opts['end'] = 10.

## Execute the solver
result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
print "CASE ONE COMPLETE"

########################################
# Case 2.
# One period with elastic labor
########################################

# Add labor to utility function
def utilityDefault(x, y):
	k1, l = x
	k0, w = y

	if l >= 1 or k0**0.33 + (0.9 * k0) - k1 + (w * l) < 0:
		return -1

	return (k0**0.33 + (0.9 * k0) - k1 + (w * l))**0.5 + (1. - l)**0.5

opts['wages'] = [5]
opts['utility'] = utilityDefault

## Execute the solver
result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
print "CASE TWO COMPLETE"

########################################
# Case 3.
# Ten periods with elastic labor
# deterministic varying wages 
########################################

## Define wages over time
# Rise in first third, constant for second third,
# and slowly declining in third third
opts['T'] = 10
wages = np.zeros(opts['T'])
period = int(np.floor(opts['T']/3.))
if period == 0:
	print "Too few periods - setting constant wages"
	wages[:] = 5
else:
	for i in xrange(period+1):
		wages[i] = i * (5. / period)# growing up to 5
	wages[period:2*period] = 5.
	for i in xrange(2*period,opts['T']):
		wages[i] = 5. - 2. * (i - 2.*period) / (opts['T'] - 2.*period)
opts['wages'] = wages

## Execute the solver
result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
print "CASE THREE COMPLETE"

########################################
# Case 4.
# Ten periods with elastic labor
# and wages which move through
# three states, high, mid, and low
# in high state, wages 10% above mean
# in low state wages 10% below mean
# movement between states described by
# transition matrix
########################################

mTrans = np.zeros((3,3))
mTrans[0,:] = [0.75,0.2,0.05]
mTrans[1,:] = [0.25,.5,.25]
mTrans[2,:] = [0.05,0.2,0.75]

# wage state
# computes the wage given the current state
def statewageDefault(state, wage):
	if state == 0:
		return wage * 0.9
	elif state == 1:
		return wage * 1
	elif state == 2:
		return wage * 1.1

opts['trans'] = mTrans
opts['states'] = [0,1,2]
opts['statewage'] = statewageDefault

## Execute the solver
result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
print "CASE FOUR COMPLETE"

########################################
# Case 5.
# Case Four, with timeline extended to
# 60 periods
########################################

opts['T'] = 60
wages = np.zeros(opts['T'])
period = int(np.floor(opts['T']/3.))
for i in xrange(period+1):
	wages[i] = i * (5. / period)# growing up to 5
wages[period:2*period] = 5.
for i in xrange(2*period,opts['T']):
	wages[i] = 5. - 2. * (i - 2.*period) / (opts['T'] - 2.*period)
opts['wages'] = wages

## Execute the solver
result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
print "CASE FIVE COMPLETE"