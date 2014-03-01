# teststochdp.py

from stochdp import *
import numpy as np
#import matplotlib.pyplot as ppt
import traceback

## Sample runs of functions provided in
## module stochdp.py

##############################################
# Setup
##############################################

# Check Numpy Version
try:
	assert(np.version.version == '1.8.0')
except:
	print "WARNING: NUMPY VERSION != 1.8.0"

# Set active cases
case1 = True
case1b = True
case2 = True
case3 = True
case4 = False
case5 = False

##############################################

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

def productionDefault(k0):
	return k0 ** (1./3.) + (1 - 0.1) * k0

# x - controls
# y - state variables
def utilityDefault(x, y):
	k1= x
	k0, w = y

	if productionDefault(k0) - k1< 0:
		return -1

	return (productionDefault(k0) - k1)**0.5

def bequestValueDefault(x):
	return x**0.5

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

# bounds for variables
opts['bounds'] = lambda s: ((0, opts['production'](s[0])),)
opts['x0'] = lambda s: np.array(sum(opts['bounds'](s)[0]) / 2.)

# Number of periods
opts['T'] = 1

## Optimization Parameters

# Degree of Chebyshev Polynomial
# and number of nodes
opts['deg'] = 10
opts['pts'] = 11

# Range of points
opts['init'] = 0.1
opts['end'] = 10.

## Execute the solver
if case1:
	try:
		result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
		print "CASE ONE COMPLETE"
	except:
		print "CASE ONE FAILED"
		print traceback.format_exc()


## Check policy function (from pen & paper calc)

# policy = k_1 = bequest
#policy = result['policy'][0][0]

# grid = k0
grid = nodes(opts['init'], opts['end'], opts['pts'])

# compute Y

#y = grid**(1./3.) + (0.9)*grid

# compute multiple

#b2y = 1. + (1. / 0.81)

#error = [y[i] - (b2y * policy)[i] for i in range(len(y))]
#error =  np.array([error[i][0] for i in range(len(error))])
#error = np.sum((error * error)**0.5)

#print error

########################################
# Case 1b.
# Like case one, but with
# modified utility function
########################################

def utilityLog(x, y):
	k1= x
	k0, w = y

	if k0**0.33 + (0.9 * k0) - k1 < 0:
		return -1

	return np.log(k0**0.33 + (0.9 * k0) - k1)

def bequestValueLog(x):
	return np.log(x)

opts['utility'] = utilityLog
opts['bequest'] = bequestValueLog

# Execute the solver
if case1b:
	try:
		result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
		print "CASE ONE PART B COMPLETE"
	except:
		print "CASE ONE PART B FAILED"
		print traceback.format_exc()

## Check policy function (from pen & paper calc)

# policy = k_1 = bequest
policy = result['policy'][0][0]

# grid = k0
grid = nodes(opts['init'], opts['end'], opts['pts'])

# compute Y

y = grid**(1./3.) + (0.9)*grid

# compute multiple

b2y = 1. + (1. / 0.9)

error = [y[i] - (b2y * policy)[i] for i in range(len(y))]
error =  np.array([error[i][0] for i in range(len(error))])

error = np.sum((error * error)**0.5)

#print "Error:", error

########################################
# Case 2.
# One period with elastic labor
########################################

# Add labor to production function
def productionDefault(k0):
	return (k0 ** .5) + 0.9 * k0
opts['production'] = productionDefault

# New bounds including labor as control
kbounds = lambda s: (0, opts['production'](s[0]) + s[1])
lbounds = lambda s: (0,1)
opts['bounds'] = lambda s: (kbounds(s), lbounds(s))
opts['x0'] = lambda s: np.array([sum(kbounds(s)) / 2.,sum(lbounds(s)) / 2.])

# Add labor to utility function
def utilityDefault(x, y):
	k1, l = x
	k0, w = y

	if l >= 1 or productionDefault(k0) + w * l - k1 < 0:
		return -1

	return np.log(productionDefault(k0) + w * l - k1) + np.log(1. - l)

def bequestValueDefault(x):
	return np.log(x)

opts['wages'] = [5]
opts['utility'] = utilityDefault
opts['bequest'] = bequestValueDefault

## Execute the solver
if case2:
	try:
		result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
		print "CASE TWO COMPLETE"
	except:
		print "CASE TWO FAILED"
		print traceback.format_exc()

# Check policy function (from pen & paper calc)
k1 = result['policy'][0][0][:,0]
l = result['policy'][0][0][:,1]

k0 = grid

e1 = 5*(1 - l) - (k0**.5 + .9 * k0 + 5*l - k1)
e2 = (10./9.) * k1 - (k0**.5 + .9 * k0 + 5*l - k1)

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
for i in xrange(period+1):
	wages[i] = 1 + i * (5. / period)# growing up to 5
wages[period:2*period] = 6.
for i in xrange(2*period,opts['T']):
	wages[i] = 1 + 5. - 2. * (i - 2.*period) / (opts['T'] - 2.*period)
opts['wages'] = wages

## Execute the solver
if case3:
	try:
		result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
		print "CASE THREE COMPLETE"
	except:
		print "CASE THREE FAILED"
		print traceback.format_exc()

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
if case4:
	try:
		result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
		print "CASE FOUR COMPLETE"
	except:
		print "CASE FOUR FAILED"
		print traceback.format_exc()

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
if case5:
	try:
		result = execute(deg = opts['deg'], pts=opts['pts'], opts=opts)
		print "CASE FIVE COMPLETE"
	except:
		print "CASE FIVE FAILED"
		print traceback.format_exc()
