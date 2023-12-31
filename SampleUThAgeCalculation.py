import sys
from itertools import starmap
import math
from math import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# some constants
L_230 = 9.17055e-06
L_234 =2.82203e-06
L_DELTA = L_230 - L_234
L_RATIO = L_230 / (L_230 - L_234)
# Newthon Raphson method
EPSILON = 1E-9
START_POINT = 0.1
MAX_TRIES_UNSTABLE_SOLUTION = 10
# output
OUTPUT_FILENAME= "new_ages.txt"

# R230_234=0.991
# R234_238=1.012
R230_234=0.892
R234_238=1.003

# Newton-Raphson
def NewtonRaphson(x0, Th230U234, U234U238, epsilon):

    # start point
    x = x0
    newx = x
    stopcriteria = False
    solutionDirection = int(0)
    switchIteration = 0

    while( not stopcriteria ):
        # Check a non zero derivative value
        # it's never going to happen with this function but still
        fder = DerAgesEq(x, U234U238)
        while ( math.isclose(fder, 0.0, abs_tol=1E-18) ) :
            x += epsilon
            fder = DerAgesEq(x, U234U238)
        newx = x - ( AgesEq(x, Th230U234, U234U238) / fder )
        stopcriteria = True if (abs(newx - x)/abs(newx) < epsilon) else False

        # Prepare next iteration
        # If the solution is changing direction it became unstable
        #  and may end up in an overflow error. Try with a different X0
        #  or otherwise declare it without solution (NaN)
        #print("--> %s | %s"%(newx,x))
        if ( solutionDirection == 0 ):
            if ( newx - x < 0.  ) : solutionDirection = int(-1)
            if ( newx - x >= 0. ) : solutionDirection = int(1)
        else:
            if ( (newx - x) * solutionDirection < -0.01*min(newx,x) ) :     # changed direction !
                newxsave = newx
                x0 = x0*10. if solutionDirection==1 else x0/10.;            # new X0 depending on the tendency
                solutionDirection = int(0)                                  # determine de direction again
                newx = x0
                switchIteration +=1
                print("[WARN] changing direction, unstable solution. Try again with x0 = %s [retry %s/%s]" % (x0, switchIteration, MAX_TRIES_UNSTABLE_SOLUTION))
                print("       --> newx = %s | x=%s"%(newxsave,x))
                if ( switchIteration >= MAX_TRIES_UNSTABLE_SOLUTION) : # give up condition
                    print("[WARN] Bad solution:")
                    print("       Th230/U234A = %s, U234/U238A = %s"%(Th230U234,U234U238))
                    print("       giving up on this entry. Returning NaN")
                    newx = nan
                    break
        x = newx # prepare next iteration

    return newx


# Age equation
def AgesEq(t, Th230U234, U234U238):

    f = ((1. - np.exp(-1. * L_230 * t)) / U234U238) - Th230U234
    f += (L_RATIO * (1. - (1. / U234U238))) * (1. - np.exp(-1. * L_DELTA * t))
    return f

# First order derivative of Age Equation. Analytical.
def DerAgesEq(t, U234U238):
    f = ((L_230 * np.exp(-1. * L_230 * t)) / U234U238)
    f += (L_RATIO * (1. - (1. / U234U238)) * (L_DELTA * np.exp(-1. * L_DELTA * t)))
    return f


#########################################
#########################################

#print( "Test AgesEq    = %s"%AgesEq(2000, 0.95, 1.2) ) 
#print( "Test DerAgesEq = %s"%DerAgesEq(2000, 1.2) )

res = NewtonRaphson(START_POINT, R230_234, R234_238, EPSILON)
print("Solution = %s"%res)

