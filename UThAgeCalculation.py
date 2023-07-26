import sys
from itertools import starmap
import math
from math import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#from myNewtonRaphson import NewtonRaphson

# some constants
L_230 = 9.15771e-06
L_234 =2.82341e-06
L_DELTA = L_230 - L_234
L_RATIO = L_230 / (L_230 - L_234)
# Newthon Raphson method
EPSILON = 1E-9
START_POINT = 0.1
MAX_TRIES_UNSTABLE_SOLUTION = 10
# output
OUTPUT_FILENAME= "new_ages.txt"

# Newton-Raphson
def NewtonRaphson(x0, thu, uu, epsilon):

    # start point
    x = x0
    newx = x
    stopcriteria = False
    solutionDirection = int(0)
    switchIteration = 0

    while( not stopcriteria ):
        # Check a non zero derivative value
        # it's never going to happen with this function but still
        fder = DerAgesEq(x, uu)
        while ( math.isclose(fder, 0.0, abs_tol=1E-18) ) :
            x += epsilon
            fder = DerAgesEq(x, uu)
        newx = x - ( AgesEq(x, thu, uu) / fder )
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
                    print("       Th230/U238A = %s, U234/U238A = %s"%(thu,uu))
                    print("       giving up on this entry. Returning NaN")
                    newx = nan
                    break
        x = newx # prepare next iteration

    return newx


# Age equation
def AgesEq(t, thu, uu):
    f  =   np.double(1. - thu - np.exp(-1. * L_230 * t))
    f +=   np.double(L_RATIO * (uu - 1) * (1. - np.exp(-1. * (L_DELTA) * t)))
    return f

# First order derivative of Age Equation. Analytical.
def DerAgesEq(t, uu):
    f =  np.double(L_230 * np.exp(-1. * L_230 * t))
    f += np.double(L_RATIO * (uu - 1) * L_DELTA * np.exp(-1. * L_DELTA * t))
    #print("--> %s | %s"%(t,f))
    return f

# Check that files have the same structure
def CheckDataConsistency(data1, data2, data3):
    if ( len(data1) != len(data2) or len(data2) != len(data3) ) : return False
    return True

# check that values to be gathered in text files are valid floats
def isFloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def GetData(fn):

    fileinput = open(fn,"r")
    data = []
    cntr_first_line = -1
    linecntr = -1
    for line in fileinput:
        linecntr += 1
        ls1 = line.split('\t')
        #data.append( ls1 )
        if ( cntr_first_line < 0 ) :
            cntr_first_line = len(ls1)
        else :
            # In case data a line is found with a different number of entries
            if ( len(ls1) != cntr_first_line ) : return 0, []

        # if the length is correct convert to float and append
        for val in ls1:
            if ( isFloat( val ) ) :
                data.append( float(val) )
            else :
                print("[FATAL] Value can not be converted to float")
                print("        Working on file %s, line %s"%(fn,linecntr))

    print("[INFO] File : %s --> %s lines read with %s entries each"%(fileinput.name, linecntr+1, cntr_first_line))

    # Convert float vector

    return cntr_first_line, data # number of columns found and the data

##############################
# Extracting parameters
n_columnsThU, dataThU = GetData("data/exampleCaveModel/Th230_U238A.txt")
n_columnsUU, dataUU = GetData("data/exampleCaveModel/U234_U238A.txt")
# ca : Control Ages
n_columnsca, data_ca = GetData("data/exampleCaveModel/control_ages.txt")

if ( not CheckDataConsistency(dataThU, dataUU, data_ca) or (n_columnsThU != n_columnsUU) or (n_columnsUU != n_columnsca) ) :
    print("[ERRO] Input data is not formatted as expected. Different number of entries or matrix dimensions ?. Giving up.")
    sys.exit(1)
else :
    print("[INFO] Data seems OK. Continue ...")

# i am sure now all three vectors have the same size
samplesize = len(data_ca)

##############################
# Find new ages using
# Newton-Raphson method.

print("[INFO] Start finding %s roots using Newton-Raphson method."%samplesize)

newControlAges = []
max = samplesize
for itr in range(0, max):
    #print("\nStart -> %s, %s, %s"%(startpoint, float(dataThU[itr]), float(dataUU[itr])))
    newControlAges.append(NewtonRaphson(START_POINT, float(dataThU[itr]), float(dataUU[itr]), EPSILON))
    #print(dataca[itr],newControlAges[itr])

print("[INFO] done with solutions")

# output file
print("[INFO] writing to file \"%s\" with %s columns" % (OUTPUT_FILENAME, n_columnsca))
outfile  = open(OUTPUT_FILENAME, 'w')
colcntr = 0
for itr in newControlAges:
    outfile.write("%s"%itr)
    if ( colcntr == n_columnsca - 1 ) :
        outfile.write("\n")
        colcntr = 0
    else :
        outfile.write('\t')
        colcntr += 1

outfile.close()


##############################
# Plot
dataca_subsample = []
for itr in range(0,max):
    dataca_subsample.append(data_ca[itr])

ages = plt.plot(dataca_subsample, newControlAges)
plt.setp(ages, color='r',linestyle = 'None', marker="o", markersize=1)
plt.xlabel('Control ages')
plt.ylabel('New calculated ages')
plt.show()

print("[INFO] Done ... besito :) !")
