
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
