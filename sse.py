import random
import math 
#Trial python stuff for SSE!

def diagonalUpdate( startConfig, bondList, operator, n, beta ):
    add = 0
    sub = 0
    time = 0
    while time < len(operator):
        if len(operator) < (4.0/3.0)*n:
            operator.extend( [0] * int(math.ceil(n/3.0)) )
        if( ( operator[ time ] > 0 ) and ( ( operator[ time ]%2 + 1 ) == 1 ) ): 
            #we have a diagonal operator at this slice 
            prob = ( 2.0 * float( len(operator) - n + 1 ) ) / ( float(len(bondList) - 1) * beta) 
            select = random.random()
            if ( (prob < 1.0 and select < prob) or prob > 1.0 ):
                operator[ time ] = 0 
                n           -= 1 
                sub += 1 
        else: #we have an identity at this slice 

            #have to propagate the spins to this point
            currentConfig = [ 0 ] * len( startConfig )
            currentConfig[ : ] = startConfig[ : ]
            for past in range(time):
                a = operator[ past ] % 2 + 1
                b = operator[ past ] / 2
                if( a == 2 ): #off diagonal
                    currentConfig[ bondList[b][0] ] = not (currentConfig[ bondList[b][0] ]) 
                    currentConfig[ bondList[b][1] ] = not (currentConfig[ bondList[b][1] ]) 
            #now we can see if a diagonal update is possible
            testBond = random.randrange( 1,len(bondList) )
            if( currentConfig[ bondList[ testBond ][ 0 ] ] ^ currentConfig[ bondList[ testBond ][ 1 ] ] ):
                prob   = beta * float(len(bondList) - 1) / (2.0 * float( len(operator) - n ) )
                select = random.random()
                if ( (prob < 1.0 and select < prob) or prob > 1.0 ):
                    operator[ time ] = 2 * testBond
                    n           += 1 
                    add += 1
        time += 1
    return n

def constructVertexList( operators, bondList, sites ): 
    vertexList  = [ -1 ]*(4 * len(operators))
    lastVertex  = [ -1 ] * sites 
    firstVertex = [ -1 ] * sites 
    for time,op in enumerate(operators):
        if (op > 0):
            v0    = 4 * time 
            b     = op/2
            sitei = bondList[ b ][ 0 ] 
            sitej = bondList[ b ][ 1 ] 
            v1    = lastVertex[ sitei ] 
            v2    = lastVertex[ sitej ]
            if( v1 > -1 ):
                vertexList[ v1 ] = v0
                vertexList[ v0 ] = v1
            else:
                firstVertex[ sitei ] = v0
            if( v2 > -1 ):
                vertexList[ v2 ] = v0 + 1
                vertexList[ v0 + 1 ] = v2
            else:
                firstVertex[ sitej ] = v0 + 1
            lastVertex[ sitei ] = v0 + 2
            lastVertex[ sitej ] = v0 + 3
    for site in range(sites):
        f = firstVertex[ site ]
        l = lastVertex[ site ]
        if( f != -1 ):
            vertexList[ f ] = l 
            vertexList[ l ] = f 
    return vertexList,firstVertex
 
def loopUpdate( operators, bondList, startConfig, sites ):

    legs = [1,0,3,2]
    Lists = constructVertexList( operators, bondList, sites )
    vertexList = Lists[ 0 ]
    firstVertex = Lists[ 1 ]
    while( (vertexList.count( -1 ) + vertexList.count( -2 ) ) < len(vertexList) ):
        #Do the loop-trace and update!
        flip = (random.random() > 0.5)
        loopStart = 0
        for v, vv in enumerate( vertexList ):
            if( vv > -1 ):
                loopStart = v
                break
        v = loopStart
        while( v < len(vertexList) ): 
            time = v/4
            a = operators[time]%2 + 1
            b = operators[time]/2
            vleg = 4*time + legs[ v%4 ] 
            vnext = vertexList[ vleg ]
            if flip: 
                if a == 2:
                    operators[ time ] = 2 * b 
                if a == 1:
                    operators[ time ] = 2 * b + 1 
                vertexList[ v ] = -2
                vertexList[ vleg ] = -2
            else: 
                vertexList[ v ] = -1
                vertexList[ vleg ] = -1
            
            v = vnext 
            if( v == loopStart ):
                break
    #update stored state
    for site in range(sites):
        v = firstVertex[ site ]
        if( v == -1 and (random.random() > 0.5 )):
            startConfig[ site ] = not( startConfig[ site ] )
        elif( vertexList[ v ] == -2 ):
            startConfig[ site ] = not( startConfig[ site ] )

# we will perform 1,000 Monte Carlo sweeps
# during each sweep, we do a diagonal update
# then a loop update, and occasionally measure


# a(p) - off diagonal vs diagonal, operatorList[ p ]%2 + 1
# b(p) - which bond index we are at

random.seed()
sites = 256 
#number of time slices in the simulation in total - can be adjusted
L = sites / 2 
bondList    = [ 0 ] * (2 * sites + 1)
startConfig = [ 0 ] * sites
beta = 16.0 
for i in range(1,sites + 1):
    bondList[ i ]         = ( i - 1, (i)%16 + 16*(( i - 1)/16) )
    bondList[ i + sites ] = ( i - 1, (i - 1 + 16)%sites )
    if( (i - 1)%2 ):
        startConfig[ i - 1 ] = True 
    else:
        startConfig[ i - 1 ] = False 

operators = [ 0 ] * L 
sweeps = 1000
energy = 0.0
n = 0 
for i in range(sweeps):

    n = diagonalUpdate( startConfig, bondList, operators, n, beta )
    loopUpdate( operators, bondList, startConfig, sites )
    energy += n/(beta*sites*sweeps) 
    print i, len(operators), n, energy

#operators.List = [4,0,9,13,6,0,0,4,13,0,9,14]
#vertexList = constructVertexList( operators, bondList, sites )

#for i,ii in enumerate(operators.List):
#    print vertexList[ 4*i ] + 1,vertexList[ 4*i + 1 ] + 1,vertexList[ 4*i + 2 ] + 1,vertexList[ 4*i + 3 ] + 1

