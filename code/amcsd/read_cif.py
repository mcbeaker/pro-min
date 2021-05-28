
from numpy import cos,sin,sqrt,matrix,pi,dot,cross,array


a2r = pi/180

def fract_to_cart(a,b,c,u,v,w):
    A = array([a,0,0])
    B = array([0,b,0])
    C = array([0,0,c])
    r = array([u,0,0])*A + array([0,v,0])*B + array([0,0,w])*C
    # print(r)
    return r

def xDist(x1,x2):
    return x2-x1
def yDist(y1,y2):
    return y2-y1
def zDist(z1,z2):
    return z2-z1


def xyz_dist(x1,y1,z1,x2,y2,z2):

    xDist = (x2-x1)*(x2-x1)
    yDist = (y2-y1)*(y2-y1)
    zDist = (z2-z1)*(z2-z1)
    sumDist = xDist + yDist + zDist
    d = sqrt(sumDist)
    return d


def diff_read_cif(filename):

    from diffpy import structure 
    # print(help(Structure))

    return(structure.loadStructure(filename))