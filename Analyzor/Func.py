import numpy as np
import math

def tbjcfit(xs,ys):
    xc,yc = xs.mean(),ys.mean()
    u,s,v = np.linalg.svd(np.array([xs-xc,ys-yc]).T)
    xi,yi = v[0]
    k = yi/(xi+1e-9*abs(xi)/xi)
    b = yc-k*xc
    return k,b

def GetLineInfo(p1,p2, L_thre = -5):
    (x1,y1),(x2,y2) = p1,p2
    if abs(x1-x2) <0.01: return 90,abs(y1-y2)
    
    L = math.hypot(x2-x1,y2-y1)

    if L<L_thre:return None, None

    return math.acos((x2-x1)/L)/math.pi*180,L
