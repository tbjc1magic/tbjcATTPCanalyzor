import numpy as np
import pylab as plt
import cv2
import seaborn as sns
import matplotlib.path as mplPath
import matplotlib.patches as patches
import math

def AveDist(x,y,k,b):
    return np.sum(np.abs(k*x+b-y)/(k*k+1))/x.shape[0]

def r2(x,y,k,b):
    return 1-np.sum(((x*k+b)-y)**2)/np.sum((y-np.mean(y))**2)

def VertexPos_(fits,y0):
    fits = [f for f in fits if f[0]>60 and f[-1]<2.5]

    num_lines = len(fits)

    if num_lines == 0: raise Exception('This is probably an unreacted beam')
    fits = fits[:1]+[f for f in fits[1:] if abs(f[1]-fits[0][1])/(abs(fits[0][1])+1e-9)>0.5]
    k1,b1 = fits[0][1:3]

    num_lines = len(fits)

    if num_lines==1:
        return (y0-b1)/k1,y0
    else:
        k2,b2 = fits[1][1:3]
        return (b2-b1)/(k1-k2), (k2*b1-k1*b2)/(k2-k1)

def VertexPos(image_,ps):
    ys,xs = np.where(image_)
    results = []

    def Error(xs,ys,xc,yc,ps):

        xc,yc = float(xc),float(yc)
        ps = ps.astype(np.float)

        def GetDist(xs,ys,xm,ym):
            Dist = np.ones(xs.shape[0])*1e5
            
            if abs(xc-xm)<0.1: 
                idx = np.where((ys-yc)*(ym-yc)>=0)
                Dist[idx] = np.abs(xs[(ys-yc)*(ym-yc)>=0]-xc)
            else:
                k = (yc-ym)/(xc-xm)
                b = (xc*ym-yc*xm)/(xc-xm)

                Di = xs-xc+k*ys-k*yc
                Di_ = xm-xc+k*ym-k*yc

                idx = np.where(Di*Di_>=0)

                Dist[idx] = np.sqrt(np.power(xs[idx]*k+b-ys[idx],2)/(k*k+1))

            return Dist

        Dist = [GetDist(xs,ys,xm,ym) for xm,ym in ps]

        return np.sum(np.min(np.stack(Dist),axis=0))


    path = mplPath.Path(ps[[0,1,2,0]])

    x = range(0,500)
    y = range(0,300)

    xv, yv = np.meshgrid(x, y)
    xv, yv = xv.reshape(-1), yv.reshape(-1)

    pix = np.stack([xv,yv]).T
    mask = path.contains_points(pix)

    res = []
    for x_,y_ in pix[mask]:
        err = Error(xs,ys,x_,y_,ps)
        res.append((x_,y_,err))

    return sorted(res,key=lambda x:x[2])[0][:-1]    
    
    
def tbjcfit(xs,ys):
    xc,yc = xs.mean(),ys.mean()
    u,s,v = np.linalg.svd(np.array([xs-xc,ys-yc]).T)
    xi,yi = v[0]
    k = yi/(xi+1e-9*abs(xi)/xi)
    b = yc-k*xc
    return k,b

def GetFit(image_, part_thresh=60, err_thresh =1.2,spread_thresh=6):

    if np.sum(image_>10) > part_thresh:
        ys,xs = np.where(image_)

        try:

            k,b = tbjcfit(xs,ys)

            if AveDist(xs,ys,k,b)<err_thresh and (np.std(xs)>spread_thresh or np.std(ys)>spread_thresh):

                return (np.sum(image_>10),k,b,AveDist(xs,ys,k,b))
            else:
                return (0,0,0,100)
        except:
            return (0,0,0,100)

    else:
        return (0,0,0,100)

def FilterBackground(image):
    hull = convexHull(image, debug_mode = False)

    ys,xs = np.where(image)
    path1 = mplPath.Path(hull[:,0,:])
    mask = path1.contains_points(zip(xs,ys))
    image_ = np.zeros(image.shape,dtype=np.uint8)
    image_[ys[mask],xs[mask]] = 255

    return image_


def GetEventPositions(pic,debug_mode=0):
    pic_ = np.copy(pic)
    points = TipFinder(pic_,debug_mode)
    xc,yc = VertexPos(pic_,points)
    
    return points, (xc,yc)

def GetEventPositions_(pic,debug_mode=0, center_width = 12, quadrant_thresh=100,
        center_thresh=300, err_thresh =12, spread_thresh=6 ):

    pic_ = np.copy(pic)

    points = TipFinder(pic_,debug_mode)
    xv,yv = points[points[:,0].argsort()][-1]

    #### get fit from each of the quadrant ###
    fits = []
    for idx in range(2):
        image_ = np.copy(pic_)
        if idx:
            image_[int(yv)-center_width:,:] = 0
        else:
            image_[:int(yv)+center_width,:] = 0

        fits.append(GetFit(cv2.erode(image_, None, iterations=1),
                    quadrant_thresh,err_thresh,spread_thresh))

    ##### get the fitting from the center line ####
    fit = fits[0] if fits[0][0]>quadrant_thresh and fits[0][-1]<err_thresh else fits[1]
    if fit[0]>quadrant_thresh and fit[-1]<err_thresh:
        y_,x_ = (150-fit[2])/(fit[1]+1e-9*(abs(fit[1])<1e-9)),150
        image_ = np.copy(pic_)
        image_[:,:int(y_)+10] = 0
        fits.append(GetFit(cv2.erode(image_, None, iterations=1),
                           center_thresh,err_thresh,spread_thresh))

    else:
        fits.append((0,0,0,100))

    ### process the fits according to each of the scenario
    xc,yc = VertexPos(fits,yv)

    if debug_mode:
        for s_,k_,b_,err_ in fits:
            if s_>60 and err_<3: plt.plot(np.linspace(0,600,2),b_+k_*np.linspace(0,600,2))
        plt.scatter([xc],[yc],marker='o',s=300,c='g')

    return points, (xc,yc)

def GetLineInfo(p1,p2, L_thre = -5):
    (x1,y1),(x2,y2) = p1,p2
    L = math.hypot(x2-x1,y2-y1)

    if L<L_thre:return None, None

    return math.acos((x2-x1)/L)/math.pi*180,L

def GetEventInfo(points,p0):
    points, p0 = np.array(points), np.array(p0)
    return [GetLineInfo(p,p0) for p in points], math.hypot(*(points[-1]-p0))

def Distance(contours,n1,n2):

    c1,c2 = contours[n1], contours[n2]

    c1,c2 = c1[:,0,:],c2[:,0,:]

    m1 = np.repeat([c1],c2.shape[0],axis=0)
    m2 = np.repeat([c2],c1.shape[0],axis=0)
    t2 = np.transpose(m2,axes=(1,0,2))
    diff = np.sqrt(np.sum(np.power(m1-t2,2),axis=2))

    return np.min(diff)

def Groups(contours):
    r,a = [],[]
    ## breadth first search
    pool = range(len(contours))
    area = 0
    while pool:
        seen = set([pool[0]])
        ans = set([pool[0]])
        sea = set([pool[0]])
        area = cv2.contourArea(contours[pool[0]])

        while sea:
            nsea = set([])

            for c in sea:
                for n in pool:
                    if n not in seen and Distance(contours,c,n)<40:
                        #print Distance(c,n),Distance(c,n)<8
                        area += cv2.contourArea(contours[n])
                        ans.add(n)
                        nsea.add(n)
                        seen.add(n)
            sea = nsea

        pool = [_ for _ in pool if _ not in ans]

        #print pool
        a.append(area)
        r.append(ans)

    r = [[contours[__] for __ in _] for _ in r]

    return zip(*[(x,y) for x,y in zip(r,a) if y>500])

def convexHull(thresh, debug_mode = 0):

    thresh = np.copy(thresh)
    if debug_mode: sns.heatmap(thresh[::-1], xticklabels=30, yticklabels=30)

    m1, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    gpc,gpa = Groups(contours)
    gpc = [np.concatenate(g,axis=0) for g in gpc]

    if debug_mode: plt.scatter(gpc[0][:,0,0],gpc[0][:,0,1])

    cnt = gpc[np.argmax(gpa)]
    hull = cv2.convexHull(cnt)

    if debug_mode: plt.plot(hull[:,0,0],hull[:,0,1],c='r')

    return hull#[:,0,:]

def MaxEnclosedTriangle(hull):
    def area(*i):
        a,b,c= hull[i,0]
        return np.abs(np.cross(b-c,c-a))

    A = 0;B = 1; C = 2; n=hull.shape[0]
    bA= A; bB= B; bC= C #The "best" triple of points
    while True: #loop A
        while True: #loop B
            while area(A, B, C) <= area(A, B, (C+1)%n): #loop C
                C = (C+1)%n
            if area(A, B, C) <= area(A, (B+1)%n, C):
                B = (B+1)%n
                continue
            else:
                break

        if area(A, B, C) > area(bA, bB, bC):
            bA = A; bB = B; bC = C

        A = (A+1)%n
        if A==B: B = (B+1)%n
        if B==C: C = (C+1)%n
        if A==0: break

    return bA,bB,bC

def TipFinder(thresh, debug_mode = 0):

    hull = convexHull(thresh, debug_mode)

    #############################
    bA,bB,bC = MaxEnclosedTriangle(hull)

    if debug_mode: plt.scatter(hull[[bA,bB,bC],0,0],hull[[bA,bB,bC],0,1],marker='*', s=500,c='g')

    return hull[[bA,bB,bC],0]
