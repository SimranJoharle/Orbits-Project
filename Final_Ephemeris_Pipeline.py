#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import get_sun
from scipy.optimize import fsolve
import math
import datetime as dt
from prettytable import PrettyTable
from astropy import units as u
from astropy.coordinates import SkyCoord
import scipy.constants as sc
from astropy.table import QTable


# In[2]:


def cnv(ra,dec):
    """ For converting ra and dec values to hms and dms respectively. """
    ra= ra.split()
    ra = ra[0]+'h'+ra[1]+'m'+ra[2]+'s'
    dec = dec.split()
    dec = dec[0]+'d'+dec[1]+'m'+dec[2]+'s'
    c = SkyCoord(ra, dec)
    r = np.radians(c.ra.value)
    d = np.radians(c.dec.value)
    return r,d


# In[3]:


def UTC2TT(tm):
    """ To convert the given time from UTC to TT. """
    times = dt.datetime.strptime(tm, "%Y-%b-%d")
    t = str(times.year)+"-"+str(times.month)+"-"+str(times.day)
    t1 = Time(t, format='isot', scale='ut1')
    return (t1.utc).tt


# In[4]:


def sun(date):
    """ To get the Sun's geocentric coordinates"""
    to = Time([date])
    ra = np.deg2rad((get_sun(to).ra.value)[0])
    dec = np.deg2rad((get_sun(to).dec.value)[0])
    r = (get_sun(to).distance.value)[0]
    x = r*np.cos(dec)*np.cos(ra)
    y = r*np.cos(dec)*np.sin(ra)
    z = r*np.sin(dec)
    return x,y,z


# In[5]:


def hd(l, m, n, r0, n0, z0, g):
    """ To calculate heliocentric distances. """
    xh = l*g - r0
    yh = m*g - n0
    zh = n*g - z0
    r = np.sqrt(xh**2 + yh**2 + zh**2)
    return r,xh,yh,zh


# In[6]:



def fapprox(da1,da2,da3,ra1,ra2,ra3,de1,de2,de3):
""" For first approximation of heliiocentric and geocentric distances. """
for i in range(0,2):
    t1 = (da3 - da2).value
    t2 = (da3 - da1).value
    t3 = (da2 - da1).value
    b1 = t1/t2
    b3 = t3/t2

    def dc(r,d):
        l = np.cos(r)*np.cos(d)
        m = np.sin(r)*np.cos(d)
        n = np.sin(d)
        return l,m,n

    l1, m1, n1 = dc(ra1, de1)
    l2, m2, n2 = dc(ra2, de2)
    l3, m3, n3 = dc(ra3, de3)
    



    r01, n01, z01 = sun(da1)
    r02, n02, z02 = sun(da2)
    r03, n03, z03 = sun(da3)



    a = np.array([[l1*b1, -l2, l3*b3], [m1*b1, -m2, m3*b3], [n1*b1, -n2, n3*b3]])
    b = np.array([(b1*r01 - r02 + b3*r03), (b1*n01 -
                                            n02 + b3*n03), (b1*z01 - z02 + b3*z03)])

    x = np.linalg.solve(a, b)

    g1,g2,g3 = x
    r1,xh1,yh1,zh1 = hd(l1, m1, n1, r01, n01, z01, g1)
    r2,xh2,yh2,zh2 = hd(l2, m2, n2, r02, n02, z02, g2)
    r3,xh3,yh3,zh3 = hd(l3, m3, n3, r03, n03, z03, g3)

    t1 = (da3 - da2).value / 58.13244087
    t2 = (da3 - da1).value / 58.13244087
    t3 = (da2 - da1).value / 58.13244087



    prec = 1
    r_old = r1
    count = 0
    while (prec > 1e-5):

        a1 = b1 + ((t1*t3)*(1+b1))/(6*(r2)**3)
        a3 = b3 + ((t1*t3)*(1+b3))/(6*(r2)**3)
        a = np.array([[l1*a1, -l2, l3*a3], [m1*a1, -m2, m3*a3], [n1*a1, -n2, n3*a3]])
        b = np.array([(a1*r01 - r02 + a3*r03), (a1*n01 - n02 + a3*n03), (a1*z01 - z02 + a3*z03)])
        try:
            Deltas = np.linalg.solve(a,b)
        except:
            Deltas = np.linalg.lstsq(a,b,rcond=None)

        g1,g2,g3 = Deltas
        r1,xh1,yh1,zh1 = hd(l1, m1, n1, r01, n01, z01, g1)
        r2,xh2,yh2,zh2 = hd(l2, m2, n2, r02, n02, z02, g2)
        r3,xh3,yh3,zh3 = hd(l3, m3, n3, r03, n03, z03, g3)


        prec = abs(r1 - r_old)/r_old
        count +=1
        r_old = r1


    imp = np.array([[r1,xh1,yh1,zh1],[r2,xh2,yh2,zh2],[r3,xh3,yh3,zh3]])


    if (i == 1):
        break

    da1 = da1 - (g1/10065.320)
    da2 = da2 - (g2/10065.320)
    da3 = da3 - (g3/10065.320)
    
return imp,l1,m1,n1,l2,m2,n2,l3,m3,n3,r01, n01, z01,r02, n02, z02,r03, n03, z03,t1,t2,t3,b1,b3


# In[7]:


def cosf(i,j,imp):
    """ To calculate cos(f). """
    r1,xh1,yh1,zh1 = imp[i]
    r2,xh2,yh2,zh2 = imp[j]
    cos2f3 = (xh1*xh2 + yh1*yh2 + zh1*zh2)/(r1*r2)
    cosf3 = np.sqrt(0.5*(cos2f3 + 1))
    return cosf3


# In[8]:


def R(i,j,t,cos,imp):
    """ To calculate the sector-area ratios. """
    r1, r2 = [imp[i,0],imp[j,0]]
    M32 = t**2 / (4 * (np.sqrt(r1*r2)*cos)**3)
    N3 = (r1+r2)/(2*np.sqrt(r1*r2)*cos)
    
    def equations(p):
        R3, g3 = p
        f = np.zeros(2)
        f[0] = (R3)**2 - (M32 / (N3 - np.cos(g3)))
        f[1] = (R3)**3 - (R3)**2 - (M32*(g3 - (np.sin(g3)*np.cos(g3)))/(np.sin(g3))**3) 
        return f

    z =  fsolve(equations, (1, np.arccos(cos)))
    
    return z[0]


# In[9]:


def finalapprox(a1,a3,l1,m1,n1,l2,m2,n2,l3,m3,n3,r01, n01,z01,r02,n02,z02,r03,n03,z03):
    """ Final approximation of the geocentric and heliocentric distances. """
    a = np.array([[l1*a1, -l2, l3*a3], [m1*a1, -m2, m3*a3], [n1*a1, -n2, n3*a3]])
    b = np.array([(a1*r01 - r02 + a3*r03), (a1*n01 -
                                            n02 + a3*n03), (a1*z01 - z02 + a3*z03)])
    try:
        x = np.linalg.solve(a, b)
    except LinAlgError:
        x = np.linalg.lstsq(a, b)[0]
    gf1,gf2,gf3 = x

    imp1 = np.array([np.array(hd(l1, m1, n1, r01, n01, z01, gf1)),np.array(hd(l2, m2, n2, r02, n02, z02, gf2)), np.array(hd(l3, m3, n3, r03, n03, z03, gf3))])
    rf1 = hd(l1, m1, n1, r01, n01, z01, gf1)[0]
    rf2 = hd(l2, m2, n2, r02, n02, z02, gf2)[0]
    rf3 = hd(l3, m3, n3, r03, n03, z03, gf3)[0]
    return gf1,gf2,gf3,rf1,rf2,rf3,imp1


# In[10]:


def latr(R1,R2,R3,t1,t2,t3,rf1,rf2,rf3,imp1):
    """ To calculate the latus rectum. """
    lr3 = (R3  *(rf1*rf2 * np.sin(2*np.arccos(cosf(0,1,imp1))))/t3)**2
    lr2 = (R2  *(rf1*rf3 * np.sin(2*np.arccos(cosf(0,2,imp1))))/t2)**2
    lr1 = (R1  *(rf2*rf3 * np.sin(2*np.arccos(cosf(1,2,imp1))))/t1)**2
    lr = (lr1+lr2+lr3)/3
    return lr,lr1,lr2,lr3


# In[11]:


def f2(i,j,imp1):
    """ To calculate 2f for further calculations. """
    r1,xh1,yh1,zh1 = imp1[i]
    r2,xh2,yh2,zh2 = imp1[j]
    cos2f3 = (xh1*xh2 + yh1*yh2 + zh1*zh2)/(r1*r2)
    return np.arccos(cos2f3)


# In[12]:


def ecc(lr,rf1,rf2,rf3,imp1): 
    """ To calculate eccentricity. """
    e = np.sqrt(((lr/rf1)-1)**2 + (((((lr/rf1)-1)*np.cos(f2(0,2,imp1)))-((lr/rf3)-1))/np.sin(f2(0,2,imp1)))**2)
    cv1 = ((lr/rf1) - 1)/e
    sv1 = ((((lr/rf1)-1)*np.cos(f2(0,2,imp1)))-((lr/rf3)-1))/(e*np.sin(f2(0,2,imp1)))
    atv1 = np.arctan2(sv1,cv1)
    if (atv1<0):
        v1 = (2*np.pi)+atv1
    else:
        v1 = atv1
        
    v2 = v1 + f2(0,1,imp1)
    v3 = v1 + f2(0,2,imp1)
    return e,v1,v2,v3    


# In[13]:


def ArgP(Pz,Py,Qz,Qy): 
    """ To calculate the Arguement of perihelion. """
    ob = np.deg2rad(23.438960)
    w_ = math.atan2(((Pz*np.cos(ob))-(Py*np.sin(ob))),((Qz*np.cos(ob))-(Qy*np.sin(ob))))
    if (np.rad2deg(w_)<0):
        w = (2*np.pi)+w_
    else :
        w  = w_
    return w    


# In[14]:


def Node(Px,Py,Qx,Qy,w): 
    """ To calculate the Longitude of Ascending Node. """
    ob = np.deg2rad(23.438960)
    cW = (Px*np.cos(w)) - (Qx*np.sin(w))
    sW =((Py*np.cos(w)) - (Qy*np.sin(w)))/np.cos(ob)
    atW = np.arctan2(sW,cW)
    if (atW<0):
        W = (2*np.pi)+atW
    else:
        W = atW
    return W


# In[15]:


def inc(Px,Qx,Qz,Qy,w,W):
    """ To calculate the angle of inclination. """
    ob = np.deg2rad(23.438960)
    ci = -((Px*np.sin(w))+(Qx*np.cos(w)))/np.sin(W)
    si = ((Qz*np.cos(ob)) - (Qy*np.sin(ob)))/np.cos(w)
    ati = np.arctan2(si,ci)

    if (ati<0):
        i = (2*np.pi)+ati
    else:
        i = ati
    return i         


# In[17]:


def deltaT(e,da1,v,t0,P):
    """ To calculate the time of Periastron passage. """
    cE = (e+np.cos(v))/(1+(e*np.cos(v)))
    sE = ((np.sin(v)*np.sqrt(1 - e**2)))/(e+np.cos(v))    
    E = np.arccos(cE)
    if((sE/cE)<0):
        E =  2*np.pi - E
                                          
    dd = UTC2TT(t0)
        
    delt = -((P*(E - e*np.sin(E)))/(2*np.pi))+P
    print('delt=',delt)
    TP = (da1 + delt)
    d = (dd-TP).value
    MA = np.rad2deg(((2*np.pi)/P)*(d))
        
    return TP,MA


# In[25]:


#To compute the Ephemeris.

def time(n,t0,t):
    """ Time Interval dates """
    t0 = dt.datetime.strptime(t0, "%Y-%m-%d")
    date = []
    for i in range(1,n+1):    
        d = t0 + dt.timedelta(days=(t*i))
        date.append(d)
    return date

def suncoords(date,n):
    """ Sun Coordinates at the required Time intervals """
    sun = []
    for i in range(0,n):
        v = get_sun(Time(date[i]))
        ra = np.deg2rad(v.ra.value)
        dec = np.deg2rad(v.dec.value)
        r = v.distance.value
        x = r*np.cos(dec)*np.cos(ra)
        y = r*np.cos(dec)*np.sin(ra)
        z = r*np.sin(dec)
        sun.append([x,y,z])
    return sun
    
def EAnom(M0,n,t,e,a,date):
    """ Calculating Eccentric Anomaly """
    p = 365.25193358350276
    P = a**1.5
    M = []
    for i in range(0,n):
        m = M0 + (2*np.pi*(date[i]-date[0]).days/(P*p))
        M.append(m)        
    
    E = np.zeros(n)    
    for i in range(n):
        def funcs(vari):
            return vari-e*np.sin(vari)-M[i]
        E[i] = fsolve(funcs, 3)
    return E

def radec(Px,Py,Pz,Qx,Qy,Qz,e,n,a,b,E,sun):
    """ Calculating radec """
    c1 = []
    for i in range(0,n):
        X = a*Px*(np.cos(E[i])-e) + b*Qx*np.sin(E[i])
        Y = a*Py*(np.cos(E[i])-e) + b*Qy*np.sin(E[i])
        Z = a*Pz*(np.cos(E[i])-e) + b*Qz*np.sin(E[i])
        c1.append([X,Y,Z])
        
    rd = []
    for j in range(0,n):
        c2 = sun[j]
        c= c1[j]
        fx = c[0]+c2[0]
        fy = c[1]+c2[1]
        fz = c[2]+c2[2]
        alpha =np.rad2deg(np.arccos(fx/np.sqrt(fx**2 + fy**2)))
        if (fy/np.sqrt(fx**2 + fy**2))<0:
            alpha = 360 - alpha
        delta =np.rad2deg(np.arcsin(fz/np.sqrt(fx**2 + fy**2 + fz**2)))
        dist =np.sqrt(fx**2 + fy**2 + fz**2)
        cc = SkyCoord(ra=alpha*u.degree, dec=delta*u.degree, frame='icrs')
        cs = cc.to_string('hmsdms')
        rd.append([cs.split()[0],cs.split()[1],dist])
    return(rd)

def makePrettyTable(table_col1, table_col2, table_col3,table_col4):
    """ For table Output """
    table = PrettyTable()
    table.add_column("Date", table_col1)
    table.add_column("RA ", table_col2)
    table.add_column("DEC ", table_col3)
    table.add_column("Distance(AU)", table_col4)
    return print(table)


def ephemeris(a,e,i,com,om,M0,ob,t0,t,n):
    """ to call the functions one by one for Ephemeris computation. """
    b = a*(1 - e**2)**0.5
    i = np.deg2rad(i)
    com = np.deg2rad(com)
    om = np.deg2rad(om)
    Mo = np.deg2rad(M0)
    ob = np.deg2rad(ob)
    
    Px = np.cos(com)*np.cos(om) - np.sin(com)*np.sin(om)*np.cos(i)
    Qx = -(np.cos(com)*np.sin(om) + np.sin(com)*np.cos(om)*np.cos(i))
    Py = (np.sin(com)*np.cos(om) + np.cos(com)*np.sin(om)*np.cos(i))*np.cos(ob) - np.sin(om)*np.sin(i)*np.sin(ob)
    Qy = (-np.sin(com)*np.sin(om) + np.cos(com)*np.cos(om)*np.cos(i))*np.cos(ob) - np.cos(om)*np.sin(i)*np.sin(ob)
    Pz = (np.sin(com)*np.cos(om) + np.cos(com)*np.sin(om)*np.cos(i))*np.sin(ob) + np.sin(om)*np.sin(i)*np.cos(ob)
    Qz = (-np.sin(com)*np.sin(om) + np.cos(com)*np.cos(om)*np.cos(i))*np.sin(ob) + np.cos(om)*np.sin(i)*np.cos(ob)
    
    date = time(n,t0,t)
    sun = suncoords(date,n)
    E = EAnom(Mo,n,t,e,a,date)
    rd = np.array(radec(Px,Py,Pz,Qx,Qy,Qz,e,n,a,b,E,sun))
    makePrettyTable(date,rd[:,0], rd[:,1], rd[:,2])
    return date,rd[:,0], rd[:,1], rd[:,2]


# In[31]:


def Orbitdetermination(csvfile,t0,t,n):
    """ to call the functions one by one for the determination of orbita parameters, followed by the ephemeris.
    csvfile = csvfile with 3 observations including date, ra and dec.
    t0 = Begining date of Ephemeris
    t = Time Interval
    n = Total Number of days."""
    in_file = csvfile
    data = np.genfromtxt(in_file,delimiter=',',dtype=str)
    
    ra1, de1 = cnv(data[0,1], data[0,2])
    ra2, de2 = cnv(data[1,1], data[1,2])
    ra3, de3 = cnv(data[2,1], data[2,2])
    
    da1 = UTC2TT(data[0,0])
    da2 = UTC2TT(data[1,0])
    da3 = UTC2TT(data[2,0])
    
    imp,l1,m1,n1,l2,m2,n2,l3,m3,n3,r01, n01, z01,r02, n02, z02,r03, n03, z03,t1,t2,t3,b1,b3 = fapprox(da1,da2,da3,ra1,ra2,ra3,de1,de2,de3)
    
    R3 = R(0,1,t3,cosf(0,1,imp),imp)
    R2 = R(0,2,t2,cosf(0,2,imp),imp)
    R1 = R(1,2,t1,cosf(1,2,imp),imp)
    
    a1 = (R2*b1)/R1
    a3 = (R2*b3)/R3
    
    gf1,gf2,gf3,rf1,rf2,rf3,imp1 = finalapprox(a1,a3,l1,m1,n1,l2,m2,n2,l3,m3,n3,r01, n01,z01,r02,n02,z02,r03,n03,z03)
    
    lr,lr1,lr2,lr3 = latr(R1,R2,R3,t1,t2,t3,rf1,rf2,rf3,imp)
    
    e,v1,v2,v3 = ecc(lr,rf1,rf2,rf3,imp1)  
    
    a = lr/(1 - e**2)
    
    P = (a**1.5)*365.25636
    
    Px = ((imp1[0,1]*rf3*np.sin(v3)) - (imp1[2,1]*rf1*np.sin(v1)))/(rf1*rf3*np.sin(f2(0,2,imp1)))
    Qx = ((imp1[2,1]*rf1*np.cos(v1)) - (imp1[0,1]*rf3*np.cos(v3)))/(rf1*rf3*np.sin(f2(0,2,imp1)))
    Py = ((imp1[0,2]*rf3*np.sin(v3)) - (imp1[2,2]*rf1*np.sin(v1)))/(rf1*rf3*np.sin(f2(0,2,imp1)))
    Qy = ((imp1[2,2]*rf1*np.cos(v1)) - (imp1[0,2]*rf3*np.cos(v3)))/(rf1*rf3*np.sin(f2(0,2,imp1)))
    Pz = ((imp1[0,3]*rf3*np.sin(v3)) - (imp1[2,3]*rf1*np.sin(v1)))/(rf1*rf3*np.sin(f2(0,2,imp1)))
    Qz = ((imp1[2,3]*rf1*np.cos(v1)) - (imp1[0,3]*rf3*np.cos(v3)))/(rf1*rf3*np.sin(f2(0,2,imp1)))
    
    w = ArgP(Pz,Py,Qz,Qy)
    W = Node(Px,Py,Qx,Qy,w)
    i = inc(Px,Qx,Qz,Qy,w,W)
    TP,M0 = deltaT(e,da1,v1,t0,P) 
    
    t_0 = (UTC2TT(t0).value).split("T")[0]
    
    print("OBITAL ELEMENTS:")
    print("-----------------")
    print(f"latus rectum = {lr}\neccentricity = {e}\na = {a}\nPeriod = {P}\nArgP = {np.rad2deg(w)}\nNode = {np.rad2deg(W)}\ninclination = {np.rad2deg(i)}\nTP = {TP.value}\nMA = {M0}")
    print("\n")
    print("EPHEMERIS:")
    date,ria, decl, dista = ephemeris(a,e,np.rad2deg(i),np.rad2deg(W),np.rad2deg(w),M0,23.438960,t_0,t,n)
    return date,ria, decl, dista


# In[ ]:




