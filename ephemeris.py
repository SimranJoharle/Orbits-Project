# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 21:44:18 2020

@author: simran
"""
# a = semi major axis
# e = eccentricity
# i = angle of inclination
# com = capital omega =  longitude of ascending node
# om = omega = arguement of perihelion
# M0 = mean anomaly for the epoch of osculation t0
# ob = epsilon = angle of obliquity
# t = time seperation
# n = number of times coordinates needed

import numpy as np
from astropy.time import Time
from astropy.coordinates import get_sun
import datetime as dt
from prettytable import PrettyTable

def time(n,t0,t):
    """ Time Interval dates """
    t0 = dt.datetime.strptime(t0, "%d/%m/%Y")
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
    
def EAnom(M0,n,t,e,a):
    """ Calculating Eccentric Anomaly """
    p = 365.25193358350276
    P = a**1.5
    E_ = []
    for i in range(1,n+1):
        M = M0 + (2*np.pi*(t*i)/(P*p))
        E = np.pi
        for j in range(0,20):
            c = (M-e*(E*np.cos(E)-np.sin(E)))/(1- e*np.cos(E))
            E = c
        E_.append(E)
    return E_

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
        delta =np.rad2deg(np.sin(fz/np.sqrt(fx**2 + fy**2 + fz**2)))
        dist =np.sqrt(fx**2 + fy**2 + fz**2)
        rd.append([alpha,delta,dist])
    return(rd)

def makePrettyTable(table_col1, table_col2, table_col3):
    """ For table Output """
    table = PrettyTable()
    table.add_column("RA (deg)", table_col1)
    table.add_column("DEC (deg)", table_col2)
    table.add_column("Distance(AU)", table_col3)
    return print(table)

def ephemeris(a,e,i,com,om,M0,ob,t0,t,n):
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
    E = EAnom(Mo,n,t,e,a)
    rd = np.array(radec(Px,Py,Pz,Qx,Qy,Qz,e,n,a,b,E,sun))
    makePrettyTable(rd[:,0], rd[:,1], rd[:,2])
    
ephemeris(2.7664122,0.0791158,10.58347,80.48632,73.98440,189.27500,23.439291,"06/05/2020",70,5)    

