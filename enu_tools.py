# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:05:39 2020

@author: Jingce
"""
import numpy as np


def days2mdh( year, days):
    """
    This function converts the day of the year, days, to the equivalent month day,
    hour, minute and second. From Vallado's original code translated to Python.

    Parameters
    ----------
    year : The year as a four digit number, e.g., 2013
    days : day of the year as a decimal number

    Returns
    -------
    mo : month as two digit integer
    day : day as two digits
    hour : hour as two digits
    minute : minute as two digits,
    sec : second as two digit integer, but can include fractional values

    Notes
    -----
    A support function in the class, currently not utilized by any methods.

    """
    # ---------------------------------------------------------------------
    #
    #                           function days2mdh
    #
    #  this function converts the day of the year, days, to the equivalent
    #  month day, hour, minute and second.
    #
    #  author        : david vallado            719-573-2600   22 jun 2002
    #
    #  revisions
    #                -
    #
    #  inputs          description              range / units
    #    year        - year                     900 .. 2100
    #    days        - julian day of the year   0.0  .. 366.0
    #
    #  outputs       :
    #    mon         - month                    1 .. 12
    #    day         - day                      1 .. 28,29,30,31
    #    hr          - hour                     0 .. 23
    #    minute      - minute                   0 .. 59
    #    sec         - second                   0.0 .. 59.999
    #
    #  locals        :
    #    dayofyr     - day of year
    #    temp        - temporary extended values
    #    inttemp     - temporary integer value
    #    i           - index
    #    lmonth(12)  - integer array containing the number of days per month
    #
    #  coupling      :
    #    none.
    #
    # [mon,day,hr,minute,sec] = days2mdh ( year,days);
    # ---------------------------------------------------------------------
    # Recoded to Python by Mark Wickert, July 2013

    lmonth = np.zeros(12)
    " --------------- set up array of days in month  --------------"
    for i in range(12):
        lmonth[i] = 31
        if i + 1 == 2:
            lmonth[i] = 28
        if i + 1 == 4 or i + 1 == 6 or i + 1 == 9 or i + 1 == 11:
            lmonth[i] = 30
    dayofyr = np.floor(days)
    " ----------------- find month and day of month ---------------"
    if np.mod(year - 1900, 4) == 0:
        lmonth[2 - 1] = 29
    i = 1 - 1
    inttemp = 0
    while (dayofyr > inttemp + lmonth[i]) and (i + 1 < 12):
        inttemp = inttemp + lmonth[i]
        i = i + 1
    mon = i + 1

    day = int(dayofyr - inttemp)
    " ----------------- find hours minutes and seconds ------------"
    temp = (days - dayofyr) * 24.0
    hr = int(temp)
    temp = (temp - hr) * 60.0
    minute = int(temp)
    sec = (temp - minute) * 60.0
    return mon, day, hr, minute, sec

def earth_model():
    """
    Define the constants from the WGS-84 ellipsoidal Earth model.


    Parameters
    ----------
    None

    Returns
    -------
    a : semi-major axis of the Earth ellipsoid model
    f : flattening

    Notes
    -----
    The World Geodetic System (WGS) is a standard for use in cartography, geodesy, and navigation.
    The latest revision is WGS-84.

    """
    a = 6378137.0  # meters
    f = 1.0 / 298.257223563
    return a, f


def ecef2llh(ecef):
    """
    Reference:RTKLIB rtkcmn.cpp function ecef2pos
    """
    x=ecef[0]
    y=ecef[1]
    z=ecef[2]
    
    llh=np.zeros(3)
    
    a, f=earth_model()
    e2=f*(2-f)
    r2=x*x+y*y

    zk=0.0
    while(np.fabs(z-zk)>=10E-4):
        zk=z
        sinp=z/np.sqrt(r2+z*z)
        N=a/np.sqrt(1.0-e2*sinp*sinp)
        z=ecef[2]+N*e2*sinp
    
    llh[0]=np.arctan(z/np.sqrt(r2))
    llh[1]=np.arctan2(np.array([y]),np.array([x]))
    #llh[1]=np.arctan(y/x)
    llh[2]=np.sqrt(r2+z*z)-N

    return llh

    
    
    
    

def llh2ecef(llh):
    """
    Convert lat,lon,hgt geographic coords to X,Y,Z Earth Centered Earth
    Fixed (ecef) or just (ecf) coords.

    Parameters
    ----------
    llh : A three element ndarray containing latitude(lat), longitude (lon), and altitude (a) or height (hgt), all in meters

    Returns
    -------
    x : The ecef x coordinate
    y : The ecef y coordinate
    z : The ecef z coordinate

    Notes
    -----
    This is a function that computes:
    N = a/sqrt( 1 - f*(2-f)*sin(lat)*sin(lat) )
    X = (N + h)*cos(lat)*cos(lon)
    Y = (N + h)*cos(lat)*sin(lon)
    Z = ((1-f)^2 * N + h)*sin(lat)
    by also calling EarthModel()

    Examples
    --------

    """

    lat = llh[0] * np.pi / 180.
    lon = llh[1] * np.pi / 180.
    hgt = llh[2]

    ecf = np.zeros(3)
    " Set up WGS-84 constants."
    a, f = earth_model()

    " Store some commonly used values."
    slat = np.sin(lat)
    N = a / np.sqrt(1 - f * (2 - f) * slat ** 2)
    Nplushgtclat = (N + hgt) * np.cos(lat)

    x = Nplushgtclat * np.cos(lon)
    y = Nplushgtclat * np.sin(lon)
    z = ((1 - f) ** 2 * N + hgt) * slat

    return np.array([x, y, z])


def ecef2enu(r_ecef, r_ref, phi_ref, lam_ref):
    """
    Convert ECEF coordinates to ENU using an ECEF reference location
    r_ref having lat = phi_ref and lon = lam_ref
    """
    # Convert lat and long angles in degress to radians
    phi_rad = phi_ref * np.pi / 180.0
    lam_rad = lam_ref * np.pi / 180.0
    # Form a 3-element column vector of the ECF (X,Y,Z) differences
    r_diff = np.array([r_ecef - r_ref]).T
    # Form the rotations transformation matrix
    A_matrix_ecef2enu = np.array([[-np.sin(lam_rad), np.cos(lam_rad), 0],
                                  [-np.sin(phi_rad) * np.cos(lam_rad), -np.sin(phi_rad) * np.sin(lam_rad),
                                   np.cos(phi_rad)],
                                  [np.cos(phi_rad) * np.cos(lam_rad), np.cos(phi_rad) * np.sin(lam_rad),
                                   np.sin(phi_rad)]])
    # Multiply the 3x3 matrix times the 3x1 column vector
    r_enu = np.dot(A_matrix_ecef2enu, r_diff)
    # Upon return flatten column vector back to a simple 1D
    # Also need to scale units of meters to what is needed
    return r_enu.flatten()

    