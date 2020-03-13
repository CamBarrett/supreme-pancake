#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:37:06 2020

@author: Cam
"""

import warnings
warnings.filterwarnings("ignore", category = RuntimeWarning)
import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt
import scipy as sc
from scipy import optimize
from scipy.stats import gaussian_kde
import astropy.units as u
import subprocess
from astropy.constants import R_sun
from astropy import constants as const
import matplotlib.mlab as mlab
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def round_down(num, div):
    return num - (num % div)

def Q_out(m_1, m_2, m_3):
    return np.nan_to_num(m_3 / (m_1 + m_2))

def P_out(m_1, m_2, a_out):
    return np.nan_to_num(2 * np.pi * np.sqrt((R_sol * a_out)**3 
                                          / (G * M_sol * (m_1 + m_2))))
def P_in(m_1, m_2, a_in):
    return np.nan_to_num(2 * np.pi * np.sqrt((R_sol * a_in)**3 
                                          / (G * M_sol * (m_1 + m_2))))
    
def crit_stability(e_out, inc, q_out):
    return np.nan_to_num((2.8 / (1 - e_out)) * ( 1 - 0.3 * inc / 
           np.pi) * ((1 + q_out) * (1 + e_out) / np.sqrt(1 - e_out))**(2/5))                   

def stability(a_out, a_in):
    return np.nan_to_num(a_out / a_in)

def t_kozai(P_out, P_in, m_1, m_2, m_3, e_out):
    return np.nan_to_num((P_out**2 / P_in) * (m_1 + m_2 + m_3) / m_3 * (1 - e_out**2)**(3/2))

def e_max(inc_0):
    return np.nan_to_num(np.sqrt(1 - 5 / 3 * np.cos(inc_0)**2))

def eps_oct(m_1, m_2, a_in, a_out, e_out):
    return np.nan_to_num((m_1 - m_2) / (m_1 + m_2) * a_in / a_out * e_out / (1 - e_out**2))
    
def orb_av_lim(q_out, a_in, a_out, e_out):
    return np.nan_to_num(f_crit * 5 * np.pi * q_out * (a_in / (a_out * (1 - e_out)))**3)

def semi_sec(e_max):
    return np.nan_to_num(np.sqrt(1 - e_max))

def round_down(num, divisor):
    return num - (num % divisor)

def save(name, array):
    np.savetxt(name, array)
    
def add_counter(file):
    with open(file, 'r') as program:
        data = program.readlines()

    with open(file, 'w') as program:
        for (number, line) in enumerate(data):
            program.write('%d  %s' % (number, line))

def s_crit(G1, G2):
    diff = np.abs(G1 - G2)
    return 10**(0.075*diff - 0.53)

def ang_sep(a, e, d, i):
    a *= R_sun.value
    d *= const.kpc.value
    apparent = a * (1+(e**2)/2) / d
    theta = apparent/2
    phi = apparent/2
    return apparent * (np.sin((i-theta)*np.pi/180) + np.sin((180-phi-i)*np.pi/180))/2

def proper_resolve(l,b,d):
    solar_x, solar_y, solar_z = 11.1, 12.4, 7.25 #km/s
    
    longitude = l
    latitude = b
    d = d
    
    phi = random.uniform(0, np.pi)
    theta = random.uniform(0, np.pi*2)
    v = np.random.normal(53,23)
    
    v_x = v * np.sin(phi) * np.cos(theta)
    v_y = v * np.sin(phi) * np.sin(theta)
    v_z = v * np.cos(phi)
    
    res_x = solar_x - v_x # in km/s, so distance here is in km
    res_y = solar_y - v_y
    res_z = solar_z - v_z
    
    
    if (np.abs(longitude) < (np.pi/2)) : longitude = np.abs(longitude)
    if ((np.pi/2) < np.abs(longitude) < np.pi) : longitude = np.abs((np.pi/2) - np.abs(longitude))
    if (np.pi < np.abs(longitude) < (3*np.pi/2)) : longitude = np.abs((3*np.pi/2) - np.abs(longitude))
    if ((3*np.pi/2) < np.abs(longitude) < (2*np.pi)) : longitude = np.abs(np.abs(longitude) - (3*np.pi/2))
    
    if (np.abs(latitude) < (np.pi/2)) : latitude = np.abs(latitude)
    if ((np.pi/2) < np.abs(latitude) < np.pi ) : latitude = np.abs((np.pi/2) - np.abs(latitude))
    
    
    
    d = d* 3.1e19 # changing into metres
    res_x = res_x * 1e3
    res_y = res_y * 1e3
    res_z = res_z * 1e3
    
    
    # last part of these is change from radians per second into milli arcsecs per year
    #    theta_x = np.abs(np.arcsin( (res_x * np.sin(longitude + (np.pi/2))) / np.sqrt(res_x**2 + d**2 - 2*res_x*d*np.cos(longitude + (np.pi/2)))) * 3.154e7 * (180/np.pi) * 1000 * 3600)
    #    theta_y = np.abs(np.arcsin( (res_y * np.sin(np.pi - longitude)) / np.sqrt(res_y**2 + d**2 - 2*res_y*d*np.cos(np.pi - longitude))) * 3.154e7 * (180/np.pi) * 1000 * 3600)
    #    theta_z = np.abs(np.arcsin( (res_z * np.sin(latitude + (np.pi/2))) / np.sqrt(res_z**2 + d**2 - 2*res_z*d*np.cos(latitude + (np.pi/2)))) * 3.154e7 * (180/np.pi) * 1000 * 3600)
    
    
    #print(theta_x)
    #print(theta_y)
    #print(theta_z)
    
    
    latitude = np.abs(np.pi - latitude)
    R_lat = np.array([[np.cos(latitude), -np.sin(latitude), 0], [np.sin(latitude), np.cos(latitude), 0], [0,0,1]])
    
    R_long = np.array([[1,0,0], [0, np.cos(longitude), -np.sin(longitude)], [0, np.sin(longitude), np.cos(longitude)]])
    
    old_coords = np.array([[res_x], [res_y], [res_z]])
    
    R = np.matmul(R_lat, R_long)
    
    new = np.matmul(R, old_coords)
    
    A = np.abs(new[0][0])
    B = np.abs(new[2][0])
    
    
    D = np.sqrt(A**2 + B**2)
    theta = np.arctan(D/d) * 3.154e7 * (180/np.pi) * 1000 * 3600
    
    return theta




def processor(name):
    data = np.loadtxt((open('output_{}.txt'.format(name),'rt').readlines()[:-1]), skiprows=4, dtype=None)
    print('===============\n\n')
    print('number of systems in {} = {}'.format(name, len(data)))
    
    G1, G2, G3 = data[:,31], data[:,34], data[:,37]
    GBP_1, GBP_2, GBP_3 = data[:,32], data[:,35], data[:,38]
    GRP_1, GRP_2, GRP_3 = data[:,33], data[:,36], data[:,39]
    semi_major, ecc, d, inc = data[:,45], data[:,47], data[:,5], data[:,49]
    star_1_id, star_2_id, star_3_id = data[:,40], data[:,41], data[:,42]
    
    index_all_obs = []
    index_Gaia_suit = []
    index_pos_trip = []
    
    for i in range(len(data)):
        if ((G1[i] < Gaia_mag_res) and (G2[i] < Gaia_mag_res) and (G3[i] < Gaia_mag_res)):
            index_all_obs.append(i)
            
#            if ((1/np.pi)*2*180*3600*np.arctan(data[:,45][i]*6.96e8 / (2*data[:,47][i]*3.1e19))) > Gaia_ang_res:
            if ang_sep(semi_major[i], ecc[i], d[i], inc[i]) > s_crit(G1[i], G2[i]):
                index_Gaia_suit.append(i)
            else:
                if (star_1_id[i] < 2) and (GBP_2[i] - GBP_1[i] < 1) and (GRP_2[i] - GRP_1[i] < 1): 
                    index_Gaia_suit.append(i)
                if (star_2_id[i] < 2) and (GBP_1[i] - GBP_2[i] < 1) and (GRP_1[i] - GRP_2[i] < 1):
                    index_Gaia_suit.append(i)
                    
                elif (np.abs(G1[i] - G2[i])<1):
                    index_Gaia_suit.append(i)
                    
                
                
                
        if (G1[i] < Gaia_mag_res and G2[i] < Gaia_mag_res and G3[i] > Gaia_mag_res) or (G1[i] < Gaia_mag_res and G2[i] > Gaia_mag_res and G3[i] < Gaia_mag_res) or (G1[i] > Gaia_mag_res and G2[i] < Gaia_mag_res and G3[i] < Gaia_mag_res):
            index_pos_trip.append(i)
    
    index_all_obs = np.array(index_all_obs)
    index_Gaia_suit = np.array(index_Gaia_suit)
    index_pos_trip = np.array(index_pos_trip)
    print("number lost due to (angular) resolution = {}".format(len(index_all_obs) - len(index_Gaia_suit)))
    Gaia_suit = data[index_Gaia_suit]
        
    index_Gaia = []
    for i in range(len(Gaia_suit)):
        theta = random.uniform(0, np.pi/2)
        if theta == 0:
            break
        else:            
            l = Gaia_suit[i][3] * np.pi/180
            b = Gaia_suit[i][4] * np.pi/180
            d = Gaia_suit[i][5] # changed into parsecs
                        
            angle = proper_resolve(l,b,d)
            if angle > prop_motion_res:
                index_Gaia.append(i)
    
    index_Gaia = np.array(index_Gaia)
    Gaia_triples_final = data[index_Gaia]
    
    
    
    Gaia_triples = data[index_Gaia_suit]
    pos_triples = data[index_pos_trip] # sort into whether it is inner star or outer that cannot see
    Gaia_id = Gaia_triples[:,44]
    print('difference after proper common motion is {}'.format(len(Gaia_triples) - len(Gaia_triples_final)))
    
    inclination = Gaia_triples[:,1]
    l = Gaia_triples[:,3] # longitude
    b = Gaia_triples[:,4] # latitude - 0 means in horizontal plane
    d = Gaia_triples[:,5] # r in kpc
    m_1 = Gaia_triples[:,7]
    m_2 = Gaia_triples[:,10]
    m_3 = Gaia_triples[:,13]
    inner_sep = Gaia_triples[:,45]
    outer_sep = Gaia_triples[:,46]
    
    
    index_inside_blind = []
    index_outside_blind = []
    for i in range(len(pos_triples)):
        if (pos_triples[:,31][i] > Gaia_mag_res or pos_triples[:,34][i] > Gaia_mag_res):
            index_inside_blind.append(i)
        if (pos_triples[:,37][i] > Gaia_mag_res):
            index_outside_blind.append(i)
    index_inside_blind = np.array(index_inside_blind)
    index_outside_blind = np.array(index_outside_blind)
    inside = pos_triples[index_inside_blind]
    outside = pos_triples[index_outside_blind]
    
    
    print('unique systems in {} = {}'.format(name, len(np.unique(Gaia_id))))
    print('number of Gaia observables is {}'.format(len(Gaia_triples)))
    print('possible missed triples = {}'.format(len(pos_triples)))
    print('of which {} are inside bind and {} are outside blind'.format(len(inside), len(outside)))
    print('{:.2f}% are the systems are gaia observable'.format(len(Gaia_triples)*100/len(data)))
    

        
    l = l * np.pi/180
    b = b * np.pi/180
    
    theta = np.pi - l   
    phi = (np.pi/2) - b # to convert into the polar angle

    x = d * np.cos(theta) * np.sin(phi)
    y = d * np.sin(theta) * np.sin(phi)
    z = d * np.cos(phi)
    
    mpl.rcParams['legend.fontsize'] = 15

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z, c='r', marker='o', s=0.09, label='Triples')
#    ax.scatter(x_pos,y_pos,z_pos, c='y', marker='o', s=0.6, label='Possible Triples')
    ax.scatter(0,0,0, c='b', marker='+', label='Earth')
    ax.plot([], [], [], ' ', label='Num Triples = {}'.format(len(Gaia_triples)))
    ax.set_xlabel('x axis (kpc)')
    ax.set_ylabel('y axis (kpc)')
    ax.set_zlabel('z axis (kpc)')
    ax.legend()
    

    plt.show()
    
    from scipy.stats import norm

    fig = plt.figure()
    num_bins = 70
## the histogram of the data
    n, bins, patches = plt.hist(d, num_bins, facecolor='b', alpha=0.4, edgecolor='black')#, num_bins, density=True, facecolor='b', alpha=0.4, edgecolor='black')


## add a 'best fit' line
#    y = norm.pdf(num_bins, np.mean(d), np.std(d))
    plt.xlabel('d (kpc)')
    plt.ylabel('Frequency')
    plt.title('{}'.format(name))
    
    fig = plt.figure()
    num_bins = 70
## the histogram of the data
    n, bins, patches = plt.hist(z, num_bins, facecolor='b', alpha=0.4, edgecolor='black')#, num_bins, density=True, facecolor='b', alpha=0.4, edgecolor='black')


## add a 'best fit' line
#    y = norm.pdf(num_bins, np.mean(d), np.std(d))
    plt.xlabel('z (kpc)')
    plt.ylabel('Frequency')
    plt.title('{}'.format(name))
#    plt.plot(num_bins, y, 'r--', label='norm pdf')
#    plt.legend()

    plt.subplots_adjust(left=0.15)
    plt.show()
    print('')
    print('')
    print('')