#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 16:15:10 2020

@author: Cam
"""

import numpy as np
from astropy.constants import R_sun
from astropy import constants as const
import matplotlib.mlab as mlab
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from Function_list import proper_resolve


#stable_data = np.loadtxt((open('output_stable.txt','rt').readlines()[:-1]), skiprows=4, dtype=None)
#LK_data = np.loadtxt((open('output_LK.txt','rt').readlines()[:-1]), skiprows=4, dtype=None)
#ecc_data = np.loadtxt((open('output_ecc.txt','rt').readlines()[:-1]), skiprows=4, dtype=None)
##file = 'fresh.txt'
#data = np.loadtxt(file, skiprows=4)

#print(len(data))
#print('')

Gaia_mag_res = 21
Gaia_ang_res = 0.4 #arcsec
v = 220 #km/s
prop_motion_res = 0.1 # in milli arcsec per year
Gaia_obs_timescale = 640 #days,, for astrometry?
Gaia_lifetime = 4*365 # 4 years in days

#name = 'stable'
#data = np.loadtxt((open('output_{}.txt'.format(name),'rt').readlines()[:-1]), skiprows=4, dtype=None)
#print('The shape is {}'.format(np.shape(data)))

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
    pos_triples = data[index_pos_trip] # sort into whether it is inner star or outer that cannot see

    
#    index_Gaia = []
#    for i in range(len(Gaia_suit)):
#        theta = random.uniform(0, np.pi/2)
#        if theta == 0:
#            break
#        else:
#            phi = (np.pi/2) - theta
#            v_h = v * np.cos(phi)
#            d = Gaia_suit[i][5]
#            angle = np.arctan( (np.abs(v_h) * 1000) / (d*3.1e19) ) * 3.154e7 * (180/np.pi) * 1000 * 3600 #in milli arcsec / yr 
#            if angle > prop_motion_res:
#                index_Gaia.append(i)
#    
#    index_Gaia = np.array(index_Gaia)
#    Gaia_triples_final = data[index_Gaia]
    
    
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
    Gaia_suit_proper = Gaia_suit[index_Gaia]
    
    print('difference after proper common motion is {}'.format(len(Gaia_suit) - len(Gaia_suit_proper)))

  

    print('Before observation period = {}'.format(len(Gaia_suit_proper)))
    
    index_period_check = []
    for i in range(len(Gaia_suit_proper)):
        if Gaia_suit_proper[i][1] > 6*Gaia_obs_timescale:
            index_period_check.append(i)
            
    index_period_check = np.array(index_period_check)
    Gaia_triples = Gaia_suit_proper[index_period_check]
    Gaia_id = Gaia_triples[:,44]  
    
    print('After observation period = {}'.format(len(Gaia_triples)))
    
    
    inclination = Gaia_triples[:,1]
    l = Gaia_triples[:,3] # longitude
    b = Gaia_triples[:,4] # latitude - 0 means in horizontal plane
    d = Gaia_triples[:,5] # r in kpc
    m_1 = Gaia_triples[:,7]
    m_2 = Gaia_triples[:,10]
    m_3 = Gaia_triples[:,13]
    inner_sep = Gaia_triples[:,45]
    outer_sep = Gaia_triples[:,46]
    
    
    
#    # Calculate the distances to each star
#    r_in = [] #distance from centre of mass to the inner binary cm
#    r_out = []
#    inn_1 = [] #distance from star 1 to the inner binary cm
#    inn_2 = []
#    cm_1 = [] #distance from star 1 to in the centre of mass
#    cm_2 = []
#    d_1 = []
#    d_2 = []
#    d_3 = []
#    for i in range(len(Gaia_triples)):
#        r_in.append(outer_sep[i] * (m_3[i] / (m_1[i] + m_2[i] + m_3[i])))
#        r_out.append(outer_sep[i] * (m_1[i] + m_2[i])/(m_1[i] + m_2[i] + m_3[i]))
#        inn_1.append(inner_sep[i] * (m_2[i])/(m_1[i] + m_2[i]))
#        inn_2.append(inner_sep[i] * (m_1[i])/(m_1[i] + m_2[i]))
#
#    r_in = np.array(r_in)    
#    r_out = np.array(r_out)
#    inn_1 = np.array(inn_1)
#    inn_2 = np.array(inn_2)
#    
#    for i in range(len(Gaia_triples)):
#        if inclination[i]==0:
#            cm_1.append(r_in[i] + inn_1[i])
#            cm_2.append(r_in[i] - inn_2[i])
#        else:
#            cm_1.append(np.sqrt(inn_1[i]**2 + r_in[i]**2 - 2*inn_1[i]*r_in[i]*np.cos((180-inclination[i])*np.pi/180)))
#            cm_2.append(np.sqrt(inn_2[i]**2 + r_in[i]**2 - 2*inn_2[i]*r_in[i]*np.cos((inclination[i])*np.pi/180)))
#        
#        d_1.append( np.sqrt((d[i]*3.1e19)**2 + (cm_1[i]*6.96e8)**2) / 3.1e19 )
#        d_2.append( np.sqrt((d[i]*3.1e19)**2 + (cm_2[i]*6.96e8)**2) / 3.1e19 )
#        d_3.append( np.sqrt((d[i]*3.1e19)**2 + (r_out[i]*6.96e8)**2) / 3.1e19 )
#    
#    cm_1 = np.array(cm_1)
#    cm_2 = np.array(cm_2)
#    d_1 = np.array(d_1)
#    d_2 = np.array(d_2)
    
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
    
    print('number of Gaia observables is {}'.format(len(Gaia_triples)))
    print('unique systems in {} = {}'.format(name, len(np.unique(Gaia_id))))
    
    print('possible missed triples = {}'.format(len(pos_triples)))
    
    print('of which {} are inside bind and {} are outside blind'.format(len(inside), len(outside)))
    
    print('{:.2f}% are the systems are gaia observable'.format(len(Gaia_triples)*100/len(data)))
    
        
#    print('closest triple = {:.3f}kpc, furthest triple = {:.3f}kpc'.format(np.min(d),np.max(d)))
    
#    theta = l*np.pi/180
#    n = 90 - b
#    phi = n * np.pi /180
#    print('phi min = {}'.format(np.min(phi)))
#    print('phi max = {}'.format(np.max(phi)))
#    phi = (np.pi/2) - (b*np.pi/180)
    
    
        
    l = l * np.pi/180
    b = b * np.pi/180
    
    theta = np.pi - l   
    phi = (np.pi/2) - b # to convert into the polar angle

    x = d * np.cos(theta) * np.sin(phi)
    y = d * np.sin(theta) * np.sin(phi)
    z = d * np.cos(phi)
    
    
#    print('max z = {}, min z = {}'.format(np.max(z), np.min(z)))
#    print('average d = {}'.format(np.average(d)))

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
    
#processor('stable')
#processor('LK')
#processor('ecc')


def magnitude_cut(data):
    
    G1, G2, G3 = data[:,31], data[:,34], data[:,37]
    GBP_1, GBP_2, GBP_3 = data[:,32], data[:,35], data[:,38]
    GRP_1, GRP_2, GRP_3 = data[:,33], data[:,36], data[:,39]
    semi_major, ecc, d, inc = data[:,45], data[:,47], data[:,5], data[:,49]
    star_1_id, star_2_id, star_3_id = data[:,40], data[:,41], data[:,42]
    
    index_all_obs = []
    index_Gaia_suit = []
    index_pos_trip_two = [] #two are observable
    index_pos_trip_one = [] #only one observable
    
    for i in range(len(data)):
        if ((G1[i] < Gaia_mag_res) and (G2[i] < Gaia_mag_res) and (G3[i] < Gaia_mag_res)):
            index_all_obs.append(i)
            
#            if ((1/np.pi)*2*180*3600*np.arctan(data[:,45][i]*6.96e8 / (2*data[:,47][i]*3.1e19))) > Gaia_ang_res:
            if ang_sep(semi_major[i], ecc[i], d[i], inc[i]) > s_crit(G1[i], G2[i]):
                index_Gaia_suit.append(i)
                # if the inner binary cannot be resolved, then theres a mag criteria that a triple can still be detected
            else:
                if (star_1_id[i] < 2) and (GBP_2[i] - GBP_1[i] < 1) and (GRP_2[i] - GRP_1[i] < 1): 
                    index_Gaia_suit.append(i)
                if (star_2_id[i] < 2) and (GBP_1[i] - GBP_2[i] < 1) and (GRP_1[i] - GRP_2[i] < 1):
                    index_Gaia_suit.append(i)
                    
                elif (np.abs(G1[i] - G2[i])<1):
                    index_Gaia_suit.append(i)
                
        if (G1[i] < Gaia_mag_res and G2[i] < Gaia_mag_res and G3[i] > Gaia_mag_res) or (G1[i] < Gaia_mag_res and G2[i] > Gaia_mag_res and G3[i] < Gaia_mag_res) or (G1[i] > Gaia_mag_res and G2[i] < Gaia_mag_res and G3[i] < Gaia_mag_res):
            index_pos_trip_two.append(i)
            
        if (G1[i] < Gaia_mag_res and G2[i] > Gaia_mag_res and G3[i] > Gaia_mag_res) or (G1[i] > Gaia_mag_res and G2[i] < Gaia_mag_res and G3[i] > Gaia_mag_res) or (G1[i] > Gaia_mag_res and G2[i] > Gaia_mag_res and G3[i] < Gaia_mag_res):
            index_pos_trip_one.append(i)
            
            
    index_all_obs = np.array(index_all_obs)
    index_Gaia_suit = np.array(index_Gaia_suit)
    index_pos_trip_two = np.array(index_pos_trip_two)
    index_pos_trip_one = np.array(index_pos_trip_one)
#    print('Systems where all stars observable = {}'.format(len(index_all_obs)))
#    print("number lost due to (angular) resolution = {}".format(len(index_all_obs) - len(index_Gaia_suit)))
#    print('Remaining after magnitude and resolution cuts = {}'.format(len(index_Gaia_suit)))
    
    Gaia_suit = data[index_Gaia_suit]
    pos_triples_two = data[index_pos_trip_two] # sort into whether it is inner star or outer that cannot see
    pos_triples_one = data[index_pos_trip_one]
    
    index_inside_blind = []
    index_outside_blind = []
    for i in range(len(pos_triples_two)):
        if (pos_triples_two[:,31][i] > Gaia_mag_res or pos_triples_two[:,34][i] > Gaia_mag_res):
            index_inside_blind.append(i)
        if (pos_triples_two[:,37][i] > Gaia_mag_res):
            index_outside_blind.append(i)
    index_inside_blind = np.array(index_inside_blind)
    index_outside_blind = np.array(index_outside_blind)
    inside_blind = pos_triples_two[index_inside_blind]
    outside_blind = pos_triples_two[index_outside_blind]

    return Gaia_suit, pos_triples_two, pos_triples_one, inside_blind, outside_blind, index_all_obs
    

def different_motions(data):
    G1, G2, G3 = data[:,31], data[:,34], data[:,37]
    GBP_1, GBP_2, GBP_3 = data[:,32], data[:,35], data[:,38]
    GRP_1, GRP_2, GRP_3 = data[:,33], data[:,36], data[:,39]
    semi_major, ecc, d, inc = data[:,45], data[:,47], data[:,5], data[:,49]
    star_1_id, star_2_id, star_3_id = data[:,40], data[:,41], data[:,42]


    index_cpm = []
    index_both_astrometric = []
    index_inner_astr_outer_cpm = []
        
    for i in range(len(data)):
        if (data[i][1] > Gaia_lifetime) and (data[i][2] > Gaia_lifetime):
            index_cpm.append(i)
        if (data[i][1] < Gaia_lifetime) and (data[i][2] < Gaia_lifetime):
            index_both_astrometric.append(i)
        if (data[i][1] < Gaia_lifetime) and (data[i][2] > Gaia_lifetime):
            index_inner_astr_outer_cpm.append(i)
    
    a=0
    index_cpm = np.array(index_cpm)
    if len(index_cpm) == 0 :
        cpm_systems = []
        a +=1
    else:
        cpm_systems = data[index_cpm]
        
    b=0
    index_both_astrometric = np.array(index_both_astrometric)
    if len(index_both_astrometric) == 0:
        both_astr_systems = []
        b+=1
    else:
        both_astr_systems = data[index_both_astrometric]
    
    c=0
    index_inner_astr_outer_cpm = np.array(index_inner_astr_outer_cpm)
    if len(index_inner_astr_outer_cpm) == 0:
        inn_astr_outer_cpm = []
        c+=1
    else:
        inn_astr_outer_cpm = data[index_inner_astr_outer_cpm]
    
    index_cpm_suit = []
    for i in range(len(cpm_systems)):
        theta = random.uniform(0, np.pi/2)
        if theta == 0:
            break
        else:
            l = cpm_systems[i][3] * np.pi/180
            b = cpm_systems[i][4] * np.pi/180
            d = cpm_systems[i][5] 
            
            angle = proper_resolve(l,b,d)
            if angle > prop_motion_res:
                index_cpm_suit.append(i)
    
    index_cpm_suit = np.array(index_cpm_suit)
    Gaia_cpm = cpm_systems[index_cpm_suit]
    
    all_systems = []
    
    if (a==0) and (b==0) and (c==0):
        all_systems = np.vstack((Gaia_cpm, both_astr_systems, inn_astr_outer_cpm))
    
    if (a==0) and (b==1) and (c==0):
        all_systems = np.vstack((Gaia_cpm, inn_astr_outer_cpm))
    
    if (a==0) and (b==0) and (c==1):
        all_systems = np.vstack((Gaia_cpm, both_astr_systems))
        
    if (a==0) and (b==1) and (c==1):
        all_systems = Gaia_cpm
    
    if (a==1) and (b==0) and (c==1):
        all_systems = both_astr_systems
    
    if (a==1) and (b==1) and (c==0):
        all_systems = inn_astr_outer_cpm
    
    if (a==1) and (b==1) and (c==1):
        print('You are a failure')
        all_systems = []
        
    return Gaia_cpm, both_astr_systems, inn_astr_outer_cpm, all_systems
    

def resulting__three_file(a,b,c):
    file = []
    if len(a)==0 and len(b)==0 and len(c)==0:
        file = []
    if len(a)!=0 and len(b)==0 and len(c)==0:
        file = a
    if len(a)==0 and len(b)!=0 and len(c)==0:
        file = b
    if len(a)==0 and len(b)==0 and len(c)!=0:
        file = c
    if len(a)!=0 and len(b)!=0 and len(c)==0:
        file = np.vstack((a,b))
    if len(a)!=0 and len(b)==0 and len(c)!=0:
        file = np.vstack((a,c))
    if len(a)==0 and len(b)!=0 and len(c)!=0:
        file = np.vstack((c,b))
    if len(a)!=0 and len(b)!=0 and len(c)!=0:
        file = np.vstack((a,b,c))
        
    return file

def resulting_two_file(a,b):
    file = []
    if len(a)==0 and len(b)==0:
        file = []
    if len(a)!=0 and len(b)!=0:
        file = np.vstack((a,b))
    if len(a)!=0 and len(b)==0:
        file = a
    if len(a)==0 and len(b)!=0:
        file = b

    return file

def gaia_observables(name):
#name = 'stable'
    data = np.loadtxt((open('output_{}.txt'.format(name),'rt').readlines()[:-1]), skiprows=4, dtype=None)
    print('===============\n\n')
    print('number of systems in {} = {}'.format(name, len(data)))
    print('')
    
    G1, G2, G3 = data[:,31], data[:,34], data[:,37]
    GBP_1, GBP_2, GBP_3 = data[:,32], data[:,35], data[:,38]
    GRP_1, GRP_2, GRP_3 = data[:,33], data[:,36], data[:,39]
    semi_major, ecc, d, inc = data[:,45], data[:,47], data[:,5], data[:,49]
    star_1_id, star_2_id, star_3_id = data[:,40], data[:,41], data[:,42]
    
    
    index_3_wd = []
    index_2_wd = []
    
    for i in range(len(data)):
        if (star_1_id[i] in {10,11,12} and star_2_id[i] in {10,11,12} and star_3_id[i] in {10,11,12}):
            index_3_wd.append(i)
        
        if (star_1_id[i] in {10,11,12} and star_2_id[i] in {10,11,12} and star_3_id[i] in {0,1,2}) or \
            (star_1_id[i] in {10,11,12} and star_3_id[i] in {10,11,12} and star_2_id[i] in {0,1,2}) or \
            (star_3_id[i] in {10,11,12} and star_2_id[i] in {10,11,12} and star_1_id[i] in {0,1,2}):
                index_2_wd.append(i)
           
    index_3_wd = np.array(index_3_wd)
    index_2_wd = np.array(index_2_wd)
    
    all_3_wd_cpm = []
    all_3_astr = [] # Upper limit as no cut made, just period sorting
    all_3_inn_astr_out_cpm = []
    
    if len(index_3_wd) == 0:
        print('All WDs systems = 0')
#        Gaia_3_wd = np.zeros(50)
        Gaia_3_wd = []
    else:
        all_wds = data[index_3_wd]
        print('All WDs systems = {}'.format(len(all_wds)))
    
        print('for all mag observable')
        all_3_wd_obs = magnitude_cut(all_wds)[0]
        print('Systems where all stars observable = {}'.format(len(magnitude_cut(all_wds)[5])))
        print("number lost due to (angular) resolution = {}".format(len(magnitude_cut(all_wds)[5]) - len(magnitude_cut(all_wds)[0])))
        print('Remaining after magnitude and resolution cuts = {}'.format(len(magnitude_cut(all_wds)[0])))
        
        print('')
        print('for a possible triple')
        all_3_wd_pos_two = magnitude_cut(all_wds)[1]
        all_3_wd_pos_one = magnitude_cut(all_wds)[2]
        all_3_inside_blind = magnitude_cut(all_wds)[3]
        all_3_outside_blind = magnitude_cut(all_wds)[4]
        
    
        print('All three wds mag observable = {}'.format(len(all_3_wd_obs)))
        print('All three wds possible triple = {}'.format(len(all_3_wd_pos_two)+len(all_3_wd_pos_one)))
        print('Of which 2 wds visible = {}, only one visible {}'.format(len(all_3_wd_pos_two), len(all_3_wd_pos_one)))
        print('For 2 visible, inside blind = {}, outside blind = {}'.format(len(all_3_inside_blind), len(all_3_outside_blind)))
        
        #Sort into different motion types
        all_3_wd_cpm = different_motions(all_3_wd_obs)[0]
        print('Number pcpm = {}'.format(len(all_3_wd_cpm)))
        all_3_astr = different_motions(all_3_wd_obs)[1] # Upper limit as no cut made, just period sorting
        print('Number both astr = {}'.format(len(all_3_astr)))
        all_3_inn_astr_out_cpm = different_motions(all_3_wd_obs)[2] # Upper limit...
        print('Number inner astr, outer cpm = {}'.format(len(all_3_inn_astr_out_cpm)))
        
#        Gaia_3_wd = np.vstack((all_3_wd_cpm, all_3_astr, all_3_inn_astr_out_cpm))
#        Gaia_3_wd = different_motions(all_3_wd_obs)[3]
        Gaia_3_wd = resulting__three_file(all_3_wd_cpm,all_3_astr,all_3_inn_astr_out_cpm)
        
        
    print('')   
    
    
    if len(index_2_wd) == 0:
        print('Two WDs, one MS = 0')
#        Gaia_2_wd = np.zeros(50)
        Gaia_2_wd = []
    else:
        two_wds = data[index_2_wd]
        print('Two WDs, one MS = {}'.format(len(two_wds)))
        
        print('for all mag observable')
        two_wd_obs = magnitude_cut(two_wds)[0]
        print('Systems where all stars observable = {}'.format(len(magnitude_cut(two_wds)[5])))
        print("number lost due to (angular) resolution = {}".format(len(magnitude_cut(two_wds)[5]) - len(magnitude_cut(two_wds)[0])))
        print('Remaining after magnitude and resolution cuts = {}'.format(len(magnitude_cut(two_wds)[0])))
        
        print('')
        print('for a possible triple')
        two_wd_pos_two = magnitude_cut(two_wds)[1]
        two_wd_pos_one = magnitude_cut(two_wds)[2]
        two_wd_inside_blind = magnitude_cut(two_wds)[3]
        two_wd_outside_blind = magnitude_cut(two_wds)[4]
    
        print('Two WDs mag observable = {}'.format(len(two_wd_obs)))
        print('Two WDs possible triple = {}'.format(len(two_wd_pos_two)+len(two_wd_pos_one)))
        print('Of which 2 are visible = {}, only one visible {}'.format(len(two_wd_pos_two), len(two_wd_pos_one)))
        print('For 2 visible, inside blind = {}, outside blind = {}'.format(len(two_wd_inside_blind), len(two_wd_outside_blind)))

        two_wd_wd_cpm = different_motions(two_wd_obs)[0]
        print('Number pcpm = {}'.format(len(two_wd_wd_cpm)))
        two_wd_astr = different_motions(two_wd_obs)[1]
        print('Number both astr = {}'.format(len(two_wd_astr)))
        two_wd_inn_astr_out_cpm = different_motions(two_wd_obs)[2]
        print('Number inner astr, outer cpm = {}'.format(len(two_wd_inn_astr_out_cpm)))
        
#        Gaia_2_wd = np.vstack((two_wd_wd_cpm, two_wd_astr, two_wd_inn_astr_out_cpm))
#        Gaia_2_wd = different_motions(two_wd_obs)[3]
        Gaia_2_wd = resulting__three_file(two_wd_wd_cpm,two_wd_astr,two_wd_inn_astr_out_cpm)
    
#    Gaia_systems = np.vstack((Gaia_3_wd, Gaia_2_wd))
    Gaia_systems = resulting_two_file(Gaia_3_wd, Gaia_2_wd)
    
    print('')
    print('number of Gaia observables is {}'.format(len(Gaia_systems)))
    print('unique systems in {} = {}'.format(name, len(np.unique(Gaia_systems[:,44]))))
    
    
    
    l = Gaia_systems[:,3] # longitude
    b = Gaia_systems[:,4] # latitude - 0 means in horizontal plane
    d = Gaia_systems[:,5] # r in kpc
    
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
    ax.plot([], [], [], ' ', label='Num Triples = {}'.format(len(Gaia_systems)))
    ax.set_xlabel('x axis (kpc)')
    ax.set_ylabel('y axis (kpc)')
    ax.set_zlabel('z axis (kpc)')
    ax.legend()
    

    plt.show()
    
    from scipy.stats import norm
    
    fig = plt.figure()
    num_bins = 70
    n, bins, patches = plt.hist(d, num_bins, facecolor='b', alpha=0.4, edgecolor='black')#, num_bins, density=True, facecolor='b', alpha=0.4, edgecolor='black')

    plt.xlabel('d (kpc)')
    plt.ylabel('Frequency')
    plt.title('{}'.format(name))
    
    fig = plt.figure()
    num_bins = 70
    n, bins, patches = plt.hist(z, num_bins, facecolor='b', alpha=0.4, edgecolor='black')#, num_bins, density=True, facecolor='b', alpha=0.4, edgecolor='black')


    plt.xlabel('z (kpc)')
    plt.ylabel('Frequency')
    plt.title('{}'.format(name))
    plt.subplots_adjust(left=0.15)
    plt.show()
    print('')
    print('')
    print('')
    
    
    return Gaia_systems

    



#gaia_observables('stable')
#gaia_observables('LK')
#gaia_observables('ecc')

All_Gaia_Observables = np.vstack((gaia_observables('stable'), gaia_observables('LK'), gaia_observables('ecc')))

l = All_Gaia_Observables[:,3] # longitude
b = All_Gaia_Observables[:,4] # latitude - 0 means in horizontal plane
d = All_Gaia_Observables[:,5] # r in kpc

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
ax.plot([], [], [], ' ', label='Num Overall Triples = {}'.format(len(All_Gaia_Observables)))
ax.set_xlabel('x axis (kpc)')
ax.set_ylabel('y axis (kpc)')
ax.set_zlabel('z axis (kpc)')
ax.legend()


plt.show()








