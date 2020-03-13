#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:34:36 2019

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

#np.set_printoptions(threshold=np.inf)

G = const.G
M_sol = const.M_sun
R_sol = const.R_sun
c = const.c
hubble_t = 13.5e9

f_crit = 1
eps_oct_lower = 0.001

Data = np.loadtxt('data.txt')
Evolutions = len(Data)
Number_of_triples = Evolutions / 5

print(Evolutions, Number_of_triples)

Triples = Data.reshape(int(Number_of_triples), 5, 31) #formatted as [[triple 0], [triple 1].......[triple n]]
             #[triple n] = [[evolution 0], ... , [evolution 4]]
             #[evolution 0] = [col 0, col 1, ... , col 30]

def round_down(num, div):
    return num - (num % div)

index_first = np.unique(round_down(np.arange(Evolutions), 5))

inc_0 = []
for i in index_first:
    non_zero = np.unique(np.nonzero(Data[i:i+4])[0])
    if 3 in non_zero:
        inc_0.append(Data[i+3, 1])
    else:
        last = np.max(non_zero)
        inc_0.append(Data[last + i, 1])
        
inc_0 = np.array(inc_0)


#inc_0 = []
#for x in range(int(Number_of_triples)):
#    for y in range(4):
#        if Triples[x][y][1] == 0:
#            inc_0.append(Triples[x][y - 1][1])
#            break
#    else: 
#        if Triples[x][3][1] > 0:
#            inc_0.append(Triples[x][3][1])
#
#inc_0 = np.array(inc_0)


a_out = Data[:,6]
a_in = Data[:,2]
e_out = Data[:,7]
e_in = Data[:,3]
inc = Data[:,1]
inc_0 =  np.repeat(inc_0, 5)   #intial mutual inclination
m_1 = Data[:,10]
m_2 = Data[:,11]
m_3 = Data[:,12]




def Q_out(m_1, m_2, m_3):
    return np.nan_to_num(m_3 / (m_1 + m_2))

def P_out(m_1, m_2, a_out):
    return np.nan_to_num(2 * np.pi * np.sqrt((R_sol * a_out)**3 
                                          / (G * M_sol * (m_1 + m_2))))

def P_in(m_1, m_2, a_in):
    return np.nan_to_num(2 * np.pi * np.sqrt((R_sol * a_in)**3 
                                          / (G * M_sol * (m_1 + m_2))))


q_out = Q_out(m_1, m_2, m_3)
P_out1 = P_out(m_1, m_2, a_out)
P_in1 = P_in(m_1, m_2, a_in)


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

crit_stab = crit_stability(e_out, inc, q_out)

stab = stability(a_out, a_in)

t_k = t_kozai(P_out1,P_in1,m_1,m_2,m_3,e_out)

e_max = e_max(inc_0)

eps_oct = eps_oct(m_1, m_2, a_in, a_out, e_out)

orb_av_lim = orb_av_lim(q_out, a_in, a_out, e_out)

semi_sec = semi_sec(e_max)


'''
==============================================================================
'''

#INDEXING and REGIME CREATION

index_stable = np.arange(len(Data))[(stab > crit_stab)]

stable_evolutions = Data[index_stable]

print('Total number in stable = {}'.format(len(index_stable)))


index_semisec = np.arange(len(Data))[(semi_sec <= orb_av_lim)]

index_stable_semisec = np.intersect1d(index_semisec, index_stable)

stable_semisec_evolutions = Data[index_stable_semisec]

print('Total number in semi sec = {}'.format(len(index_stable_semisec)))


index_eccentric = np.arange(len(Data))[(np.abs(eps_oct) > eps_oct_lower)]

index_stable_eccentric = np.intersect1d(index_eccentric, index_stable)

index_eccentric_notsec = np.setdiff1d(index_stable_eccentric, index_stable_semisec)

stable_eccentric_evolutions = Data[index_eccentric_notsec]

print('Total number in eccentric = {}'.format(len(index_eccentric_notsec)))

    
index_LK = np.arange(len(Data))[(39.2 * np.pi / 180 < inc_0) & (inc_0 < 140.8 * np.pi / 180)]

index_stable_LK = np.intersect1d(index_LK, index_stable)

index_LK_notsec = np.setdiff1d(index_stable_LK, index_stable_semisec)

index_LK_notsec_notecc = np.setdiff1d(index_LK_notsec, index_stable_eccentric)

stable_LK_evolutions = Data[index_LK_notsec_notecc]

print('Total bumber in LK = {}'.format(len(index_LK_notsec_notecc)))

'''
==============================================================================
'''

index_first = np.unique(round_down(np.arange(Evolutions), 5))

index_last = []
for i in index_first:
    non_zero = np.unique(np.nonzero(Data[i:i+4])[0]) #goes to 3rd one
    last = np.max(non_zero)
    if Data[i + last][22] in {10,11,12} and Data[i + last][23] in {10,11,12} and Data[i + last][24] in ({0,1,2,10,11,12}) or \
        Data[i + last][22] in {10,11,12} and Data[i + last][24] in {10,11,12} and Data[i + last][23] in ({0,1,2,10,11,12})or \
        Data[i + last][23] in {10,11,12} and Data[i + last][24] in {10,11,12} and Data[i + last][22] in ({0,1,2,10,11,12}):
            index_last.append(i + last)

index_last = np.array(index_last)

index_last_stable = index_last[(stab[index_last] > crit_stab[index_last])]

index_last_semi = index_last_stable[(semi_sec[index_last_stable] <= orb_av_lim[index_last_stable])]

index_last_ecc = index_last_stable[(np.abs(eps_oct[index_last_stable]) > eps_oct_lower)]
index_eccentric_notsec = np.setdiff1d(index_last_ecc, index_last_semi)

index_last_LK = index_last_stable[(39.2 * np.pi / 180 < inc_0[index_last_stable]) & (inc_0[index_last_stable] < 140.8 * np.pi / 180)]
index_LK_notsec = np.setdiff1d(index_last_LK, index_last_semi)
index_LK_notsec_not_ecc = np.setdiff1d(index_LK_notsec, index_eccentric_notsec)


print('Number in stable = {}'.format(len(index_last_stable)))
print('Number in semi = {}'.format(len(index_last_semi)))
print('Number in ecc = {}'.format(len(index_eccentric_notsec)))
print('Number in LK = {}'.format(len(index_LK_notsec_not_ecc)))


evolutions_stable = Data[index_last_stable]

times_stable = []
for i in round_down(index_last_stable, 5):
    row = []
    for j in range(1,4):
        if len(row) == 2:
            break
        
        if Data[i+j, 28] != 0:
            row.append(Data[i+j, 28])
    times_stable.append(row)
    
times_stable = np.array(times_stable)

stable_final = []
for i in range(len(index_last_stable)):
    stable_final.append(np.append(evolutions_stable[i], times_stable[i]))

stable_final = np.array(stable_final)



evolutions_semi = Data[index_last_semi]

times_semi = []
for i in round_down(index_last_semi, 5):
    row = []
    for j in range(1,4):
        if len(row) == 2:
            break
        
        if Data[i+j, 28] != 0:
            row.append(Data[i+j, 28])
    times_semi.append(row)
    
times_semi = np.array(times_semi)

semi_final = []
for i in range(len(index_last_semi)):
    semi_final.append(np.append(evolutions_semi[i], times_semi[i]))

semi_final = np.array(semi_final)


evolutions_ecc = Data[index_eccentric_notsec]

times_ecc = []
for i in round_down(index_eccentric_notsec, 5):
    row = []
    for j in range(1,4):
        if len(row) == 2:
            break
        
        if Data[i+j, 28] != 0:
            row.append(Data[i+j, 28])
    times_ecc.append(row)
    
times_ecc = np.array(times_ecc)

ecc_final = []
for i in range(len(index_eccentric_notsec)):
    ecc_final.append(np.append(evolutions_ecc[i], times_ecc[i]))

ecc_final = np.array(ecc_final)


evolutions_LK = Data[index_LK_notsec_not_ecc]

times_LK = []
for i in round_down(index_LK_notsec_not_ecc, 5):
    row = []
    for j in range(1,4):
        if len(row) == 2:
            break
        
        if Data[i+j, 28] != 0:
            row.append(Data[i+j, 28])
    times_LK.append(row)
    
times_LK = np.array(times_LK)

LK_final = []
for i in range(len(index_LK_notsec_not_ecc)):
    LK_final.append(np.append(evolutions_LK[i], times_LK[i]))

LK_final = np.array(LK_final)

'''
==============================================================================
'''








def save(name, array):
    np.savetxt(name, array)
    
#stable_seismic_evolutions = Data[index_last_semisec]
#stable_evolutions = Data[index_last_stable]
#stable_LK_evolutions = Data[index_last_LK]
#stable_eccentric_evolutions = Data[index_last_eccentric]

#save('semi_sec.txt', stable_semisec_evolutions)
save('stable.txt', stable_final)
save('LK.txt', LK_final)
save('eccentric.txt', ecc_final)

stab = np.loadtxt('stable.txt')
LK = np.loadtxt('LK.txt')
ecc = np.loadtxt('eccentric.txt')

print('shapes are {}, {}, {}'.format(np.shape(stab), np.shape(LK), np.shape(ecc)))

with open('stable.txt', 'r') as program:
    data = program.readlines()

with open('stable.txt', 'w') as program:
    for (number, line) in enumerate(data):
        program.write('%d  %s' % (number, line))
        
with open('LK.txt', 'r') as program:
    data = program.readlines()

with open('LK.txt', 'w') as program:
    for (number, line) in enumerate(data):
        program.write('%d  %s' % (number, line))

with open('eccentric.txt', 'r') as program:
    data = program.readlines()

with open('eccentric.txt', 'w') as program:
    for (number, line) in enumerate(data):
        program.write('%d  %s' % (number, line))
        
#stab = np.loadtxt('stable.txt')
#LK = np.loadtxt('LK.txt')
#ecc = np.loadtxt('eccentric.txt')
#
#print('shapes are {}, {}, {}'.format(np.shape(stab), np.shape(LK), np.shape(ecc)))



#subprocess.call(['./get_outputs.sh'])



'''
==============================================================================
'''
#
##Make graph of semisec first evolution
#
#m_1_semi = np.average(Data[index_first_semisec,10])
#m_2_semi = np.average(Data[index_first_semisec,11])
#m_3_semi = np.average(Data[index_first_semisec,12])
#
#m_1_ecc = np.average(Data[index_first_eccentric,10])
#m_2_ecc = np.average(Data[index_first_eccentric,11])
#m_3_ecc = np.average(Data[index_first_eccentric,12])
#
#m_1_stab = np.average(Data[index_first_stable,10])
#m_2_stab = np.average(Data[index_first_stable,11])
#m_3_stab = np.average(Data[index_first_stable,12])
#
#
#inc_g = 0
#e_ing = 0.99
#eps_octg = 0.001
#
#q_out_semi = m_3_semi / (m_1_semi + m_2_semi)
#q_out_stab = m_3_stab / (m_1_stab + m_2_stab)
#
#e_outg = np.linspace(0, 1, 1001)
#
#stab_line = (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / np.sqrt(1 - e_outg))**(2/5)
#
#semisecular_line = (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3)
#
#eccentric_line = (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#def f(e_outg):
#    return (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / 
#                    np.sqrt(1 - e_outg))**(2/5) - (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi /
#                           np.sqrt(1 - e_ing))**(1/3)
#
#def g(e_outg):
#    return (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3) - (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#
#root1 = optimize.brentq(f, 0, 0.999, maxiter=100)
#root2 = optimize.brentq(g, 0, 0.999, maxiter=2000)
#
#y = np.nan_to_num(Data[index_first_semisec,6] / Data[index_first_semisec,2])
#x = Data[index_first_semisec,7]
#y = y[x != 0]
#x = x[x != 0]
#y = y[x <= 1]
#x = x[x <= 1]
#
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(e_outg, stab_line)
#ax.plot(e_outg[e_outg < root1], semisecular_line[e_outg < root1])
#ax.plot(e_outg[e_outg > root2], eccentric_line[e_outg > root2])
#
#ax.scatter(x, y, c=z, cmap='Greys', s=50, edgecolor='')
#
#ax.set_yscale('log')
#ax.set_title('Semisec First Evolution')
#ax.text(0.2, 0.95, 'Evolutions = {}'.format(len(index_first_semisec)), horizontalalignment='center', 
#     verticalalignment='center', transform=ax.transAxes)
#plt.show()
#
#'''
#==============================================================================
#'''
#
##Make graph of semisec last evolution
#
#m_1_semi = np.median(Data[index_last_semisec,10])
#m_2_semi = np.median(Data[index_last_semisec,11])
#m_3_semi = np.median(Data[index_last_semisec,12])
#
#m_1_ecc = np.median(Data[index_last_eccentric,10])
#m_2_ecc = np.median(Data[index_last_eccentric,11])
#m_3_ecc = np.median(Data[index_last_eccentric,12])
#
#m_1_stab = np.median(Data[index_last_stable,10])
#m_2_stab = np.median(Data[index_last_stable,11])
#m_3_stab = np.median(Data[index_last_stable,12])
#
#
#inc_g = 0
#e_ing = 0.99
#eps_octg = 0.001
#
#q_out_semi = m_3_semi / (m_1_semi + m_2_semi)
#q_out_stab = m_3_stab / (m_1_stab + m_2_stab)
#
#e_outg = np.linspace(0, 1, 1001)
#
#stab_line = (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / np.sqrt(1 - e_outg))**(2/5)
#
#semisecular_line = (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3)
#
#eccentric_line = (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#def f(e_outg):
#    return (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / 
#                    np.sqrt(1 - e_outg))**(2/5) - (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi /
#                           np.sqrt(1 - e_ing))**(1/3)
#
#def g(e_outg):
#    return (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3) - (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#
#root1 = optimize.brentq(f, 0, 0.999, maxiter=100)
#root2 = optimize.brentq(g, 0, 0.999, maxiter=2000)
#
#y = np.nan_to_num(Data[index_last_semisec,6] / Data[index_last_semisec,2])
#x = Data[index_last_semisec,7]
#y = y[x != 0]
#x = x[x != 0]
#y = y[x <= 1]
#x = x[x <= 1]
#
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(e_outg, stab_line)
#ax.plot(e_outg[e_outg < root1], semisecular_line[e_outg < root1])
#ax.plot(e_outg[e_outg > root2], eccentric_line[e_outg > root2])
#
#ax.scatter(x, y, c=z, cmap='Greys', s=50, edgecolor='')
#
#ax.set_yscale('log')
#ax.set_title('Semisec Last Evolution')
#ax.text(0.2, 0.95, 'Evolutions = {}'.format(len(index_last_semisec)), horizontalalignment='center', 
#     verticalalignment='center', transform=ax.transAxes)
#plt.show()
#
#'''
#==============================================================================
#'''
#
##Make graph of LK first evolution
#
#m_1_semi = np.average(Data[index_first_semisec,10])
#m_2_semi = np.average(Data[index_first_semisec,11])
#m_3_semi = np.average(Data[index_first_semisec,12])
#
#m_1_ecc = np.average(Data[index_first_eccentric,10])
#m_2_ecc = np.average(Data[index_first_eccentric,11])
#m_3_ecc = np.average(Data[index_first_eccentric,12])
#
#m_1_stab = np.average(Data[index_first_stable,10])
#m_2_stab = np.average(Data[index_first_stable,11])
#m_3_stab = np.average(Data[index_first_stable,12])
#
#
#inc_g = 0
#e_ing = 0.99
#eps_octg = 0.001
#
#q_out_semi = m_3_semi / (m_1_semi + m_2_semi)
#q_out_stab = m_3_stab / (m_1_stab + m_2_stab)
#
#e_outg = np.linspace(0, 1, 1001)
#
#stab_line = (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / np.sqrt(1 - e_outg))**(2/5)
#
#semisecular_line = (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3)
#
#eccentric_line = (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#def f(e_outg):
#    return (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / 
#                    np.sqrt(1 - e_outg))**(2/5) - (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi /
#                           np.sqrt(1 - e_ing))**(1/3)
#
#def g(e_outg):
#    return (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3) - (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#
#root1 = optimize.brentq(f, 0, 0.999, maxiter=100)
#root2 = optimize.brentq(g, 0, 0.999, maxiter=2000)
#
#y = np.nan_to_num(Data[index_first_LK,6] / Data[index_first_LK,2])
#x = Data[index_first_LK,7]
#y = y[x != 0]
#x = x[x != 0]
#y = y[x <= 1]
#x = x[x <= 1]
#
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(e_outg, stab_line)
#ax.plot(e_outg[e_outg < root1], semisecular_line[e_outg < root1])
#ax.plot(e_outg[e_outg > root2], eccentric_line[e_outg > root2])
#
#ax.scatter(x, y, c=z, cmap='Greys', s=50, edgecolor='')
#
#ax.set_yscale('log')
#ax.set_title('LK First Evolution')
#ax.text(0.2, 0.95, 'Evolutions = {}'.format(len(index_first_LK)), horizontalalignment='center', 
#     verticalalignment='center', transform=ax.transAxes)
#plt.show()
#
#'''
#==============================================================================
#'''
#
##Make graph of LK last evolution
#
#m_1_semi = np.median(Data[index_last_semisec,10])
#m_2_semi = np.median(Data[index_last_semisec,11])
#m_3_semi = np.median(Data[index_last_semisec,12])
#
#m_1_ecc = np.median(Data[index_last_eccentric,10])
#m_2_ecc = np.median(Data[index_last_eccentric,11])
#m_3_ecc = np.median(Data[index_last_eccentric,12])
#
#m_1_stab = np.median(Data[index_last_stable,10])
#m_2_stab = np.median(Data[index_last_stable,11])
#m_3_stab = np.median(Data[index_last_stable,12])
#
#
#inc_g = 0
#e_ing = 0.99
#eps_octg = 0.001
#
#q_out_semi = m_3_semi / (m_1_semi + m_2_semi)
#q_out_stab = m_3_stab / (m_1_stab + m_2_stab)
#
#e_outg = np.linspace(0, 1, 1001)
#
#stab_line = (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / np.sqrt(1 - e_outg))**(2/5)
#
#semisecular_line = (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3)
#
#eccentric_line = (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#def f(e_outg):
#    return (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / 
#                    np.sqrt(1 - e_outg))**(2/5) - (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi /
#                           np.sqrt(1 - e_ing))**(1/3)
#
#def g(e_outg):
#    return (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3) - (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#
#root1 = optimize.brentq(f, 0, 0.999, maxiter=100)
#root2 = optimize.brentq(g, 0, 0.999, maxiter=2000)
#
#y = np.nan_to_num(Data[index_last_LK,6] / Data[index_last_LK,2])
#x = Data[index_last_LK,7]
#y = y[x != 0]
#x = x[x != 0]
#y = y[x <= 1]
#x = x[x <= 1]
#
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(e_outg, stab_line)
#ax.plot(e_outg[e_outg < root1], semisecular_line[e_outg < root1])
#ax.plot(e_outg[e_outg > root2], eccentric_line[e_outg > root2])
#
#ax.scatter(x, y, c=z, cmap='Greys', s=50, edgecolor='')
#
#ax.set_yscale('log')
#ax.set_title('LK Last Evolution')
#ax.text(0.2, 0.95, 'Evolutions = {}'.format(len(index_last_LK)), horizontalalignment='center', 
#     verticalalignment='center', transform=ax.transAxes)
#plt.show()
#
#'''
#==============================================================================
#'''
#
##Make graph of eccentric first evolution
#
#m_1_semi = np.average(Data[index_first_semisec,10])
#m_2_semi = np.average(Data[index_first_semisec,11])
#m_3_semi = np.average(Data[index_first_semisec,12])
#
#m_1_ecc = np.average(Data[index_first_eccentric,10])
#m_2_ecc = np.average(Data[index_first_eccentric,11])
#m_3_ecc = np.average(Data[index_first_eccentric,12])
#
#m_1_stab = np.average(Data[index_first_stable,10])
#m_2_stab = np.average(Data[index_first_stable,11])
#m_3_stab = np.average(Data[index_first_stable,12])
#
#
#inc_g = 0
#e_ing = 0.99
#eps_octg = 0.001
#
#q_out_semi = m_3_semi / (m_1_semi + m_2_semi)
#q_out_stab = m_3_stab / (m_1_stab + m_2_stab)
#
#e_outg = np.linspace(0, 1, 1001)
#
#stab_line = (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / np.sqrt(1 - e_outg))**(2/5)
#
#semisecular_line = (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3)
#
#eccentric_line = (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#def f(e_outg):
#    return (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / 
#                    np.sqrt(1 - e_outg))**(2/5) - (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi /
#                           np.sqrt(1 - e_ing))**(1/3)
#
#def g(e_outg):
#    return (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3) - (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#
#root1 = optimize.brentq(f, 0, 0.999, maxiter=100)
#root2 = optimize.brentq(g, 0, 0.999, maxiter=2000)
#
#y = np.nan_to_num(Data[index_first_eccentric,6] / Data[index_first_eccentric,2])
#x = Data[index_first_eccentric,7]
#y = y[x != 0]
#x = x[x != 0]
#y = y[x <= 1]
#x = x[x <= 1]
#
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(e_outg, stab_line)
#ax.plot(e_outg[e_outg < root1], semisecular_line[e_outg < root1])
#ax.plot(e_outg[e_outg > root2], eccentric_line[e_outg > root2])
#
#ax.scatter(x, y, c=z, cmap='Greys', s=50, edgecolor='')
#
#ax.set_yscale('log')
#ax.set_title('Eccentric First Evolution')
#ax.text(0.2, 0.95, 'Evolutions = {}'.format(len(index_first_eccentric)), horizontalalignment='center', 
#     verticalalignment='center', transform=ax.transAxes)
#plt.show()
#
#'''
#==============================================================================
#'''
#
##Make graph of LK last evolution
#
#m_1_semi = np.median(Data[index_last_semisec,10])
#m_2_semi = np.median(Data[index_last_semisec,11])
#m_3_semi = np.median(Data[index_last_semisec,12])
#
#m_1_ecc = np.median(Data[index_last_eccentric,10])
#m_2_ecc = np.median(Data[index_last_eccentric,11])
#m_3_ecc = np.median(Data[index_last_eccentric,12])
#
#m_1_stab = np.median(Data[index_last_stable,10])
#m_2_stab = np.median(Data[index_last_stable,11])
#m_3_stab = np.median(Data[index_last_stable,12])
#
#
#inc_g = 0
#e_ing = 0.99
#eps_octg = 0.001
#
#q_out_semi = m_3_semi / (m_1_semi + m_2_semi)
#q_out_stab = m_3_stab / (m_1_stab + m_2_stab)
#
#e_outg = np.linspace(0, 1, 1001)
#
#stab_line = (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / np.sqrt(1 - e_outg))**(2/5)
#
#semisecular_line = (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3)
#
#eccentric_line = (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#def f(e_outg):
#    return (2.8 / (1 - e_outg)) * ( 1 - 0.3 * inc_g / 
#           np.pi) * ((1 + q_out_stab) * (1 + e_outg) / 
#                    np.sqrt(1 - e_outg))**(2/5) - (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi /
#                           np.sqrt(1 - e_ing))**(1/3)
#
#def g(e_outg):
#    return (1 / (1 - e_outg)) * (f_crit * 5 * np.pi * q_out_semi / 
#                   np.sqrt(1 - e_ing))**(1/3) - (m_1_ecc - m_2_ecc) / (m_1_ecc + m_2_ecc) * (1 / eps_octg) * e_outg / (1 - e_outg**2)
#
#
#root1 = optimize.brentq(f, 0, 0.999, maxiter=100)
#root2 = optimize.brentq(g, 0, 0.999, maxiter=2000)
#
#y = np.nan_to_num(Data[index_last_eccentric,6] / Data[index_last_eccentric,2])
#x = Data[index_last_eccentric,7]
#y = y[x != 0]
#x = x[x != 0]
#y = y[x <= 1]
#x = x[x <= 1]
#
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(e_outg, stab_line)
#ax.plot(e_outg[e_outg < root1], semisecular_line[e_outg < root1])
#ax.plot(e_outg[e_outg > root2], eccentric_line[e_outg > root2])
#
#ax.scatter(x, y, c=z, cmap='Greys', s=50, edgecolor='')
#
#ax.set_yscale('log')
#ax.set_title('Eccentric Last Evolution')
#ax.text(0.2, 0.95, 'Evolutions = {}'.format(len(index_last_eccentric)), horizontalalignment='center', 
#     verticalalignment='center', transform=ax.transAxes)
#plt.show()
#
#'''
#==============================================================================
#'''
#
##Making histograms for the masses to choose between medians or means
#
#
#bins = np.linspace(0.05, 8, 35)
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_title('First semisec masses')
#ax.hist(Data[index_first_semisec,10], bins, alpha=0.5, edgecolor='k', label='m1')
#ax.hist(Data[index_first_semisec,11], bins, alpha=0.5, label='m2')
#ax.hist(Data[index_first_semisec,12], bins, alpha=0.5, label='m3')
#ax.legend()
#plt.show()
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_title('First eccentric masses')
#ax.hist(Data[index_first_eccentric,10], bins, alpha=0.5, edgecolor='k', label='m1')
#ax.hist(Data[index_first_eccentric,11], bins, alpha=0.5, label='m2')
#ax.hist(Data[index_first_eccentric,12], bins, alpha=0.5, label='m3')
#ax.legend()
#plt.show()
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_title('First LK masses')
#ax.hist(Data[index_first_LK,10], bins, alpha=0.5, edgecolor='k', label='m1')
#ax.hist(Data[index_first_LK,11], bins, alpha=0.5, label='m2')
#ax.hist(Data[index_first_LK,12], bins, alpha=0.5, label='m3')
#ax.legend()
#plt.show()
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_title('Last semisec masses')
#ax.hist(Data[index_last_semisec,10], bins, alpha=0.5, edgecolor='k', label='m1')
#ax.hist(Data[index_last_semisec,11], bins, alpha=0.5, label='m2')
#ax.hist(Data[index_last_semisec,12], bins, alpha=0.5, label='m3')
#ax.legend()
#plt.show()
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_title('Last eccentric masses')
#ax.hist(Data[index_last_eccentric,10], bins, alpha=0.5, edgecolor='k', label='m1')
#ax.hist(Data[index_last_eccentric,11], bins, alpha=0.5, label='m2')
#ax.hist(Data[index_last_eccentric,12], bins, alpha=0.5, label='m3')
#ax.legend()
#plt.show()
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_title('Last LK masses')
#ax.hist(Data[index_last_LK,10], bins, alpha=0.5, edgecolor='k', label='m1')
#ax.hist(Data[index_last_LK,11], bins, alpha=0.5, label='m2')
#ax.hist(Data[index_last_LK,12], bins, alpha=0.5, label='m3')
#ax.legend()
#plt.show()
