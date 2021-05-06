# Importation des packages
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import csv
import sys
csv.field_size_limit(sys.maxsize)


exp = 'twomey-'

# lire data
with open(exp+'C-3500-k-09-data.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile)
    i = 0
    for temp in spamreader:
        row = []
        for item in temp:
            row.append(float(item))
        if i == 0:
            dt_list_0 = row
        elif i == 1:
            pres_vapsat_0 = row
        elif i == 2:
            S_0 = row
        elif i == 3:
            pres_vap_0 = row
        elif i == 4:
            masse_condensation_0 = row
        elif i == 5:
            masse_activation_0 = row
        elif i == 6:
            masse_precipitation_0 = row
        elif i == 7:
            N_tot_0 = row
        elif i == 8:
            rayon_0 = row

        i = i+1


with open(exp+'C-3500-k-08-data.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile)
    i = 0
    for temp in spamreader:
        row = []
        for item in temp:
            row.append(float(item))
        if i == 0:
            dt_list_1 = row
        elif i == 1:
            pres_vapsat_1 = row
        elif i == 2:
            S_1 = row
        elif i == 3:
            pres_vap_1 = row
        elif i == 4:
            masse_condensation_1 = row
        elif i == 5:
            masse_activation_1 = row
        elif i == 6:
            masse_precipitation_1 = row
        elif i == 7:
            N_tot_1 = row
        elif i == 8:
            rayon_1 = row

        i = i+1

with open(exp+'C-600-k-05-data.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile)
    i = 0
    for temp in spamreader:
        row = []
        for item in temp:
            row.append(float(item))
        if i == 0:
            dt_list_2 = row
        elif i == 1:
            pres_vapsat_2 = row
        elif i == 2:
            S_2 = row
        elif i == 3:
            pres_vap_2 = row
        elif i == 4:
            masse_condensation_2 = row
        elif i == 5:
            masse_activation_2 = row
        elif i == 6:
            masse_precipitation_2 = row
        elif i == 7:
            N_tot_2 = row
        elif i == 8:
            rayon_2 = row

        i = i+1


with open(exp+'C-250-k-05-data.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile)
    i = 0
    for temp in spamreader:
        row = []
        for item in temp:
            row.append(float(item))
        if i == 0:
            dt_list_3 = row
        elif i == 1:
            pres_vapsat_3 = row
        elif i == 2:
            S_3 = row
        elif i == 3:
            pres_vap_3 = row
        elif i == 4:
            masse_condensation_3 = row
        elif i == 5:
            masse_activation_3 = row
        elif i == 6:
            masse_precipitation_3 = row
        elif i == 7:
            N_tot_3 = row
        elif i == 8:
            rayon_3 = row

        i = i+1

with open(exp+'C-53-k-05-data.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile)
    i = 0
    for temp in spamreader:
        row = []
        for item in temp:
            row.append(float(item))
        if i == 0:
            dt_list_4 = row
        elif i == 1:
            pres_vapsat_4 = row
        elif i == 2:
            S_4 = row
        elif i == 3:
            pres_vap_4 = row
        elif i == 4:
            masse_condensation_4 = row
        elif i == 5:
            masse_activation_4 = row
        elif i == 6:
            masse_precipitation_4 = row
        elif i == 7:
            N_tot_4 = row
        elif i == 8:
            rayon_4 = row

        i = i+1

# Graphiques
#
plt.figure(1)
plt.plot(np.array(dt_list_0)/3600,pres_vapsat_0)
plt.plot(np.array(dt_list_1)/3600,pres_vapsat_1)
plt.plot(np.array(dt_list_2)/3600,pres_vapsat_2)
plt.plot(np.array(dt_list_3)/3600,pres_vapsat_3)
plt.plot(np.array(dt_list_4)/3600,pres_vapsat_4)
plt.title('Pression de vapeur saturante')
plt.xlabel('Temps (h)')
plt.ylabel(r'$e_{sw}$ (Pa)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'esw.png')

#
plt.figure(2)
plt.plot(np.array(dt_list_0)/3600, S_0)
plt.plot(np.array(dt_list_1)/3600, S_1)
plt.plot(np.array(dt_list_2)/3600, S_2)
plt.plot(np.array(dt_list_3)/3600, S_3)
plt.plot(np.array(dt_list_4)/3600, S_4)
plt.title('Rapport saturant')
plt.xlabel('Temps (h)')
plt.ylabel(r'$S$')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'S.png')

#
plt.figure(3)
plt.plot(np.array(dt_list_0)/3600,pres_vap_0)
plt.plot(np.array(dt_list_1)/3600,pres_vap_1)
plt.plot(np.array(dt_list_2)/3600,pres_vap_2)
plt.plot(np.array(dt_list_3)/3600,pres_vap_3)
plt.plot(np.array(dt_list_4)/3600,pres_vap_4)
plt.title('Pression de vapeur')
plt.xlabel('Temps (h)')
plt.ylabel(r'$e$ (Pa)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'e.png')

#
plt.figure(4)
plt.plot(np.array(dt_list_0)/3600, masse_condensation_0)
plt.plot(np.array(dt_list_1)/3600, masse_condensation_1)
plt.plot(np.array(dt_list_2)/3600, masse_condensation_2)
plt.plot(np.array(dt_list_3)/3600, masse_condensation_3)
plt.plot(np.array(dt_list_4)/3600, masse_condensation_4)
plt.title('Masse condensée')
plt.xlabel('Temps (h)')
plt.ylabel(r'$(q_w)_{condensation}$ (kg/kg)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'qw_condensation.png')

#
plt.figure(5)
plt.plot(np.array(dt_list_0)/3600, masse_activation_0)
plt.plot(np.array(dt_list_1)/3600, masse_activation_1)
plt.plot(np.array(dt_list_2)/3600, masse_activation_2)
plt.plot(np.array(dt_list_3)/3600, masse_activation_3)
plt.plot(np.array(dt_list_4)/3600, masse_activation_4)
plt.title('Masse activée')
plt.xlabel('Temps (h)')
plt.ylabel(r'$(q_w)_{activation}$ (kg/kg)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'qw_activation.png')

#
plt.figure(6)
plt.plot(np.array(dt_list_0)/3600, masse_precipitation_0)
plt.plot(np.array(dt_list_1)/3600, masse_precipitation_1)
plt.plot(np.array(dt_list_2)/3600, masse_precipitation_2)
plt.plot(np.array(dt_list_3)/3600, masse_precipitation_3)
plt.plot(np.array(dt_list_4)/3600, masse_precipitation_4)
plt.title('Masse précipitée')
plt.xlabel('Temps (h)')
plt.ylabel(r'$(q_w)_{précipitation}$ (kg/kg)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'qw_precipitation.png')

#
plt.figure(7)
plt.plot(np.array(dt_list_0)/3600, N_tot_0)
plt.plot(np.array(dt_list_1)/3600, N_tot_1)
plt.plot(np.array(dt_list_2)/3600, N_tot_2)
plt.plot(np.array(dt_list_3)/3600, N_tot_3)
plt.plot(np.array(dt_list_4)/3600, N_tot_4)
plt.title('Concentration en nombre')
plt.xlabel('Temps (h)')
plt.ylabel(r'$N$ (m$^-3$)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'N_tot.png')

#
plt.figure(8)
plt.plot(np.array(dt_list_0)/3600, np.array(rayon_0)*1e6)
plt.plot(np.array(dt_list_1)/3600, np.array(rayon_1)*1e6)
plt.plot(np.array(dt_list_2)/3600, np.array(rayon_2)*1e6)
plt.plot(np.array(dt_list_3)/3600, np.array(rayon_3)*1e6)
plt.plot(np.array(dt_list_4)/3600, np.array(rayon_4)*1e6)
plt.title('Rayon')
plt.xlabel('Temps (h)')
plt.ylabel(r'$r$ ($\mu$m)')
plt.legend(['$C_{twomey} = 3500$, $k_{twomey} = 0.9$','$C_{twomey} = 3500$, $k_{twomey} = 0.8$', '$C_{twomey} = 600$, $k_{twomey} = 0.5$', '$C_{twomey} = 250$, $k_{twomey} = 0.5$', '$C_{twomey} = 53$, $k_{twomey} = 0.5$'])
plt.savefig('figures/'+exp+'rayon.png')

plt.show()
                             








