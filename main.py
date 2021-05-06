# Projet de modelisation d'un brouillard chaud
# Dans le cadre du cours SCA7051 - Physique et dynamique des nuages
# Département des sciences de la Terre et de l'atmosphère
# Université du Québec à Montréal
# Auteurs: Audrey Lauer, Émile Cardinal
# Date de création : 2021-04-13
# Date de remise : 2021-05-06

# Importation des packages
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import csv

# bools pour tester des choses
test = False
version_prelim = False

# nom de l'expérience pour l'enregistrement des données
exp = 'p-1000-'

# Constantes d'entrées
dt = 1 # sec
pres = 1000e2 # Pa
Rd = 287.04 # J/kgK
Rv = 461.51 # J/kgK
Lv = 2.4656e6 # J/kg
cp = 1005 # J/kg
cv = cp - Rd
rho_w = 997 
H = 250

# Constantes Twomey
C_twomey = 3500
k_twomey = 0.9

# Variables initiales
S = [1]
temperature_initial = 15 + 273.15 # degC
pres_vapsat_initial = 1704.2 # Pa H2O, le esw
pres_vap_initial    = pres_vapsat_initial * S[0] # Pa H2O, le e
e_prime = [pres_vap_initial]
qv = (Rd/Rv)*(pres_vap_initial/pres) # kg H2O/kg air

if version_prelim:
    rayon_initial = 0.2e-6 # m
    concentration_initial = 1000
    masse_initial = qv
else:
    rayon_initial = 0.2e-6
    concentration_initial = 1000
    masse_initial = 0

# Boucle temporelle
# temperature et sa variation
taux_constant = True

# juste du refroidissement pour 8h
if taux_constant:
    timerange = (8*60*60) # 8 heures
    dt_list = np.arange(0, timerange+1, dt) # Pas de temps dans une liste

    taux_refroidissement = np.zeros(len(dt_list)+1)
    temperature = np.zeros(len(dt_list)+1)
    temperature[0] = temperature_initial
    for i in range(1,len(dt_list)):
        taux_refroidissement[i] = -5/(24*60*60) # K/s
        temperature[i] = temperature[i-1] + dt*taux_refroidissement[i]

# refroidissement pour 8h et réchauffement pour 4h
elif not taux_constant:
    timerange = (12*60*60)
    dt_list = np.arange(0, timerange+1, dt) # Pas de temps dans une liste

    taux_refroidissement = np.zeros(timerange+1)
    temperature = np.zeros(timerange+1)
    temperature[0] = temperature_initial
    taux_1 = -5/(24*60*60)
    temps_1 = 8*60*60
    taux_2 = 12/(24*60*60)
    temps_2 = 12*60*60

    for i in range(len(dt_list)):
        t = dt_list[i]
        if (t <= temps_1):
            taux_refroidissement[i] = taux_1 
        elif (t > temps_1) and (i <= temps_2):
            taux_refroidissement[i] = taux_2

        if i == 0:
            temperature[i] = temperature_initial
        else:
            temperature[i] = temperature[i-1] + taux_refroidissement[i]*dt

# Valeurs initiales
N_CCN = [0]
pres_vap = [pres_vap_initial]
pres_vapsat = [pres_vapsat_initial]
e_prime = [pres_vap_initial]
C_prime = [0]
C_double_prime = [0]
P = [0]
S_double_prime = [1]
rayon = [rayon_initial]
masse_activation = [masse_initial]
masse_condensation = [0]
masse_precipitation = [0]
N_precipitation = [0]
masse_totale = [0]
N_total = [0]

# Début de la boucle temporelle
for i in range(1,np.size(dt_list)): 
    # Ajout d'un element dans les variables enregistrees
    N_CCN.append([])
    pres_vap.append([])
    pres_vapsat.append([])
    S.append([])
    e_prime.append([])
    C_prime.append([])
    C_double_prime.append([])
    P.append([])
    S_double_prime.append([])
    rayon.append([])
    masse_activation.append([])
    masse_condensation.append([])
    masse_precipitation.append([])
    N_precipitation.append([])
    masse_totale.append([])
    N_total.append([])

    # Pression de vapeur saturante (Clausius Clapeyron)
    pres_vapsat[i] = pres_vapsat[i-1] + (Lv/Rv)*pres_vapsat[i-1]*taux_refroidissement[i]*dt/temperature[i]**2

    # Calcul de S prime (C = 0)
    P[i] = -S[i-1] * (Lv/Rv)/temperature[i-1]**2 * taux_refroidissement[i]
    S_prime = S[i-1] + P[i]*dt
    e_prime[i] = pres_vapsat[i] * S_prime

    # Condensation
    C_prime[i] = (1/pres_vapsat[i-1]) * (e_prime[i] - e_prime[i-1])/dt
    S_double_prime[i] = S[i-1] + (P[i] - C_prime[i])*dt

    # Activation des aerosols (regarder ce qui reste qui n'a pas condensé)
    if version_prelim: # test sans jamais activer les aerosols
        S[i] = S_double_prime[i]
        pres_vap[i] = pres_vapsat[i] * S[i]
        delta_masse_condensation = -1*(e_prime[i] - e_prime[i-1])/dt * Rd / (Rv*pres)
        masse_condensation[i] = masse_condensation[i-1] + delta_masse_condensation*dt

        rayon[i] = ( (3 / (4*np.pi*rho_w)) * (masse_condensation[i])/(1000) )**(1/3)

    if not version_prelim: # vrai code
        if (S_double_prime[i] < 1): # pas d'activation (pas assez de vapeur)

            # calculer la masse_condensation avec le C_prime
            masse_condensation[i] = 0.
            delta_masse_condensation = 0.
 
            #C_ajuste = P[i] - (1 - S[i-1])/dt
            #S[i] = S[i-1] + (P[i] - C_ajuste)*dt
            S[i] = S_double_prime[i]
            pres_vap[i] = pres_vapsat[i] * S[i]

            N_CCN[i] = 0.
            masse_activation[i] = 0.
            delta_masse_activation = 0.

        elif (S_double_prime[i] >= 1): # activation (assez de vapeur)

            # calculer la masse_condensation avec le C_prime
            delta_masse_condensation = -1*(e_prime[i] - e_prime[i-1])/dt * Rd / (Rv*pres)
            masse_condensation[i] = masse_condensation[i-1] + delta_masse_condensation*dt

            # Relation de Twomey
            N_CCN[i] = C_twomey * (S_double_prime[i] - 1)**k_twomey
        
            # (qw)_activation
            masse_activation[i] = rayon[i-1]**3 * 4*np.pi*rho_w*N_CCN[i] / 3 
            delta_masse_activation = ( masse_activation[i] - masse_activation[i-1] )/dt
            C_double_prime[i] = - delta_masse_activation * (Rv*pres) / Rd / pres_vapsat[i]

            S_double_prime[i] = S[i-1] + (P[i] - C_prime[i] - C_double_prime[i])*dt

            # Re-regarder s'il reste assez de vapeur (?)
            if (S_double_prime[i] < 1):

                #C_ajuste = P[i] - (1 - S[i-1])/dt
                #pres_vap[i] = C_ajuste*pres_vapsat[i]*dt + pres_vap[i-1]
                #S[i] = 1
                S[i] = S_double_prime[i]
                pres_vap[i] = pres_vapsat[i] * S[i]

            else:

                S[i] = S_double_prime[i]
                pres_vap[i] = pres_vapsat[i] * S[i]
  
    # Précipitation
    vitesse = 1.19e8 * rayon[i-1]**2
    delta_masse_precipitation = vitesse * masse_totale[i-1] / H
    masse_precipitation[i] = masse_precipitation[i-1] + delta_masse_precipitation*dt
    delta_N_precipitation = vitesse * N_total[i-1] / H
    N_precipitation[i] = N_precipitation[i-1] + delta_N_precipitation*dt

    # Calcul de la masse à la fin du pas de temps
    masse_totale[i] = masse_totale[i-1] + (delta_masse_condensation + delta_masse_activation + delta_masse_precipitation)*dt
    N_total[i] = N_total[i-1] + (delta_N_precipitation + (N_CCN[i] - N_CCN[i-1])/dt)*dt

    # rayon final
    rayon[i] = ( (3 / (4*np.pi*rho_w)) * (masse_totale[i])/(N_total[i]) )**(1/3)


# Fin boucle temporelle


# Écriture des données dans un fichier csv
with open(exp+'data.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(dt_list)
     wr.writerow(pres_vapsat)
     wr.writerow(S)
     wr.writerow(pres_vap)
     wr.writerow(masse_condensation)
     wr.writerow(masse_activation)
     wr.writerow(masse_precipitation)
     wr.writerow(np.array(N_CCN) + np.array(N_precipitation))
     wr.writerow(rayon)

# Graphiques
#
plt.figure(1)
plt.plot(dt_list/3600,pres_vapsat)
plt.title('Pression de vapeur saturante')
plt.xlabel('Temps (h)')
plt.ylabel(r'$e_{sw}$ (Pa)')
plt.savefig('figures/'+exp+'esw.png')

#
plt.figure(2)
plt.plot(dt_list/3600, S)
plt.title('Rapport saturant')
plt.xlabel('Temps (h)')
plt.ylabel(r'$S$')
plt.savefig('figures/'+exp+'S.png')

#
plt.figure(3)
plt.plot(dt_list/3600,pres_vap)
plt.title('Pression de vapeur')
plt.xlabel('Temps (h)')
plt.ylabel(r'$e$ (Pa)')
plt.savefig('figures/'+exp+'e.png')

#
plt.figure(4)
plt.plot(dt_list/3600, masse_condensation)
plt.title('Masse condensée')
plt.xlabel('Temps (h)')
plt.ylabel(r'$(q_w)_{condensation}$ (kg/kg)')
plt.savefig('figures/'+exp+'qw_condensation.png')

#
plt.figure(5)
plt.plot(dt_list/3600, masse_activation)
plt.title('Masse activée')
plt.xlabel('Temps (h)')
plt.ylabel(r'$(q_w)_{activation}$ (kg/kg)')
plt.savefig('figures/'+exp+'qw_activation.png')

#
plt.figure(6)
plt.plot(dt_list/3600, masse_precipitation)
plt.title('Masse précipitée')
plt.xlabel('Temps (h)')
plt.ylabel(r'$(q_w)_{précipitation}$ (kg/kg)')
plt.savefig('figures/'+exp+'qw_precipitation.png')

#
plt.figure(7)
plt.plot(dt_list/3600, np.array(N_CCN) + np.array(N_precipitation))
plt.title('Concentration en nombre')
plt.xlabel('Temps (h)')
plt.ylabel(r'$N$ (m$^-3$)')
plt.savefig('figures/'+exp+'N_tot.png')

#
plt.figure(8)
plt.plot(dt_list/3600, np.array(rayon)*1e6)
plt.title('Rayon')
plt.xlabel('Temps (h)')
plt.ylabel(r'$r$ ($\mu$m)')
plt.savefig('figures/'+exp+'rayon.png')
