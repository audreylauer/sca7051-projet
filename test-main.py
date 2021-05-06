# Projet de modelisation d'un brouillard chaud
# Dans le cadre du cours SCA7051 - Physique et dynamique des nuages
# Département des sciences de la Terre et de l'atmosphère
# Université du Québec à Montréal
# Auteurs: Audrey Lauer, Émile Cardinal
# Date de création : 2021-04-13
# Date de remise : 2021-05-04

# Importation des packages
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

test = False
version_prelim = False

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
C_twomey = 2000
k_twomey = 0.7

# Variables initiales
S = [1]
temperature_initial = 15 + 273.15 # degC
pres_vapsat_initial = 1704.2 # Pa H2O, le esw
pres_vap_initial    = 1704.2 # Pa H2O, le e
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
timerange = (12*60*60)*dt # 8 heures
dt_list = np.arange(0, timerange+1, dt) # Pas de temps dans une liste
taux_constant = False

if taux_constant:
    taux_refroidissement = -10/(24*60*60) # K/s
    temperature = temperature_initial + dt_list*taux_refroidissement # Températures pour un pas de temps de 8h

elif not taux_constant:
    taux_refroidissement = np.zeros(timerange+1)
    temperature = np.zeros(timerange+1)
    temperature[0] = temperature_initial
    taux_1 = -5/(24*60*60)
    temps_1 = 8*60*60
    taux_2 = -5/(24*60*60)
    temps_2 = 12*60*60

    for i in range(timerange+1):
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

is_it_activated = False

# Début de la boucle temporelle
for i in range(1,timerange+1): 
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

    # Variables pronostiques
    pres_vapsat[i] = pres_vapsat[i-1] + (Lv/Rv)*pres_vapsat[i-1]*taux_refroidissement[i]*dt/temperature[i]**2

    # delta masses
    # masse totale du temps précédent
    masse_totale = masse_activation[i-1] + masse_condensation[i-1] + masse_precipitation[i-1]
    N_total = N_CCN[i-1] + N_precipitation[i-1]

    # Précipitation
    vitesse = 1.19e8 * rayon[i-1]**2
    delta_masse_precipitation = vitesse * masse_totale / H
    masse_precipitation[i] = masse_precipitation[i-1] + delta_masse_precipitation*dt
    delta_N_precipitation = vitesse * N_total / H
    N_precipitation[i] = N_precipitation[i-1] + delta_N_precipitation*dt

    # Calcul de S prime
    P[i] = -S[i-1] * (Lv/Rv)/temperature[i-1]**2 * taux_refroidissement[i]
    S_prime = S[i-1] + P[i]*dt
    e_prime[i] = pres_vapsat[i] * S_prime

    # Condensation
    FW_plus_FD = 1/7e-10
    delta_masse_condensation = (4*np.pi*rho_w*Rd)/(pres*FW_plus_FD) * temperature[i-1] * N_total * rayon[i-1] * (S[i-1] - 1)
    masse_condensation[i] = masse_condensation[i-1] + delta_masse_condensation*dt
    C_prime[i] = (Rv/Rd) * (pres/pres_vapsat[i-1]) * delta_masse_condensation

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
 
            print('<1')

            # calculer la masse_condensation avec le C_prime
            #masse_condensation[i] = 0.
 
            #C_ajuste = P[i] - (1 - S[i-1])/dt
            #S[i] = S[i-1] + (P[i] - C_ajuste)*dt
            S[i] = S_double_prime[i]
            pres_vap[i] = pres_vapsat[i] * S[i]

            N_CCN[i] = 0.
            masse_activation[i] = 0.

        else:

            is_it_activated = True

            # Relation de Twomey
            N_CCN[i] = C_twomey * (S_double_prime[i] - 1)**k_twomey
            #if not is_it_activated:
            #    N_CCN[i] = C_twomey * (S_double_prime[i] - 1)**k_twomey
            #elif is_it_activated:
            #    N_CCN[i] = C_twomey * S_double_prime[i] - N_CCN[i-1]
            #    N_CCN[i] = max(0., N_CCN[i])
        
            # (qw)_activation
            #delta_masse_activation = (4*np.pi*rho_w*Rd)/(pres*FW_plus_FD) * temperature[i-1] * N_CCN[i] * rayon[i-1] * (S_double_prime[i] - 1)
            #masse_activation[i] = masse_condensation[i-1] + delta_masse_condensation*dt
            masse_activation[i] = rayon[i-1]**3 * 4*np.pi*rho_w*N_CCN[i] / 3 
            delta_masse_activation = ( masse_activation[i] - masse_activation[i-1] )/dt
            C_double_prime[i] = masse_activation[i] * (Rv*pres) / Rd / pres_vapsat[i]

            S_double_prime[i] = S[i-1] + (P[i] - C_prime[i] - C_double_prime[i])*dt

            # Re-regarder s'il reste assez de vapeur (?)
            if (S_double_prime[i] < 1):

                print('>1 puis <1')
                print(i)
                C_ajuste = P[i] - (1 - S[i-1])/dt
                #pres_vap[i] = C_ajuste*pres_vapsat[i]*dt + pres_vap[i-1]
                #S[i] = 1
                S[i] = S_double_prime[i]
                pres_vap[i] = pres_vapsat[i] * S[i]

            else:

                S[i] = S_double_prime[i]
                pres_vap[i] = pres_vapsat[i] * S[i]
  

    # rayon final
    rayon[i] = ( (3 / (4*np.pi*rho_w)) * (masse_activation[i] + masse_condensation[i] + masse_precipitation[i])/(N_CCN[i] + N_precipitation[i]) )**(1/3)

# Fin boucle temporelle

# Tests
if test:
    print(temperature[0]) # Température initiale
    print(temperature[4320]) # T après 1.2 heures
    print(temperature[timerange]) # T après 8 heures
    print(pres_vapsat[0]) # esw initial
    print(pres_vapsat[4320]) # esw après 1.2 heures
    print(pres_vapsat[timerange]) # esw après 8 heures

# Graphiques
 
plt.figure(1)
plt.plot(dt_list[1:-1],S_double_prime[1:-1])
plt.title('S_double_prime')

#
plt.figure(2)
plt.plot(dt_list,S)
plt.title('S')

#
plt.figure(3)
plt.plot(dt_list,pres_vapsat)
plt.title('e_sw')

#
plt.figure(4)
plt.plot(dt_list,pres_vap)
plt.title('e')


plt.show()

# Je me surprends moi même d'avoir réussi à faire ça. J'ai bel et bien checker pour esw et c'est pas une droite, mais une très légerte courbe.
# J'ai fait le test en changeant dT/dt pour -20/86400, j'ai été récupéré la valeur de esw que j'obtenais pour T=14C, et je l'ai comparée avec la table des constantes, et c'est pareil (1597.6), donc ça marche!
