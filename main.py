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
version_prelim = True

# Constantes d'entrées
dt = 1 # sec
timerange = (8*60*60)*dt # 8 heures
pres = 1000e2 # Pa
taux_refroidissement = -10/(24*60*60) # K/s
Rd = 287.04 # J/kgK
Rv = 461.51 # J/kgK
Lv = 2.4656e6 # J/kg
cp = 1005 # J/kg
cv = cp - Rd

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

# Boucle temporelle
dt_list = np.arange(0, timerange+1, dt) # Pas de temps dans une liste
temperature = temperature_initial + dt_list*taux_refroidissement # Températures pour un pas de temps de 8h

C_prime = [0]
P = [0]
# Pression de vapeur saturante pour un pas de 8h (fonction de T)
pres_vapsat = []
pres_vapsat.append(pres_vapsat_initial)
for i in range(1,timerange+1): # Début boucle temporelle
    # Variables pronostiques
    pres_vapsat.append([])
    pres_vapsat[i] = pres_vapsat[i-1] + (Lv/Rv)*pres_vapsat[i-1]*taux_refroidissement*dt/temperature[i]**2

    # Calcul de S (avec approximation)
    S.append([])
    e_prime.append([])
    C_prime.append([])
    P.append([])
#    P = -1 * ( S[i-1] / pres_vapsat[i-1] ) * (pres_vapsat[i] - pres_vapsat[i-1])/dt
    P[i] = -S[i-1] * (Lv/Rv)/temperature[i-1]**2 * taux_refroidissement
    S_prime = S[i-1] + P[i]*dt
    e_prime[i] = pres_vapsat[i] * S_prime
    C_prime[i] = 1/pres_vapsat[i] * (e_prime[i] - e_prime[i-1])/dt
    S_double_prime = S[i-1] + (P[i] - C_prime[i])*dt
    S[i] = S_double_prime

    # Variables diagnostiques
#    rayon[i] = 3*qw[i] / (4*pi*rho_w*N[i])

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
# Pour temperature
plot_temperature = plt.figure(1)
plt.plot(dt_list,temperature)

# Pour presvapsat
plot_presvapsat = plt.figure(2)
plt.plot(dt_list,pres_vapsat)

plt.show()

# Je me surprends moi même d'avoir réussi à faire ça. J'ai bel et bien checker pour esw et c'est pas une droite, mais une très légerte courbe.
# J'ai fait le test en changeant dT/dt pour -20/86400, j'ai été récupéré la valeur de esw que j'obtenais pour T=14C, et je l'ai comparée avec la table des constantes, et c'est pareil (1597.6), donc ça marche!
