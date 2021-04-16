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

# Constantes d'entrées
dt = 1 # sec
pres = 1000e2 # Pa
taux_refroidissement = -2/(24*60*60) # K/s
Rd = 287.04 # J/kgK
Rv = 461.51 # J/kgK
Lv = 2.4656e6 # J/kg
cp = 1005 # J/kg
cv = cp - Rd

# Variables initiales
s = 1
temperature_initiale = 15 + 273.15 # degC
pres_vapsat_initiale = 1704.2 # Pa H2O, le esw
pres_vap_initiale    = 1704.2 # Pa H2O, le e
qv = (Rd/Rv)*(pres_vap_initiale/pres) # kg H2O/kg air

# Boucles temporelles
timerange = (8*60*60)*dt # 8 heures
dt_list = np.arange(0, timerange+1, dt)
    
# Températures pour un pas de temps de 8h
temperature = temperature_initiale + dt_list*taux_refroidissement

# Pression de vapeur saturante pour un pas de 8h (fonction de T)
pres_vapsat = []
pres_vapsat.append(pres_vapsat_initiale)
for i in range(1,timerange+1):
    pres_vapsat.append([])
    pres_vapsat[i] = pres_vapsat[i-1] + (Lv/Rv)*pres_vapsat[i-1]*taux_refroidissement*dt/temperature[i]**2

# Tests
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
