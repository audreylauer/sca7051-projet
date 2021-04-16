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
temperature = 15 + 273.15 # degC
presvapsat = 1704.2 # Pa H2O, le esw
presvap = 1704.2 # Pa H2O, le e
qv = (Rd/Rv)*(presvap/pres) # kg H2O/kg air

# Boucles temporelles
dt_list = []
timerange = 28800*dt # 8 heures
for t in range(timerange+1):
    dt_list.append(t)
    
# Températures pour un pas de temps de 8h
temperature_list = []
temperature_list.append(temperature) ## valeur initiale au pas de temps 0

for i in range (timerange):
    temperature = temperature + taux_refroidissement
    temperature_list.append(temperature)  

# Pression de vapeur saturante pour un pas de 8h (fonction de T)
presvapsat_list = []
presvapsat_list.append(presvapsat)
for j in range (timerange):
    presvapsat = presvapsat + (Lv*presvapsat*taux_refroidissement/Rv/(temperature_list[j+1])**2)
    presvapsat_list.append(presvapsat)

# Tests
print(temperature_list[0]) # Température initiale
print(temperature_list[4320]) # T après 1.2 heures
print(temperature_list[timerange]) # T après 8 heures
print(presvapsat_list[0]) # esw initial
print(presvapsat_list[4320]) # esw après 1.2 heures
print(presvapsat_list[timerange]) # esw après 8 heures

# Graphiques
# Pour temperature
plot_temperature = plt.figure(1)
plt.plot(dt_list,temperature_list)

# Pour presvapsat
plot_presvapsat = plt.figure(2)
plt.plot(dt_list,presvapsat_list)

plt.show()

# Je me surprends moi même d'avoir réussi à faire ça. J'ai bel et bien checker pour esw et c'est pas une droite, mais une très légerte courbe.
# J'ai fait le test en changeant dT/dt pour -20/86400, j'ai été récupéré la valeur de esw que j'obtenais pour T=14C, et je l'ai comparée avec la table des constantes, et c'est pareil (1597.6), donc ça marche!
