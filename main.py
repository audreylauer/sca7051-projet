# Projet de modelisation d'un brouillard chaud
# Dans le cadre du cours SCA7051 - Physique et dynamique des nuages
# Département des sciences de la Terre et de l'atmosphère
# Université du Québec à Montréal
# Auteurs: Audrey Lauer, Émile Cardinal
# Date de création : 2021-04-13
# Date de remise : 2021-05-04
#
#########################
#
# Importation des packages
import numpy as np
import math
import matplotlib as mpl
import pandas as pd

#
########## Constantes d'entrées
#
dt = 1 # sec
pres = 1000e2 # Pa
taux_refroidissement = -2/(24*60*60) # K/s
Rd = 287.04 # J/kgK
Rv = 461.51 # J/kgK
Lv = 2.4656e6 # J/kg
cp = 1005 # J/kg
cv = cp-Rd

#
########## Variables initiales
#

s = 1
temperature = 15 + 273.15 # degC
presvapsat = 1704.2 # Pa H2O, le esw
presvap = 1704.2 # Pa H2O, le e
qv = (Rd/Rv)*(presvap/pres) # kg H2O/kg air
print(qv)


#
########## Boucle temporelle
#


