# sca7051-projet

projet de modélisation d'un brouillard chaud.

- dT/dt fixe
- dp/dt = 0 car brouillard mince

## Étape 1 - calculer de_sw/dt

$$\frac{d e_{sw}}{dt} = \frac{L_v e_{sw}}{R_vT^2} \frac{dT}{dt}$$

## Étape 2 - S

Il faut expliciter les termes de production  P et de consommation C.

$$\frac{dS}{dt} = P - C$$

P dépend uniquement du refroidissement/réchauffement. S'il y a refroidissement, il y a production ($P > 0$).

$$P = - \frac{S}{e_{sw}} \frac{d e_{sw}}{dt} = - S \frac{L_v}{R_v T^2} \frac{dT}{dt}$$

$$C = - \frac{1}{e_{sw}} \frac{R_v}{R_d} p \frac{dq_w}{dt}$$

