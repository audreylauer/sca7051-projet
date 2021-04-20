# sca7051-projet

projet de modélisation d'un brouillard chaud.

- dT/dt fixe
- dp/dt = 0 car brouillard mince

## Étape 1 - calculer de_sw/dt

![equation](https://latex.codecogs.com/gif.latex?\frac{d&space;e_{sw}}{dt}&space;=&space;\frac{L_v&space;e_{sw}}{R_vT^2}&space;\frac{dT}{dt})

[comment]: <> ($$\frac{d e_{sw}}{dt} = \frac{L_v e_{sw}}{R_vT^2} \frac{dT}{dt}$$)

## Étape 2 - S

Il faut expliciter les termes de production  P et de consommation C.

![equation](https://latex.codecogs.com/gif.latex?\frac{dS}{dt}&space;=&space;P&space;-&space;C)

[comment]: <> ($$\frac{dS}{dt} = P - C$$)


P dépend uniquement du refroidissement/réchauffement. S'il y a refroidissement, il y a production.

![equation](https://latex.codecogs.com/gif.latex?P&space;=&space;-&space;\frac{S}{e_{sw}}&space;\frac{d&space;e_{sw}}{dt}&space;=&space;-&space;S&space;\frac{L_v}{R_v&space;T^2}&space;\frac{dT}{dt})

![equation](https://latex.codecogs.com/gif.latex?C&space;=&space;\frac{1}{e_{sw}}&space;\frac{de}{dt}&space;=&space;-&space;\frac{1}{e_{sw}}&space;\frac{R_v}{R_d}&space;p&space;\frac{dq_w}{dt})

[comment]: <> ($$P = - \frac{S}{e_{sw}} \frac{d e_{sw}}{dt} = - S \frac{L_v}{R_v T^2} \frac{dT}{dt}$$)
[comment]: <> ($$C = \frac{1}{e_{sw}} \frac{de}{dt}         = - \frac{1}{e_{sw}} \frac{R_v}{R_d} p \frac{dq_w}{dt}$$)

