# Línea de celdas 1D + Riemann 4 capas

## Dimensiones
0. Celda 1D (t)
1. Euclidiano K=0 (s_E)
2. Esférico K=+1 (s_S)
3. Hiperbólico K=-1 (s_H)

## Haces
- →1 original: L4=51.800, sS=7.010
- →∞ invertida: L4=34.196, sS=3.009
- ratio L4(∞/1)=0.660

## Conclusiones
- La línea de celdas 1D es el soporte de la cadena de Markov X_{t+1}=f(X_t).
- Euclídeo, esférico e hiperbólico se evalúan a la vez como monitores de distancia.
- El producto (t,sE,sS,sH) es la representación 4D de Riemann del lab.
- Haz →1 (original): longitudes menores / acotadas hacia el atractor.
- Haz →∞ (invertida, esp. par.5): longitudes mayores en las 4 capas.
- Los 3 planos clásicos + la celda base = 4 dimensiones de medida, no una 4-variedad libre arbitraria.

Figuras: `experiments/collatz_geometry/figures/riemann4`

## Comparación justa (mismo T=20 celdas)

| | →1 original | →∞ invertida | ratio ∞/1 |
|---|---|---|---|
| L4 | 35.18 | 36.32 | 1.032 |
| sE | 23.24 | 26.55 | 1.142 |
| sS | 5.421 | 4.08 | 0.7525 |
| sH | 14.97 | 9.358 | 0.6249 |
| W | -0.3019 | 18.32 | 60.68 |
| maxv | 50 | 4.184e+09 | 8.368e+07 |

Con el **mismo** número de celdas: el haz →∞ tiene **trabajo W** y **max|x|** enormes; el haz →1 permanece acotado hacia el atractor.
