# Cadena de Markov dependiente — Collatz Invertida

## Definición

$$X_{n+1}=f(X_n),\quad P=1$$

$$f(x)=\begin{cases}3x+1 & \text{paridad/dígito par}\\ x/2 & \text{impar}\end{cases}$$

## Propiedad de Markov

$$P(X_{n+1}\mid X_n,\ldots,X_0)=P(X_{n+1}\mid X_n)$$

## Verificación

- Determinista OK (invertida): **True**
- Estados visitados: 1055
- Desviación de paridades vs 0.5 bajo f: **0.269** (i.i.d. ~0.007)

## Conclusión

Dependiente del **presente**. No independiente. El ~50% no es moneda paso a paso.

Figuras en `/home/ashpokemon/Proyectos/06_Trading_Simulacion/collatz_ensemble_bot/experiments/collatz_geometry/figures/markov`