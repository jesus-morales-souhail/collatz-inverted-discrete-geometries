# Deterministic Markov chain for the inverted map

## Definition

\[
X_{n+1}=g(X_n),\qquad
g(n)=\begin{cases}3n+1 & n\text{ even}\\ n/2 & n\text{ odd}\end{cases}
\]

with probability one.

## Markov property

\[
P(X_{n+1}\mid X_n,\ldots,X_0)=P(X_{n+1}\mid X_n).
\]

If the next value were independent of the present state, the map \(g\) would play no role and there would be no trajectories to study.

## Check

On finite samples of seeds, each visited state has a unique successor under \(g\). Conditional on the present parity, the next parity is far from a fair coin (in particular even is always followed by odd).

Figures: `figures/markov/`.
