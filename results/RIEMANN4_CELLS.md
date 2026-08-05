# One-dimensional cells and four path monitors

## Coordinates

0. Cell index $(t\$) along the orbit.
1. Cumulative Euclidean path length $(s_E\$).
2. Cumulative spherical path length $(s_S\$).
3. Cumulative hyperbolic path length $(s_H\$).

The tuple $((t,s_E,s_S,s_H)\$) is a product of monitors, not an independent smooth 4-manifold with a built-in parity.

## Bundles

- Toward 1: iterate $(f\$).
- Toward large values: iterate $(g\$) (including the even half-integer family where divergence is proved).

With a fixed number of steps, $(g\$)-orbits typically show much larger peak values and larger work $(\$sum F_t$,d_t\$) with $(F=\$Delta$log(1+x)\$), while $(f\$)-orbits remain closer to the attractor region.

Figures: `figures/riemann4/`.
