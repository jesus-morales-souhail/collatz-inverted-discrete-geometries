# Dos ecuaciones Collatz — Euclídeo y Esférico

## Ecuaciones

| Mapa | Par | Impar | Comportamiento |
|------|-----|-------|----------------|
| Original | n/2 | 3n+1 | Atractor 4→2→1 |
| Invertida | 3n+1 | n/2 | Crecimiento ~50% a h corto; ciclos frecuentes |

## Resultados clave

- ORIGINAL: 100.00% de semillas en 0..19999 alcanzan el atractor {1,2,4} (pasos med. 79.0).
- INVERTIDA: tras 10 pasos, 49.96% cumple x_10 > x_0 (métrica empírica ~50–52%).
- INVERTIDA ciclos: 100.0% de muestra corta cae en ciclo ≤200 pasos → no es escape universal a ∞.
- EUCLÍDEO: expansión media original=0.559, invertida=0.505 (fracción expand>1: orig 0.11, inv 0.08).
- ESFÉRICO: longitud geodésica media original=6.785, invertida=9.306; distancia fin-inicio media orig=0.233, inv=0.328.
- EXPANSIÓN LOCAL: original mean_ratio=1.750 (par 0.500 / impar 3.001); invertida mean_ratio=1.750 (par 3.001 / impar 0.499).
- MATEMÁTICA APLICABLE: ambos mapas son funciones N→N (o Z) medibles; se embuten isométricamente en R^d (fase) y se proyectan a S^{d-1}. La diferencia no es 'euclídeo vs esférico' sino la dinámica del mapa: el original colapsa al ciclo; la invertida explora más la esfera y el plano log antes de ciclar o crecer a horizonte finito.

## Hardware
- torch_cuda / NVIDIA GeForce RTX 5070 Ti

JSON completo: `two_maps_euclidean_spherical.json`

## Figuras

- `orbits_log.png`
- `euclidean_log_phase.png`
- `spherical_pca_clouds.png`
- `euclidean_expansion_hist.png`