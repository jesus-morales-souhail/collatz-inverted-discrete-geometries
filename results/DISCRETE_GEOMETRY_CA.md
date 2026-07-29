# Geometrías discretas + Collatz + CA

## Principio
No se define par/impar en una variedad continua sin discretizar.

## Modelos
- E: Z+
- S: Z/NZ
- H: (valor, nivel árbol)
- 4D: (E, S, H_val, H_level)
- CA: reglas 30, 90, 110, 184

## Conclusiones
- Euclidiano Z: f empuja a 4-2-1 (conjetura); g mezcla divergencia y ciclos (~no siempre ∞).
- Esférico Z/N: no existe ∞; f visita 1 más que g; g tiene más exploración en el anillo.
- Hiperbólico: g sube de nivel (capacidad exp) cuando el valor crece; f baja de nivel hacia el atractor.
- Producto 4D: modo f → componentes E/H acotadas hacia 1; modo g → E/H crecen; S siempre finita.
- CA elementales: Markov espacial local; misma lógica que Collatz (dependencia del presente).
- Grieta: paridad Collatz = cinta temporal; CA = cinta espacial; juntas muestran dependencia ≠ independencia.

Figuras: `/home/ashpokemon/Proyectos/06_Trading_Simulacion/collatz_ensemble_bot/experiments/collatz_geometry/figures/discrete_geo`