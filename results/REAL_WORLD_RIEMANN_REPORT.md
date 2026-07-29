# Collatz en geometrías (incl. Riemann) × mundo real

## Espacios

| Espacio | Curvatura K | Qué mide |
|---------|-------------|----------|
| Euclidiano R²/R³ | 0 | longitud / fin-inicio en fase log |
| Esfera S^{n} | +1 | geodésica esférica |
| Hiperbólico H² | −1 | geodésica Poincaré |
| Riemann pullback | (según g) | log-scale, conformal, Fisher |

## Por clase de activo (semillas = dígitos de precio real)

### commodity
- crece original=0.177, invertida=0.369
- ratio path inv/orig: euclid=1.2779505617851623, sphere=2.633183675731083, hyper=1.578063673722175, riemann_log=3.205939419098326
- assets: PAXGUSDT, GLD, USO

### crypto
- crece original=0.145, invertida=0.298
- ratio path inv/orig: euclid=1.3080002957493613, sphere=2.6432951200972625, hyper=1.6024258291210085, riemann_log=3.2152106261317854
- assets: BTCUSDT, ETHUSDT, SOLUSDT, DOGEUSDT

### equity
- crece original=0.199, invertida=0.382
- ratio path inv/orig: euclid=1.2186890687195748, sphere=2.4818801323268036, hyper=1.4927078021344695, riemann_log=3.019728928762649
- assets: SPY, QQQ, AAPL, NVDA

### rare_earth
- crece original=0.189, invertida=0.395
- ratio path inv/orig: euclid=1.23405877056443, sphere=2.4536377469346315, hyper=1.5047105171161366, riemann_log=2.9810071934619997
- assets: REMX, MP, UUUU

## Conclusiones

- GEOMETRÍAS USADAS: R^d (K=0), S^{n} (K=+1), H^2 Poincaré (K=-1), y longitudes riemannianas (euclídea, log-scale, conformal, Fisher) sobre el pullback de la órbita en el plano log-fase.
- MUNDO REAL: semillas = último dígito (y opcionalmente entero escalado) de precios reales crypto/bolsa/commodities/tierras raras.
- CRECIMIENTO GLOBAL: original crece en 17.7% de semillas-reales; invertida en 35.8% (comparar con ~50% abstracto a h=10).
- ESFERA: longitud geodésica inv/orig media = 2.548 (>1 ⇒ la invertida recorre más la cúpula sobre datos reales).
- EUCLÍDEO: path inv/orig media = 1.257.
- HIPERBÓLICO: path inv/orig media = 1.540.
- RIEMANN: las longitudes con g log-scale / conformal / Fisher miden la misma curva con distinta noción de distancia (escala-invariante, conformal, simplex). No hace falta curvatura variable a lo largo del camino para comparar mapas: los espacios modelo K∈{-1,0,+1} ya separan contracción vs exploración.
- APLICACIÓN REAL: si en un activo la invertida muestra ratio>>1 en esfera/H^2 y frac_grew≈0.5, hay 'grieta geométrica' de exploración; eso NO implica edge de trading hasta superar costes (como vimos en el lab paper).