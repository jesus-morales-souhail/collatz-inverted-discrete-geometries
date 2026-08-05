# Integer seeds from external numerical series

Some checks use last digits (or scaled integers) taken from public price series as **sources of integer seeds**. The maps $f$ and $g$ are still pure arithmetic. The point is only that the exploration difference between $g$ and $f$ is not an artefact of the artificial range $\{0,\ldots,N\}$.

Across the series stored in the corresponding JSON:

- Growth fractions under $g$ exceed those under $f$.
- Path-length ratios (inverted over normal) in Euclidean, spherical and hyperbolic monitors are typically larger than one.

No statement is made about financial use of these series.

Full table: `real_world_riemann.json`.
