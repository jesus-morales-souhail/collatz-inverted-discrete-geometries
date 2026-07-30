# Claims

| Claim | Status |
|-------|--------|
| $f$ is the classical Collatz map | Definition |
| $g$ swaps the even and odd branches of $f$ | Definition |
| $X_{n+1}=g(X_n)$ is a deterministic Markov chain | True by construction |
| Successive parities under $g$ are i.i.d. fair coins | False (even maps to odd under $3n+1$) |
| Every positive integer under $f$ reaches 4-2-1 | Open conjecture |
| Every integer orbit under $g$ diverges to $\infty$ | False (cycles exist) |
| On $\mathbb{Z}/N\mathbb{Z}$ orbits cannot diverge to $\infty$ | True (finite state space) |
| $x_0=2k+\frac{1}{2}$ under floor-even $\to 3x+1$ diverges | Theorem (this repo) |
| $g$ explores more than $f$ on the discrete E/S/H/product monitors here | Empirical |
| Continuous manifolds carry a natural parity for Collatz | Not claimed |
| Discrete flux $=-\kappa\nabla\varphi$ factors under DFT for constant $\kappa$ (circular) | Standard convolution theorem |
| Same flux form on cell fields $\varphi^E,\varphi^S,\varphi^H$ | Construction (this repo) |
| $f$ or $g$ is a Fourier/Ohm law | Not claimed |
| Bitstreams from $f$/$g$ on lattice planes pass a self-contained statistical battery | Empirical: mostly fail in the sample here |
| $f$ or $g$ is cryptographically secure as a PRNG | Not claimed |
| Full NIST / TestU01 / PractRand certification | Not claimed |
| Log-orbit $Y_t=\log(1+|X_t|)$ has negative drift under $f$ on average | Heuristic (Lagarias-type); empirical diagnostics only |
| Fitted OU SDE, skew-product RDS, PDMP or MSM for these maps | Not claimed |
