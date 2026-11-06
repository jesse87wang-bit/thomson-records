# A lower-energy configuration for the Thomson problem, N = 2502

| | Energy (Coulomb, unit sphere) |
|---|---|
| This configuration | **3060842.2207570** |
| Cambridge Cluster Database (s = 1), as listed on 2026-11-06 | 3060842.2403807 |
| Improvement | **0.0196237** |

- Points: N = 2502, coordinates in `thomson_N2502.xyz` (one point per line, x y z, 16 decimals).
- Symmetry elements detected in the new configuration (tolerance 1e-6; the detector may miss some): no rotation axes.
- Largest tangential gradient component: 3.5e-09.
- Smallest distance between points: 0.0683.

## Verify it yourself

```
python3 verify.py --n 2502 --candidate thomson_N2502.xyz --record <the database file points_2502_1.dat>
```

Run this inside the `N2502/` directory. `verify.py` needs NumPy and SciPy (SciPy is used only for nearest-neighbour
checks): pairwise distances in float64 and exact summation with `math.fsum`. It does not share any code with the
search program.

## How it was found

Two lower-energy candidates were first found by replica-exchange basin hopping in the full 3N-dimensional space, starting from the structure listed in the Cambridge Cluster Database (random displacements and rotations of small patches of points, each followed by L-BFGS relaxation on the sphere). Their displacement fields relative to that structure were then superposed and the result relaxed with L-BFGS. Finally the configuration was refined with a gradient-only Barzilai-Borwein iteration until the largest tangential gradient was below 1e-8.

The search code was written with AI assistance (Claude Code). Energies were evaluated in C during the search; this
configuration was re-checked with the independent script above.

## Credit and licence

Found by Tao Wang. Coordinates and scripts are released under CC BY 4.0. The previous best structure is from
the Cambridge Cluster Database (D. J. Wales et al.); see also P. Amore and D. J. Wales, arXiv:2610.08122 (2026).
