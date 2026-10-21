# Changelog
 
## 0.9.1 (2026-10-21)
 
- Fixed: with `vclass(pd)`, the solver stopped with "PD cell SDP failed: optimal_inaccurate" whenever the treated unit lies inside the donors' convex hull. The attainability SDP in `interior_pd.py` was solved at tolerance 1e-9, where Clarabel reports a solution as inaccurate, and the code accepted only an optimal status. It now solves at tolerance 1e-6 and also accepts an inaccurate solution whose slack is clearly away from zero (|slack| > 1e-4). In the California placebos this affected Iowa, Nebraska and South Dakota. Results are otherwise unchanged. The diagonal class and the exterior cases were not affected.
 
## 0.9.0 (2026-10-19)
 
- First public release.