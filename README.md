# GLIMPSE
 
GLIMPSE is a collection of MATLAB implementations of resonance-based low-regularity integrators and related numerical methods developed in the context of the MSCA project [GLIMPSE](https://doi.org/10.3030/101064261) and subsequent research projects.
 
The code is organised as a small MATLAB library (under `src/`) together with reproducibility scripts (under `examples/`) corresponding to the related publications listed below.
 
## Repository structure
 
- `src/`
  - `low_regularity_integrators/`
    - `KdV/`
    - `cubic_NLS/`
    - `quadratic_NLS/`
    - `SchroedingerMap/`
    - low-regularity integrators for a range of systems
  - `classical_numerical_methods/`
    - reference / baseline methods used in comparisons
  - `auxilliary_modules/`
    - common helpers (1D/2D utilities, vector-field operations, etc.)
  - `set_paths.m`
    - adds all required subfolders to the MATLAB path
- `examples/`
  - `introductory_example/` (recommended starting point)
  - `[1]_...` to `[5]_...` folders containing scripts used for the experiments in the publications listed below
 
## Requirements
 
- MATLAB (a recent version)
 
## Quickstart
 
1. Start MATLAB with the repository root as the working directory.
2. Add the library to your MATLAB path by running:
 
```matlab
run('src/set_paths.m')
```
 
3. Run the introductory example:
 
```matlab
cd('examples/introductory_example')
run('convergence_experiments_nls.m')
run('plotting.m')
```
 
The introductory example writes/reads data files under `examples/introductory_example/data/` and writes figures under `examples/introductory_example/images/`.
 
## Paper-specific examples
 
The folder `examples/` contains experiments grouped by paper:
 
- `examples/[1]_symplectic_RK_resonance-based_schemes/`
- `examples/[2]_symmetric_resonance-based_schemes/`
- `examples/[3]_schroedinger_maps/`
- `examples/[4]_long_time_error_low-reg_integrators/`
- `examples/[5]_nlse_with_noise_dispersion/`
 
Each folder contains one or more top-level driver scripts (`*.m`) and a `data/` subfolder containing cached results / reference solutions.
 
## Citation
 
If you use this software, please cite the Zenodo record:
 
- https://doi.org/10.5281/zenodo.15733738
 
If you use this code in an academic paper, please also cite the related publications:
 
```bibtex
@article{maierhofer2022symplectic,
     title={Bridging the gap: symplecticity and low regularity in {R}unge--{K}utta resonance-based schemes},
     author={Georg Maierhofer and Katharina Schratz},
     journal={Mathematics of Computation},
     year={2025}
  }
 
@article{alamabronsard2023symmetric,
     title={Symmetric resonance based integrators and forest formulae},
     author={Yvonne Alama Bronsard and Yvain Bruned and Georg Maierhofer and Katharina Schratz},
     journal={Foundations of Computational Mathematics},
     year={2026}
  }
 
@article{banica2024schroedingermaps,
     title={Numerical integration of {S}chr\"odinger maps via the {H}asimoto transform},
     author={Valeria Banica and Georg Maierhofer and Katharina Schratz},
     journal={SIAM Journal on Numerical Analysis},
     volume={62},
     number={1},
     pages={322--352},
     year={2024},
     publisher={SIAM}
  }
 
@article{feng2024longtime,
     title={Long-time error bounds of low-regularity integrators for nonlinear {S}chr\"odinger equations},
     author={Yue Feng and Georg Maierhofer and Katharina Schratz},
     journal={Mathematics of Computation},
     volume={93},
     number={348},
     pages={1569--1598},
     year={2024}
  }
 
@misc{cui2025wongzakai,
     title={A {W}ong--{Z}akai resonance-based integrator for the nonlinear {S}chr\"odinger equation with white noise dispersion},
     author={Jianbo Cui and Georg Maierhofer},
     year={2025},
     eprint={2503.19346},
     archivePrefix={arXiv},
     primaryClass={math.NA}
  }
```
 
1. Maierhofer, G. and Schratz, K., “Bridging the gap: symplecticity and low regularity in Runge-Kutta resonance-based schemes”, Math. Comp., 2025.
2. Alama Bronsard, Y., Bruned, Y., Maierhofer, G. and Schratz, K., “Symmetric resonance based integrators and forest formulae”, Found. Comp. Math., 2026.
3. Banica, V., Maierhofer, G. and Schratz, K., “Numerical integration of Schrödinger maps via the Hasimoto transform”, SIAM J. Numer. Anal. 62(1), pp.322-352, 2024.
4. Feng, Y., Maierhofer, G. and Schratz, K., “Long-time error bounds of low-regularity integrators for nonlinear Schrödinger equations”, Math. Comp. 93(348), pp.1569-1598, 2024.
5. Cui, J. and Maierhofer, G., “A Wong--Zakai resonance-based integrator for the nonlinear Schrödinger equation with white noise dispersion”, 2025. arXiv:2503.19346.
 
## License
 
This project is licensed under the MIT license (see `LICENSE`).
