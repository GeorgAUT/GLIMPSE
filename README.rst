GLIMPSE
=======

This repository contains a collection of Matlab implementations of resonance-based low-regularity integrators and related methods that were studied over the course of the MSCA project `GLIMPSE <https://doi.org/10.3030/101064261>`_ and subsequent research projects. The main routines can be found in ``src/resonance_based_schemes`` and examples illustrating the application of these methods as used in the below mentioned papers are provided in ``examples/``. To run the examples please first run ``src/set_paths.m``.

A good starting point to familiarise yourself with the code is given in ``examples/introductory_example/`` where you can start by running ``convergence_experiments_nls.m`` followed by ``plotting.m``.



Acknowledgements
================

If you use this code in an academic paper, please cite [1]_, [2]_, [3]_, [4]_, [5]_:

 @misc{maierhofer2022symplectic,
      title={Bridging the gap: symplecticity and low regularity in {R}unge--{K}utta resonance-based schemes}, 
      author={Georg Maierhofer and Katharina Schratz},
      year={2022},
      eprint={2205.05024},
      archivePrefix={arXiv},
      primaryClass={math.NA}
   }

 @misc{alamabronsard2023symmetric,
      title={Symmetric resonance based integrators and forest formulae}, 
      author={Yvonne Alama Bronsard and Yvain Bruned and Georg Maierhofer and Katharina Schratz},
      year={2023},
      eprint={2305.16737},
      archivePrefix={arXiv},
      primaryClass={math.NA}
   }

 @misc{banica2024schroedingermaps,
      title={Numerical integration of {S}chr\"odinger maps via the {H}asimoto transform}, 
      author={Valeria Banica and Georg Maierhofer and Katharina Schratz},
      journal={SIAM Journal on Numerical Analysis},
      volume={62},
      number={1},
      pages={322--352},
      year={2024},
      publisher={SIAM}
   }

 @misc{feng2024longtime,
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



.. [1] Maierhofer, G. and Schratz, K., “Bridging the gap: symplecticity and low regularity in Runge-Kutta resonance-based schemes”, 2022. arXiv.2205.05024, to appear in Math. Comp.

.. [2] Alama Bronsard, Y., Bruned, Y., Maierhofer, G. and Schratz, K., “Symmetric resonance based integrators and forest formulae”, 2023. arXiv.2305.16737.

.. [3] Banica, V., Maierhofer, G. and Schratz, K., “Numerical integration of Schrödinger maps via the Hasimoto transform”, SIAM J. Numer. Anal. 62(1), pp.322-352, 2024.

.. [4] Feng, Y., Maierhofer, G. and Schratz, K., “Long-time error bounds of low-regularity integrators for nonlinear Schrödinger equations”, Math. Comp. 93(348), pp.1569-1598, 2024.

.. [5] Cui, J. and Maierhofer, G., “A Wong--Zakai resonance-based integrator for the nonlinear Schrödinger equation with white noise dispersion”, 2025. arXiv:2503.19346.
