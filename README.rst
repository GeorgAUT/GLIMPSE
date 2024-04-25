GLIMPSE
=======

This repository contains a collection of Matlab implementations of resonance-based low-regularity integrators and related methods that were studied over the course of the MSCA project `GLIMPSE <https://doi.org/10.3030/101064261>`_. The main routines can be found in ``src/resonance_based_schemes`` and examples illustrating the application of these methods as used in the below mentioned papers are provided in ``examples/``. To run the examples please first run ``src/set_paths.m``.

A good starting point to familiarise yourself with the code is given in ``examples/introductory_example/``, where you can start by running ``convergence_experiments_nls.m`` followed by ``plotting.m``.



Acknowledgements
================

If you use this code in an academic paper, please cite [1]_, [2]_, [3]_, [4]_::

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

 @misc{banica2022schroedingermaps,
      title={Numerical integration of {S}chr\"odinger maps via the {H}asimoto transform}, 
      author={Valeria Banica and Georg Maierhofer and Katharina Schratz},
      year={2022},
      eprint={2211.01282},
      archivePrefix={arXiv},
      primaryClass={math.NA}
   }

 @misc{feng2023longtime,
      title={Long-time error bounds of low-regularity integrators for nonlinear {S}chr\"odinger equations}, 
      author={Yue Feng and Georg Maierhofer and Katharina Schratz},
      year={2023},
      eprint={2302.00383},
      archivePrefix={arXiv},
      primaryClass={math.NA}
   }


.. [1] Maierhofer, G. and Schratz, K., “Bridging the gap: symplecticity and low regularity in Runge-Kutta resonance-based schemes”, 2022. arXiv.2205.05024.

.. [2] Alama Bronsard, Y., Bruned, Y., Maierhofer, G. and Schratz, K., “Symmetric resonance based integrators and forest formulae”, 2023. arXiv.2305.16737.

.. [3] Banica, V., Maierhofer, G. and Schratz, K., “Numerical integration of Schrödinger maps via the {H}asimoto transform”, 2022. arXiv.2211.01282, to appear in SIAM J. Numer. Anal.

.. [4] Feng, Y., Maierhofer, G. and Schratz, K., “Long-time error bounds of low-regularity integrators for nonlinear Schrödinger equations”, 2023. arXiv.2302.00383, to appear in Math. Comput.
