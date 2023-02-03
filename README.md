# Parametric Level Set Framework

A comprehensive implementation of the Parametric Level Set approach for inverse problems, aimed at improving reconstruction accuracy.

## Overview

ParametricLevelSet is a MATLAB-based toolkit that employs the Parametric Level Set approach in addressing inverse problems. The framework is designed to be user-friendly and accessible to a wide range of users.

The general least-squares formulation of inverse problems can be expressed as follows:

$$\min_{x \in \mathbb{R}^N} \| F(A x) - y \|^2 $$

where $A$ is a linear operator and $F$ is a proper, convex and lower-semicontinuous function. This problem represents the so-called primal formulation of the minimization problem, and $x$ is the primal variable of interest that is being recovered.


## License

For distribution purposes, commercial or otherwise, please contact [Ajinkya Kadu](https://ajinkyakadu.github.io) for further information.

## Requirements

The framework has been tested on Matlab 2016b.

## Usage Instructions

The framework includes examples of usage in the test folder.

## Reference

When using this code, please cite the following publication:
```
@article{Kadu2017,
  doi = {10.1109/tci.2016.2640761},
  url = {https://doi.org/10.1109/tci.2016.2640761},
  year = {2017},
  month = {jun},
  publisher = {Institute of Electrical and Electronics Engineers ({IEEE})},
  volume = {3},
  number = {2},
  pages = {305--315},
  author = {Ajinkya Kadu and Tristan van Leeuwen and Wim A. Mulder},
  title = {Salt Reconstruction in Full-Waveform Inversion With a Parametric Level-Set Method},
  journal = {{IEEE} Transactions on Computational Imaging}
}
```
A preprint of the article can be found [here](https://arxiv.org/pdf/1610.00251.pdf)


## Bug Reports

Should you encounter any issues, please do not hesitate to contact [Ajinkya Kadu](mailto:ajinkyakadu125@gmail.com).
