# ParametricLevelSet  

An easy-to-use framework of Parametric Level Set approach for inverse problem to improve the reconstructions.

## Introduction
**ParametricLevelSet** is a MATLAB toolbox for Parametric Level Set approach in inverse problems.

The classic least-squares formulation of a general inverse problem in written as


![Problem formulation][generalPrimalFormulation] ,

where ![A](https://latex.codecogs.com/svg.latex?A) denotes a linear operator and
![F](https://latex.codecogs.com/svg.latex?F) is a proper, convex and lower-semicontinuous function. This Problem refers to the so-called _primal_ formulation of the minimization problem and ![x in R^N](https://latex.codecogs.com/svg.latex?x\in\mathbb{R}^N) is known as the primal variable we are interested in recovering.


[generalPrimalFormulation]: https://latex.codecogs.com/svg.latex?\min_{x}&space;F(Ax) "Problem formulation"
[matA]: https://latex.codecogs.com/svg.latex?A "A"
[funcF]: https://latex.codecogs.com/svg.latex?F" "F"

## Authors
* Ajinkya Kadu ([a.a.kadu@uu.nl](mailto:a.a.kadu@uu.nl))*
* Tristan van Leeuwen ([T.vanLeeuwen@uu.nl](mailto:T.vanLeeuwen@uu.nl))*

\*Mathematical Institute, Utrecht University, The Netherlands

## License
If you plan to distribute the software (commercially or not), please contact [Ajinkya Kadu](https://ajinkyakadu125.github.io) for more information.

## Dependencies
This framework has been tested on Matlab 2016b.


## Usage
The examples can be found in test folder.

## Citation

If you use this code please use the following citation
```
@article{Kadu2017,
  doi = {10.1109/tci.2016.2640761},
  url = {https://doi.org/10.1109/tci.2016.2640761},
  year  = {2017},
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

## Reporting Bugs
In case you experience any problems, please contact [Ajinkya Kadu](mailto:a.a.kadu@uu.nl)
