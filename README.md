# PySDTest

### A Stata/Python Package for Stochastic Dominance Tests

**PySDTest** is a Python 3 implementation of routines for Stochastic Dominance tests. This package was created by Kyungho Lee in collaboration with Yoon-Jae Whang. 

Supported testing procedures include, but are not limited to, tests through the least favorable case approach of [Barrett and Donald (2003, BD)](https://doi.org/10.1111/1468-0262.00390), subsampling approach of [Linton, Maasoumi and Whang (2005, LMW)](https://ideas.repec.org/a/oup/restud/v72y2005i3p735-765.html), contact set approach of [Linton, Song and Whang (2010, LSW)](https://econpapers.repec.org/article/eeeeconom/v_3a154_3ay_3a2010_3ai_3a2_3ap_3a186-202.htm) and sele recentering method of [Donald and Hsu (2014, DH)](https://www.tandfonline.com/doi/full/10.1080/07474938.2013.833813). In addition, the numerical delta method by [Hong and Li (2018, HL)](https://doi.org/10.1016/j.jeconom.2018.06.007
) can be used to compute the critical value. 

**PySDTest** also supports the combination of various resampling methods, test statistics and procedures for approximating the limiting distribution. Furthremore, a user can flexibly conduct a statistical test for advanced/complex hypotheses such as stochastic maximality among distributions.

We also introduce the Stata command **pysdtest** based on the python package for Stata users.

## Paper

https://arxiv.org/abs/2307.10694

The [paper](https://arxiv.org/abs/2307.10694) contains practical guidance for using **PySDTest**. In addition, we briefly give an overview about concepts of stochastic dominance and testing methods. 

If you use our package, please cite the paper:

Lee, Kyungho, and Yoon-Jae Whang. "PySDTest: a Python Package for Stochastic Dominance Tests." arXiv preprint arXiv:2307.10694 (2023).

## Installation

### Python package

**PySDTest** is listed on the Python Package Index (PyPI) which is a repository of software for the Python programming language. If Python is installed on your computer, you can install **PySDTest** by entering:

```python
pip install PySDTest
```

in Windows cmd or Mac (or Linux) terminal. For detailed information about installing a Python package, see [Python Package User Guide in PyPI](https://packaging.python.org/tutorials/installing-packages/).

### Stata command

We also provide  the stata command **pysdtest** that is based on the Python package. To use the command, Stata with version higher than 16.0 and installation of the package **PySDTest** are required. The **pysdtest** module (.ado and .sthlp files) can be installed in Stata by the following command:

```{stata}
 net install pysdtest, from("https://raw.githubusercontent.com/lee-kyungho/pysdtest/main/Stata") replace
```
