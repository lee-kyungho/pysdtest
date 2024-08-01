# PySDTest

 A Stata/Python Package for Stochastic Dominance Tests

PySDTest is a Python 3 implementation of routines for Stochastic Dominance tests. This package was created by Kyungho Lee in collaboration with Yoon-Jae Whang.

This package implements stochastic dominance tests, including, but not limited to, [Barrett and Donald (2003, BD)](https://doi.org/10.1111/1468-0262.00390), [Linton, Maasoumi and Whang (2005, LMW)](https://ideas.repec.org/a/oup/restud/v72y2005i3p735-765.html), [Linton, Song and Whang (2010, LSW)](https://econpapers.repec.org/article/eeeeconom/v_3a154_3ay_3a2010_3ai_3a2_3ap_3a186-202.htm), [Donald and Hsu (2014, DH)](https://www.tandfonline.com/doi/full/10.1080/07474938.2013.833813). PySDTest also contains stochastic dominance tests by applying the numerical delta method [(Hong and Li (2018, HL))](https://www.sciencedirect.com/science/article/abs/pii/S0304407618300988.html).

## Paper

https://arxiv.org/abs/2307.10694

The [paper](https://arxiv.org/abs/2307.10694) contains practical guidance for using **PySDTest**. In addition, we briefly give an overview about concepts of stochastic dominance and testing methods. 

If you use our package, please cite the paper:

Lee, Kyungho, and Yoon-Jae Whang. "PySDTest: a Python Package for Stochastic Dominance Tests." arXiv preprint arXiv:2307.10694 (2023).

## Features

- Stochastic dominance (SD) testing via imposing the least favorable case (BD, LMW)
- SD testing via Contact-set approach (LSW)
- SD testing via Selective recentering approach (DH)
- SD testing via Numerical delta method (HL)
- Advanced hypotheses testing (Joint hypotheses testing)
- Plotting (s-th order) CDFs
- Resampling: bootstrap, subsampling, and paired-bootstrap


## Installation

### Python package

**PySDTest** is listed on the Python Package Index (PyPI) which is a repository of software for the Python programming language. If Python is installed on your computer, you can install **PySDTest** by entering:

```python
pip install PySDTest
```

in Windows cmd or Mac (or Linux) terminal. For detailed information about installing a Python package, see [Python Package User Guide in PyPI](https://packaging.python.org/tutorials/installing-packages/).

### Stata command

**pysdtest** is the stata command based on the Python package. To use the command, Stata with version higher than 16.0 and installation of the package **PySDTest** are required. The **pysdtest** module (.ado and .sthlp files) can be installed in Stata by the following command:

```{stata}
 net install pysdtest, from("https://raw.githubusercontent.com/lee-kyungho/pysdtest/main/Stata") replace
```
