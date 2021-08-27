# PySDTest

 A Python Package for Stochastic Dominance Tests

PySDTest is a Python 3 implementation of routines for Stochastic Dominance tests. This package was created by Kyungho Lee in collaboration with Yoon-Jae Whang.

This package implements stochastic dominance tests proposed by [Barrett and Donald (2003, BD)](https://doi.org/10.1111/1468-0262.00390), [Linton, Maasoumi and Whang (2005, LMW)](https://ideas.repec.org/a/oup/restud/v72y2005i3p735-765.html), [Linton, Song and Whang (2010, LSW)](https://econpapers.repec.org/article/eeeeconom/v_3a154_3ay_3a2010_3ai_3a2_3ap_3a186-202.htm), [Donald and Hsu (2014, DH)](https://www.tandfonline.com/doi/full/10.1080/07474938.2013.833813). PySDTest also contains stochastic dominance tests by applying the numerical delta method [(Hong and Li (2018, HL))](https://www.sciencedirect.com/science/article/abs/pii/S0304407618300988.html). For details of each model, refer to the corresponding paper. 

## Installation

**PySDTest** is listed on the Python Package Index (PyPI) which is a repository of software for the Python programming language. If Python is installed on your computer, you can install **PySDTest** by entering:

```python
pip install PySDTesT
```

in Windows cmd or Mac (or Linux) terminal. For detailed information about installing a Python package, see [Python Package User Guide in PyPI](https://packaging.python.org/tutorials/installing-packages/).

## Importation

After installing **PySDTest**, an user can import it through the following code:

```python
import pysdtest
```
## Features

- Stochastic dominance (SD) testing via imposing the least favorable case (BD, LMW)
- SD testing via Contact-set approach (LSW)
- SD testing via Selective recentering approach (DH)
- SD testing via Numerical delta method (HL)
- Advanced hypothesis testing (Joint hypotheses testing)
- Plotting (s-th order) CDFs
- Resampling: bootstrap, subsampling, and paired-bootstrap

## Integration with Other Software

Python provides significant advantage regarding interaction with other programming languages and statistical software. Our package can be used in MATLAB, R, Julia, and Stata. 

In MATLAB, an user can implement **PySDTest** as:

```matlab
% Calling PySDTest in MATLAB
py.pysdtest
```

The reticulate library in R and the PyCall package in Julia serve similar feature, which calls Python modules in R and Julia, respectively:

```r
# Calling PySDTest in R

library(reticulate)
pysdtest <- import("pysdtest")
```

```julia
# Calling PySDTest in Julia

using PyCall
pysdtest = pyimport("pysdtest")
```

In Stata, an user can call **PySDTest** by plugging python codes in between Python: and end. For example, the following codes work for importing **PySDTest**:

```python
# /* Calling PySDTest in Stata */
python:
import pysdtest
end
```

In addition, a user can use python command in Stata to activate Python environment in Stata. For more detailed information, please refer to the [Using MATLAB with Python](https://www.mathworks.com/products/matlab/matlab-and-python.html) in MATLAB, [reticulate](https://rstudio.github.io/reticulate/articles/calling_python.html) library in R, [PyCall](https://github.com/JuliaPy/PyCall.jl) package in Julia and [Python integration feature](https://www.stata.com/stata16/python-integration/) in Stata.
