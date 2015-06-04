# PyHLA (v1.0.0)
Python for HLA analysis: summary, association analysis, zygosity test and interaction test

# Requirement

* [Python 2.7](https://www.python.org/)
* [pandas](http://pandas.pydata.org/)
* [SciPy](http://www.scipy.org/)
* [StatsModels](http://statsmodels.sourceforge.net/)
* [PyQt4](https://wiki.python.org/moin/PyQt4) (If you want to use the GUI)

# Installation

The easiest way to install Python and the required packages: install **FREE** scientific python distributions such as [Anaconda](http://continuum.io/downloads) and [Enthought Canopy](https://www.enthought.com/products/canopy/) which are already integrated the core scientific analytic and scientific Python packages such as `SciPy`, `pandas`, `StatsModels` and `PyQt4`.

In case you want to install all package by yourself, you can try the following commands:
```
sudo pip install pandas
sudo pip install git+http://github.com/scipy/scipy/
sudo pip install statsmodels
```
Please follow this [guild](http://pyqt.sourceforge.net/Docs/PyQt4/installation.html) to install PyQt4, if you want to use the GUI. 

# Getting Started

The latest PyHLA is available [here](https://github.com/felixfan/PyHLA/archive/v1.0.0.tar.gz).

or, you can clone this repository via the command

```
git clone https://github.com/felixfan/PyHLA.git
```

Once you have installed PyHLA as well as the required packages, typing

```
$ python PyHLA.py -h
```

will print a list of all command-line options. Tutorials of PyHLA can be found in the [wiki](https://github.com/felixfan/PyHLA/wiki). 
