# SAConvergenceAnalysis
This library attempts to collect the convergence analyses for ranking purposes from the literature, and the ranking is based explicitly on the sensitivity analysis results.

The library currently collects Position Factor [1, 2], TDCC (top-down coefficient of concordance) with Savage score [3, 4, 7], and Reliability [5].

Sep 18th: Now the library also supports Stat<sub>ranking</sub> created by Sarrazin et al. [8]. The code for using Stat<sub>ranking</sub> is similar as TDCC.

The "Tutorial" Jupyter Notebook shows how to use each of the ranking measurements with the help of SALib Python Library [6].
## Requirements
[Numpy](https://numpy.org/), [Scipy](https://www.scipy.org/), Python 3
#### Optional (In order for the Jupyter Notebook "Tutorial" to work)
[SALib](https://github.com/SALib/SALib)

## Reference
[1]	M. V. Ruano, J. Ribes, A. Seco, and J. Ferrer, “An improved sampling strategy based on trajectory design for application of the Morris method to systems with many input factors,” Environ. Model. Softw., vol. 37, pp. 103–109, 2012, doi: 10.1016/j.envsoft.2012.03.008.

[2]	A. Cosenza, G. Mannina, P. A. Vanrolleghem, and M. B. Neumann, “Global sensitivity analysis in wastewater applications: A comprehensive comparison of different methods,” Environ. Model. Softw., vol. 49, pp. 40–52, 2013, doi: 10.1016/j.envsoft.2013.07.009.

[3]	R. L. Iman and W. J. Conover, “A Measure of Top-Down Correlation,” Technometrics, vol. 29, no. 3, p. 351, Aug. 1987, doi: 10.2307/1269344.

[4]	R. Confalonieri, G. Bellocchi, S. Bregaglio, M. Donatelli, and M. Acutis, “Comparison of sensitivity analysis techniques: A case study with the rice model WARM,” Ecol. Modell., vol. 221, no. 16, pp. 1897–1906, 2010, doi: 10.1016/j.ecolmodel.2010.04.021.

[5]	S. Razavi and H. V. Gupta, “A new framework for comprehensive, robust, and efficient global sensitivity analysis: 2. Application,” Water Resour. Res., vol. 52, no. 1, pp. 440–455, Jan. 2016, doi: 10.1002/2015WR017559.

[6]	J. Herman and W. Usher, “SALib: An open-source Python library for Sensitivity Analysis,” J. Open Source Softw., vol. 2, no. 9, p. 97, Jan. 2017, doi: 10.21105/joss.00097.

[7]	I. R. Savage, “Contributions to the Theory of Rank Order Statistics-the Two-Sample Case,” Ann. Math. Stat., vol. 27, no. 3, pp. 590–615, 1956, [Online]. Available: https://www.jstor.org/stable/2237370.

[8] F. Sarrazin, F. Pianosi, and T. Wagener, “Global Sensitivity Analysis of environmental models: Convergence and validation,” Environ. Model. Softw., vol. 79, pp. 135–152, 2016, doi: 10.1016/j.envsoft.2016.02.005.
