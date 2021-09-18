import numpy as np
from .utils import ranking, ranking_matrix

'''
This library attempts to collect convergence analyses for ranking purposes from the literature.
The ranking is based explicitly on the sensitivity analysis results.
For the universal setting, the most important model parameter, which has the highest sensitivity, has the rank of 1.
The least important model parameter, which has the lowest sensitivity, has the rank of k.
k is the total number of model parameters in the targeted model.

Reference

[1] M. V. Ruano, J. Ribes, A. Seco, and J. Ferrer, “An improved sampling strategy based on trajectory design for 
    application of the Morris method to systems with many input factors,” Environ. Model. Softw., vol. 37, 
    pp. 103–109, 2012, doi: 10.1016/j.envsoft.2012.03.008.
[2] A. Cosenza, G. Mannina, P. A. Vanrolleghem, and M. B. Neumann, “Global sensitivity analysis in wastewater 
    applications: A comprehensive comparison of different methods,” Environ. Model. Softw., vol. 49, pp. 40–52,
    2013, doi: 10.1016/j.envsoft.2013.07.009.
[3] R. L. Iman and W. J. Conover, “A Measure of Top-Down Correlation,” Technometrics, vol. 29, no. 3, p. 351, 
    Aug. 1987, doi: 10.2307/1269344.
[4] R. Confalonieri, G. Bellocchi, S. Bregaglio, M. Donatelli, and M. Acutis, “Comparison of sensitivity analysis
    techniques: A case study with the rice model WARM,” Ecol. Modell., vol. 221, no. 16, pp. 1897–1906, 2010,
    doi: 10.1016/j.ecolmodel.2010.04.021.
[5] S. Razavi and H. V. Gupta, “A new framework for comprehensive, robust, and efficient global sensitivity
    analysis: 2. Application,” Water Resour. Res., vol. 52, no. 1, pp. 440–455, Jan. 2016, doi: 10.1002/2015WR017559.
[6] I. R. Savage, “Contributions to the Theory of Rank Order Statistics-the Two-Sample Case,” Ann. Math. Stat.,
    vol. 27, no. 3, pp. 590–615, 1956, [Online]. Available: https://www.jstor.org/stable/2237370.
[7] F. Sarrazin, F. Pianosi, and T. Wagener, “Global Sensitivity Analysis of environmental models: Convergence and
    validation,” Environ. Model. Softw., vol. 79, pp. 135–152, 2016, doi: 10.1016/j.envsoft.2016.02.005.
   
'''

def position_factor(S1, S2):
    '''

    Parameters
    ----------
    S1 : numpy.array
        The sensitivity measures of model parameters obtained by a specific sample set    
    
    S2 : numpy.array
        The sensitivity measures of model parameters obtained by another sample set

    Returns
    -------
    float
        A scalar that indicates the agreement on the ranking between S1 and S2
        if this scalar is zero, the ranking results are the same 
        
    
    Position Factor was initially proposed by Ruano et al. [1], and it was modified by Cosenza et al. [2]
    to take the absolute value of the difference between two ranking results.

    '''
    
    rank1 = ranking(S1)
    rank2 = ranking(S2)
    
    PF = np.abs(rank1 - rank2) / ((rank1 + rank2) /2.)
    
    return np.sum(PF)


def SSTDCC(S_resample):
    '''
    
    Parameters
    ----------
    S_resample : numpy.array
        An array with a size of R by k, where R is the number of bootstrap resamples,
        and k is the number of model parameters.
        This array consists of sensitivity measures obtained from every single bootstrap resample.

    Returns
    -------
    TDCC : float
        A scalar that indicates the agreement of ranking between every possible pair of bootstrap resamples
        if this scalar is close to 1, the experiment has a high reproducibility under the same setting.
        
        
    The top-down coefficient of concordance (TDCC) is capable of measuring the level of agreement between the
    ranking from either multiple sensitivity analysis methods or multiple bootstrap resamples [3].
    The savage score [6] is later employed with TDCC by Confalonieri et al. [4].

    ''' 
    
    num_resamples, k = np.shape(S_resample)
    ss_matrix = np.zeros((num_resamples, k))
    for j in range(num_resamples):
        ss_matrix[j, :] = ranking(S_resample[j, :])
        for i in range(k):
            ss_matrix[j, i] = ss(ss_matrix[j, i], k)
            
    TDCC = np.sum((np.sum(ss_matrix, axis = 0) ** 2)) - (num_resamples**2 * k)
    
    TDCC = TDCC / (num_resamples**2 * (k - ss(1, k)))
    
    return TDCC

def ss(rank, k):
    ss_value = 0.
    for i in range(int(rank), k + 1):
        ss_value += 1./i
    return ss_value


def Reliability(S1, S_resample):
    '''

    Parameters
    ----------
    S1 : numpy.array
        The sensitivity measures of model parameters obtained by the original sample set (before bootstrap resample)
    S_resample : numpy.array
        An array with the size of R by k, where R is the number of bootstrap resamples,
        and k is the number of model parameters.
        This array consists of sensitivity measures obtained from every single bootstrap resample.

    Returns
    -------
    numpy.array
        An array with the size of 1 by k, where k is the number of model parameters.
        Each index of the array is a scalar indicating the reliability of the corresponding model parameter. 
        
    The measurement of reliability was invented by Razavi and Gupta [5]. It provides a percentage for every single model
    parameter to show how many resamples give the same ranking as the original sample set.
    '''
    S1_rank = ranking(S1)
    S_resamples_rank = ranking_matrix(S_resample)
    
    m, n = np.shape(S_resamples_rank)
    count = np.zeros(n)
    for i in range(n):
        count[i] = (S_resamples_rank[:, i] == S1_rank[i]).sum()
    
    return count/float(m)

def Sarrazin_stat_ranking(S_resample, conf_level = 95):
    '''

    Parameters
    ----------
    S_resample : numpy.array
        An array with the size of R by k, where R is the number of bootstrap resamples,
        and k is the number of model parameters.
        This array consists of sensitivity measures obtained from every single bootstrap resample.
    conf_level : Int, optional
        The desired confidence interval. The default is 95.

    Returns
    -------
    rho_stat : float
        A scalar that indicates the convergence status of the ranking results by checking
        the agreement between the resamples. The ranking result is considered to be converged
        when the returned value is less than 1. If the value is 1, the difference of ranking in
        resamples for the most sensitive model parameters is less or equal to one position on
        average.
        
    This measurement is based on the work of Sarrazin et al [7], and it is a modified version of 
    the Spearman's rank correlation coefficient. 

    '''
    num_resamples, k = np.shape(S_resample)
    item = num_resamples * (num_resamples - 1) / 2
    rho_arr = np.zeros(int(item))

    l = 0
    
    for i in range(num_resamples - 1):
        for j in range(i + 1, num_resamples):
            rho_arr[l] = Sarrazin_ranking(S_resample[i, :], S_resample[j, :], k)
            l += 1

    rho_stat = np.percentile(rho_arr, conf_level, interpolation = 'linear')
    
    return rho_stat

def Sarrazin_ranking(S_j, S_k, k):
    
    R_j = ranking(S_j)
    R_k = ranking(S_k)
    
    Y1 = 0.
    Y2 = 0.
    for i in range(k):
        Y1 = Y1 + np.abs(R_j[i] - R_k[i]) * (np.max([S_j[i], S_k[i]])**2)
        Y2 = Y2 + np.max([S_j[i], S_k[i]])**2
    rho = Y1 / Y2
    
    return rho
