from scipy.stats import rankdata
import numpy as np


def ranking(arr, _type = 'Descending'):
    
    rank = len(arr) - rankdata(arr, method = 'max') + 1
    
    if _type == 'Ascending':
        rank = rankdata(arr, method = 'max')
    
    return rank

def ranking_matrix(arr, _type = 'Descending'):
    

    m, n = np.shape(arr)
    arr_rank = np.zeros((m, n))
    for i in range(m):
        arr_rank[i, :] = ranking(arr[i, :], _type)
            
    return arr_rank

