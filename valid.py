# import numpy as np

# coef = [3.02385, 1.91351, 0.288506, 8.9231]
# coef = [- 3.27011, -0.661767, -2.32141]

# ans = np.roots(coef)

# print(ans)


import numpy as np


def learnnb(X, Y):
    maxfeatval = X.max().max()
    m = len(X)
    n = len(X[0])

    # priorp
    y0_count = np.sum(Y[Y==0])
    y1_count = m - y0_count
    priorp = np.array([y0_count/m, y1_count/m])

    # condp
    condp = np.zeros(np.zeros((n, maxfeatval, 2)))
    for i in range(m):
        for j in range(n):
            label = Y[i]
            condp[j][X[i][j]][label] += 1
    condp[:][:][1] /= y1_count
    condp[:][:][0] /= y0_count
    