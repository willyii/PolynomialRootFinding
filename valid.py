import numpy as np
from numpy.random import randint, uniform

test_file = open("test/validation.test", "w")

for _ in range(10000):
    d = randint(1, 7)
    coef = []
    for _ in range(d+1):
        coef.append(uniform(-10, 10))
    roots = np.roots(coef)
    roots = sorted([ x.real for x in roots if x.imag == 0])
    roots = ["{:.6f}".format(x) for x in roots]
    coef_str = " ".join(map(str, coef))
    root_str = " ".join(map(str, roots))

    test_file.writelines(coef_str + "\n")
    if(len(roots) == 0):
        test_file.writelines("#" + "\n")
    else:
        test_file.writelines(root_str + "\n")
test_file.close()
