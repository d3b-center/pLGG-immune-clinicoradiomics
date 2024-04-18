import time

import scipy.stats as st
import numpy as np


def timeit(s):
    e = time.time()
    return time.strftime("%H:%M:%S", time.gmtime(e - s))


def ci(scores, c=0.95):
    u = scores.mean()
    l, h = st.norm.interval(c, loc=u, scale=st.sem(scores))
    return c, u, l, h


def bins(x):
    return round(1 + (3.322 * np.log10(x)))
    