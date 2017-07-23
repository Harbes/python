import numpy as np


class mle:
    def __init__(self, arr):
        if type(arr) is not np.ndarray:
            self.arr = np.array(arr)
        else:
            self.arr = arr
    def newton_raphson(self):











