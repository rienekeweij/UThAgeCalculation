import sys
from itertools import starmap
import math
from math import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# some constants
L_230 = 9.15771e-06
L_234 =2.82341e-06
L_DELTA = L_230 - L_234
L_RATIO = L_230 / (L_230 - L_234)
# Newthon Raphson method
EPSILON = 1E-9
START_POINT = 0.1
MAX_TRIES_UNSTABLE_SOLUTION = 10
# output
OUTPUT_FILENAME= "new_ages.txt"

