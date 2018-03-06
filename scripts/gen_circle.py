import csv
import math
import random

print("px py c")
for i in range(10):
    for j in range(10):
        d2center = math.sqrt( (i-4.5)**2 + (j-4.5)**2 )
        if d2center <= 4.5 and d2center >= 2.5:
            print(i, j, 1)

