import csv
import math
import random

#base vector
px = 1.
py = 0.
qx = 0.5
qy = math.sqrt(3)/2.

l = 15  # the length of the qudrilateral
r = l*math.sqrt(3)/4.  # the circle inside the mesh quadrilateral
# center of the circle
cx = l*3./4.
cy = r

points = []
for i in range(l):
    for j in range(l):
        tx = px*i+qx*j
        ty = py*i+qy*j

        if (tx-cx)**2 + (ty-cy)**2 < r**2:
            points.append([tx, ty, 1])

# randomly choose 100 points to be category 2
for i in range(int(0.2*len(points))):
    index = random.randint(0, len(points)-1)
    points[index][2] = 2

# choose an inner circle to be category 3
# choose an intercetion circle to be category 4
for point in points:
    if (point[0]-0.8*cx)**2 + (point[1]-cy)**2 < (r/3.)**2:
        point[2] = 3

    if (point[0]-1.1*l)**2 + (point[1]-cy)**2 < (r/3.)**2:
        point[2] = 4

# write to csv file
with open('mesh.csv', "wb") as file:
    writer = csv.writer(file)
    writer.writerow(['px','py','c'])
    for point in points:
        writer.writerow(point)
