import csv
import sys
import math
import numpy as np

fname = sys.argv[1]

# calculate the (x,y,z) of the center (155, 28)
centerRA_radian = math.radians(155)
centerDEC_radian = math.radians(28)

center_cx = math.cos(centerRA_radian)*math.cos(centerDEC_radian)
center_cy = math.sin(centerRA_radian)*math.cos(centerDEC_radian)
center_cz = math.sin(centerDEC_radian)

def valid(row):
    abs_mag = float(row['absMagR_5logh'])-float(row['kcorrR01'])
    synthR = float(row['synthR'])
    z_photoz = float(row['z_photoz'])
    if abs_mag > -21.1 and abs_mag < -19.6 and  synthR > 10.2 and synthR < 17.72 \
       and z_photoz >= 0.08 and z_photoz < 0.09:
        # distance to (155,28) <= 70 Mpc
        comove_dist_Mpc = float(row['comove_dist_Mpc'])
        x = comove_dist_Mpc*float(row['cx'])
        y = comove_dist_Mpc*float(row['cy'])
        z = comove_dist_Mpc*float(row['cz'])
        center_x = comove_dist_Mpc*center_cx
        center_y = comove_dist_Mpc*center_cy
        center_z = comove_dist_Mpc*center_cz
        a = np.array([x,y,z])
        center_a = np.array([center_x,center_y,center_z])
        dist = np.linalg.norm(a-center_a)
        if dist <= 70.0:
            return True
        else:
            return False
    else:
        return False

validrows = []

with open(fname, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if valid(row) is True:
            validrows.append(row)

    print("{} valid rows.".format(len(validrows)))

with open('filtered_'+fname, 'w', newline='') as csvfile:
    fieldnames = ['ra', 'dec', 'z_photoz', 'x', 'y', 'z']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for r in validrows:
        comove_dist_Mpc = float(r['comove_dist_Mpc'])
        x = comove_dist_Mpc*float(r['cx'])
        y = comove_dist_Mpc*float(r['cy'])
        z = comove_dist_Mpc*float(r['cz'])
        writer.writerow({
            'ra':r['ra'], 'dec':r['dec'], 'z_photoz': r['z_photoz'],
            'x':x, 'y':y, 'z':z
            })
