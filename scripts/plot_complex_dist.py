import argparse
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Plot complex histogram')
    parser.add_argument('filename', help='Input complex file')

    args = parser.parse_args()

    x = []
    with open(args.filename, 'r') as f:
        for line in f:
            row = line.split()
            if row[0] == 'number':
                continue
            else:
                x.append(float(row[-1]))

    plt.figure(figsize=(16,4))
    # the histogram of the data
    plt.subplot(121)
    n, bins, patches = plt.hist(x, 50, facecolor='green', cumulative=False)
    plt.xlabel('Distance')
    plt.ylabel('Number of Simplices')
    plt.title('Histogram')
    plt.grid(True)

    plt.subplot(122)
    n, bins, patches = plt.hist(x, 50, facecolor='green', cumulative=True)
    plt.xlabel('Distance')
    # plt.ylabel('Number of Simplices')
    plt.title('Cumulative Histogram')
    plt.grid(True)

    plt.savefig(args.filename+'_hist.png', dpi=300)
