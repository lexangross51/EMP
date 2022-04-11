import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import json

def main():
    elems = []
    points = []
    r = []
    z = []

    with open("mesh\mesh.json") as file:
        data = json.load(file)
        points = data['points']
        elems = data['elems']

    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_subplot()

        
    for elem in elems:
        rk = points[elem[0]][0]
        zk = points[elem[0]][1]
        rk1 = points[elem[3]][0]
        zk1 = points[elem[3]][1]
        rect = Rectangle((rk, zk), rk1 - rk, zk1 - zk, edgecolor = 'black', facecolor = 'white')
        ax.add_patch(rect)

        r.append(rk)
        r.append(rk)
        r.append(rk1)
        r.append(rk1)

        z.append(zk)
        z.append(zk1)
        z.append(zk1)
        z.append(zk)


    plt.scatter(r, z, s = 5, color = "black")
    plt.show()


if __name__ == "__main__":
    main()