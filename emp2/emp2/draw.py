import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal as dcm

# Внутренние узлы
fe = []

with open ("fe.txt") as file:
	for x in file:
		fe.append(dcm(x))

def main():
	plt.grid()

	plt.scatter(fe, [0 for x in fe], s = 20, color = "black")

	plt.show()

if __name__ == "__main__":
	main()