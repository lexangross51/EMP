import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal as dcm

# Внутренние узлы
x_int = []
y_int = []

# Граничные узлы
x_bord = []
y_bord = []

# Фиктивные узлы
x_fict = []
y_fict = []

with open ("internal.txt") as file:
	for coord in file:
		x, y = coord.split(" ")
		x_int.append(dcm(x))
		y_int.append(dcm(y))

with open ("border.txt") as file:
	for coord in file:
		x, y = coord.split(" ")
		x_bord.append(dcm(x))
		y_bord.append(dcm(y))		

with open ("fictitious.txt") as file:
	for coord in file:
		x, y = coord.split(" ")
		x_fict.append(dcm(x))
		y_fict.append(dcm(y))		

def main():
	plt.grid()

	plt.scatter(x_int, y_int, s = 20, color = "c")

	plt.scatter(x_bord, y_bord, s = 20, color = "k")

	plt.scatter(x_fict, y_fict, s = 10, color = "gray")	

	plt.show()

if __name__ == "__main__":
	main()