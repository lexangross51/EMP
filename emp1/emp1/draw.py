import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal as dcm

# Внутренние узлы
x_int = []
y_int = []

# Граничные узлы
x_bord = []
y_bord = []

# Краевые условия
x_bc1 = []
y_bc1 = []

x_bc2 = []
y_bc2 = []

x_bc3 = []
y_bc3 = []

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

with open ("first.txt") as file:
	for coord in file:
		x, y = coord.split(" ")
		x_bc1.append(dcm(x))
		y_bc1.append(dcm(y))				

with open ("second.txt") as file:
	for coord in file:
		x, y = coord.split(" ")
		x_bc2.append(dcm(x))
		y_bc2.append(dcm(y))						

with open ("third.txt") as file:
	for coord in file:
		x, y = coord.split(" ")
		x_bc3.append(dcm(x))
		y_bc3.append(dcm(y))						

def main():
	plt.grid()

	# Рисуем внутренние точки
	plt.scatter(x_int, y_int, s = 20, color = "c")

	# Рисуем граничные точки
	plt.scatter(x_bord, y_bord, s = 20, color = "black")

	# Рисуем краевые условия
	plt.scatter(x_bc1, y_bc1, s = 20, color = "red")
	plt.scatter(x_bc2, y_bc2, s = 20, color = "green")
	plt.scatter(x_bc3, y_bc3, s = 20, color = "blue")

	plt.show()

if __name__ == "__main__":
	main()