import random


x0 = 0
y0 = 15
xf = 20
yf = 10
h = 0.01
xa = x0

while xa < xf-h:
    xa  = round(xa + h, 2)
    array.append((xa, random.randint(1, 10)))

array.append((xf, yf))
for i in array:
    print(i)