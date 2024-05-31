import numpy as np
from numpy import *
import matplotlib as plt

## 1.1 Scalar Variables
"""
a = 3
b = 5.0
print(a)
print("a is: ", a, "\n", "b is: ", b, "\n", "a*b is: ", a*b)

print(a+b, a-b, a*b, a/b, a%b, a**b, a//b)
"""

## 1.2 Array Variables
"""
list_a = ['a', 'b', 'c', 1, 2, 3]
print(
    list_a[0], "\n",
    list_a[3:], "\n",
    list_a[:3], "\n",
    len(list_a), "\n",
    list_a + list_a, "\n",
    list_a * 3
)


x = np.array([0,1,2,3,4,5])
y = np.array([1,2,3,4,5,6])
z = x*y

print("x is: ", x, "\n", "y is: ", y, "\n", "z is: ", z)
print("the dot prod. with x. is: ", x.dot(y), "\n", "the dot prod. with np. is: ", np.dot(x,y),  "\n")
"""

## 1.3 Matrix variables
"""
v = np.array([[1,2,3], [7,11,13]])
w = np.array([[2,3,5], [-1,-2,-3]])
print(
    v.shape,  "\n",
    np.shape(v),  "\n",
    w.shape,  "\n",
    np.shape(w), "\n",
    v.dot(w.T)
)

u = np.array([7,8,9])

print(
    v*u, "\n",
    v.dot(u), "\n",
    u.T
)
"""

## 2 Plotting using Matplotlib
"""
x = np.arange(0,10,0.01)
y = np.sin(x)

figure = plt.plot(x,y)
plt.xlabel('distance [m]')
plt.ylabel('Amplitude')
plt.title('y = sin(x)')
plt.show()
"""

## 3 Task

#loadtxt(RecentIndices.txt)
test = genfromtxt('C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/REU2024_Tutorials/Python_Tutorial_1/RecentIndices.txt')
print(test, "\n", test[0,0])
