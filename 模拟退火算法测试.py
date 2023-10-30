# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 22:24:30 2023

@author: admin
"""
import numpy as np
import matplotlib.pyplot as plt
def xfunction(x):
    return (x**2-5*x)*np.sin(x**2)

#初始参数
T = 100 #初始温度
T_end = 0.000001 #终止温度
coldrate = 0.999 #冷却速率：T =100 ,T1 = T*0.999 = 99.9 ，T2=T1*0.999=... Tn > T_end
#生成区间0~5的服从均匀分布的随机数
x= np.random.uniform(0,5)
 
while T>T_end:
    y = xfunction(x)
    
    x_new = x + np.random.uniform(-1,1) #添加扰动
    
    if 0 <= x_new <= 5:
        y_new = xfunction(x_new)
        if y_new <y: #min说明这是个好解
            x = x_new
        else:
            p = np.exp(-(y_new-y)/ T)
            r = np.random.uniform(0,1)
            if p > r:
                x =x_new
    T = T * coldrate
x1 = [i for i in np.linspace(0,5)]
y1 = map(xfunction,x1)
plt.plot(x1, list(y1))
plt.plot(x,xfunction(x),'om')
plt.show()
print(x,xfunction(x))



