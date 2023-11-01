# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 21:58:43 2023
参考github代码：https://github.com/kellenf/TSP_collection
@author: 张启元
"""
## 20231031 
#修改了源代码，将其变为给定初始点（终点）求最优的路径，目前可以得出一条规划好的路线，但可能不是最优的。
#尤其是在产生新解的方法上有待进一步考虑与改进
##

import random
import math
import numpy as np
import matplotlib.pyplot as plt


class SA(object):
    def __init__(self, num_city, data):
        self.T0 = 1000
        self.Tend = 1e-10
        self.rate = 0.9995
        self.num_city = num_city
        self.scores = []
        self.location = data
        # fires中存每一个个体是下标的list
        self.fires = []
        self.dis_mat = self.compute_dis_mat(num_city, data) #城市距离矩阵
        self.fire = self.greedy_init(self.dis_mat,num_city) #贪心算法返回的初始路径标号
        # 显示初始化后的路径
        init_pathlen = 1. / self.compute_pathlen(self.fire, self.dis_mat)   #贪婪算法求得的路径长度
        init_best = self.location[self.fire]        #贪婪算法求得的初始路径位置
        # 存储存储每个温度下的最终路径，画出收敛图
        self.iter_x = [0]
        self.iter_y = [1. / init_pathlen]
    #贪婪策略求初始路径
    def greedy_init(self,dis_mat, num_city):
        #生成一个列表0,1,...,num_city-1
        rest = [x for x in range(0, num_city-1)]    #剔除最后一个点，即原点
        current = 0
        #去掉列表rest中值为current
        rest.remove(current)
        # 找到一条最近邻路径
        result_one = [] #存储路径点标号
        result_one.append(current)
        while len(rest) != 0:
            #暂存最小值
            tmp_min = math.inf
            tmp_choose = -1
            for x in rest:
                if dis_mat[current][x] < tmp_min:
                    tmp_min = dis_mat[current][x]
                    tmp_choose = x
            #找到距离当前点最近的点，并更新之成为最近的点
            current = tmp_choose
            result_one.append(tmp_choose)   #在路径下标数组不断更新最新找到的最近点
            rest.remove(tmp_choose)         #将已经找到的点剔除
        result_one.append(101)  #默认最后一个点为原点，添加到路径末尾
        return result_one

    # 初始化一条随机路径
    # def random_init(self, num_city):
    #     tmp = [x for x in range(num_city)]
    #     random.shuffle(tmp)
    #     return tmp

    # 计算不同城市之间的距离
    def compute_dis_mat(self, num_city, location):
        dis_mat = np.zeros((num_city, num_city)) 
        for i in range(num_city):
            for j in range(num_city):
                if i == j:
                    dis_mat[i][j] = np.inf
                    continue
                a = location[i]
                b = location[j]
                #tmp = np.sqrt(sum([(x[0] - x[1]) ** 2 for x in zip(a, b)])) #Zip函数 ,计算欧拉距离
                #计算经纬度距离
                tmp = 6370*math.acos(math.cos(a[0]-b[0])*math.cos(a[1])*math.cos(b[1])+math.sin(a[1])*math.sin(b[1]))
                dis_mat[i][j] = tmp
        return dis_mat

    # 计算路径长度
    def compute_pathlen(self, path, dis_mat):
        a = path[-2]
        b = path[-1]
        result = dis_mat[a][b]  #先计算回到原点的距离
        for i in range(len(path) - 2):
            a = path[i]
            b = path[i + 1]
            result += dis_mat[a][b]
        return result

    # 计算一个温度下产生的一个群体的长度
    def compute_paths(self, paths):
        result = []
        #计算100条按照贪婪算法得到的路径长度，返回结果
        for one in paths:
            length = self.compute_pathlen(one, self.dis_mat)
            result.append(length)
        return result

    # 产生一个新的解：随机交换两个元素的位置
    def get_new_fire(self, fire):
        fire = fire.copy()
        #生成1到100间的2个随机数
        # temp = fire[a]
        # fire[a] = fire[b]
        # fire[b] = temp
       #a, b = np.random.choice(t, 2)
        t = [x for x in range(1,len(fire)-1)]
        a, b = np.random.choice(t, 2)
        fire[a:b] = fire[a:b][::-1]
        return fire

    # 退火策略，根据温度变化有一定概率接受差的解
    def eval_fire(self, raw, get, temp):
        len1 = self.compute_pathlen(raw, self.dis_mat)
        len2 = self.compute_pathlen(get, self.dis_mat)
        dc = len2 - len1
        p = max(1e-1, np.exp(-dc / temp))
        if len2 < len1:
            return get, len2
        elif np.random.rand() <= p:
            return get, len2
        else:
            return raw, len1

    # 模拟退火总流程
    def sa(self):
        count = 0
        # 记录最优解
        best_path = self.fire   #
        best_length = self.compute_pathlen(self.fire, self.dis_mat)     #初始最优路径长度

        while self.T0 > self.Tend:
            count += 1
            # 产生在这个温度下的随机解
            tmp_new = self.get_new_fire(self.fire.copy())
            # 根据温度判断是否选择这个解
            self.fire, file_len = self.eval_fire(best_path, tmp_new, self.T0)
            # 更新最优解
            if file_len < best_length:
                best_length = file_len
                best_path = self.fire
            # 降低温度
            self.T0 *= self.rate
            # 记录路径收敛曲线
            self.iter_x.append(count)
            self.iter_y.append(best_length)
            print(count, best_length)
        return best_length, best_path

    def run(self):
        best_length, best_path = self.sa()
        return self.location[best_path], best_length


# 读取数据
# def read_tsp(path):
#     lines = open(path, 'r').readlines()
#     assert 'NODE_COORD_SECTION\n' in lines
#     index = lines.index('NODE_COORD_SECTION\n')
#     data = lines[index + 1:-1]
#     tmp = []
#     for line in data:
#         line = line.strip().split(' ')
#         if line[0] == 'EOF':
#             continue
#         tmpline = []
#         for x in line:
#             if x == '':
#                 continue
#             else:
#                 tmpline.append(float(x))
#         if tmpline == []:
#             continue
#         tmp.append(tmpline)
#     data = tmp
#     return data

#计算距离矩阵测试
# def compute_dis_mat(num_city, location):
#     dis_mat = np.zeros((num_city, num_city)) 
#     for i in range(num_city):
#         for j in range(num_city):
#             if i == j:
#                 dis_mat[i][j] = np.inf
#                 continue
#             a = location[i]
#             b = location[j]
#             #tmp = np.sqrt(sum([(x[0] - x[1]) ** 2 for x in zip(a, b)])) #Zip函数 ,计算欧拉距离
#             #计算经纬度距离
#             tmp = 6370*math.acos(math.cos(a[0]-b[0])*math.cos(a[1])*math.cos(b[1])+math.sin(a[1])*math.sin(b[1]))
#             dis_mat[i][j] = tmp
#     return dis_mat

#贪婪算法求初始路径测试
# def greedy_init(dis_mat, num_total, num_city):
#     start_index = 0
#     #生成一个列表0,1,...,num_city-1
#     rest = [x for x in range(0, num_city-1)]
#     current = start_index
#     #去掉列表rest中值为current
#     rest.remove(current)
#     # 找到一条最近邻路径
#     result_one = [] #存储路径点标号
#     result_one.append(current)
#     while len(rest) != 0:
#         #暂存最小值
#         tmp_min = math.inf
#         tmp_choose = -1
#         for x in rest:
#             if dis_mat[current][x] < tmp_min:
#                 tmp_min = dis_mat[current][x]
#                 tmp_choose = x
#         #找到距离当前点最近的点，并更新之成为最近的点
#         current = tmp_choose
#         result_one.append(tmp_choose)   #在路径下标数组不断更新最新找到的最近点
#         rest.remove(tmp_choose)         #将已经找到的点剔除
#     result_one.append(102)
#     return result_one

f=open(r"TSP_Data.txt")
line = f.readline()
data_list = []
while line :
    num = list(map(float,line.split()))
    data_list.append(num)
    line = f.readline()
f.close()
temp_data = np.array(data_list)
temp_data = np.vstack([temp_data[:,0:2], temp_data[:,2:4],temp_data[:,4:6], temp_data[:,6:8]])    
#不添加索引
data = np.empty((102,2))
data[0,:],data[101,:]= [70,40],[70,40]
data[1:101,0:2] =temp_data


#在第一列添加索引
index_num = np.array(list(range(1, 101)))
data_index = np.empty((102,3))
data_index[1:101,0] = index_num
data_index[1:101,1:3] =temp_data
data_index[0,:]= [0,70,40]
data_index[101,:]= [101,70,40]


#经纬度数据，将角度转换成弧度
data = data*math.pi/180

#show_data = np.vstack([data, data[0]])
Best_path = np.empty((102,3))   #第一列为路径规划后目标点索引，第二、三列为经纬度坐标
Best,Fly_time = math.inf,None
#代入退火算法模型求解
model = SA(num_city=data.shape[0], data=data.copy())
#返回路径坐标与路径
path, path_len = model.run()
Fly_time = path_len/1000     #路径路程/飞行时速1000km/h
print("path length is ",path_len,"fly time is ",Fly_time)
if path_len < Best:
    Best = path_len
    Best_path[:,1:3] = path


#绘制航行路线
#还原出经纬度数据，将弧度转换成角度
Best_path = Best_path*180/(math.pi)
#找到对应的索引下标
for i in range(Best_path.shape[0]):
    for j in range(Best_path.shape[0]):
        if((round(Best_path[i,1],4)==round(data_index[j,1],4)) and (round(Best_path[i,2],4)==round(data_index[j,2],4))):
           Best_path[i,0] = data_index[j,0]
           break 
Best_path[101,0] = 101              
# 加上一行因为会回到起点
#Best_path = np.vstack([Best_path, Best_path[0]])

plt.scatter(Best_path[:, 1], Best_path[:,2])
plt.plot(Best_path[:, 1], Best_path[:, 2])
plt.title('Optimized route results')
#标出起始点
plt.annotate(' Start 1 ', (70, 40), xytext=(56, 39),
              arrowprops=dict(arrowstyle='->')) 
#标出最后一个目标点
plt.annotate(f" Last target {int(Best_path[0,0])}" , (Best_path[0,1],Best_path[0,2]), xytext=(Best_path[0,1]+2, Best_path[0,2]-3),
              arrowprops=dict(arrowstyle='->'))

plt.savefig('Imag/路线_1000_e10_9995.jpg', dpi=600)
plt.show()

iterations = model.iter_x
best_record = model.iter_y
plt.plot(iterations, best_record)
plt.title('Convergence curves')
plt.show()
