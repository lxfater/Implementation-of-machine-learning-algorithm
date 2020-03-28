from copy import deepcopy
import numpy as np
import math
# 读入文件
readin = open('iris.txt', 'r')
X = []
try:
    readLine = readin.readline()
    while readLine:
        readLine = readLine.strip('\n')
        l = readLine.split(',')
        l.pop()
        l = [  float(x) for x in l ]
        X.append(l)
        readLine = readin.readline()
finally:
    readin.close()

X = np.array(X)
# 距离
def dist(a, b, ax=1):
    return np.linalg.norm(a - b, axis=ax)

def sse(X,clusters,k):
    sum = 0
    for i in range(k):
    # 生成对应中心点的数组
        points = [X[j] for j in range(len(X)) if clusters[j] == i]
    # 对于这些点的集合中每一个点
        temp = 0
        for j in range(len(points)):
    # 根据每一个点到中心点的距离产生距离数组
            temp = temp+math.pow(dist(points[j], C[i],ax=0), 2)
        sum += temp
    return sum

# # 设定参数k
k = 3

# X coordinates of random centroids
C_x = np.random.randint(3, 8, size=k)
# Y coordinates of random centroids
C_y = np.random.randint(2, 4, size=k)
# Z coordinates of random centroids
C_z = np.random.randint(3, 8, size=k)
# K coordinates of random centroids
C_k = np.random.randint(0, 3, size=k)

C = np.array(list(zip(C_x, C_y, C_z, C_k)), dtype=np.float32)
# # 随机选择中心点



# 初始化中心点数组为0
C_old = np.zeros(C.shape)
# 给数据上标签的方法分配数据
clusters = np.zeros(len(X))
# Error func. - Distance between new centroids and old centroids
error = dist(C, C_old, None)

# e = 0.0001
# sse_old = 0
# sse_deta = 1

# Loop will run till the error becomes zero
while error != 0:

    for i in range(len(X)):
        # 根据每一个点到中心点的距离产生距离数组
        distances = dist(X[i], C)
        # 获取距离最小的点的下标
        cluster = np.argmin(distances)
        # 上标签
        clusters[i] = cluster
    # 储存以前的中心点
    C_old = deepcopy(C)
    # 重新计算中心点
    for i in range(k):
        # 生成对应中心点的数组
        points = [X[j] for j in range(len(X)) if clusters[j] == i]
        # 生成平均值
        if len(points) == 0:
            # 空集合任意选一点
            C[i] = np.random.randint(0, np.max(X), size=4)
        else:
            C[i] = np.mean(points, axis=0)
    # 再一次计算中心点偏移量
    error = dist(C, C_old, None)


#计算sse




print(sse(X,clusters,k))
print(C)

    # 再一次计算中心点偏移量