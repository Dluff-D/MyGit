
# 基本参数
L = 80 # RVE尺寸
fiber = 1 # 纤维种类
void = 3  # 孔隙种类

import numpy as np
R_f = [5] # 纤维半径，几种就填几个，三种就填一个1×3的数组
print(f"纤维半径为：{R_f}")

np.random.seed(23) # 设置随机种子，相同的种子每次运行都会生成相同的随机值，用于调试。完成调试后，将此行代码注释。
R_v = np.random.uniform(0, 2, size=void) # 生成一个随机数组，范围为 [0，2)，区间左开右闭，大小为 1×void，类型为浮点数
R_v = np.round(R_v, 1) # 孔隙半径，浮点数保留一位小数
print(f"孔隙半径为：{R_v}")

R = np.concatenate((R_f, R_v), axis=0) # 按行拼接数组R_f和R_v
print(f"半径为：{R}")

# 体积分数
Num_f=[50] # 纤维的数量
Num_v=[20,30,40] # 孔隙的数量
Num = np.concatenate((Num_f, Num_v), axis=0) # 按行拼接数组Num_f和Num_v
print(f"数量为：{Num}")

area = np.pi * (R ** 2) # 计算圆的面积：pi × R^2
Volume = Num * area / (L ** 2) # 逐元素相乘，并除以 L 的平方
Volume_f = np.sum(Volume[:fiber]) # 前fiber个元素之和，即纤维体积分数
print(f"纤维体积分数为：{Volume_f}")
Volume_v = np.sum(Volume[-void:]) # 后void个元素之和，即孔隙体积分数
print(f"孔隙体积分数为：{Volume_v}")

# 随机位置刷新，确定三个位置的圆的数量的配比
Number = np.sum(Num) # 所有圆的数量之和
print(f"所有圆的数量为：{Number}")
corner = np.random.choice([0, 4]) # 从 [0, 4] 中随机选择一个元素，角点圆数量=coner/4
print(f"角点圆数量为：{corner/4}")

import random
r_min = np.min(R) * 2 # 最小圆的直径
print(f"最小圆的直径为：{r_min}")
cir_max = np.floor( L / r_min * 2 -2) # 两对边界最多能放下的圆的数量，向下取整，数量往小了取。
cir_max = int(cir_max) # 将 cir_max 从浮点数转换为整数
print(f"两对边界最多能放下的圆的数量为：{cir_max}")

random.seed(42)  # 设置 random 的随机种子，之前那个是 numpy 的。完成调试后，将此行代码注释。
edge= random.randrange(0, cir_max, 2) # 从 [0 , cir_max) 中随机选择一个偶数，类型为整数，边界圆数量=edge/2
print(f"边界圆数量为：{edge/2}")

inner = Number-corner/4-edge/2 # 除了角点圆和边界圆，剩下的都是内部圆，内部圆数量=inner
print(f"内部圆数量为：{inner}")

# 画角点圆
import math
Cir = Num # 圆的计数器，每一种半径的圆还剩下多少个没画
non_zero_circle = [i for i, x in enumerate(Cir) if x != 0] # 检查非零元素，我们要画的圆只能是从数量不为0的圆里选
position = np.ones((int(corner + edge + inner), 3)) # 用于存放圆的位置信息，每一行存放一个圆的位置信息：圆心坐标和半径。初始化为(corner+edge+inner)行3列，每画一个圆就修改一行。由于角点圆和边界圆会有衍生圆，所以圆的总数量不等于Number，而是等于corner+edge+inner。
Circle = -1 # 画的圆的序号，从-1开始，这样序号就会从0开始
if corner == 4: # 如果有角点圆，从数量不为0的圆里随机选一个画。逻辑：选确定半径，再确定圆心，最后再补全周期性。
    corner_x = L/2
    corner_y = L/2 # 从右上角开始画

    index = random.choice(non_zero_circle) # 从数量不为0的圆里选到的圆，的索引
    pos_r = R[index] # 从数量不为0的圆里选到的圆，的半径

    # 以右上角为圆心，画一个半径为 pos_r 的圆，在这个圆里（不包括圆的边界）随机选一个点作为角点圆的圆心
    theta = random.uniform(0, 2 * math.pi) # 随机选择一个角度 θ (0 到 2π)，注意区分 random.uniform 和 np.random.uniform，两者用法不同
    radii = random.uniform(0, pos_r) # 随机选择一个半径 radii (0 到 pos_r)
    pos_x = np.round(corner_x + radii * math.cos(theta), 1) # 随机选择的圆心 x 坐标，保留一位小数
    pos_y = np.round(corner_y + radii * math.sin(theta), 1) # 随机选择的圆心 y 坐标，保留一位小数

    Circle = Circle + 1 # 画的圆的序号加1
    position[Circle] = np.array([pos_x, pos_y, pos_r]) # 记录角点圆的位置信息
    Cir[index] = Cir[index] - 1 # 选到的圆的数量减1

    print(f"第{Circle}个圆的种类索引：{index}")
    print(f"第{Circle}个圆的位置信息：圆心在({position[Circle,0]},{position[Circle,1]})，半径为{position[Circle,2]}")
    print(f"第{index}种圆还剩多少个：{Cir[index]}")

    Circle = Circle + 1 # 左上角补一个圆
    position[Circle] = np.array([pos_x - L, pos_y, pos_r]) # 记录左上角的角点圆的位置信息
    Circle = Circle + 1 # 左下角补一个圆
    position[Circle] = np.array([pos_x - L, pos_y - L, pos_r]) # 记录左下角的角点圆的位置信息
    Circle = Circle + 1 # 右下角补一个圆
    position[Circle] = np.array([pos_x, pos_y - L, pos_r]) # 记录右下角的角点圆的位置信息


# 画边界圆
edge_circle = edge // 2 # 边界圆的数量
for i in range(edge_circle): # 画 edge_circle 个边界圆
    non_zero_circle = [i for i, x in enumerate(Cir) if x != 0] # 检查非零元素，我们要画的圆只能是从数量不为0的圆里选
    index = random.choice(non_zero_circle) # 从数量不为0的圆里选到的圆，的索引
    pos_r = R[index] # 从数量不为0的圆里选到的圆，的半径

    if i < edge_circle // 2: # 左右边界画一半数量的边界圆
        pos_x = np.round(random.uniform(-L / 2 - pos_r, -L / 2 + pos_r)) # 随机选择圆心的 x 坐标，保留一位小数
        pos_y = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r)) # 随机选择圆心的 y 坐标，保留一位小数

        # 检查是否与已画的圆重叠
        while True:
            overlap = False
            for j in range(Circle + 1):
                if (pos_x - position[j, 0]) ** 2 + (pos_y - position[j, 1]) ** 2 < (pos_r + position[j, 2]) ** 2: # 判断圆心之间的距离是否小于两个圆的半径之和
                    overlap = True
                    break   
            if overlap: # 如果重叠，重新选择圆心
                pos_x = np.round(random.uniform(-L / 2 - pos_r, -L / 2 + pos_r), 1) # 随机选择圆心的 x 坐标，保留一位小数
                pos_y = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 y 坐标，保留一位小数
            else:
                break
        
        # 不重叠，就在左边界画一个圆
        Circle = Circle + 1 # 已经画好的圆的序号加1
        position[Circle] = np.array([pos_x, pos_y, pos_r]) # 记录边界圆的位置信息
        # 右边界补一个圆
        Circle = Circle + 1 # 已经画好的圆的序号加1
        position[Circle] = [pos_x + L, pos_y, pos_r] # 记录边界圆的位置信息
        Cir[index] = Cir[index] - 1 # 选到的圆的数量减1



    else: # 上下边界再画另一半数量的边界圆（数量不一定要各一半？可以随机，后续还可以调整，先这么写）
        pos_x = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 x 坐标，保留一位小数
        pos_y = np.round(random.uniform( L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 y 坐标，保留一位小数

        # 检查是否有圆与已画的圆重叠
        while True:
            overlap = False
            for j in range(Circle + 1):
                if (pos_x - position[j, 0]) ** 2 + (pos_y - position[j, 1]) ** 2 < (pos_r + position[j, 2]) ** 2:
                    overlap = True
                    break   
            if overlap:
                pos_x = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 x 坐标，保留一位小数
                pos_y = np.round(random.uniform( L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 y 坐标，保留一位小数
            else:
                break
        
        # 不重叠，就在上边界画一个圆
        Circle = Circle + 1 # 已经画好的圆的序号加1
        position[Circle] = [pos_x, pos_y, pos_r] # 记录边界圆的位置信息
        # 下边界补一个圆
        Circle = Circle + 1 # 已经画好的圆的序号加1
        position[Circle] = [pos_x, pos_y - L, pos_r] # 记录边界圆的位置信息
        Cir[index] = Cir[index] - 1

    i = i + 1

'''
    print(f"第{Circle}个圆的种类索引：{index}")
    print(f"第{Circle}个圆的位置信息：圆心在({position[Circle,0]},{position[Circle,1]})，半径为{position[Circle,2]}")
    print(f"第{index}种圆还剩多少个：{Cir[index]}")
'''


# 画内部圆
# 纤维采用RSE算法
planet = Cir[:fiber]
print(f"内部圆要画的纤维个数：{planet}")
Plan = np.sum(planet)
Inner = np.ones((int(Plan),3)) # 行星圆的位置信息，初始化为Plan行3列，每画一个纤维就修改一行

non_zero_planet = [i for i, x in enumerate(planet) if x != 0] # 检查非零元素，我们要画的纤维只能是从数量不为0的纤维里选
index = random.choice(non_zero_planet) # 从数量不为0的纤维里选到的纤维，的索引
center_r = R[index] # 从数量不为0的纤维里选到的纤维，的半径

inner_num = 0 # 行星圆的计数器
center_x = np.round(random.uniform(-L / 10,  L / 10), 1) # 随机选择圆心的 x 坐标，保留一位小数
center_y = np.round(random.uniform(-L / 10,  L / 10), 1) # 随机选择圆心的 y 坐标，保留一位小数

Circle = Circle + 1 # 已经画好的圆的序号加1
position[Circle] = np.array([center_x, center_y, center_r]) # 记录内部圆的位置信息
Cir[index] = Cir[index] - 1 # 选到的圆的数量减1，即纤维的数量减1，因为纤维在前面，与Cir的索引重合，所以索引不用变
print(f"画了第一个纤维，还剩多少个：{Cir[:fiber]}")

l_min = 0.1 # 纤维之间的最小间距
attempt = 0 # 尝试次数
max_attempt = 2000 # 最大尝试次数
i = 0 # 纤维圆的计数器
while i < int(Plan) - 1: # 画 Plan-1 个内部圆，每一个纤维的圆心就会变成下一个行星圆的center_x和center_y。Plan是一个浮点数，需要转换为整数。
    non_zero_planet = [i for i, x in enumerate(planet) if x != 0] # 检查非零元素，我们要画的圆只能是从数量不为0的圆里选
    index = random.choice(non_zero_planet) # 从数量不为0的纤维里选到的纤维，的索引
    pos_r = R[index] # 从数量不为0的纤维里选到的纤维，的半径
    l_max = (center_r + pos_r) / 10 - 4 * l_min # 控制体积分数的关键参数
    
    # 以第一个纤维的圆心为圆心，画一个半径在 lmin+r1+r2 和 lmax+r1+r2 之间的圆环，在这个圆环里（不包括圆环的边界？）随机选一个点作为下一个内部圆的圆心
    theta = random.uniform(0, 2 * math.pi) # 随机选择一个角度 θ (0 到 2π)
    radii = random.uniform(l_min + center_r + pos_r, l_max + center_r +pos_r) # 随机选择一个半径 radii (lmin+r1+r2 到 lmax+r1+r2 之间)
    pos_x = np.round(center_x + radii * math.cos(theta), 1) # 随机选择的圆心 x 坐标，保留一位小数
    pos_y = np.round(center_y + radii * math.sin(theta), 1) # 随机选择的圆心 y 坐标，保留一位小数

    # 检查是否是内部圆，以及是否与已画的圆重叠
    for attempt in range(max_attempt):
        overlap = False
        for j in range(Circle + 1):
            if (pos_x - position[j, 0]) ** 2 + (pos_y - position[j, 1]) ** 2 < (pos_r + position[j, 2]) ** 2:
                overlap = True
                break   # break 只能跳出最内层的循环（for 和 while，if不是循环），所以这里只能跳出 for j in range(Circle + 1) 这个循环
        if overlap or pos_x < -L / 2 + pos_r or pos_x > L / 2 - pos_r or pos_y < -L / 2 + pos_r or pos_y > L / 2 - pos_r: # 如果重叠或者超出范围，重新选择圆心和半径
            non_zero_planet = [i for i, x in enumerate(planet) if x != 0] # 检查非零元素，我们要画的圆只能是从数量不为0的圆里选
            index = random.choice(non_zero_planet) # 从数量不为0的圆里选到的圆，的索引
            pos_r = R[index] # 从数量不为0的圆里选到的圆，的半径
            l_max = (center_r + pos_r) / 10 - 4 * l_min # 重选了半径，所以lmax也要重新计算

            theta = random.uniform(0, 2 * math.pi) # 随机选择一个角度 θ (0 到 2π)
            radii = random.uniform(l_min + center_r + pos_r, l_max + center_r +pos_r) # 随机选择一个半径 radii (lmin+r1+r2 到 lmax+r1+r2 之间)
            pos_x = np.round(center_x + radii * math.cos(theta), 1) # 随机选择的圆心 x 坐标，保留一位小数
            pos_y = np.round(center_y + radii * math.sin(theta), 1) # 随机选择的圆心 y 坐标，保留一位小数
            attempt += 1
        else: 
            Circle = Circle + 1 # 已经画好的圆的序号加1
            position[Circle] = np.array([pos_x, pos_y, pos_r]) # 记录内部圆的位置信息
            Cir[index] = Cir[index] - 1 # 选到的圆的数量减1
            # 更新下一个纤维的圆心
            i = i + 1 # 纤维的计数器加1
            print(f"纤维还剩多少个：{Cir[:fiber]}")
            Inner[i] = np.array([pos_x, pos_y, pos_r]) # 记录纤维的位置信息
            break # 跳出 for attempt in range(max_attempt) 这个循环
    if attempt == max_attempt: # 如果尝试次数达到最大尝试次数，就换下一个行星圆
        print("更换行星圆")
        inner_num = inner_num + 1
        if inner_num > Plan -1 : # 放下纤维的总体空间还是有的，就是被纤维之间分割成小块了，如何突破这一限制？
                                # 或者让边界上的纤维稍微多一点，我没有采用论文中提到的在生成的同时处理边界的周期性，而是先生成了边界的圆，再生成内部的圆。
                                # 是不是完全采用论文的方法就可以解决这个问题？
                                # 可以让纤维移动起来，把空间挪出来。比如让纤维像小球一样受到重力，向下移动，同时保持不重叠。
                                # 如何紧密排列，再随机碰撞？
            '''       
            print("纤维实在画不下了，不信你看")
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            fig, ax = plt.subplots()
            ax.set_xlim(-L/2, L/2)
            ax.set_ylim(-L/2, L/2)
            ax.set_aspect('equal')
            for i in range(len(position)):
                circle = plt.Circle((position[i, 0], position[i, 1]), position[i, 2], edgecolor='blue', fill=False)
                ax.add_artist(circle)
            plt.show()
            ''' 
            break # 跳出 while i < int(Plan) - 1 这个循环        
        center_x = Inner[inner_num, 0]
        center_y = Inner[inner_num, 1]
        center_r = Inner[inner_num, 2]


print(f"RSE算法数量结算——纤维还剩多少个：{Cir[:fiber]}")
zero = np.sum(Cir[:fiber]) # RSE算法下纤维还剩多少个没画完
Cir[:fiber] = 0 # 纤维的数量归零
position = position[:-zero] # 去掉没有画的纤维的位置信息

# 孔隙采用Hard-Core算法
satellite = Cir[-void:]
print(f"孔隙还要画多少个：{satellite}")
Sate = np.sum(satellite)

attempt = 0 # 尝试次数
max_attempt = 1000 # 最大尝试次数

void_num = 0 # 孔隙的计数器
while void_num < Sate: # 画 Sate 个孔隙
    non_zero_circle = [i for i, x in enumerate(Cir) if x != 0] # 检查非零元素，我们要画的圆只能是从数量不为0的圆里选。
    index = max(non_zero_circle, key=lambda idx: R[idx]) # 从数量不为0的圆里选到的半径最大的圆，的索引。
                                                        # 和前面不同，这里是选最大的半径的圆，合理的顺序应该是先放大的圆，再放小的圆。
                                                        # 否则小圆把空间分割得很碎，大圆就没地方放了。
                                                        # 和我们先用RSE算法放纤维，再用Hard-Core算法放孔隙的逻辑是一样的，也是先放大的，再放小的。
    pos_r = R[index] # 从数量不为0的圆里选到的圆，的半径
    pos_x = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 x 坐标，保留一位小数。这个范围是不考虑与边界圆重叠的内部圆的理论范围。
    pos_y = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 y 坐标，保留一位小数
    # 检查是否与已画的圆重叠
    for attempt in range(max_attempt):
        overlap = False
        for j in range(Circle + 1):
            if (pos_x - position[j, 0]) ** 2 + (pos_y - position[j, 1]) ** 2 < (pos_r + position[j, 2]) ** 2: # 判断圆心之间的距离是否小于两个圆的半径之和
                overlap = True
                break   
        if overlap: # 如果重叠，重新选择圆心
            pos_x = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 x 坐标，保留一位小数
            pos_y = np.round(random.uniform(-L / 2 + pos_r,  L / 2 - pos_r), 1) # 随机选择圆心的 y 坐标，保留一位小数
            attempt += 1
        else:
            Circle = Circle + 1 # 已经画好的圆的序号加1
            position[Circle] = np.array([pos_x, pos_y, pos_r])
            Cir[index] = Cir[index] - 1
            void_num = void_num + 1
            print(f"孔隙还剩多少个：{Cir[-void:]}")
            break

print(f"Hard-Core算法数量结算——孔隙还剩多少个：{satellite}")
print(f"两个算法数量结算——所有的圆还剩多少个：{Cir}")

# 画图
import matplotlib.pyplot as plt
import matplotlib.patches as patches
fig, ax = plt.subplots()
ax.set_xlim(-L/2, L/2)
ax.set_ylim(-L/2, L/2)
ax.set_aspect('equal')
for i in range(len(position)):
    circle = plt.Circle((position[i, 0], position[i, 1]), position[i, 2], edgecolor='blue', fill=False)
    ax.add_artist(circle)
plt.show()

