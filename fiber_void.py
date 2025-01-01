
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

# 随机位置刷新
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
if corner == 4: # 如果有角点圆




