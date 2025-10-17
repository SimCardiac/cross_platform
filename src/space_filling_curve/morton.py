import numpy as np
import matplotlib.pyplot as plt

def morton_encode(x, y):
    """
    将2D坐标(x,y)编码为Morton码(Z-order)
    通过交错x和y的二进制位
    """
    z = 0
    for i in range(32):  # 假设32位整数
        z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1))
    return z

def morton_decode(z):
    """
    将Morton码解码为2D坐标(x,y)
    """
    x, y = 0, 0
    for i in range(32):
        x |= (z & (1 << (2 * i))) >> i
        y |= (z & (1 << (2 * i + 1))) >> (i + 1)
    return x, y

def morton_curve(n):
    """
    生成n阶Morton曲线(Z-order curve)
    n: 阶数，生成 2^n x 2^n 的网格
    """
    size = 2 ** n
    total_points = size * size
    points = np.zeros((total_points, 2), dtype=int)
    
    # 遍历所有网格点，按Morton顺序排列
    for x in range(size):
        for y in range(size):
            z = morton_encode(x, y)
            if z < total_points:
                points[z] = [x, y]
    
    return points

def xy2morton(x, y):
    """
    直接将坐标(x,y)转换为Morton码
    优化版本：使用位操作
    """
    def part1by1(n):
        """将整数的位分散开，每隔一位插入0"""
        n &= 0x0000ffff
        n = (n | (n << 8)) & 0x00FF00FF
        n = (n | (n << 4)) & 0x0F0F0F0F
        n = (n | (n << 2)) & 0x33333333
        n = (n | (n << 1)) & 0x55555555
        return n
    
    return part1by1(x) | (part1by1(y) << 1)

# 使用示例
if __name__ == "__main__":
    n = 4
    order = 2 ** n
    points = morton_curve(n)
    
    plt.figure(figsize=(8, 8))
    plt.plot(points[:, 0], points[:, 1], 'b-', linewidth=1.5, alpha=0.7)
    
    # 标记起点和终点
    plt.scatter(points[0, 0], points[0, 1], color='green', s=100, label='Start', zorder=5)
    plt.scatter(points[-1, 0], points[-1, 1], color='red', s=100, label='End', zorder=5)
    
    plt.title(f'Morton Curve / Z-Order Curve (Order {n})')
    plt.legend()
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()