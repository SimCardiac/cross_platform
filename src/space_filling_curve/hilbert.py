import numpy as np
import matplotlib.pyplot as plt

def d2xy(n, d):
    """
    将Hilbert曲线上的距离d转换为2D坐标(x,y)
    n: 网格大小为 2^n x 2^n
    d: Hilbert曲线上的距离 (0 到 4^n - 1)
    """
    x = y = 0
    s = 1
    while s < n:
        rx = 1 & (d // 2)
        ry = 1 & (d ^ rx)
        x, y = rot(s, x, y, rx, ry)
        x += s * rx
        y += s * ry
        d //= 4
        s *= 2
    return x, y

def xy2d(n, x, y):
    """
    将2D坐标(x,y)转换为Hilbert曲线上的距离
    n: 网格大小为 2^n x 2^n
    """
    d = 0
    s = n // 2
    while s > 0:
        rx = 1 if (x & s) > 0 else 0
        ry = 1 if (y & s) > 0 else 0
        d += s * s * ((3 * rx) ^ ry)
        x, y = rot(s, x, y, rx, ry)
        s //= 2
    return d

def rot(n, x, y, rx, ry):
    """
    旋转和翻转四象限
    """
    if ry == 0:
        if rx == 1:
            x = n - 1 - x
            y = n - 1 - y
        x, y = y, x
    return x, y

def hilbert_curve(n):
    """
    生成n阶Hilbert曲线
    n: 阶数，生成 2^n x 2^n 的网格
    返回按Hilbert顺序排列的坐标点数组
    """
    size = 2 ** n
    total_points = size * size
    points = np.zeros((total_points, 2), dtype=int)
    
    for d in range(total_points):
        x, y = d2xy(size, d)
        points[d] = [x, y]
    
    return points

# 使用示例
if __name__ == "__main__":
    nn = 4
    order = 2 ** nn
    points = hilbert_curve(nn)
    
    plt.figure(figsize=(8, 8))
    plt.plot(points[:, 0], points[:, 1], 'r-', linewidth=1.5, alpha=0.7)
    
    # 标记起点和终点
    plt.scatter(points[0, 0], points[0, 1], color='green', s=100, label='Start', zorder=5)
    plt.scatter(points[-1, 0], points[-1, 1], color='red', s=100, label='End', zorder=5)
    
    # 添加一些示例点的Morton索引标注
    for i in [0, 1, 2, 3, len(points)-1]:
        plt.annotate(f'{i}', xy=(points[i, 0], points[i, 1]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.title(f'Hilbert Curve (Order {nn})')
    plt.legend()
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()