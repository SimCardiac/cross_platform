# 泊松测试的数据模型

## 实体

- TestCase（测试用例）
  - id: string（字符串）
  - description: string（描述）
  - grid: { nx: int, ny: int }（网格）
  - domain: { x0: double, x1: double, y0: double, y1: double }（域）
  - pair: string（对）# (u_exact, f) 的标识符
  - bc: { west: string, east: string, south: string, north: string }（边界条件）# 例如，Dirichlet/Neumann
  - tol: { l2: double, linf: double }（容差）

- Pair（解对）
  - id: string（标识符）
  - u_exact(x,y): function（精确解函数）
  - f(x,y): function（右侧函数）
  - derivs: optional function pointers for boundary derivatives if needed（可选的边界导数函数指针）

- Boundary（边界）
  - type: enum { Dirichlet, Neumann }（类型）
  - value(x,y): function（值函数）# 对于 Dirichlet u_exact
  - flux(x,y): function（通量函数）# 对于 Neumann ∂u/∂n

## 注意事项

- 域默认值：单位正方形 [0,1] × [0,1]
- 网格默认值：为方便起见使用 2 的幂；测试可能使用 (16,16)、(32,32)、(64,64)
- 容差默认值：l2=1e-10、linf=1e-8；可以通过 CLI 为每个 TestCase 覆盖
- 对标识符：poly2、sinpi、cospi
