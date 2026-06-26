# 文档构建与发布指南

本项目使用 [MkDocs](https://www.mkdocs.org/) + [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) 构建文档站点，通过 GitHub Pages 发布。

## 环境要求

- Python 3.8+
- 项目虚拟环境已配置（`.venv/`）

## 安装

```bash
# 进入项目根目录
cd cross_platform

# 激活虚拟环境并安装依赖
.venv/bin/pip install mkdocs mkdocs-material
```

## 本地预览

```bash
# 启动开发服务器（支持热重载）
.venv/bin/mkdocs serve
```

浏览器打开 `http://127.0.0.1:8000`，修改文档后自动刷新。

## 发布到 GitHub Pages

```bash
# 构建并推送到 gh-pages 分支
.venv/bin/mkdocs gh-deploy
```

执行后文档将在 1-2 分钟内上线：

🔗 **https://SimCardiac.github.io/cross_platform/**

## 配置文件

站点配置位于项目根目录的 `mkdocs.yml`：

```yaml
site_name: Cross-Platform PDE Solver Framework
theme:
  name: material            # Material for MkDocs 主题
  features:
    - navigation.instant    # 即时导航
    - navigation.tabs       # 顶部标签导航
    - content.code.copy     # 代码块复制按钮
  palette:
    - scheme: default       # 亮色模式
    - scheme: slate         # 暗色模式

markdown_extensions:
  - pymdownx.arithmatex     # LaTeX 数学公式
  - pymdownx.highlight      # 代码高亮
  - pymdownx.superfences    # 增强围栏代码块

extra_javascript:
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js  # MathJax 渲染
```

## 文档结构

```
docs/
├── index.md                          # 首页
├── project_overview.md               # 项目概述
├── installation.md                   # 安装指南
├── documentation_guide.md            # 本文档
├── 01_poisson_vertex_centered.md     # 顶点 Poisson
├── 02_poisson_cell_centered.md       # 单元 Poisson
├── 03_poisson_staggered.md           # 交错 Poisson
├── 04_heat_vertex_centered.md        # 顶点热传导
├── 05_heat_cell_centered.md          # 单元热传导
├── 06_heat_staggered.md              # 交错热传导
├── 07_stokes_projection.md           # Stokes 投影法
├── 08_stokes_monolithic.md           # Stokes 单片法
├── poisson_neumann_bc.md             # Neumann 边界条件
├── poisson_compatibility_condition.md # 相容性条件
├── poisson_convergence_results.md    # Poisson 收敛结果
└── heat_convergence_results.md       # 热传导收敛结果
```

## 添加新页面

1. 在 `docs/` 目录下创建 `.md` 文件
2. 在 `mkdocs.yml` 的 `nav` 部分添加条目：

```yaml
nav:
  - 新章节:
      - 新页面: new_page.md
```

## 数学公式

使用标准 LaTeX 语法，由 MathJax 渲染：

```latex
行内公式: $-\Delta u = f$

块级公式:
$$
u(x,y) = \sin(\pi x) \sin(\pi y)
$$

带编号的公式:
$$
\begin{equation}
\frac{\partial u}{\partial t} = \alpha \Delta u
\label{eq:heat}
\end{equation}
$$
```

## 代码块

支持语法高亮和复制按钮：

````markdown
```cpp
PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt,
                               const ManufacturedSolution &mms,
                               const BoundaryCondition &bc) {
  // ...
}
```
````

## 主题特性

Material for MkDocs 支持：

| 特性 | 说明 |
|------|------|
| 亮色/暗色切换 | 右上角按钮 |
| 代码复制 | 悬停代码块出现复制按钮 |
| 即时搜索 | `Ctrl+K` 或 `/` |
| 导航标签 | 顶部标签页 |
| 响应式布局 | 自适应移动端 |
