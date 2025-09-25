# DFTatom.jl

这是一个用于原子Hartree-Fock和密度泛函理论（DFT）计算的Julia程序包。

## 项目目标

本项目旨在实现一个能够对孤立原子进行**自旋极化** Hartree-Fock 和 DFT (使用**局域自旋密度近似 LSDA**) 计算的工具。最终目标是能够计算并可视化 **H 和 C 原子** 的低能单粒子能级 (所有占据态和第一个未占据态) 和相应的波函数。

## 项目结构

```
DFTatom.jl/
├── src/
│   └── DFTatom.jl      # 主模块，包含核心函数
├── test/
│   └── runtests.jl     # 测试脚本
├── docs/
│   ├── make.jl
│   └── src/
│       └── index.md    # 项目文档（当前文件）
├── Project.toml        # 项目依赖和配置
└── README.md           # 项目简介
```

## 实现计划 (基于基组和 Roothaan-Hall 方程)

我们将采用在现代量子化学中标准的**基组法**来解决此问题。此方法将微分方程问题转换为矩阵代数问题，更易于在计算机上实现。我们将使用 `GaussianBasis.jl` 这个包来加载预定义的高斯基组，例如 STO-3G。

核心是求解 **Roothaan-Hall 方程**: **F C = S C E**

1.  **加载基组**: 使用 `GaussianBasis.jl` 为给定的原子（例如 H 或 C）加载一个基组（例如 `"STO-3G"`）。这将为我们提供一组基函数 `χ`。

2.  **计算积分**: 在 SCF 循环开始前，计算所有必需的积分。这些积分只依赖于基函数，因此只需计算一次。
    *   **S**: 重叠积分矩阵
    *   **T**: 动能积分矩阵
    *   **V**: 核吸引积分矩阵
    *   **ERI**: 双电子排斥积分张量

3.  **实现 SCF 循环**:
    a. **猜测密度矩阵 (P)**: 对电子密度进行初始猜测，通常使用核心哈密顿顿（`H_core = T + V`）的结果来构建。
    b. **构建 Fock 矩阵 (F)**: Fock 矩阵是系统的有效哈密顿顿。它是核心哈密顿顿与根据当前密度矩阵 `P` 计算出的双电子排斥部分的组合。
    c. **求解 Roothaan-Hall 方程**: 求解广义特征值问题 `F C = S C E`。这将得到一组新的轨道系数 `C` 和轨道能 `E`。
    d. **构建新的密度矩阵**: 使用新的轨道系数 `C` 和电子排布来构建新的密度矩阵 `P`。
    e. **检查收敛**: 如果新的密度矩阵与旧的足够接近（或者能量变化足够小），则计算收敛。否则，返回步骤 (b) 并使用新的密度矩阵重复该过程。

4.  **分离 HF 和 DFT**: HF 和 DFT 的区别在于如何构建 Fock 矩阵的第 (b) 步。HF 使用精确交换，而 DFT 使用交换相关泛函。我们将实现两个不同版本的 `build_fock_matrix` 函数。

## 函数和参数定义 (基组法)

```julia
using GaussianBasis

"""
定义原子（比以前更简单）
"""
struct Atom
    symbol::String # 例如 "H", "C"
    charge::Int
    multiplicity::Int # 自旋多重度 (2S+1)
end

"""
计算所有必需的分子积分
"""
function compute_integrals(basis::BasisSet)
    # S = ...
    # T = ...
    # V = ...
    # ERI = ...
    # return (S=S, T=T, V=V, ERI=ERI)
end

"""
主 SCF 循环函数
"""
function run_scf(atom::Atom, basis_name::String; max_iter=100, tol=1e-7)
    # 1. 加载基组
    # basis = BasisSet(basis_name, atom.symbol)
    
    # 2. 计算积分
    # integrals = compute_integrals(basis)
    
    # 3. 初始猜测
    # H_core = integrals.T + integrals.V
    # F = H_core
    # C, E = solve_eigenproblem(F, integrals.S)
    # P = make_density_matrix(C, n_electrons_up, n_electrons_down)
    
    # 4. 迭代循环
    # for iter in 1:max_iter
    #     F = build_fock_matrix(H_core, integrals.ERI, P) # HF 或 DFT 版本
    #     C, E = solve_eigenproblem(F, integrals.S)
    #     P_new = make_density_matrix(C, ...)
    #     
    #     # 检查收敛
    #     if converged
    #         break
    #     end
    #     P = P_new
    # end
end

"""
从密度矩阵构建 Fock 矩阵
"""
function build_fock_matrix(H_core, ERI, P_up, P_down)
    # J = ... # 从总密度 P_up + P_down 和 ERI 构建 Coulomb 部分
    # K_up = ... # 从 P_up 和 ERI 构建 HF 交换部分
    # K_down = ...
    # F_up = H_core + J - K_up
    # F_down = H_core + J - K_down
    # return F_up, F_down
end

"""
从轨道系数构建密度矩阵
"""
function make_density_matrix(C, n_occupied_up, n_occupied_down)
    # P_up = C_up[:, 1:n_occupied_up] * C_up[:, 1:n_occupied_up]'
    # P_down = ...
    # return P_up, P_down
end

"""
求解广义特征值问题 F C = S C E
"""
function solve_eigenproblem(F, S)
    # ... 使用标准线性代数方法 ...
    # return C, E
end
```