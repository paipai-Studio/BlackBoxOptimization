# BlackBoxOptimization

基于 MoonBit 实现的黑盒优化算法框架，支持多种经典优化算法，提供统一的接口和基准测试框架。

## 支持的算法

| 类别 | 算法 | 说明 |
|------|------|------|
| 单点搜索 | 随机搜索 (Random Search) | 最简单的基线算法 |
| 单点搜索 | 爬山法 (Hill Climbing) | 基于邻域的局部搜索 |
| 单点搜索 | 模拟退火 (Simulated Annealing) | 接受劣解跳出局部最优 |
| 群体智能 | 粒子群优化 (PSO) | 模仿鸟群觅食行为 |
| 进化算法 | 差分进化 (DE) | 基于种群的进化优化 |
| 进化算法 | 遗传算法 (GA) | 选择、交叉、变异 |
| 确定性方法 | Nelder-Mead 单纯形法 | 无导数的几何变形方法 |
| 群体智能 | 萤火虫算法 (Firefly) | 基于萤火虫吸引力模型 |

## 项目结构

```
BlackBoxOptimization/
├── lib/                          # 核心库包
│   ├── moon.pkg                  # 包配置
│   └── lib.mbt                   # 所有优化算法实现
├── main/                         # 主程序入口
│   ├── moon.pkg                  # 包配置
│   └── main.mbt                  # 基准测试入口
├── test/                         # 测试包
│   ├── moon.pkg                  # 包配置
│   └── lib_test.mbt              # 单元测试
├── example/                      # 示例包
│   ├── moon.pkg                  # 包配置
│   └── simple_usage.mbt          # 使用示例
├── .github/workflows/            # CI 配置
│   └── ci.yml                    # GitHub Actions
├── moon.mod.json                 # 模块配置
└── README.md
```

## 快速开始

### 环境要求

- [MoonBit](https://moonbitlang.com/) 0.14.0+

### 安装

```bash
# 克隆项目
git clone https://github.com/your-username/BlackBoxOptimization.git
cd BlackBoxOptimization
```

### 运行示例

```bash
# 运行简单示例
moon run example

# 运行基准测试
moon run main
```

### 运行测试

```bash
# 运行所有测试
moon test

# 代码检查
moon check --all
```

## 使用方法

### 基本用法

```moonbit
import lib

fn main {
  // 1. 定义问题
  let dim = 3
  let bounds = @lib.make_bounds(-5.0, 5.0)
  let rng = @lib.make_rng(123.0)

  // 2. 选择目标函数（Sphere 函数）
  // f(x) = sum(x_i^2)，最优值为 0

  // 3. 配置并运行算法
  let config = @lib.make_pso_config(
    1000,      // max_evaluations
    dim,       // dimension
    bounds,    // bounds
    30,        // population_size
    0.7,       // inertia_weight
    1.5,       // cognitive_coeff
    1.5        // social_coeff
  )
  let result = @lib.pso(config, @lib.sphere_function, rng)

  // 4. 获取结果
  println("最优值: \{result.best_value}")
  println("评估次数: \{result.evaluations}")
}
```

### 使用不同算法

```moonbit
// 随机搜索
let rs_config = @lib.make_random_search_config(1000, dim, bounds)
let rs_result = @lib.random_search(rs_config, f, rng)

// 爬山法
let hc_config = @lib.make_hill_climbing_config(1000, dim, bounds, 0.5, 10)
let hc_result = @lib.hill_climbing(hc_config, f, rng)

// 模拟退火
let sa_config = @lib.make_simulated_annealing_config(1000, dim, bounds, 100.0, 0.995, 0.5)
let sa_result = @lib.simulated_annealing(sa_config, f, rng)

// 粒子群优化
let pso_config = @lib.make_pso_config(1000, dim, bounds, 30, 0.7, 1.5, 1.5)
let pso_result = @lib.pso(pso_config, f, rng)

// 差分进化
let de_config = @lib.make_de_config(1000, dim, bounds, 30, 0.8, 0.9)
let de_result = @lib.differential_evolution(de_config, f, rng)

// 遗传算法
let ga_config = @lib.make_ga_config(1000, dim, bounds, 30, 0.8, 0.05)
let ga_result = @lib.genetic_algorithm(ga_config, f, rng)

// Nelder-Mead
let nm_config = @lib.make_nm_config(1000, dim, bounds, 1.0, 2.0, 0.5, 0.5)
let nm_result = @lib.nelder_mead(nm_config, f, rng)

// 萤火虫算法
let ff_config = @lib.make_firefly_config(1000, dim, bounds, 30, 1.0, 1.0, 0.2)
let ff_result = @lib.firefly_algorithm(ff_config, f, rng)
```

### 内置测试函数

```moonbit
// Sphere 函数: f(x) = sum(x_i^2)，最优值: 0
let f1 = @lib.sphere_function

// Rastrigin 函数: f(x) = 10n + sum(x_i^2 - 10*cos(2*pi*x_i))，最优值: 0
let f2 = @lib.rastrigin_function

// Rosenbrock 函数: f(x) = sum(100*(x_{i+1} - x_i^2)^2 + (1 - x_i)^2)，最优值: 0
let f3 = @lib.rosenbrock_function
```

### 基准测试

```moonbit
// 运行单轮基准测试
let results = @lib.run_benchmark_at_evals(f, "Sphere", dim, bounds, 0.0, 5000, rng)
@lib.print_benchmark_table(results, "Sphere", 5000)

// 运行多次试验统计分析
@lib.run_multiple_trials(f, "Sphere", dim, bounds, 0.0, 10000, 20, 100.0)
```

## 算法参数说明

### PSO (粒子群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 粒子数量 | 30 |
| inertia_weight | 惯性权重 | 0.7 |
| cognitive_coeff | 认知系数 | 1.5 |
| social_coeff | 社会系数 | 1.5 |

### DE (差分进化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 30 |
| mutation_factor | 变异因子 F | 0.8 |
| crossover_probability | 交叉概率 CR | 0.9 |

### GA (遗传算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 30 |
| crossover_rate | 交叉率 | 0.8 |
| mutation_rate | 变异率 | 0.05 |

### Nelder-Mead
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| reflection_coeff | 反射系数 α | 1.0 |
| expansion_coeff | 扩展系数 γ | 2.0 |
| contraction_coeff | 收缩系数 ρ | 0.5 |
| shrink_coeff | 缩小系数 σ | 0.5 |

### Firefly (萤火虫算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 萤火虫数量 | 30 |
| initial_attractiveness | 初始吸引力 β₀ | 1.0 |
| light_absorption | 光吸收系数 γ | 1.0 |
| step_size | 随机步长 α | 0.2 |

## 性能对比（Sphere 函数，5维）

| 算法 | 1000 evals | 5000 evals | 10000 evals |
|------|------------|------------|-------------|
| 随机搜索 | 3.17 | 0.92 | 0.92 |
| 爬山法 | 0.0094 | 0.0094 | 0.0094 |
| 模拟退火 | 0.57 | 0.0059 | 0.0059 |
| **PSO** | 1.05e-3 | **2.3e-14** | **4.3e-28** |
| **DE** | 0.20 | **2.3e-8** | **7.8e-17** |
| GA | 1.00 | 0.24 | 0.05 |
| **Nelder-Mead** | 5.3e-32 | **1.6e-199** | **0** |
| Firefly | 3.39 | 0.11 | 0.02 |

## CI/CD

项目使用 GitHub Actions 进行持续集成，配置文件位于 [.github/workflows/ci.yml](file:///home/developer/BlackBoxOptimization/.github/workflows/ci.yml)。

CI 流程包括：
1. **Check**: `moon check --all` - 代码质量检查
2. **Test**: `moon test` - 运行所有单元测试
3. **Build**: `moon build --release` - 构建发布版本

## 许可证

本项目基于 [Apache-2.0 license](LICENSE) 开源。
