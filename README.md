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
| 群体智能 | 蚁群算法 (ACO) | 基于信息素的群智能优化 |
| 群体智能 | 蜂群算法 (Bee Colony) | 模拟蜜蜂觅食行为 |
| 局部搜索 | 禁忌搜索 (Tabu Search) | 使用禁忌表避免循环 |
| 局部搜索 | 引导局部搜索 (Guided LS) | 使用惩罚项逃出局部最优 |
| 进化策略 | CMA-ES | 协方差矩阵自适应进化策略 |
| 估计分布 | 交叉熵方法 (CEM) | 基于高斯分布的迭代估计 |
| 群体智能 | 布谷鸟搜索 (Cuckoo Search) | 基于巢寄生和Lévy飞行 |
| 群体智能 | 灰狼优化 (GWO) | 模拟灰狼领导狩猎机制 |
| 确定性方法 | Hooke-Jeeves 模式搜索 | 无导数的直接搜索方法 |
| 群体智能 | 鲸群优化 (WOA) | 模拟座头鲸气泡网捕食行为 |
| 群体智能 | 蝙蝠算法 (Bat Algorithm) | 模拟蝙蝠回声定位行为 |
| 群体智能 | 和声搜索 (Harmony Search) | 模拟音乐即兴创作过程 |
| 群体智能 | 正弦余弦算法 (SCA) | 利用正弦余弦函数振荡搜索 |
| 群体智能 | 飞蛾火焰优化 (MFO) | 模拟飞蛾螺旋飞行路径 |
| 群体智能 | 蚁狮优化 (ALO) | 模拟蚁狮陷阱捕食行为 |

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
moon check
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

// 蚁群算法 (ACO)
let aco_config = @lib.make_aco_config(1000, dim, bounds, 20, 0.1, 1.0, 10.0)
let aco_result = @lib.aco(aco_config, f, rng)

// 蜂群算法 (Bee Colony)
let bc_config = @lib.make_bee_colony_config(1000, dim, bounds, 25, 100)
let bc_result = @lib.bee_colony(bc_config, f, rng)

// 禁忌搜索 (Tabu Search)
let ts_config = @lib.make_tabu_search_config(1000, dim, bounds, 50, 30, 0.5)
let ts_result = @lib.tabu_search(ts_config, f, rng)

// 引导局部搜索 (Guided LS)
let ls_config = @lib.make_local_search_config(1000, dim, bounds, 30, 0.5, 0.5)
let ls_result = @lib.local_search(ls_config, f, rng)

// CMA-ES
let cma_config = @lib.make_cmaes_config(1000, dim, bounds, 0.5, 20, 10)
let cma_result = @lib.cma_es(cma_config, f, rng)

// 交叉熵方法 (CEM)
let cem_config = @lib.make_cem_config(1000, dim, bounds, 100, 0.1, 0.7)
let cem_result = @lib.cross_entropy_method(cem_config, f, rng)

// 布谷鸟搜索 (Cuckoo Search)
let cs_config = @lib.make_cuckoo_config(1000, dim, bounds, 25, 0.25, 1.5)
let cs_result = @lib.cuckoo_search(cs_config, f, rng)

// 灰狼优化 (GWO)
let gwo_config = @lib.make_gwo_config(1000, dim, bounds, 25)
let gwo_result = @lib.grey_wolf_optimizer(gwo_config, f, rng)

// Hooke-Jeeves 模式搜索
let hj_config = @lib.make_hooke_jeeves_config(1000, dim, bounds, 0.5, 0.5, 1.0e-6)
let hj_result = @lib.hooke_jeeves(hj_config, f, rng)

// 鲸群优化 (WOA)
let woa_config = @lib.make_woa_config(1000, dim, bounds, 25, 1.0)
let woa_result = @lib.whale_optimization(woa_config, f, rng)

// 蝙蝠算法 (Bat Algorithm)
let bat_config = @lib.make_bat_config(1000, dim, bounds, 25, 0.0, 1.0, 0.5, 1.0)
let bat_result = @lib.bat_algorithm(bat_config, f, rng)

// 和声搜索 (Harmony Search)
let hs_config = @lib.make_hs_config(1000, dim, bounds, 20, 0.9, 0.3, 0.05)
let hs_result = @lib.harmony_search(hs_config, f, rng)

// 正弦余弦算法 (SCA)
let sca_config = @lib.make_sca_config(1000, dim, bounds, 25)
let sca_result = @lib.sine_cosine_algorithm(sca_config, f, rng)

// 飞蛾火焰优化 (MFO)
let mfo_config = @lib.make_mfo_config(1000, dim, bounds, 25, 1.0)
let mfo_result = @lib.moth_flame_optimization(mfo_config, f, rng)

// 蚁狮优化 (ALO)
let alo_config = @lib.make_alo_config(1000, dim, bounds, 25)
let alo_result = @lib.antlion_optimizer(alo_config, f, rng)
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

### ACO (蚁群算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| ant_count | 蚂蚁数量 | 20 |
| evaporation_rate | 信息素蒸发率 ρ | 0.1 |
| pheromone_weight | 信息素权重 α | 1.0 |
| heuristic_weight | 启发式权重 β | 10.0 |

### Bee Colony (蜂群算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| scout_bee_count | 侦察蜂数量 | 25 |
| max_trials | 最大尝试次数 | 100 |

### Tabu Search (禁忌搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| tabu_list_size | 禁忌表大小 | 50 |
| neighborhood_size | 邻域大小 | 30 |
| step_size | 步长 | 0.5 |

### Guided LS (引导局部搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| neighborhood_size | 邻域大小 | 30 |
| step_size | 步长 | 0.5 |
| lambda | 惩罚系数 λ | 0.5 |

### CMA-ES
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| sigma | 初始标准差 σ | 0.5 |
| lambda | 种群大小 λ | 20 |
| mu | 父代数量 μ | 10 |

### CEM (交叉熵方法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| sample_size | 采样数量 | 100 |
| elite_fraction | 精英比例 | 0.1 |
| smoothing | 平滑系数 | 0.7 |

### Cuckoo Search (布谷鸟搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| nest_count | 巢穴数量 | 25 |
| abandon_probability | 弃巢概率 p_a | 0.25 |
| levy_parameter | Lévy飞行参数 β | 1.5 |

### GWO (灰狼优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| wolf_count | 灰狼数量 | 25 |

### Hooke-Jeeves (模式搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| initial_step | 初始步长 | 0.5 |
| step_reduction | 步长缩减因子 | 0.5 |
| min_step | 最小步长 | 1e-6 |

## 性能对比

### Sphere 函数（5维）

| 算法 | 1000 evals | 5000 evals | 10000 evals |
|------|------------|------------|-------------|
| 随机搜索 | 3.17 | 0.92 | 0.92 |
| 爬山法 | 0.0094 | 0.0094 | 0.0094 |
| 模拟退火 | 0.57 | 0.0059 | 0.0059 |
| **PSO** | **8.56e-4** | **2.3e-15** | **3.4e-29** |
| **DE** | 0.16 | **1.0e-8** | **8.7e-18** |
| GA | 0.70 | 0.056 | 0.011 |
| Nelder-Mead | 49.38 | 49.38 | 49.38 |
| 萤火虫 | 4.25 | 0.92 | 0.92 |
| 蚁群算法 | 2.43 | 2.43 | 2.12 |
| 蜂群算法 | 0.27 | **3.3e-4** | **6.2e-10** |
| 禁忌搜索 | 0.017 | 0.015 | 0.015 |
| 引导局部搜索 | 0.017 | 0.017 | 0.017 |
| CMA-ES | 0.024 | 3.6e-4 | 9.2e-7 |
| CEM | 0.31 | 3.5e-7 | 1.0e-7 |
| 布谷鸟搜索 | 0.74 | 0.064 | 0.064 |
| GWO | 0.0016 | 1.0e-9 | 6.6e-16 |
| **Hooke-Jeeves** | **1.7e-12** | **1.7e-12** | **1.7e-12** |
| WOA | 0.0020 | 9.7e-10 | 9.0e-16 |
| 蝙蝠算法 | 0.0062 | 6.5e-5 | 6.5e-5 |
| 和声搜索 | 0.0028 | 3.4e-5 | 2.3e-5 |
| **SCA** | **1.1e-14** | **9.1e-47** | **1.6e-99** |
| MFO | 13.47 | 13.47 | 0.81 |
| ALO | 2.63 | 2.49 | 2.47 |

### Rastrigin 函数（5维，10000 evals）

| 算法 | 最优值 | 相对误差 |
|------|--------|----------|
| 随机搜索 | 17.38 | 17.38 |
| 爬山法 | 52.75 | 52.75 |
| 模拟退火 | 23.41 | 23.41 |
| PSO | 2.98 | 2.98 |
| DE | 1.35 | 1.35 |
| GA | 1.66 | 1.66 |
| Nelder-Mead | 61.39 | 61.39 |
| 萤火虫 | 17.78 | 17.78 |
| 蚁群算法 | 9.49 | 9.49 |
| 蜂群算法 | 1.11 | 1.11 |
| 禁忌搜索 | 6.91 | 6.91 |
| 引导局部搜索 | 15.64 | 15.64 |
| CMA-ES | 22.88 | 22.88 |
| **CEM** | **2.3e-5** | **2.3e-5** |
| 布谷鸟搜索 | 11.74 | 11.74 |
| GWO | 0.0025 | 0.0025 |
| Hooke-Jeeves | 45.77 | 45.77 |
| **WOA** | **1.3e-11** | **1.3e-11** |
| 蝙蝠算法 | 13.20 | 13.20 |
| 和声搜索 | 0.61 | 0.61 |
| **SCA** | **0** | **0** |
| MFO | 31.25 | 31.25 |
| ALO | 14.39 | 14.39 |

### Rosenbrock 函数（5维，10000 evals）

| 算法 | 最优值 | 相对误差 |
|------|--------|----------|
| 随机搜索 | 22.78 | 22.78 |
| 爬山法 | 5.63 | 5.63 |
| 模拟退火 | 2.74 | 2.74 |
| PSO | 4.09 | 4.09 |
| **DE** | **4.0e-8** | **0** |
| GA | 17.26 | 17.26 |
| Nelder-Mead | 43504.88 | 43504.88 |
| 萤火虫 | 91.21 | 91.21 |
| 蚁群算法 | 7.47 | 7.47 |
| 蜂群算法 | 0.26 | 0.26 |
| 禁忌搜索 | 1.95 | 1.95 |
| 引导局部搜索 | 3.67 | 3.67 |
| CMA-ES | 1.25 | 1.25 |
| CEM | 2.35 | 2.35 |
| 布谷鸟搜索 | 15.13 | 15.13 |
| GWO | 0.039 | 0.039 |
| **Hooke-Jeeves** | **0.29** | **0.29** |
| **WOA** | **0.0022** | **0.0022** |
| 蝙蝠算法 | 130.18 | 130.18 |
| 和声搜索 | 0.68 | 0.68 |
| SCA | 0.0046 | 0.0046 |
| MFO | 6.92 | 6.92 |
| ALO | 183.73 | 183.73 |

## 算法分类特点

| 类别 | 代表算法 | 特点 |
|------|----------|------|
| **单点搜索** | 随机搜索、爬山法、模拟退火 | 简单、易于实现，适合低维问题 |
| **群体智能** | PSO、萤火虫、蚁群、蜂群、布谷鸟、灰狼、鲸群、蝙蝠、和声、正弦余弦、飞蛾火焰、蚁狮 | 基于群体行为，全局搜索能力强 |
| **进化算法** | DE、GA | 基于种群进化，鲁棒性好 |
| **估计分布** | CEM | 基于概率分布估计，适合连续优化 |
| **确定性直接搜索** | Nelder-Mead、Hooke-Jeeves | 无需梯度，适合不可导函数 |
| **局部搜索** | 禁忌搜索、引导局部搜索 | 强化局部开发，避免陷入局部最优 |
| **进化策略** | CMA-ES | 协方差自适应，适合病态问题 |

## 许可证

本项目基于 [Apache-2.0 license](LICENSE) 开源。
