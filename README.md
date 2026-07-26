# BlackBoxOptimization

基于 MoonBit 实现的黑盒优化算法框架，支持 75 种经典优化算法，提供统一的接口和基准测试框架。

## 算法分类

本框架实现了 **75 种**黑盒优化算法，按设计灵感和机制分为以下类别：

### 🔍 单点搜索算法
最简单的优化策略，从单个初始解出发，通过迭代改进逐步逼近最优值。

| 算法 | 简称 | 说明 |
|------|------|------|
| 随机搜索 | Random Search | 随机采样搜索空间中的点作为基线对比 |
| 爬山法 | Hill Climbing | 从当前解的邻域中选择更优解进行迭代 |
| 模拟退火 | Simulated Annealing | 以概率接受劣解，跳出局部最优陷阱 |

> **通用参数**: 所有算法均包含 `max_evals`（最大评估次数）、`dim`（维度）、`bounds`（搜索边界）三个基础参数。

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| Random Search | 无需额外参数 | - |
| Hill Climbing | 无需额外参数 | - |
| Simulated Annealing | `initial_temp`（初始温度）, `cooling_rate`（降温速率）, `step_size`（步长） | 100.0, 0.995, 0.5 |

### 🔧 局部搜索算法
强化局部开发能力，通过记忆或惩罚机制避免陷入循环。

| 算法 | 简称 | 说明 |
|------|------|------|
| 禁忌搜索 | Tabu Search | 使用禁忌表记录近期操作，避免循环 |
| 引导局部搜索 | Guided LS | 添加惩罚项引导搜索逃出局部最优 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| Tabu Search | `tabu_list_size`（禁忌表大小）, `neighborhood_size`（邻域大小）, `step_size`（步长） | 50, 30, 0.5 |
| Guided LS | `neighborhood_size`（邻域大小）, `step_size`（步长）, `lambda`（惩罚系数） | 30, 0.5, 0.5 |

### 🧬 进化算法
模拟自然选择与遗传过程，基于种群的进化搜索机制。

| 算法 | 简称 | 说明 |
|------|------|------|
| 差分进化 | DE | 通过差分变异和交叉实现种群进化 |
| 遗传算法 | GA | 选择、交叉、变异三阶段模拟自然选择 |
| 帝国主义竞争算法 | ICA | 模拟帝国殖民、同化和竞争机制 |
| 生物地理学优化 | BBO | 基于栖息地适宜度与物种迁移模型 |
| 回溯搜索算法 | BSA | 基于历史种群记忆的定向交叉变异 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| DE | `population_size`, `mutation_factor` (F), `crossover_probability` (CR) | 30, 0.8, 0.9 |
| GA | `population_size`, `crossover_rate`, `mutation_rate` | 30, 0.8, 0.05 |
| ICA | `population_size`, `imperialist_count`, `assimilation_rate`, `revolution_rate` | 25, 3, 0.3, 0.3 |
| BBO | `population_size`, `mutation_rate` | 25, 0.01 |
| BSA | `population_size` | 25 |

### 🐝 群体智能算法
本框架最大的算法类别，模拟自然界生物群体行为进行协同搜索。

| 算法 | 简称 | 说明 |
|------|------|------|
| 粒子群优化 | PSO | 模仿鸟群觅食的群体协同行为 |
| 萤火虫算法 | Firefly | 基于萤火虫发光吸引力模型 |
| 蚁群算法 | ACO | 基于信息素的蚂蚁路径搜索 |
| 蜂群算法 | Bee Colony | 模拟蜜蜂觅食与舞蹈交流行为 |
| 布谷鸟搜索 | Cuckoo Search | 基于巢寄生和Lévy飞行的搜索 |
| 灰狼优化 | GWO | 模拟灰狼领导的三层狩猎机制 |
| 鲸群优化 | WOA | 模拟座头鲸气泡网捕食行为 |
| 蝙蝠算法 | Bat Algorithm | 模拟蝙蝠回声定位行为 |
| 和声搜索 | Harmony Search | 模拟音乐即兴创作过程 |
| 正弦余弦算法 | SCA | 利用正弦余弦函数进行振荡搜索 |
| 飞蛾火焰优化 | MFO | 模拟飞蛾螺旋飞行路径 |
| 蚁狮优化 | ALO | 模拟蚁狮陷阱捕食行为 |
| 麻雀搜索算法 | SSA | 模拟麻雀觅食和反捕食行为 |
| 哈里斯鹰优化 | HHO | 模拟哈里斯鹰协同捕食行为 |
| 花粉传播算法 | FPA | 模拟花粉全局和局部传播过程 |
| 蜻蜓算法 | DA | 模拟蜻蜓分离、对齐、聚集行为 |
| 乌鸦搜索算法 | CSA | 模拟乌鸦记忆和追逐行为 |
| 樽海鞘群算法 | SalpSA | 模拟樽海鞘链状群体结构 |
| Monarch蝴蝶优化 | MBO | 模拟蝴蝶迁徙和调整行为 |
| 人工生态系统优化 | AEO | 模拟生态系统生产者-消费者-分解者 |
| 入侵杂草优化 | IWO | 模拟杂草生长扩散与竞争 |
| 细菌觅食优化 | BFO | 模拟大肠杆菌觅食行为 |
| 混合蛙跳算法 | SFLA | 模拟青蛙群体觅食与信息共享 |
| 萤火虫群优化 | GSO | 模拟萤火虫发光吸引行为 |
| 动物迁徙优化 | AMO | 模拟动物群体迁徙跟随行为 |
| 磷虾群算法 | KH | 模拟磷虾诱导运动、觅食和扩散 |
| 共生生物搜索 | SOS | 模拟生物间互利共生/偏利共生/寄生关系 |
| 大象牧群优化 | EHO | 模拟大象族群分离与母系领导机制 |
| 高斯骨干粒子群优化 | GBPSO | 基于高斯分布的简化粒子群算法 |
| 飞蛾搜索算法 | MSA | 基于飞蛾 Lévy 飞行和螺旋搜索 |
| 路径finder算法 | PFA | 基于领导者-跟随者层级移动 |
| 菌落捕食算法 | CPA | 模拟捕食-逃跑概率行为 |
| 细菌菌落优化 | BCO | 基于细菌趋化聚集行为 |
| 飞狐优化算法 | FFOA | 带惯性权重的群体飞行搜索 |
| 人工鱼群算法 | AFSA | 模拟鱼群觅食、聚群和追尾行为 |
| 烟花算法 | FWA | 模拟烟花爆炸及火花分布机制 |
| 鸽群优化 | PIO | 模拟信鸽导航的地图指南针与地标两阶段 |
| 蜜獾算法 | HBA | 模拟蜜獾挖洞与觅食的密度因子搜索 |
| 鹈鹕优化算法 | POA | 模拟鹈鹕捕猎的两阶段策略 |
| 金枪鱼群优化 | TSA | 基于螺旋运动和觅食行为的群体搜索 |
| 金豺优化 | GJO | 模拟雄雌金豺协同捕猎的双引导搜索 |
| 蝗虫优化算法 | GOA | 模拟蝗虫群体觅食的吸引-排斥机制 |
| 海鸥优化算法 | SOA | 模拟海鸥迁徙和螺旋攻击行为 |
| 猫群优化 | CSO | 模拟猫的追踪和搜寻两种模式切换 |
| 大猩猩部队优化 | GTO | 模拟大猩猩群体社会等级与领导机制 |
| 黏菌算法 | SMA | 模拟黏菌觅食行为的群体路径搜索 |
| 苍鹰优化器 | AO | 模拟苍鹰捕猎的四个阶段策略 |
| 黑猩猩优化算法 | ChOA | 模拟黑猩猩狩猎的驱动-追逐-竞争行为 |
| 蒲公英优化器 | DO | 模拟蒲公英种子风力传播的上升-下降搜索 |
| 斑马优化算法 | ZOA | 模拟斑马群体探索-开发机制 |

> **通用群体智能参数**: 大部分算法默认 `population_size = 25`。下表列出各算法的特有参数。

| 算法 | 特有参数 | 推荐值 |
|------|----------|--------|
| PSO | `inertia_weight`, `cognitive_coeff`, `social_coeff` | 0.7, 1.5, 1.5 |
| Firefly | `initial_attractiveness` (β₀), `light_absorption` (γ), `step_size` (α) | 1.0, 1.0, 0.2 |
| ACO | `ant_count`, `evaporation_rate` (ρ), `pheromone_weight` (α), `heuristic_weight` (β) | 20, 0.1, 1.0, 10.0 |
| Bee Colony | `scout_bee_count`, `max_trials` | 25, 100 |
| Cuckoo Search | `nest_count`, `abandon_probability` (p_a), `levy_parameter` (β) | 25, 0.25, 1.5 |
| Bat Algorithm | `freq_min`, `freq_max`, `loudness`, `pulse_rate` | 0.0, 1.0, 0.5, 1.0 |
| Harmony Search | `memory_size`, `hmcr`, `par`, `fretwidth` | 20, 0.9, 0.3, 0.05 |
| SSA | `discoverer_ratio`, `warning_threshold` | 0.2, 0.8 |
| FPA | `switch_probability` | 0.8 |
| DA | `separation_weight`, `alignment_weight`, `cohesion_weight` | 1.0, 1.0, 1.0 |
| CSA | `flight_length`, `awareness_prob` | 2.0, 0.1 |
| MBO | `migration_period`, `adjustment_rate` | 5, 1.0 |
| IWO | `initial_pop`, `max_pop`, `min_seeds`, `max_seeds`, `initial_sigma`, `final_sigma` | 10, 40, 1, 5, 2.0, 0.1 |
| BFO | `chemotactic_steps`, `reproduction_steps`, `step_size` | 5, 3, 0.2 |
| SFLA | `memeplex_count`, `iterations_per_group` | 5, 10 |
| GSO | `luciferin_decay`, `luciferin_enhancement`, `neighbor_count`, `sensing_radius`, `step_size` | 0.4, 0.6, 5, 5.0, 1.0 |
| AMO | `migration_interval` | 5 |
| KH | `max_speed`, `foraging_speed` | 0.01, 0.02 |
| EHO | `clan_count`, `alpha` | 5, 0.5 |
| AFSA | `visual`, `step`, `crowd_factor` | 2.0, 0.5, 0.5 |
| CSO | `mixture_ratio` | 0.3 |
| WCA | `num_rivers` | 4 |
| MFO | `spiral_constant` | 1.0 |
| WOA | `spiral_constant` | 1.0 |
| 其余算法 (GWO, SCA, ALO, HHO, SalpSA, AEO, SOS, GBPSO, MSA, PFA, CPA, BCO, FFOA, FWA, PIO, HBA, POA, TSA, GJO, GOA, SOA, GTO, SMA, AO, ChOA, DO, ZOA) | 仅需 `population_size` | 25 |

### 📐 确定性直接搜索
无需梯度信息，通过几何模式进行直接搜索，适合不可导函数。

| 算法 | 简称 | 说明 |
|------|------|------|
| Nelder-Mead 单纯形法 | NM | 无导数的单纯形几何变形方法 |
| Hooke-Jeeves 模式搜索 | HJ | 无导数的直接搜索方法 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| Nelder-Mead | `reflection_coeff` (α), `expansion_coeff` (γ), `contraction_coeff` (ρ), `shrink_coeff` (σ) | 1.0, 2.0, 0.5, 0.5 |
| Hooke-Jeeves | `initial_step`, `step_reduction`, `min_step` | 0.5, 0.5, 1e-6 |

### 🧮 进化策略
基于概率分布的参数自适应优化，适合复杂病态问题。

| 算法 | 简称 | 说明 |
|------|------|------|
| 协方差矩阵自适应进化策略 | CMA-ES | 协方差矩阵自适应调整，适合病态问题 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| CMA-ES | `sigma` (初始标准差), `lambda` (种群大小), `mu` (父代数量) | 0.5, 20, 10 |

### 📊 估计分布算法
基于概率分布估计的迭代优化方法。

| 算法 | 简称 | 说明 |
|------|------|------|
| 交叉熵方法 | CEM | 基于高斯分布的迭代估计 |
| 随机分形搜索 | SFS | 基于分形发现和学习两阶段搜索 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| CEM | `sample_size`, `elite_fraction`, `smoothing` | 100, 0.1, 0.7 |
| SFS | `population_size`, `discovery_iterations`, `learning_iterations` | 25, 3, 2 |

### ⚛️ 物理启发算法
基于物理定律和自然现象的搜索机制。

| 算法 | 简称 | 说明 |
|------|------|------|
| 引力搜索算法 | GSA | 基于牛顿万有引力定律的搜索 |
| 多宇宙优化 | MVO | 基于白洞/黑洞/虫洞的宇宙机制 |
| 均衡优化 | EO | 基于控制理论均衡池与指数衰减策略 |
| 风驱动优化 | WDO | 模拟大气风速更新与科里奥利力 |
| 水循环算法 | WCA | 基于水循环和河流流动过程 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| GSA | `gravity_constant` (G₀), `decay_alpha` (α) | 100.0, 0.99 |
| MVO | `wormhole_exploitation`, `travelling_distance_rate` | 1.0, 0.6 |
| EO | 仅需 `population_size` | 25 |
| WDO | 仅需 `population_size` | 25 |
| WCA | `population_size`, `num_rivers` | 25, 4 |

### 🎓 教学模型算法
模拟教学与学习过程的优化方法。

| 算法 | 简称 | 说明 |
|------|------|------|
| 教与学优化 | TLBO | 模拟班级教学和学习过程，无需额外参数 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| TLBO | 仅需 `population_size` | 30 |

### 🔢 数学启发算法
基于数学原理和运算规则的搜索机制。

| 算法 | 简称 | 说明 |
|------|------|------|
| Jaya算法 | Jaya | 基于最优-最差个体更新策略 |
| 算术优化算法 | AOA | 基于算术运算符(除/乘/减/加)的搜索 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| Jaya | 仅需 `population_size` | 25 |
| AOA | 仅需 `population_size` | 25 |

### 👥 人类行为算法
模拟人类社会活动和思维过程的优化方法。

| 算法 | 简称 | 说明 |
|------|------|------|
| 头脑风暴优化 | BSO | 模拟人类头脑风暴聚类的创意生成过程 |
| 足球联赛竞争算法 | SLC | 模拟足球联赛团队竞争与球员交换 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| BSO | 仅需 `population_size` | 25 |
| SLC | 仅需 `population_size` | 25 |

## 项目结构

```
BlackBoxOptimization/
├── lib/                          # 核心库包
│   ├── moon.pkg.json             # 包配置
│   └── lib.mbt                   # 所有优化算法实现
├── main/                         # 主程序入口
│   ├── moon.pkg.json             # 包配置
│   └── main.mbt                  # 基准测试入口
├── test/                         # 测试包
│   ├── moon.pkg.json             # 包配置
│   └── lib_test.mbt              # 单元测试（140 个测试用例）
├── example/                      # 示例包
│   ├── moon.pkg.json             # 包配置
│   └── simple_usage.mbt          # 使用示例（61 个示例）
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
# 运行所有测试（140 个测试用例）
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

// 麻雀搜索算法 (SSA)
let ssa_config = @lib.make_ssa_config(1000, dim, bounds, 25)
let ssa_result = @lib.sparrow_search_algorithm(ssa_config, f, rng)

// 哈里斯鹰优化 (HHO)
let hho_config = @lib.make_hho_config(1000, dim, bounds, 25)
let hho_result = @lib.harris_hawks_optimization(hho_config, f, rng)

// 花粉传播算法 (FPA)
let fpa_config = @lib.make_fpa_config(1000, dim, bounds, 25)
let fpa_result = @lib.flower_pollination_algorithm(fpa_config, f, rng)

// 蜻蜓算法 (DA)
let da_config = @lib.make_da_config(1000, dim, bounds, 25)
let da_result = @lib.dragonfly_algorithm(da_config, f, rng)

// 乌鸦搜索算法 (CSA)
let csa_config = @lib.make_csa_config(1000, dim, bounds, 25)
let csa_result = @lib.crow_search_algorithm(csa_config, f, rng)

// 樽海鞘群算法 (SalpSA)
let salp_config = @lib.make_salp_sa_config(1000, dim, bounds, 25)
let salp_result = @lib.salp_swarm_algorithm(salp_config, f, rng)

// Monarch蝴蝶优化 (MBO)
let mbo_config = @lib.make_mbo_config(1000, dim, bounds, 25)
let mbo_result = @lib.monarch_butterfly_optimization(mbo_config, f, rng)

// 人工生态系统优化 (AEO)
let aeo_config = @lib.make_aeo_config(1000, dim, bounds, 25)
let aeo_result = @lib.artificial_ecosystem_optimization(aeo_config, f, rng)

// 入侵杂草优化 (IWO)
let iwo_config = @lib.make_iwo_config(1000, dim, bounds, 10, 40, 1, 5, 2.0, 0.1, 3.0)
let iwo_result = @lib.invasive_weed_optimization(iwo_config, f, rng)

// 细菌觅食优化 (BFO)
let bfo_config = @lib.make_bfo_config(1000, dim, bounds, 25, 5, 3, 0.2)
let bfo_result = @lib.bacteria_foraging_optimization(bfo_config, f, rng)

// 混合蛙跳算法 (SFLA)
let sfla_config = @lib.make_sfla_config(1000, dim, bounds, 30, 5, 10)
let sfla_result = @lib.shuffled_frog_leaping_algorithm(sfla_config, f, rng)

// 萤火虫群优化 (GSO)
let gso_config = @lib.make_gso_config(1000, dim, bounds, 25, 0.4, 0.6, 5, 5.0, 1.0)
let gso_result = @lib.glowworm_swarm_optimization(gso_config, f, rng)

// 引力搜索算法 (GSA)
let gsa_config = @lib.make_gsa_config(1000, dim, bounds, 25, 100.0, 0.99)
let gsa_result = @lib.gravitational_search_algorithm(gsa_config, f, rng)

// 教与学优化 (TLBO)
let tlbo_config = @lib.make_tlbo_config(1000, dim, bounds, 30)
let tlbo_result = @lib.teaching_learning_based_optimization(tlbo_config, f, rng)

// 动物迁徙优化 (AMO)
let amo_config = @lib.make_amo_config(1000, dim, bounds, 25, 5)
let amo_result = @lib.animal_migration_optimization(amo_config, f, rng)

// 多宇宙优化 (MVO)
let mvo_config = @lib.make_mvo_config(1000, dim, bounds, 25)
let mvo_result = @lib.multi_verse_optimizer(mvo_config, f, rng)

// 磷虾群算法 (KH)
let kh_config = @lib.make_kh_config(1000, dim, bounds, 25)
let kh_result = @lib.krill_herd(kh_config, f, rng)

// 帝国主义竞争算法 (ICA)
let ica_config = @lib.make_ica_config(1000, dim, bounds, 25)
let ica_result = @lib.imperialist_competitive_algorithm(ica_config, f, rng)

// 随机分形搜索 (SFS)
let sfs_config = @lib.make_sfs_config(1000, dim, bounds, 25)
let sfs_result = @lib.stochastic_fractal_search(sfs_config, f, rng)

// 共生生物搜索 (SOS)
let sos_config = @lib.make_sos_config(1000, dim, bounds, 25)
let sos_result = @lib.symbiotic_organisms_search(sos_config, f, rng)

// 生物地理学优化 (BBO)
let bbo_config = @lib.make_bbo_config(1000, dim, bounds, 25, 0.01)
let bbo_result = @lib.biogeography_based_optimization(bbo_config, f, rng)

// 大象牧群优化 (EHO)
let eho_config = @lib.make_eho_config(1000, dim, bounds, 25, 5, 0.5)
let eho_result = @lib.elephant_herding_optimization(eho_config, f, rng)

// 均衡优化 (EO)
let eo_config = @lib.make_eo_config(1000, dim, bounds, 25)
let eo_result = @lib.equilibrium_optimizer(eo_config, f, rng)

// 风驱动优化 (WDO)
let wdo_config = @lib.make_wdo_config(1000, dim, bounds, 25)
let wdo_result = @lib.wind_driven_optimization(wdo_config, f, rng)

// Jaya算法
let jaya_config = @lib.make_jaya_config(1000, dim, bounds, 25)
let jaya_result = @lib.jaya_optimization(jaya_config, f, rng)

// 足球联赛竞争算法 (SLC)
let slc_config = @lib.make_slc_config(1000, dim, bounds, 25)
let slc_result = @lib.soccer_league_competition(slc_config, f, rng)

// 高斯骨干粒子群优化 (GBPSO)
let gbpso_config = @lib.make_gbpso_config(1000, dim, bounds, 25)
let gbpso_result = @lib.gaussian_bare_bones_pso(gbpso_config, f, rng)

// 飞蛾搜索算法 (MSA)
let msa_config = @lib.make_msa_config(1000, dim, bounds, 25)
let msa_result = @lib.moth_search_algorithm(msa_config, f, rng)

// 路径finder算法 (PFA)
let pfa_config = @lib.make_pfa_config(1000, dim, bounds, 25)
let pfa_result = @lib.pathfinder_algorithm(pfa_config, f, rng)

// 菌落捕食算法 (CPA)
let cpa_config = @lib.make_cpa_config(1000, dim, bounds, 25)
let cpa_result = @lib.colony_predation_algorithm(cpa_config, f, rng)

// 细菌菌落优化 (BCO)
let bco_config = @lib.make_bco_config(1000, dim, bounds, 25)
let bco_result = @lib.bacteria_colony_optimization(bco_config, f, rng)

// 飞狐优化算法 (FFOA)
let ffoa_config = @lib.make_ffoa_config(1000, dim, bounds, 25)
let ffoa_result = @lib.flying_fox_optimization(ffoa_config, f, rng)

// 人工鱼群算法 (AFSA)
let afsa_config = @lib.make_afsa_config(1000, dim, bounds, 25, 2.0, 0.5, 0.5)
let afsa_result = @lib.artificial_fish_swarm(afsa_config, f, rng)

// 烟花算法 (FWA)
let fwa_config = @lib.make_fwa_config(1000, dim, bounds, 25)
let fwa_result = @lib.fireworks_algorithm(fwa_config, f, rng)

// 头脑风暴优化 (BSO)
let bso_config = @lib.make_bso_config(1000, dim, bounds, 25)
let bso_result = @lib.brain_storm_optimization(bso_config, f, rng)

// 回溯搜索算法 (BSA)
let bsa_config = @lib.make_bsa_config(1000, dim, bounds, 25)
let bsa_result = @lib.backtracking_search_algorithm(bsa_config, f, rng)

// 鸽群优化 (PIO)
let pio_config = @lib.make_pio_config(1000, dim, bounds, 25)
let pio_result = @lib.pigeon_inspired_optimization(pio_config, f, rng)

// 算术优化算法 (AOA)
let aoa_config = @lib.make_aoa_config(1000, dim, bounds, 25)
let aoa_result = @lib.arithmetic_optimization_algorithm(aoa_config, f, rng)

// 蜜獾算法 (HBA)
let hba_config = @lib.make_hba_config(1000, dim, bounds, 25)
let hba_result = @lib.honey_badger_algorithm(hba_config, f, rng)

// 鹈鹕优化算法 (POA)
let poa_config = @lib.make_poa_config(1000, dim, bounds, 25)
let poa_result = @lib.pelican_optimization_algorithm(poa_config, f, rng)

// 金枪鱼群优化 (TSA)
let tsa_config = @lib.make_tsa_config(1000, dim, bounds, 25)
let tsa_result = @lib.tuna_swarm_optimization(tsa_config, f, rng)

// 金豺优化 (GJO)
let gjo_config = @lib.make_gjo_config(1000, dim, bounds, 25)
let gjo_result = @lib.golden_jackal_optimization(gjo_config, f, rng)

// 蝗虫优化算法 (GOA)
let goa_config = @lib.make_goa_config(1000, dim, bounds, 25)
let goa_result = @lib.grasshopper_optimization_algorithm(goa_config, f, rng)

// 海鸥优化算法 (SOA)
let soa_config = @lib.make_soa_config(1000, dim, bounds, 25)
let soa_result = @lib.seagull_optimization_algorithm(soa_config, f, rng)

// 猫群优化 (CSO)
let cso_config = @lib.make_cso_config(1000, dim, bounds, 25)
let cso_result = @lib.cat_swarm_optimization(cso_config, f, rng)

// 大猩猩部队优化 (GTO)
let gto_config = @lib.make_gto_config(1000, dim, bounds, 25)
let gto_result = @lib.gorilla_troops_optimizer(gto_config, f, rng)

// 水循环算法 (WCA)
let wca_config = @lib.make_wca_config(1000, dim, bounds, 25, 4)
let wca_result = @lib.water_cycle_algorithm(wca_config, f, rng)
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

## 性能对比

### Sphere 函数（5维）

| 算法 | 1000 evals | 5000 evals | 10000 evals |
|------|------------|------------|-------------|
| 随机搜索 | 3.17 | 0.92 | 0.92 |
| 爬山法 | 0.0094 | 0.0094 | 0.0094 |
| 模拟退火 | 0.57 | 0.0059 | 0.0059 |
| **PSO** | **8.56e-4** | **2.28e-15** | **3.38e-29** |
| **DE** | 0.156 | **1.02e-8** | **8.65e-18** |
| GA | 0.698 | 0.0559 | 0.0111 |
| Nelder-Mead | 49.38 | 49.38 | 49.38 |
| 萤火虫 | 4.25 | 0.92 | 0.92 |
| 蚁群算法 | 2.43 | 2.43 | 2.12 |
| 蜂群算法 | 0.273 | **3.33e-4** | **6.23e-10** |
| 禁忌搜索 | 0.0172 | 0.0151 | 0.0151 |
| 引导局部搜索 | 0.0172 | 0.0172 | 0.0172 |
| CMA-ES | 0.0236 | 3.56e-4 | 9.22e-7 |
| CEM | 0.310 | 3.52e-7 | 1.04e-7 |
| 布谷鸟搜索 | 0.744 | 0.0636 | 0.0636 |
| GWO | 0.00163 | 1.00e-9 | 6.61e-16 |
| **Hooke-Jeeves** | **1.71e-12** | **1.71e-12** | **1.71e-12** |
| WOA | 0.0020 | 9.70e-10 | 9.03e-16 |
| 蝙蝠算法 | 0.00625 | 6.49e-5 | 6.49e-5 |
| 和声搜索 | 0.00284 | 3.39e-5 | 2.31e-5 |
| **SCA** | **1.14e-14** | **9.11e-47** | **1.60e-99** |
| MFO | 0.140 | 1.63e-8 | 2.02e-17 |
| ALO | 1.45 | 1.41 | 0.385 |
| SSA | 1.24e-4 | 6.28e-13 | 1.95e-17 |
| **HHO** | **~0** | **0** | **0** |
| FPA | 0.596 | 0.0728 | 0.0542 |
| DA | 1.589 | 2.667 | 2.635 |
| CSA | 0.00365 | 1.70e-7 | 1.38e-7 |
| SalpSA | 0.731 | 0.0591 | 0.00506 |
| MBO | 5.73e-5 | 2.53e-5 | 2.47e-5 |
| AEO | 0.0355 | 0.00394 | 0.0105 |
| IWO | 3.58 | 0.00242 | 0.00106 |
| BFO | 1.68 | 0.00735 | 0.00735 |
| SFLA | 9.86e-4 | 3.48e-5 | 7.36e-9 |
| GSO | 7.61 | 6.39 | 6.39 |
| GSA | 18.11 | 18.11 | 18.11 |
| **TLBO** | **1.41e-4** | **4.26e-24** | **9.02e-47** |
| **AMO** | **0.0171** | **5.02e-7** | **8.41e-12** |
| MVO | 13.47 | 13.47 | 13.47 |
| KH | 0.193 | 0.193 | 0.193 |
| ICA | 0.906 | 0.906 | 0.906 |
| **SFS** | **1.82e-4** | **5.63e-5** | **3.36e-6** |
| **SOS** | **3.17e-4** | **~0** | **~0** |
| BBO | 13.46 | 8.47 | 4.10 |
| EHO | 0.0112 | 7.69e-4 | 2.37e-4 |
| EO | 6.88 | 0.598 | 0.00655 |
| WDO | 3.96 | 10.78 | 0.142 |
| **Jaya** | **0.0133** | **~0** | **~0** |
| **SLC** | **2.87e-5** | **2.70e-5** | **2.70e-5** |
| **GBPSO** | **~0** | **~0** | **~0** |
| MSA | 13.47 | 2.77 | 0.413 |
| PFA | 0.226 | 0.109 | 0.353 |
| CPA | 1.90 | 0.592 | 2.73 |
| **BCO** | **0.245** | **0.00157** | **0.00270** |
| FFOA | 0.701 | 0.759 | 0.346 |
| AFSA | 0.134 | 0.0183 | 0.0183 |
| FWA | 7.28 | 4.21 | 4.21 |
| BSO | 0.00321 | 0.0000290 | 0.0000208 |
| BSA | 1.57 | 0.182 | 0.00170 |
| PIO | 0.00676 | 1.83e-10 | 1.88e-18 |
| AOA | 0.00280 | 5.27e-5 | 2.86e-4 |
| HBA | 0.00949 | 0.0183 | 6.31e-4 |
| POA | 0.00888 | 7.24e-10 | 6.18e-17 |
| TSA | 0 | 0 | 0 |
| GJO | 1.87e-7 | 1.07e-26 | 3.75e-56 |
| GOA | 0.0553 | 0.00414 | 1.65e-4 |
| **SOA** | **4.89e-5** | **~0** | **~0** |
| CSO | 0.0218 | 2.57e-5 | ~0 |
| **GTO** | **0.0573** | **7.99e-6** | **9.28e-6** |
| WCA | 0.0343 | 0.0129 | 0.0117 |
| SMA | 0.0390 | 0.0360 | 0.0360 |
| AO | 0.117 | 0.00331 | 0.000536 |
| ChOA | 1.18e-10 | ~0 | ~0 |
| DO | 0.279 | 0.396 | 0.366 |
| ZOA | 0.0113 | 0.000251 | 0.0000328 |

### Rastrigin 函数（5维）

| 算法 | 5000 evals | 10000 evals |
|------|------------|-------------|
| 随机搜索 | 21.84 | 17.38 |
| 爬山法 | 53.53 | 52.75 |
| 模拟退火 | 23.41 | 23.41 |
| PSO | 2.99 | 2.98 |
| DE | 8.29 | 1.35 |
| GA | 9.76 | 1.66 |
| Nelder-Mead | 61.39 | 61.39 |
| 萤火虫 | 22.99 | 17.78 |
| 蚁群算法 | 9.49 | 9.49 |
| 蜂群算法 | 3.92 | 1.11 |
| 禁忌搜索 | 6.91 | 6.91 |
| 引导局部搜索 | 23.31 | 15.64 |
| CMA-ES | 22.99 | 22.88 |
| CEM | 0.593 | 2.30e-5 |
| 布谷鸟搜索 | 11.74 | 11.74 |
| GWO | 5.02 | 0.0025 |
| Hooke-Jeeves | 45.77 | 45.77 |
| **WOA** | **0.995** | **1.27e-11** |
| 蝙蝠算法 | 13.37 | 13.20 |
| 和声搜索 | 2.03 | 0.607 |
| **SCA** | **0** | **0** |
| MFO | 3.99 | 1.99 |
| ALO | 11.45 | 11.45 |
| SSA | 8.95 | 8.95 |
| **HHO** | **0** | **0** |
| FPA | 9.43 | 5.09 |
| DA | 49.12 | 49.75 |
| CSA | 9.06e-5 | 8.92e-5 |
| SalpSA | 25.56 | 25.58 |
| MBO | 0.00525 | 0.00515 |
| AEO | 2.82 | 5.38 |
| IWO | 2.69 | 1.23 |
| BFO | 3.02 | 3.02 |
| SFLA | 1.99 | 1.99 |
| GSO | 22.23 | 22.23 |
| GSA | 65.32 | 65.32 |
| **TLBO** | **0.0333** | **0** |
| AMO | 6.97 | 6.96 |
| MVO | 48.86 | 48.86 |
| KH | 6.30 | 7.86 |
| ICA | 17.97 | 17.97 |
| **SFS** | **0.454** | **0.0113** |
| SOS | 2.98 | 2.98 |
| BBO | 37.32 | 37.32 |
| EHO | 5.11 | 5.04 |
| EO | 50.30 | 45.77 |
| WDO | 125.66 | 101.08 |
| Jaya | 7.82 | 6.71 |
| **SLC** | **2.01** | **2.01** |
| GBPSO | 9.95 | 9.95 |
| MSA | 31.65 | 35.67 |
| PFA | 21.95 | 24.16 |
| CPA | 25.99 | 23.13 |
| **BCO** | **2.29** | **1.58** |
| FFOA | 24.17 | 25.06 |
| AFSA | 14.88 | 7.53 |
| FWA | 27.79 | 27.79 |
| BSO | 2.998 | 1.993 |
| BSA | 29.99 | 22.36 |
| PIO | 5.97 | 5.97 |
| AOA | 0.552 | 0.0680 |
| HBA | 7.46 | 6.36 |
| POA | 1.12 | 0.00254 |
| TSA | 0 | 0 |
| GJO | 1.99 | 14.92 |
| GOA | 4.25 | 7.04 |
| SOA | 2.98 | 1.99 |
| CSO | 6.11 | 5.97 |
| **GTO** | **1.30** | **1.17** |
| WCA | 3.19 | 4.32 |
| SMA | 17.04 | 17.04 |
| AO | 8.48 | 4.82 |
| ChOA | 9.65 | 9.65 |
| DO | 14.99 | 13.64 |
| ZOA | 8.42 | 5.07 |

### Rosenbrock 函数（5维）

| 算法 | 5000 evals | 10000 evals |
|------|------------|-------------|
| 随机搜索 | 108.84 | 22.78 |
| 爬山法 | 7.06 | 5.63 |
| 模拟退火 | 4.98 | 2.74 |
| PSO | 4.35 | 4.09 |
| **DE** | **0.0676** | **3.98e-8** |
| GA | 28.92 | 17.26 |
| Nelder-Mead | 43504.88 | 43504.88 |
| 萤火虫 | 91.21 | 91.21 |
| 蚁群算法 | 14.26 | 7.47 |
| 蜂群算法 | 1.36 | 0.265 |
| 禁忌搜索 | 1.95 | 1.95 |
| 引导局部搜索 | 3.67 | 3.67 |
| CMA-ES | 3.22 | 1.25 |
| CEM | 2.38 | 2.35 |
| 布谷鸟搜索 | 17.80 | 15.13 |
| GWO | 5.90 | 0.0392 |
| **Hooke-Jeeves** | **0.440** | **0.288** |
| **WOA** | **2.10** | **0.00221** |
| 蝙蝠算法 | 130.18 | 130.18 |
| 和声搜索 | 0.699 | 0.679 |
| SCA | 0.0193 | 0.00461 |
| MFO | 4.34 | 2.49 |
| ALO | 34.06 | 34.06 |
| SSA | 2.09 | 1.31 |
| **HHO** | **0.00227** | **0.00227** |
| FPA | 3.06 | 0.346 |
| DA | 1025.52 | 976.07 |
| **CSA** | **1.06e-4** | **5.25e-5** |
| SalpSA | 9.62 | 2.27 |
| MBO | 0.0239 | 0.00623 |
| AEO | 5.62 | 2.95 |
| IWO | 1.42 | 1.48 |
| BFO | 2.41 | 1.87 |
| SFLA | 17.78 | 4.69 |
| GSO | 1856.04 | 451.42 |
| GSA | 4343.62 | 4343.62 |
| TLBO | 1.72 | 1.57 |
| AMO | 3.91 | 3.60 |
| MVO | 4875.10 | 4875.10 |
| KH | 9.59 | 13.22 |
| ICA | 9.57 | 32.58 |
| **SFS** | **0.0123** | **0.0128** |
| SOS | 4.51 | 4.16 |
| BBO | 2373.73 | 636.06 |
| EHO | 1.54 | 0.394 |
| EO | 4.90 | 0.160 |
| WDO | 55358.40 | 559.77 |
| **Jaya** | **0.706** | **0.0925** |
| SLC | 62.59 | 62.59 |
| **GBPSO** | **0.00110** | **0.000123** |
| MSA | 544.40 | 173.63 |
| PFA | 7.19 | 11.27 |
| CPA | 16.30 | 82.19 |
| **BCO** | **3.22** | **1.89** |
| FFOA | 86.24 | 316.22 |
| AFSA | 7.87 | 6.79 |
| FWA | 522.59 | 121.52 |
| BSO | 2.39 | 0.407 |
| BSA | 31.70 | 6.59 |
| PIO | 3.47 | 1.72 |
| AOA | 0.0424 | 0.0327 |
| HBA | 4.38 | 0.263 |
| POA | 0.999 | 0.966 |
| TSA | 3.77 | 3.77 |
| GJO | 0.898 | 0.775 |
| GOA | 0.922 | 3.65 |
| SOA | 1.49 | 0.0162 |
| CSO | 1.58 | 0.618 |
| GTO | 3.00 | 0.0833 |
| WCA | 3.55 | 2.74 |
| SMA | 38.67 | 38.24 |
| AO | 4.35 | 3.44 |
| ChOA | 2.07 | 1.99 |
| DO | 128.57 | 127.91 |
| ZOA | 1.37 | 0.821 |

## 算法分类特点总结

| 类别 | 代表算法 | 特点 |
|------|----------|------|
| **🔍 单点搜索** | 随机搜索、爬山法、模拟退火 | 简单、易于实现，适合低维问题 |
| **🔧 局部搜索** | 禁忌搜索、引导局部搜索 | 强化局部开发，避免陷入局部最优 |
| **🧬 进化算法** | DE、GA、ICA、BBO、BSA | 基于种群进化，鲁棒性好 |
| **🐝 群体智能** | PSO、萤火虫、蚁群、蜂群、灰狼、鲸群、蝙蝠、和声、麻雀搜索、哈里斯鹰、GOA、SOA、CSO、GTO、SMA、AO、ChOA、DO、ZOA 等 | 基于群体行为，全局搜索能力强 |
| **📐 确定性直接搜索** | Nelder-Mead、Hooke-Jeeves | 无需梯度，适合不可导函数 |
| **🧮 进化策略** | CMA-ES | 协方差自适应，适合病态问题 |
| **📊 估计分布** | CEM、SFS | 基于概率分布估计，适合连续优化 |
| **⚛️ 物理启发** | GSA、MVO、EO、WDO、WCA | 基于物理规律，自适应搜索 |
| **🎓 教学模型** | TLBO | 无需算法特定参数，收敛速度快 |
| **🔢 数学启发** | Jaya、AOA | 基于数学原理的搜索机制 |
| **👥 人类行为** | BSO、SLC | 模拟人类社会行为过程 |

## 许可证

本项目基于 [Apache-2.0 license](LICENSE) 开源。