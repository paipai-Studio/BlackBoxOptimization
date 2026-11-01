# BlackBoxOptimization

基于 MoonBit 实现的黑盒优化算法框架，支持 92 种经典优化算法，提供统一的接口和基准测试框架。

## 算法分类

本框架实现了 **92 种**黑盒优化算法，按设计灵感和机制分为以下类别：

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
| 自适应差分进化 | jDE | 每个个体自学习 F/CR，成功才保留 |
| 历史记忆自适应差分进化 | SHADE | 成功参数记忆 + 外部归档 + p-best |
| 进化规划 | EP | 对数正态自适应步长 + (μ+μ) 选择 |
| 遗传算法 | GA | 选择、交叉、变异三阶段模拟自然选择 |
| 帝国主义竞争算法 | ICA | 模拟帝国殖民、同化和竞争机制 |
| 生物地理学优化 | BBO | 基于栖息地适宜度与物种迁移模型 |
| 回溯搜索算法 | BSA | 基于历史种群记忆的定向交叉变异 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| DE | `population_size`, `mutation_factor` (F), `crossover_probability` (CR) | 30, 0.8, 0.9 |
| jDE | 仅需 `population_size` | 25 |
| SHADE | 仅需 `population_size` | 25 |
| EP | 仅需 `population_size` | 25 |
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
| 蝠鲼觅食优化 | MRFO | 链式/螺旋/翻滚三种觅食机制 (2020) |
| 饥饿游戏搜索 | HGS | 饥饿权重驱动的自适应搜索 (2021) |
| 爬行动物搜索算法 | RSA | 包围/狩猎/协调多阶段搜索 (2022) |
| 海洋捕食者算法 | MPA | 布朗/Levy 三阶段速度比捕食机制 (2020) |
| 金鹰优化器 | GEO | 攻击-巡航双倾向的正交向量搜索 (2021) |
| 非洲秃鹫优化 | AVOA | 饥饿因子驱动的探索-开发切换 (2021) |
| 白鲨优化器 | WSO | 速度惯性+猎物跟随的白鲨觅食 (2022) |
| 长鼻浣熊优化 | COA | 捕猎鬣蜥与逃避捕食者两阶段 (2023) |
| 鱼鹰优化算法 | OOA | 俯冲捕鱼 + 携带到合适位置 (2023) |
| 塔斯马尼亚恶魔优化 | TDO | 食腐/捕猎双策略 + 拖拽食物 (2022) |
| 草原犬鼠优化 | PDO | 觅食/挖洞/警戒/平静四阶段 (2022) |
| 北方苍鹰优化 | NGO | 俯冲猎物 + 尾部追击双阶段 (2021) |
| 白鲸优化 | BWO | 游动/捕食/鲸落三阶段 (2022) |
| 梯度优化器 | GBO | 梯度搜索规则 + 局部逃逸算子 (2020) |
| 生长优化器 | GO | 学习阶段 + 反思阶段双阶段 (2023) |
| 海象优化器 | WaOA | 觅食/迁徙/躲避天敌三行为 (2023) |
| 克隆选择算法 | CLONALG | 克隆-变异-选择的人工免疫机制 (2002) |
| 人工蜂群算法 | ABC | 雇佣/观察/侦察蜂三阶段觅食 (2007) |
| 电磁机制算法 | EM | 电荷吸引/排斥的库仑类搜索 (2003) |
| 蜣螂优化器 | DBO | 滚球/繁殖/觅食/偷窃四角色 (2022) |
| 被囊群算法 | TSO | 喷流推进 + 重力/水流被囊群搜索 (2020) |
| 笑翠鸟优化算法 | KOA | 栖枝扫描 + 俯冲扑食两相 (2024) |
| 火鹰优化算法 | FHO | 投掷火棒惊猎/捕获两相 (2022) |
| 霜冰优化算法 | RIME | 软霜搜索 + 硬霜穿刺 + 附着因子 (2023) |
| 向量加权平均优化 | INFO | best/better/worst 加权均值驱动 (2022) |
| 蛇优化器 | SO | 温度驱动的雄雌蛇觅食/争斗/探索 (2022) |

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
| 其余算法 (GWO, SCA, ALO, HHO, SalpSA, AEO, SOS, GBPSO, MSA, PFA, CPA, BCO, FFOA, FWA, PIO, HBA, POA, TSA, GJO, GOA, SOA, GTO, SMA, AO, ChOA, DO, ZOA, MRFO, HGS, RSA, MPA, GEO, AVOA, WSO, COA, OOA, TDO, PDO, NGO, BWO, GBO, GO, WaOA, CLONALG, ABC, EM, DBO, TSO, KOA, FHO, RIME, INFO, SO) | 仅需 `population_size` | 25 |

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
| (μ+λ) 进化策略 | ES | μ 父代产生 λ 高斯变异子代，择优存活 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| CMA-ES | `sigma` (初始标准差), `lambda` (种群大小), `mu` (父代数量) | 0.5, 20, 10 |
| ES | `mu` (父代数量), `lambda` (子代数量), `sigma_init` (初始步长) | 30, 100, 0.3 |

### 📊 估计分布算法
基于概率分布估计的迭代优化方法。

| 算法 | 简称 | 说明 |
|------|------|------|
| 交叉熵方法 | CEM | 基于高斯分布的迭代估计 |
| 随机分形搜索 | SFS | 基于分形发现和学习两阶段搜索 |
| 基于群体的增量学习 | PBIL | 连续分布估计与精英统计更新 |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| CEM | `sample_size`, `elite_fraction`, `smoothing` | 100, 0.1, 0.7 |
| SFS | `population_size`, `discovery_iterations`, `learning_iterations` | 25, 3, 2 |
| PBIL | `sample_size`, `elite_size`, `learning_rate` | 100, 10, 0.3 |

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
| Runge Kutta 优化器 | RUN | 基于 RK4 斜率加权组合的位置更新 (2021) |

| 算法 | 关键参数 | 推荐值 |
|------|----------|--------|
| Jaya | 仅需 `population_size` | 25 |
| AOA | 仅需 `population_size` | 25 |
| RUN | 仅需 `population_size` | 25 |

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
│   ├── moon.pkg                  # 包配置
│   ├── core.mbt                  # RNG / 数学工具 / 核心类型
│   ├── classic_algorithms.mbt    # 随机搜索/爬山/SA/PSO/DE/GA/NM
│   ├── swarm_algorithms_1.mbt    # 群体智能算法（Firefly..SSA）
│   ├── swarm_algorithms_2.mbt    # 群体智能算法（HHO..KH）
│   ├── meta_algorithms.mbt       # 物理/人文启发算法（ICA..GBPSO）
│   ├── recent_algorithms.mbt     # 新一代算法（MSA..WCA）
│   └── benchmark.mbt             # 测试函数 + 基准框架 + 算法注册表
├── main/                         # 主程序入口
│   ├── moon.pkg                  # 包配置
│   └── main.mbt                  # 基准测试入口
├── test/                         # 测试包
│   ├── moon.pkg                  # 包配置
│   ├── lib_test.mbt              # 单元测试
│   └── extra_test.mbt            # 补充测试（收敛性 / RNG / 注册表）
├── example/                      # 示例包
│   ├── moon.pkg                  # 包配置
│   └── simple_usage.mbt          # 使用示例（注册表驱动，覆盖全部算法）
├── .github/workflows/            # CI 配置
│   └── ci.yml                    # GitHub Actions
├── moon.mod                      # 模块配置
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
# 运行所有测试（172 个测试用例）
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
// Sphere 函数: f(x) = sum(x_i^2)，[-5,5]，最优值: 0
let f1 = @lib.sphere_function

// Rastrigin 函数: f(x) = 10n + sum(x_i^2 - 10*cos(2*pi*x_i))，[-5.12,5.12]，最优值: 0
let f2 = @lib.rastrigin_function

// Rosenbrock 函数: f(x) = sum(100*(x_{i+1} - x_i^2)^2 + (1 - x_i)^2)，[-5,5]，最优值: 0
let f3 = @lib.rosenbrock_function

// Ackley 函数，[-32.768,32.768]，最优值: 0
let f4 = @lib.ackley_function

// Griewank 函数，[-600,600]，最优值: 0
let f5 = @lib.griewank_function

// Schwefel 函数: 418.9829*n - sum(x_i*sin(sqrt(|x_i|)))，[-500,500]，最优值: 0
let f6 = @lib.schwefel_function

// Levy 函数，[-10,10]，最优值: 0
let f7 = @lib.levy_function

// Michalewicz 函数 (m=10)，[0,pi]，5维最优值: -4.687658
let f8 = @lib.michalewicz_function

// Zakharov 函数，[-10,10]，最优值: 0
let f9 = @lib.zakharov_function

// Booth 函数 (2维: 仅使用 x_0, x_1)，[-10,10]，最优值: 0（在 (1,3)）
let f10 = @lib.booth_function

// Matyas 函数 (2维: 仅使用 x_0, x_1)，[-10,10]，最优值: 0
let f11 = @lib.matyas_function
```

### 基准测试

```moonbit
// 运行单轮基准测试
let results = @lib.run_benchmark_at_evals(f, "Sphere", dim, bounds, 0.0, 5000, rng)
@lib.print_benchmark_table(results, "Sphere", 5000)

// 运行多次试验统计分析
@lib.run_multiple_trials(f, "Sphere", dim, bounds, 0.0, 10000, 20, 100.0)
```

### 基准函数性能对比（5维，20次试验平均值，10000 次评估）

| 算法 | Sphere | Rastrigin | Rosenbrock |
|------|--------|-----------|------------|
| Random Search | 1.23097 | 17.3331 | 81.8003 |
| Hill Climbing | 0.0117986 | 43.5007 | 3.24077 |
| Simulated Annealing | 0.0119307 | 21.1379 | 2.36552 |
| PSO | 8.570e-28 | 0.895466 | 0.481973 |
| Differential Evolution | 3.893e-17 | 2.43324 | 5.023e-07 |
| Genetic Algorithm | 0.0258742 | 1.81865 | 10.8283 |
| Nelder-Mead | 0 | 36.6178 | 160.92 |
| Firefly | 0.00151743 | 2.49215 | 3.34032 |
| ACO | 8.994e-08 | 9.30287 | 1.03295 |
| Bee Colony | 3.874e-23 | 2.92059 | 1.32541 |
| Tabu Search | 0.0105786 | 10.0375 | 2.61421 |
| Local Search | 0.0515364 | 10.052 | 2.5924 |
| CMA-ES | 8.402e-17 | 21.5905 | 2.22265 |
| CEM | 8.536e-08 | 0.199023 | 2.61812 |
| Cuckoo Search | 0.155951 | 13.2493 | 18.4164 |
| GWO | 2.377e-08 | 4.97527 | 2.81724 |
| Hooke-Jeeves | 1.479e-12 | 41.191 | 0.687194 |
| WOA | 1.930e-13 | 0.0994959 | 0.824381 |
| Bat Algorithm | 0.0303282 | 13.5217 | 72.9674 |
| Harmony Search | 1.803e-05 | 0.0611436 | 1.13162 |
| SCA | 3.570e-87 | 0 | 0.00265951 |
| MFO | 1.906e-16 | 2.34164 | 2.55408 |
| ALO | 0.336594 | 12.9665 | 33.0866 |
| SSA | 4.865e-17 | 5.92334 | 2.65036 |
| HHO | 0 | 0 | 0.00108404 |
| FPA | 0.0217742 | 4.11448 | 2.58233 |
| DA | 0.641755 | 26.799 | 347.421 |
| CSA | 4.162e-07 | 0.000298703 | 5.380e-05 |
| SalpSA | 0.00220677 | 9.37451 | 0.140312 |
| MBO | 1.410e-05 | 0.00269612 | 0.00125645 |
| AEO | 0.00220394 | 4.12605 | 2.75217 |
| IWO | 0.00155442 | 1.69164 | 1.38689 |
| BFO | 0.00434427 | 7.02834 | 1.4009 |
| SFLA | 5.733e-05 | 2.17137 | 2.7386 |
| GSO | 0.000127826 | 11.3255 | 2.91445 |
| GSA | 1.01661 | 25.6955 | 81.6524 |
| TLBO | 5.846e-46 | 0.222673 | 1.47914 |
| AMO | 1.913e-11 | 5.23707 | 2.82492 |
| MVO | 0.461486 | 6.89621 | 64.359 |
| KH | 0.282429 | 10.6311 | 17.2275 |
| ICA | 0.224188 | 8.8171 | 27.1254 |
| SFS | 3.677e-05 | 0.0117503 | 1.19748 |
| SOS | 2.349e-34 | 2.63664 | 0.519568 |
| BBO | 0.0748988 | 3.77072 | 25.4932 |
| EHO | 0.000152811 | 5.40249 | 2.15261 |
| EO | 0.00310971 | 20.4471 | 1.54658 |
| WDO | 1.81616 | 102.334 | 8693.35 |
| Jaya | 3.257e-29 | 6.19249 | 0.103739 |
| SLC | 9.336e-05 | 9.18207 | 5.74533 |
| GBPSO | 2.337e-105 | 13.0339 | 0.786705 |
| MSA | 1.54279 | 32.5245 | 215.551 |
| PFA | 3.229e-06 | 5.52268 | 0.281237 |
| CPA | 0.000360239 | 2.31977 | 2.533 |
| BCO | 0.00281821 | 3.95452 | 3.11246 |
| FFOA | 7.182e-05 | 6.28177 | 1.24135 |
| AFSA | 0.0607567 | 10.8102 | 6.47047 |
| FWA | 0.000474598 | 6.0104 | 3.66727 |
| BSO | 2.218e-05 | 2.989 | 0.513668 |
| BSA | 0.00273547 | 13.3233 | 7.72253 |
| PIO | 1.047e-09 | 10.646 | 2.44087 |
| AOA | 0.000188293 | 0.071434 | 0.0523264 |
| HBA | 0.00205751 | 3.21057 | 0.236659 |
| POA | 9.903e-14 | 1.17062 | 1.09585 |
| TSA | 0 | 0 | 3.69103 |
| GJO | 1.400e-16 | 5.67975 | 36.9972 |
| GOA | 0.00103025 | 9.41585 | 8.54423 |
| SOA | 8.909e-44 | 1.84074 | 7.82782 |
| CSO | 3.691e-07 | 8.71592 | 0.572322 |
| GTO | 0.000902703 | 1.41789 | 2.83956 |
| WCA | 0.0126508 | 5.88644 | 2.95402 |
| SMA | 0.315557 | 10.4473 | 41.8348 |
| AO | 0.00430357 | 4.47815 | 3.97172 |
| ChOA | 7.760e-62 | 7.16724 | 25.3332 |
| DO | 0.821756 | 12.0033 | 127.478 |
| ZOA | 0.00074102 | 5.06592 | 0.513748 |

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