# BlackBoxOptimization

基于 MoonBit 实现的黑盒优化算法框架，支持 65 种经典优化算法，提供统一的接口和基准测试框架。

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
| 群体智能 | 麻雀搜索算法 (SSA) | 模拟麻雀觅食和反捕食行为 |
| 群体智能 | 哈里斯鹰优化 (HHO) | 模拟哈里斯鹰协同捕食行为 |
| 群体智能 | 花粉传播算法 (FPA) | 模拟花粉全局和局部传播过程 |
| 群体智能 | 蜻蜓算法 (DA) | 模拟蜻蜓分离、对齐、聚集行为 |
| 群体智能 | 乌鸦搜索算法 (CSA) | 模拟乌鸦记忆和追逐行为 |
| 群体智能 | 樽海鞘群算法 (SalpSA) | 模拟樽海鞘链状群体结构 |
| 群体智能 | Monarch蝴蝶优化 (MBO) | 模拟蝴蝶迁徙和调整行为 |
| 群体智能 | 人工生态系统优化 (AEO) | 模拟生态系统生产者-消费者-分解者 |
| 群体智能 | 入侵杂草优化 (IWO) | 模拟杂草生长扩散与竞争 |
| 群体智能 | 细菌觅食优化 (BFO) | 模拟大肠杆菌觅食行为 |
| 群体智能 | 混合蛙跳算法 (SFLA) | 模拟青蛙群体觅食与信息共享 |
| 群体智能 | 萤火虫群优化 (GSO) | 模拟萤火虫发光吸引行为 |
| 物理启发 | 引力搜索算法 (GSA) | 基于牛顿万有引力定律的搜索 |
| 教学模型 | 教与学优化 (TLBO) | 模拟班级教学和学习过程 |
| 群体智能 | 动物迁徙优化 (AMO) | 模拟动物群体迁徙跟随行为 |
| 物理启发 | 多宇宙优化 (MVO) | 基于白洞/黑洞/虫洞的宇宙机制 |
| 群体智能 | 磷虾群算法 (KH) | 模拟磷虾诱导运动、觅食和扩散 |
| 进化算法 | 帝国主义竞争算法 (ICA) | 模拟帝国殖民、同化和竞争 |
| 估计分布 | 随机分形搜索 (SFS) | 基于分形发现和学习两阶段搜索 |
| 群体智能 | 共生生物搜索 (SOS) | 模拟生物间互利共生/偏利共生/寄生关系 |
| 进化算法 | 生物地理学优化 (BBO) | 基于栖息地适宜度与物种迁移模型 |
| 群体智能 | 大象牧群优化 (EHO) | 模拟大象族群分离与母系领导机制 |
| 物理启发 | 均衡优化 (EO) | 基于控制理论均衡池与指数衰减策略 |
| 物理启发 | 风驱动优化 (WDO) | 模拟大气风速更新与科里奥利力 |
| 数学启发 | Jaya算法 | 基于最优-最差个体更新策略 |
| 人类行为 | 足球联赛竞争算法 (SLC) | 模拟足球联赛团队竞争与球员交换 |
| 群体智能 | 高斯骨干粒子群优化 (GBPSO) | 基于高斯分布的简化粒子群算法 |
| 群体智能 | 飞蛾搜索算法 (MSA) | 基于飞蛾 Lévy 飞行和螺旋搜索 |
| 群体智能 | 路径finder算法 (PFA) | 基于领导者-跟随者层级移动 |
| 群体智能 | 菌落捕食算法 (CPA) | 模拟捕食-逃跑概率行为 |
| 群体智能 | 细菌菌落优化 (BCO) | 基于细菌趋化聚集行为 |
| 群体智能 | 飞狐优化算法 (FFOA) | 带惯性权重的群体飞行搜索 |
| 群体智能 | 人工鱼群算法 (AFSA) | 模拟鱼群觅食、聚群和追尾行为 |
| 群体智能 | 烟花算法 (FWA) | 模拟烟花爆炸及火花分布机制 |
| 人类行为 | 头脑风暴优化 (BSO) | 模拟人类头脑风暴聚类的创意生成过程 |
| 进化算法 | 回溯搜索算法 (BSA) | 基于历史种群记忆的交叉变异 |
| 群体智能 | 鸽群优化 (PIO) | 模拟信鸽导航的地图指南针与地标两阶段 |
| 数学启发 | 算术优化算法 (AOA) | 基于算术运算符(除/乘/减/加)的搜索 |
| 群体智能 | 蜜獾算法 (HBA) | 模拟蜜獾挖洞与觅食的密度因子搜索 |
| 群体智能 | 鹈鹕优化算法 (POA) | 模拟鹈鹕捕猎的两阶段策略 |
| 群体智能 | 金枪鱼群优化 (TSA) | 基于螺旋运动和觅食行为的群体搜索 |
| 群体智能 | 金豺优化 (GJO) | 模拟雄雌金豺协同捕猎的双引导搜索 |

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
│   └── lib_test.mbt              # 单元测试（120 个测试用例）
├── example/                      # 示例包
│   ├── moon.pkg.json             # 包配置
│   └── simple_usage.mbt          # 使用示例（51 个示例）
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
# 运行所有测试（120 个测试用例）
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

### WOA (鲸群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 鲸鱼数量 | 25 |
| spiral_constant | 螺旋常数 b | 1.0 |

### Bat Algorithm (蝙蝠算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 蝙蝠数量 | 25 |
| freq_min | 最小频率 | 0.0 |
| freq_max | 最大频率 | 1.0 |
| loudness | 响度 | 0.5 |
| pulse_rate | 脉冲率 | 1.0 |

### Harmony Search (和声搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| memory_size | 和声记忆库大小 | 20 |
| hmcr | 和声记忆考虑率 | 0.9 |
| par | 音调调整率 | 0.3 |
| fretwidth | 频宽 | 0.05 |

### SCA (正弦余弦算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### MFO (飞蛾火焰优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 飞蛾数量 | 25 |
| spiral_constant | 螺旋常数 b | 1.0 |

### ALO (蚁狮优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 蚁狮数量 | 25 |

### SSA (麻雀搜索算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |
| discoverer_ratio | 发现者比例 | 0.2 |
| warning_threshold | 警戒阈值 | 0.8 |

### HHO (哈里斯鹰优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### FPA (花粉传播算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |
| switch_probability | 全局/本地切换概率 | 0.8 |

### DA (蜻蜓算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |
| separation_weight | 分离权重 | 1.0 |
| alignment_weight | 对齐权重 | 1.0 |
| cohesion_weight | 聚集权重 | 1.0 |

### CSA (乌鸦搜索算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 乌鸦数量 | 25 |
| flight_length | 飞行长度 | 2.0 |
| awareness_prob | 感知概率 | 0.1 |

### SalpSA (樽海鞘群算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 樽海鞘数量 | 25 |

### MBO (Monarch蝴蝶优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 蝴蝶数量 | 25 |
| migration_period | 迁徙周期 | 5 |
| adjustment_rate | 调整率 | 1.0 |

### AEO (人工生态系统优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 生态系统大小 | 25 |

### IWO (入侵杂草优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| initial_pop | 初始种群 | 10 |
| max_pop | 最大种群 | 40 |
| min_seeds | 最小种子数 | 1 |
| max_seeds | 最大种子数 | 5 |
| initial_sigma | 初始标准差 | 2.0 |
| final_sigma | 最终标准差 | 0.1 |
| power_exponent | 调制指数 | 3.0 |

### BFO (细菌觅食优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 细菌数量 | 25 |
| chemotactic_steps | 趋化步数 | 5 |
| reproduction_steps | 繁殖代数 | 3 |
| step_size | 步长 | 0.2 |

### SFLA (混合蛙跳算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 青蛙总数 | 30 |
| memeplex_count | 模因复合体数 | 5 |
| iterations_per_group | 组内迭代次数 | 10 |

### GSO (萤火虫群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 萤火虫数量 | 25 |
| luciferin_decay | 荧光素衰减率 | 0.4 |
| luciferin_enhancement | 荧光素增强率 | 0.6 |
| neighbor_count | 邻居数 | 5 |
| sensing_radius | 感知半径 | 5.0 |
| step_size | 步长 | 1.0 |

### GSA (引力搜索算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 智能体数量 | 25 |
| gravity_constant | 初始引力常数 G₀ | 100.0 |
| decay_alpha | 衰减系数 α | 0.99 |

### TLBO (教与学优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 学生数量 | 30 |

### AMO (动物迁徙优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 动物群大小 | 25 |
| migration_interval | 迁徙间隔 | 5 |

### MVO (多宇宙优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 宇宙数量 | 25 |
| wormhole_exploitation | 虫洞开发率 | 1.0 |
| travelling_distance_rate | 距离旅行率 | 0.6 |

### KH (磷虾群算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 磷虾数量 | 25 |
| max_speed | 最大速度 | 0.01 |
| foraging_speed | 觅食速度 | 0.02 |

### ICA (帝国主义竞争算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 国家数量 | 25 |
| imperialist_count | 帝国主义者数量 | 3 |
| assimilation_rate | 同化率 | 0.3 |
| revolution_rate | 革命率 | 0.3 |

### SFS (随机分形搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 搜索点数量 | 25 |
| discovery_iterations | 发现阶段迭代 | 3 |
| learning_iterations | 学习阶段迭代 | 2 |

### SOS (共生生物搜索)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### BBO (生物地理学优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 栖息地数量 | 25 |
| mutation_rate | 变异率 | 0.01 |

### EHO (大象牧群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 大象数量 | 25 |
| clan_count | 族群数量 | 5 |
| alpha | 更新权重 α | 0.5 |

### EO (均衡优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### WDO (风驱动优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 空气质点数量 | 25 |

### Jaya算法
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### SLC (足球联赛竞争算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 球员数量 | 25 |

### GBPSO (高斯骨干粒子群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 粒子数量 | 25 |

### MSA (飞蛾搜索算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 飞蛾数量 | 25 |

### PFA (路径finder算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 动物数量 | 25 |

### CPA (菌落捕食算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 个体数量 | 25 |

### BCO (细菌菌落优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 细菌数量 | 25 |

### FFOA (飞狐优化算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 飞狐数量 | 25 |

### AFSA (人工鱼群算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 人工鱼数量 | 25 |
| visual | 视野范围 | 2.0 |
| step | 移动步长 | 0.5 |
| crowd_factor | 拥挤度因子 | 0.5 |

### FWA (烟花算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 烟花数量 | 25 |

### BSO (头脑风暴优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 创意(个体)数量 | 25 |

### BSA (回溯搜索算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### PIO (鸽群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 鸽群大小 | 25 |

### AOA (算术优化算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 种群大小 | 25 |

### HBA (蜜獾算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 蜜獾数量 | 25 |

### POA (鹈鹕优化算法)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 鹈鹕数量 | 25 |

### TSA (金枪鱼群优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 金枪鱼数量 | 25 |

### GJO (金豺优化)
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| population_size | 金豺数量 | 25 |

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

## 算法分类特点

| 类别 | 代表算法 | 特点 |
|------|----------|------|
| **单点搜索** | 随机搜索、爬山法、模拟退火 | 简单、易于实现，适合低维问题 |
| **群体智能** | PSO、萤火虫、蚁群、蜂群、布谷鸟、灰狼、鲸群、蝙蝠、和声、正弦余弦、飞蛾火焰、蚁狮、麻雀搜索、哈里斯鹰、花粉传播、蜻蜓、乌鸦、樽海鞘、Monarch蝴蝶、人工生态系统、入侵杂草、细菌觅食、混合蛙跳、萤火虫群、动物迁徙、磷虾群、共生生物搜索、大象牧群优化、飞蛾搜索、路径finder、菌落捕食、细菌菌落、飞狐优化、人工鱼群、烟花、鸽群、蜜獾、鹈鹕、金枪鱼群、金豺、高斯骨干粒子群 | 基于群体行为，全局搜索能力强 |
| **进化算法** | DE、GA、帝国主义竞争算法、生物地理学优化、回溯搜索算法、Jaya算法、足球联赛竞争算法 | 基于种群进化，鲁棒性好 |
| **人类行为** | 头脑风暴优化、足球联赛竞争算法 | 模拟人类社会行为过程 |
| **数学启发** | 算术优化算法、Jaya算法 | 基于数学原理的搜索机制 |
| **估计分布** | CEM、随机分形搜索 | 基于概率分布估计，适合连续优化 |
| **确定性直接搜索** | Nelder-Mead、Hooke-Jeeves | 无需梯度，适合不可导函数 |
| **局部搜索** | 禁忌搜索、引导局部搜索 | 强化局部开发，避免陷入局部最优 |
| **进化策略** | CMA-ES | 协方差自适应，适合病态问题 |
| **物理启发** | GSA、MVO、均衡优化、风驱动优化 | 基于物理规律，自适应搜索 |
| **教学模型** | TLBO | 无需算法特定参数，收敛速度快 |

## 许可证

本项目基于 [Apache-2.0 license](LICENSE) 开源。