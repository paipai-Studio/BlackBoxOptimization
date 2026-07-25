// 简单使用示例

fn main {
  println("=== Black-Box Optimization 使用示例 ===")
  println("")

  let dim = 3
  let bounds = @lib.make_bounds(-5.0, 5.0)
  let rng = @lib.make_rng(123.0)

  println("1. 运行随机搜索算法:")
  let rs_config = @lib.make_random_search_config(1000, dim, bounds)
  let rs_result = @lib.random_search(rs_config, @lib.sphere_function, rng)
  println("   最优值: \{rs_result.best_value}")
  println("   评估次数: \{rs_result.evaluations}")
  println("")

  println("2. 运行爬山法算法:")
  let hc_config = @lib.make_hill_climbing_config(1000, dim, bounds, 0.5, 10)
  let hc_result = @lib.hill_climbing(hc_config, @lib.sphere_function, rng)
  println("   最优值: \{hc_result.best_value}")
  println("   评估次数: \{hc_result.evaluations}")
  println("")

  println("3. 运行模拟退火算法:")
  let sa_config = @lib.make_simulated_annealing_config(1000, dim, bounds, 100.0, 0.995, 0.5)
  let sa_result = @lib.simulated_annealing(sa_config, @lib.sphere_function, rng)
  println("   最优值: \{sa_result.best_value}")
  println("   评估次数: \{sa_result.evaluations}")
  println("")

  println("4. 运行交叉熵方法 (CEM):")
  let cem_config = @lib.make_cem_config(1000, dim, bounds, 100, 0.1, 0.7)
  let cem_result = @lib.cross_entropy_method(cem_config, @lib.sphere_function, rng)
  println("   最优值: \{cem_result.best_value}")
  println("   评估次数: \{cem_result.evaluations}")
  println("")

  println("5. 运行灰狼优化 (GWO):")
  let gwo_config = @lib.make_gwo_config(1000, dim, bounds, 25)
  let gwo_result = @lib.grey_wolf_optimizer(gwo_config, @lib.sphere_function, rng)
  println("   最优值: \{gwo_result.best_value}")
  println("   评估次数: \{gwo_result.evaluations}")
  println("")

  println("6. 运行 Hooke-Jeeves 模式搜索:")
  let hj_config = @lib.make_hooke_jeeves_config(1000, dim, bounds, 0.5, 0.5, 1.0e-6)
  let hj_result = @lib.hooke_jeeves(hj_config, @lib.sphere_function, rng)
  println("   最优值: \{hj_result.best_value}")
  println("   评估次数: \{hj_result.evaluations}")
  println("")

  println("7. 运行鲸群优化 (WOA):")
  let woa_config = @lib.make_woa_config(1000, dim, bounds, 25, 1.0)
  let woa_result = @lib.whale_optimization(woa_config, @lib.sphere_function, rng)
  println("   最优值: \{woa_result.best_value}")
  println("   评估次数: \{woa_result.evaluations}")
  println("")

  println("8. 运行蝙蝠算法 (Bat Algorithm):")
  let bat_config = @lib.make_bat_config(1000, dim, bounds, 25, 0.0, 1.0, 0.5, 1.0)
  let bat_result = @lib.bat_algorithm(bat_config, @lib.sphere_function, rng)
  println("   最优值: \{bat_result.best_value}")
  println("   评估次数: \{bat_result.evaluations}")
  println("")

  println("9. 运行和声搜索 (Harmony Search):")
  let hs_config = @lib.make_hs_config(1000, dim, bounds, 20, 0.9, 0.3, 0.05)
  let hs_result = @lib.harmony_search(hs_config, @lib.sphere_function, rng)
  println("   最优值: \{hs_result.best_value}")
  println("   评估次数: \{hs_result.evaluations}")
  println("")

  println("10. 运行正弦余弦算法 (SCA):")
  let sca_config = @lib.make_sca_config(1000, dim, bounds, 25)
  let sca_result = @lib.sine_cosine_algorithm(sca_config, @lib.sphere_function, rng)
  println("   最优值: \{sca_result.best_value}")
  println("   评估次数: \{sca_result.evaluations}")
  println("")

  println("11. 运行飞蛾火焰优化 (MFO):")
  let mfo_config = @lib.make_mfo_config(1000, dim, bounds, 25, 1.0)
  let mfo_result = @lib.moth_flame_optimization(mfo_config, @lib.sphere_function, rng)
  println("   最优值: \{mfo_result.best_value}")
  println("   评估次数: \{mfo_result.evaluations}")
  println("")

  println("12. 运行蚁狮优化 (ALO):")
  let alo_config = @lib.make_alo_config(1000, dim, bounds, 25)
  let alo_result = @lib.antlion_optimizer(alo_config, @lib.sphere_function, rng)
  println("   最优值: \{alo_result.best_value}")
  println("   评估次数: \{alo_result.evaluations}")
  println("")

  println("13. 运行麻雀搜索算法 (SSA):")
  let ssa_config = @lib.make_ssa_config(1000, dim, bounds, 25)
  let ssa_result = @lib.sparrow_search_algorithm(ssa_config, @lib.sphere_function, rng)
  println("   最优值: \{ssa_result.best_value}")
  println("   评估次数: \{ssa_result.evaluations}")
  println("")

  println("14. 运行哈里斯鹰优化 (HHO):")
  let hho_config = @lib.make_hho_config(1000, dim, bounds, 25)
  let hho_result = @lib.harris_hawks_optimization(hho_config, @lib.sphere_function, rng)
  println("   最优值: \{hho_result.best_value}")
  println("   评估次数: \{hho_result.evaluations}")
  println("")

  println("15. 运行花粉传播算法 (FPA):")
  let fpa_config = @lib.make_fpa_config(1000, dim, bounds, 25)
  let fpa_result = @lib.flower_pollination_algorithm(fpa_config, @lib.sphere_function, rng)
  println("   最优值: \{fpa_result.best_value}")
  println("   评估次数: \{fpa_result.evaluations}")
  println("")

  println("16. 运行蜻蜓算法 (DA):")
  let da_config = @lib.make_da_config(1000, dim, bounds, 25)
  let da_result = @lib.dragonfly_algorithm(da_config, @lib.sphere_function, rng)
  println("   最优值: \{da_result.best_value}")
  println("   评估次数: \{da_result.evaluations}")
  println("")

  println("17. 比较不同评估预算下的性能:")
  let results = @lib.run_benchmark_at_evals(
    @lib.sphere_function,
    "Sphere Function",
    dim,
    bounds,
    0.0,
    5000,
    rng
  )
  @lib.print_benchmark_table(results, "Sphere Function", 5000)
}