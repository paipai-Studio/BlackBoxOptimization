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

  println("7. 比较不同评估预算下的性能:")
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