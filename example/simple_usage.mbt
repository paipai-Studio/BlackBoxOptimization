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

  println("4. 比较不同评估预算下的性能:")
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
