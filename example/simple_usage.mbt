// 简单使用示例：通过算法注册表运行全部算法

///|
fn main {
  println("=== Black-Box Optimization 使用示例 ===")
  println("")

  let dim = 3
  let bounds = @lib.make_bounds(-5.0, 5.0)
  let algos = @lib.default_algorithms()

  println("1. 使用注册表运行全部 \{algos.length()} 个算法")
  println("   测试函数: Sphere, 维度: \{dim}, 评估预算: 1000")
  println("")

  for i in 0..<algos.length() {
    let algo = algos[i]
    let result = (algo.run)(
      1000,
      dim,
      bounds,
      @lib.make_rng(123.0),
      @lib.sphere_function,
    )
    println("\{i + 1}. \{algo.name}")
    println("   最优值: \{result.best_value}")
    println("   评估次数: \{result.evaluations}")
    println("")
  }

  println("2. 比较不同评估预算下的性能:")
  let results = @lib.run_benchmark_at_evals(
    @lib.sphere_function,
    "Sphere Function",
    dim,
    bounds,
    0.0,
    5000,
    @lib.make_rng(123.0),
  )
  @lib.print_benchmark_table(results, "Sphere Function", 5000)
}
