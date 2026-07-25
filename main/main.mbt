// Black-Box Optimization Framework - 主入口

struct FuncConfig {
  f: (Array[Double]) -> Double
  name: String
  bounds: @lib.Bounds
  optimal: Double
}

fn main {
  println("╔══════════════════════════════════════════════════════════╗")
  println("║      黑盒优化算法 - 基准测试与准确率对比                ║")
  println("╚══════════════════════════════════════════════════════════╝")
  println("")

  let dim = 5

  let sphere_config = FuncConfig::{
    f: @lib.sphere_function,
    name: "Sphere Function",
    bounds: @lib.make_bounds(-5.0, 5.0),
    optimal: 0.0
  }

  let rastrigin_config = FuncConfig::{
    f: @lib.rastrigin_function,
    name: "Rastrigin Function",
    bounds: @lib.make_bounds(-5.12, 5.12),
    optimal: 0.0
  }

  let rosenbrock_config = FuncConfig::{
    f: @lib.rosenbrock_function,
    name: "Rosenbrock Function",
    bounds: @lib.make_bounds(-5.0, 5.0),
    optimal: 0.0
  }

  // ===== 第一部分：修复验证 =====
  println("┌──────────────────────────────────────────────────────────┐")
  println("│ 第一部分：数学函数验证                                    │")
  println("└──────────────────────────────────────────────────────────┘")
  println("")

  println("  cos(0)     = \{@lib.cos_approx(0.0)} (expected: 1.0)")
  println("  cos(pi)    = \{@lib.cos_approx(@lib.PI)} (expected: -1.0)")
  println("  cos(pi/2)  = \{@lib.cos_approx(@lib.PI / 2.0)} (expected: 0.0)")
  println("  cos(3pi)   = \{@lib.cos_approx(3.0 * @lib.PI)} (expected: -1.0)")
  println("  cos(16.0)  = \{@lib.cos_approx(16.0)} (test large input)")
  println("")

  let zeros: Array[Double] = [0.0, 0.0, 0.0, 0.0, 0.0]
  let rast_at_optimal = @lib.rastrigin_function(zeros)
  println("  Rastrigin at (0,0,0,0,0) = \{rast_at_optimal} (expected: 0.0)")
  println("")

  // ===== 第二部分：评估次数对比 =====
  println("┌──────────────────────────────────────────────────────────┐")
  println("│ 第二部分：不同评估次数的性能对比                          │")
  println("└──────────────────────────────────────────────────────────┘")
  println("")

  let base_rng = @lib.make_rng(42.0)

  println("──────────────────────────────────────────────────────────")
  @lib.run_evals_comparison(sphere_config.f, sphere_config.name, dim, sphere_config.bounds, sphere_config.optimal, base_rng)

  println("──────────────────────────────────────────────────────────")
  @lib.run_evals_comparison(rastrigin_config.f, rastrigin_config.name, dim, rastrigin_config.bounds, rastrigin_config.optimal, base_rng)

  println("──────────────────────────────────────────────────────────")
  @lib.run_evals_comparison(rosenbrock_config.f, rosenbrock_config.name, dim, rosenbrock_config.bounds, rosenbrock_config.optimal, base_rng)

  // ===== 第三部分：多次运行统计 =====
  println("┌──────────────────────────────────────────────────────────┐")
  println("│ 第三部分：多次运行统计 (n=20, eval_budget=10000)          │")
  println("└──────────────────────────────────────────────────────────┘")
  println("")

  let num_trials = 20
  let eval_budget = 10000

  println("=== Sphere Function ===")
  @lib.run_multiple_trials(sphere_config.f, sphere_config.name, dim, sphere_config.bounds, sphere_config.optimal, eval_budget, num_trials, 100.0)
  println("")

  println("=== Rastrigin Function ===")
  @lib.run_multiple_trials(rastrigin_config.f, rastrigin_config.name, dim, rastrigin_config.bounds, rastrigin_config.optimal, eval_budget, num_trials, 200.0)
  println("")

  println("=== Rosenbrock Function ===")
  @lib.run_multiple_trials(rosenbrock_config.f, rosenbrock_config.name, dim, rosenbrock_config.bounds, rosenbrock_config.optimal, eval_budget, num_trials, 300.0)

  println("")
  println("╔══════════════════════════════════════════════════════════╗")
  println("║                    测试完成                              ║")
  println("╚══════════════════════════════════════════════════════════╝")
}