// Black-Box Optimization Framework - 主入口

///|
struct FuncConfig {
  f : (Array[Double]) -> Double
  name : String
  bounds : @lib.Bounds
  optimal : Double
}

///|
fn make_func_configs() -> Array[FuncConfig] {
  [
    FuncConfig::{
      f: @lib.sphere_function,
      name: "Sphere",
      bounds: @lib.make_bounds(-5.0, 5.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.rastrigin_function,
      name: "Rastrigin",
      bounds: @lib.make_bounds(-5.12, 5.12),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.rosenbrock_function,
      name: "Rosenbrock",
      bounds: @lib.make_bounds(-5.0, 5.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.ackley_function,
      name: "Ackley",
      bounds: @lib.make_bounds(-32.768, 32.768),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.griewank_function,
      name: "Griewank",
      bounds: @lib.make_bounds(-600.0, 600.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.schwefel_function,
      name: "Schwefel",
      bounds: @lib.make_bounds(-500.0, 500.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.levy_function,
      name: "Levy",
      bounds: @lib.make_bounds(-10.0, 10.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.michalewicz_function,
      name: "Michalewicz",
      bounds: @lib.make_bounds(0.0, 3.141592653589793),
      optimal: -4.687658,
    },
    FuncConfig::{
      f: @lib.zakharov_function,
      name: "Zakharov",
      bounds: @lib.make_bounds(-10.0, 10.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.booth_function,
      name: "Booth",
      bounds: @lib.make_bounds(-10.0, 10.0),
      optimal: 0.0,
    },
    FuncConfig::{
      f: @lib.matyas_function,
      name: "Matyas",
      bounds: @lib.make_bounds(-10.0, 10.0),
      optimal: 0.0,
    },
  ]
}

///|
fn main {
  println(
    "╔══════════════════════════════════════════════════════════╗",
  )
  println(
    "║      黑盒优化算法 - 基准测试与准确率对比                ║",
  )
  println(
    "╚══════════════════════════════════════════════════════════╝",
  )
  println("")

  let dim = 5
  let funcs = make_func_configs()

  // ===== 第一部分：数学函数验证 =====
  println(
    "┌──────────────────────────────────────────────────────────┐",
  )
  println(
    "│ 第一部分：数学函数验证                                    │",
  )
  println(
    "└──────────────────────────────────────────────────────────┘",
  )
  println("")

  println("  cos(0)     = \{@lib.cos_approx(0.0)} (expected: 1.0)")
  println("  cos(pi)    = \{@lib.cos_approx(@lib.PI)} (expected: -1.0)")
  println("  cos(pi/2)  = \{@lib.cos_approx(@lib.PI / 2.0)} (expected: 0.0)")
  println("  cos(3pi)   = \{@lib.cos_approx(3.0 * @lib.PI)} (expected: -1.0)")
  println("")

  let zeros : Array[Double] = [0.0, 0.0, 0.0, 0.0, 0.0]
  for i in 0..<funcs.length() {
    let cfg = funcs[i]
    println(
      "  \{cfg.name} at zeros = \{(cfg.f)(zeros)} (known optimum: \{cfg.optimal})",
    )
  }
  println("")

  // ===== 第二部分：评估次数对比 =====
  println(
    "┌──────────────────────────────────────────────────────────┐",
  )
  println(
    "│ 第二部分：不同评估次数的性能对比                          │",
  )
  println(
    "└──────────────────────────────────────────────────────────┘",
  )
  println("")

  let base_rng = @lib.make_rng(42.0)
  for i in 0..<funcs.length() {
    let cfg = funcs[i]
    println(
      "──────────────────────────────────────────────────────────",
    )
    @lib.run_evals_comparison(
      cfg.f,
      cfg.name,
      dim,
      cfg.bounds,
      cfg.optimal,
      base_rng,
    )
  }

  // ===== 第三部分：多次运行统计 =====
  println(
    "┌──────────────────────────────────────────────────────────┐",
  )
  println(
    "│ 第三部分：多次运行统计 (n=20, eval_budget=10000)          │",
  )
  println(
    "└──────────────────────────────────────────────────────────┘",
  )
  println("")

  let num_trials = 20
  let eval_budget = 10000

  // Keep the multi-trial statistics to the three classical functions used
  // by the README benchmark table to bound the total run time.
  let stat_seeds : Array[Double] = [100.0, 200.0, 300.0]
  for i in 0..<3 {
    let cfg = funcs[i]
    println("=== \{cfg.name} Function ===")
    @lib.run_multiple_trials(
      cfg.f,
      cfg.name,
      dim,
      cfg.bounds,
      cfg.optimal,
      eval_budget,
      num_trials,
      stat_seeds[i],
    )
    println("")
  }

  println("")
  println(
    "╔══════════════════════════════════════════════════════════╗",
  )
  println(
    "║                    测试完成                              ║",
  )
  println(
    "╚══════════════════════════════════════════════════════════╝",
  )
}
