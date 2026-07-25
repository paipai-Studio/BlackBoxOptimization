// Black-Box Optimization Framework
// 黑盒优化算法演示

// ===== 数学常量 =====

let PI: Double = 3.14159265358979323846
let TWO_PI: Double = 6.28318530717958647692

// ===== 工具函数 =====

// 伪随机数生成器
struct RNG {
  seed: Double
}

fn RNG::next(self: RNG) -> (RNG, Double) {
  let new_seed = (self.seed * 9301.0 + 49297.0) % 233280.0
  ({ seed: new_seed }, new_seed / 233280.0)
}

// 将角度规约到 [-PI, PI] 范围
fn reduce_to_pi(x: Double) -> Double {
  let k = (x / TWO_PI).round().to_double()
  x - k * TWO_PI
}

// 使用泰勒展开近似cos（已添加参数约简）
fn cos_approx(x: Double) -> Double {
  let rx = reduce_to_pi(x)
  let x2 = rx * rx
  let x4 = x2 * x2
  let x6 = x4 * x2
  let x8 = x4 * x4
  1.0 - x2 / 2.0 + x4 / 24.0 - x6 / 720.0 + x8 / 40320.0 - x6 * x2 / 3628800.0 + x8 * x2 / 479001600.0
}

// 使用泰勒展开近似exp（支持负数输入）
fn exp_approx(x: Double) -> Double {
  // 对于负数，利用 exp(-x) = 1/exp(x)
  if x < 0.0 {
    return 1.0 / exp_approx(-x)
  }
  let x2 = x * x
  let x3 = x2 * x
  let x4 = x3 * x
  let x5 = x4 * x
  let x6 = x5 * x
  1.0 + x + x2 / 2.0 + x3 / 6.0 + x4 / 24.0 + x5 / 120.0 + x6 / 720.0
}

// 自定义max
fn max(a: Double, b: Double) -> Double {
  if a > b { a } else { b }
}

// 自定义min
fn min(a: Double, b: Double) -> Double {
  if a < b { a } else { b }
}

// 计算绝对值
fn abs_val(x: Double) -> Double {
  if x < 0.0 { -x } else { x }
}

// ===== 核心类型定义 =====

// 边界约束
struct Bounds {
  min: Double
  max: Double
}

// 优化结果
struct OptimizationResult {
  best_solution: Array[Double]
  best_value: Double
  evaluations: Int
}

// 基准测试结果
struct BenchmarkResult {
  algorithm_name: String
  best_value: Double
  evaluations: Int
  relative_error: Double
}

// 统计结果（多次运行）
struct Stats {
  min_value: Double
  max_value: Double
  avg_value: Double
}

// ===== 辅助函数 =====

// 在给定边界内生成随机解
fn random_solution(dim: Int, bounds: Bounds, rng: RNG) -> (RNG, Array[Double]) {
  let solution: Array[Double] = []
  let mut current_rng = rng
  for i in 0..<dim {
    let (new_rng, r) = RNG::next(current_rng)
    current_rng = new_rng
    let rand_val: Double = bounds.min + (r * (bounds.max - bounds.min))
    solution.push(rand_val)
  }
  (current_rng, solution)
}

// 评估目标函数
fn evaluate(f: (Array[Double]) -> Double, solution: Array[Double]) -> Double {
  f(solution)
}

// ===== 随机搜索算法 =====

struct RandomSearchConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
}

fn random_search(config: RandomSearchConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations

  let mut current_rng = rng
  let (rng1, initial_solution) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  
  let mut best_solution = initial_solution
  let mut best_value = evaluate(f, best_solution)
  let mut evaluations = 1

  for i in 1..<max_evals {
    let (rng2, candidate) = random_solution(dim, bounds, current_rng)
    current_rng = rng2
    let value = evaluate(f, candidate)
    evaluations = evaluations + 1

    if value < best_value {
      best_value = value
      best_solution = candidate
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 爬山法算法 =====

struct HillClimbingConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  step_size: Double
  restarts: Int
}

fn generate_neighbor(solution: Array[Double], step_size: Double, bounds: Bounds, rng: RNG) -> (RNG, Array[Double]) {
  let neighbor: Array[Double] = []
  let mut current_rng = rng
  for i in 0..<solution.length() {
    let (new_rng, r) = RNG::next(current_rng)
    current_rng = new_rng
    let perturbation: Double = (r * 2.0 - 1.0) * step_size
    let new_val: Double = max(bounds.min, min(bounds.max, solution[i] + perturbation))
    neighbor.push(new_val)
  }
  (current_rng, neighbor)
}

fn hill_climb(f: (Array[Double]) -> Double, config: HillClimbingConfig, start_solution: Array[Double], start_value: Double, remaining_evals: Int, rng: RNG) -> (Array[Double], Double, Int, RNG) {
  let mut current = start_solution
  let mut current_value = start_value
  let mut evaluations = 0
  let mut current_rng = rng

  while evaluations < remaining_evals {
    let (rng1, neighbor) = generate_neighbor(current, config.step_size, config.bounds, current_rng)
    current_rng = rng1
    let neighbor_value = evaluate(f, neighbor)
    evaluations = evaluations + 1

    if neighbor_value < current_value {
      current = neighbor
      current_value = neighbor_value
    }
  }

  (current, current_value, evaluations, current_rng)
}

fn hill_climbing(config: HillClimbingConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations

  let mut best_solution: Array[Double] = []
  let mut best_value = 1.0e308
  let mut total_evaluations = 0
  let mut current_rng = rng

  for restart in 0..<config.restarts {
    if total_evaluations >= max_evals {
      break
    }

    let (rng1, solution) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    let value = evaluate(f, solution)
    total_evaluations = total_evaluations + 1

    let remaining = max_evals - total_evaluations
    let (final_sol, final_val, evals, rng2) = hill_climb(f, config, solution, value, remaining, current_rng)
    current_rng = rng2
    total_evaluations = total_evaluations + evals

    if final_val < best_value {
      best_value = final_val
      best_solution = final_sol
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: total_evaluations
  }
}

// ===== 模拟退火算法 =====

struct SimulatedAnnealingConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  initial_temperature: Double
  cooling_rate: Double
  step_size: Double
}

fn acceptance_probability(current_value: Double, new_value: Double, temperature: Double) -> Double {
  if new_value < current_value {
    return 1.0
  }
  let delta = new_value - current_value
  exp_approx(-delta / temperature)
}

fn generate_neighbor_sa(solution: Array[Double], step_size: Double, bounds: Bounds, rng: RNG) -> (RNG, Array[Double]) {
  let neighbor: Array[Double] = []
  let mut current_rng = rng
  for i in 0..<solution.length() {
    let (new_rng, r) = RNG::next(current_rng)
    current_rng = new_rng
    let perturbation: Double = (r * 2.0 - 1.0) * step_size
    let new_val: Double = max(bounds.min, min(bounds.max, solution[i] + perturbation))
    neighbor.push(new_val)
  }
  (current_rng, neighbor)
}

fn simulated_annealing(config: SimulatedAnnealingConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let initial_temp = config.initial_temperature
  let cooling_rate = config.cooling_rate
  let step_size = config.step_size

  let mut current_rng = rng
  let (rng1, initial_solution) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  
  let mut current_solution = initial_solution
  let mut current_value = evaluate(f, current_solution)
  let mut best_solution = current_solution
  let mut best_value = current_value
  let mut temperature = initial_temp
  let mut evaluations = 1

  while evaluations < max_evals {
    let (rng2, neighbor) = generate_neighbor_sa(current_solution, step_size, bounds, current_rng)
    current_rng = rng2
    let neighbor_value = evaluate(f, neighbor)
    evaluations = evaluations + 1

    let prob = acceptance_probability(current_value, neighbor_value, temperature)
    let (rng3, rand_val) = RNG::next(current_rng)
    current_rng = rng3

    if rand_val < prob {
      current_solution = neighbor
      current_value = neighbor_value

      if current_value < best_value {
        best_value = current_value
        best_solution = current_solution
      }
    }

    temperature = temperature * cooling_rate
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 测试函数 =====

// Sphere函数：f(x) = sum(x_i^2)，最优值为0
fn sphere_function(x: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<x.length() {
    sum = sum + x[i] * x[i]
  }
  sum
}

// Rastrigin函数，最优值为0
fn rastrigin_function(x: Array[Double]) -> Double {
  let n: Double = x.length().to_double()
  let mut sum = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum = sum + xi * xi - 10.0 * cos_approx(2.0 * PI * xi)
  }
  sum + 10.0 * n
}

// Rosenbrock函数，最优值为0
fn rosenbrock_function(x: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<(x.length() - 1) {
    let xi = x[i]
    let xi1 = x[i + 1]
    let term1 = 100.0 * (xi1 - xi * xi) * (xi1 - xi * xi)
    let term2 = (1.0 - xi) * (1.0 - xi)
    sum = sum + term1 + term2
  }
  sum
}

// ===== 基准测试框架 =====

// 计算相对误差
fn calc_relative_error(best_value: Double, known_optimal: Double) -> Double {
  let diff = abs_val(best_value - known_optimal)
  if known_optimal == 0.0 {
    if best_value == 0.0 { 0.0 } else { diff }
  } else {
    diff / abs_val(known_optimal)
  }
}

// 数组求和
fn array_sum(arr: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<arr.length() {
    sum = sum + arr[i]
  }
  sum
}

// 数组最小值
fn array_min(arr: Array[Double]) -> Double {
  let mut min_val = arr[0]
  for i in 1..<arr.length() {
    if arr[i] < min_val {
      min_val = arr[i]
    }
  }
  min_val
}

// 数组最大值
fn array_max(arr: Array[Double]) -> Double {
  let mut max_val = arr[0]
  for i in 1..<arr.length() {
    if arr[i] > max_val {
      max_val = arr[i]
    }
  }
  max_val
}

// 数组平均值
fn array_avg(arr: Array[Double]) -> Double {
  array_sum(arr) / arr.length().to_double()
}

// 运行单个算法的基准测试
fn run_single_benchmark(
  algorithm: String,
  f: (Array[Double]) -> Double,
  config: RandomSearchConfig,
  bounds: Bounds,
  known_optimal: Double,
  rng: RNG
) -> BenchmarkResult {
  let result = random_search(config, f, rng)
  let error = calc_relative_error(result.best_value, known_optimal)
  BenchmarkResult::{
    algorithm_name: algorithm,
    best_value: result.best_value,
    evaluations: result.evaluations,
    relative_error: error
  }
}

// 运行所有算法在特定评估次数下的基准测试
fn run_benchmark_at_evals(
  f: (Array[Double]) -> Double,
  func_name: String,
  dim: Int,
  bounds: Bounds,
  known_optimal: Double,
  max_evals: Int,
  rng: RNG
) -> Array[BenchmarkResult] {
  let results: Array[BenchmarkResult] = []
  
  // 1. 随机搜索
  let rs_config = RandomSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds
  }
  let rs_result = random_search(rs_config, f, rng)
  let rs_error = calc_relative_error(rs_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Random Search",
    best_value: rs_result.best_value,
    evaluations: rs_result.evaluations,
    relative_error: rs_error
  })
  
  // 2. 爬山法
  let hc_config = HillClimbingConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    step_size: 0.5,
    restarts: 10
  }
  let hc_result = hill_climbing(hc_config, f, rng)
  let hc_error = calc_relative_error(hc_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Hill Climbing",
    best_value: hc_result.best_value,
    evaluations: hc_result.evaluations,
    relative_error: hc_error
  })
  
  // 3. 模拟退火
  let sa_config = SimulatedAnnealingConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    initial_temperature: 100.0,
    cooling_rate: 0.995,
    step_size: 0.5
  }
  let sa_result = simulated_annealing(sa_config, f, rng)
  let sa_error = calc_relative_error(sa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Simulated Annealing",
    best_value: sa_result.best_value,
    evaluations: sa_result.evaluations,
    relative_error: sa_error
  })
  
  results
}

// 打印基准测试结果表格
fn print_benchmark_table(results: Array[BenchmarkResult], func_name: String, max_evals: Int) -> Unit {
  println("  \{func_name} (evals = \{max_evals})")
  println("  +--------------------+-------------------+------------+-------------------+")
  println("  | Algorithm          | Best Value        | Evals      | Relative Error    |")
  println("  +--------------------+-------------------+------------+-------------------+")
  for i in 0..<results.length() {
    let r = results[i]
    let err_str = "\{r.relative_error:.6f}"
    println("  | \{r.algorithm_name} | \{r.best_value} | \{r.evaluations} | \{err_str} |")
  }
  println("  +--------------------+-------------------+------------+-------------------+")
  println("")
}

// 统计多次运行的结果
fn run_multiple_trials(
  f: (Array[Double]) -> Double,
  func_name: String,
  dim: Int,
  bounds: Bounds,
  known_optimal: Double,
  max_evals: Int,
  num_trials: Int,
  base_seed: Double
) -> Unit {
  println("=== 多次运行统计 (num_trials = \{num_trials}) ===")
  println("")
  
  let rs_values: Array[Double] = []
  let hc_values: Array[Double] = []
  let sa_values: Array[Double] = []
  
  for trial in 0..<num_trials {
    let seed = base_seed + trial.to_double() * 1000.0
    let rng = RNG::{ seed: seed }
    
    // 随机搜索
    let rs_config = RandomSearchConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds
    }
    let rs_result = random_search(rs_config, f, rng)
    rs_values.push(rs_result.best_value)
    
    // 爬山法
    let hc_config = HillClimbingConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      step_size: 0.5,
      restarts: 10
    }
    let hc_result = hill_climbing(hc_config, f, rng)
    hc_values.push(hc_result.best_value)
    
    // 模拟退火
    let sa_config = SimulatedAnnealingConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      initial_temperature: 100.0,
      cooling_rate: 0.995,
      step_size: 0.5
    }
    let sa_result = simulated_annealing(sa_config, f, rng)
    sa_values.push(sa_result.best_value)
  }
  
  // 打印统计结果
  println("  Algorithm          | Min Value           | Max Value           | Avg Value           | Avg Error")
  println("  +-------------------+---------------------+---------------------+---------------------+------------------+")
  
  // 随机搜索统计
  let rs_min = array_min(rs_values)
  let rs_max = array_max(rs_values)
  let rs_avg = array_avg(rs_values)
  let rs_err = calc_relative_error(rs_avg, known_optimal)
  println("  | Random Search     | \{rs_min} | \{rs_max} | \{rs_avg} | \{rs_err} |")
  
  // 爬山法统计
  let hc_min = array_min(hc_values)
  let hc_max = array_max(hc_values)
  let hc_avg = array_avg(hc_values)
  let hc_err = calc_relative_error(hc_avg, known_optimal)
  println("  | Hill Climbing     | \{hc_min} | \{hc_max} | \{hc_avg} | \{hc_err} |")
  
  // 模拟退火统计
  let sa_min = array_min(sa_values)
  let sa_max = array_max(sa_values)
  let sa_avg = array_avg(sa_values)
  let sa_err = calc_relative_error(sa_avg, known_optimal)
  println("  | Simulated Annealing | \{sa_min} | \{sa_max} | \{sa_avg} | \{sa_err} |")
  
  println("")
}

// 运行不同评估次数的对比
fn run_evals_comparison(
  f: (Array[Double]) -> Double,
  func_name: String,
  dim: Int,
  bounds: Bounds,
  known_optimal: Double,
  rng: RNG
) -> Unit {
  println("=== 评估次数对比 ===")
  println("  函数: \{func_name}")
  println("  维度: \{dim}")
  println("  已知最优值: \{known_optimal}")
  println("")
  
  let eval_budgets: Array[Int] = [100, 500, 1000, 5000, 10000]
  
  for i in 0..<eval_budgets.length() {
    let max_evals = eval_budgets[i]
    let results = run_benchmark_at_evals(f, func_name, dim, bounds, known_optimal, max_evals, rng)
    print_benchmark_table(results, func_name, max_evals)
  }
}

// ===== 主程序 =====

fn main {
  println("╔══════════════════════════════════════════════════════════╗")
  println("║      黑盒优化算法 - 基准测试与准确率对比                ║")
  println("╚══════════════════════════════════════════════════════════╝")
  println("")

  let dim = 5
  
  // 测试函数配置
  let sphere_config = {
    f: sphere_function,
    name: "Sphere Function",
    bounds: Bounds::{ min: -5.0, max: 5.0 },
    optimal: 0.0
  }
  
  let rastrigin_config = {
    f: rastrigin_function,
    name: "Rastrigin Function",
    bounds: Bounds::{ min: -5.12, max: 5.12 },
    optimal: 0.0
  }
  
  let rosenbrock_config = {
    f: rosenbrock_function,
    name: "Rosenbrock Function",
    bounds: Bounds::{ min: -5.0, max: 5.0 },
    optimal: 0.0
  }

  // ===== 第一部分：修复验证 =====
  println("┌──────────────────────────────────────────────────────────┐")
  println("│ 第一部分：数学函数验证                                    │")
  println("└──────────────────────────────────────────────────────────┘")
  println("")
  
  // 验证cos函数
  println("  cos(0)     = \{cos_approx(0.0)} (expected: 1.0)")
  println("  cos(pi)    = \{cos_approx(PI)} (expected: -1.0)")
  println("  cos(pi/2)  = \{cos_approx(PI / 2.0)} (expected: 0.0)")
  println("  cos(3pi)   = \{cos_approx(3.0 * PI)} (expected: -1.0)")
  println("  cos(16.0)  = \{cos_approx(16.0)} (test large input)")
  println("")
  
  // 验证Rastrigin函数在最优点的值
  let zeros: Array[Double] = [0.0, 0.0, 0.0, 0.0, 0.0]
  let rast_at_optimal = rastrigin_function(zeros)
  println("  Rastrigin at (0,0,0,0,0) = \{rast_at_optimal} (expected: 0.0)")
  println("")

  // ===== 第二部分：评估次数对比 =====
  println("┌──────────────────────────────────────────────────────────┐")
  println("│ 第二部分：不同评估次数的性能对比                          │")
  println("└──────────────────────────────────────────────────────────┘")
  println("")
  
  let base_rng = RNG::{ seed: 42.0 }
  
  // Sphere函数测试
  println("──────────────────────────────────────────────────────────")
  run_evals_comparison(sphere_config.f, sphere_config.name, dim, sphere_config.bounds, sphere_config.optimal, base_rng)
  
  // Rastrigin函数测试
  println("──────────────────────────────────────────────────────────")
  run_evals_comparison(rastrigin_config.f, rastrigin_config.name, dim, rastrigin_config.bounds, rastrigin_config.optimal, base_rng)
  
  // Rosenbrock函数测试
  println("──────────────────────────────────────────────────────────")
  run_evals_comparison(rosenbrock_config.f, rosenbrock_config.name, dim, rosenbrock_config.bounds, rosenbrock_config.optimal, base_rng)

  // ===== 第三部分：多次运行统计 =====
  println("┌──────────────────────────────────────────────────────────┐")
  println("│ 第三部分：多次运行统计 (n=20, eval_budget=10000)          │")
  println("└──────────────────────────────────────────────────────────┘")
  println("")
  
  let num_trials = 20
  let eval_budget = 10000
  
  // Sphere函数
  println("=== Sphere Function ===")
  run_multiple_trials(sphere_config.f, sphere_config.name, dim, sphere_config.bounds, sphere_config.optimal, eval_budget, num_trials, 100.0)
  println("")
  
  // Rastrigin函数
  println("=== Rastrigin Function ===")
  run_multiple_trials(rastrigin_config.f, rastrigin_config.name, dim, rastrigin_config.bounds, rastrigin_config.optimal, eval_budget, num_trials, 200.0)
  println("")
  
  // Rosenbrock函数
  println("=== Rosenbrock Function ===")
  run_multiple_trials(rosenbrock_config.f, rosenbrock_config.name, dim, rosenbrock_config.bounds, rosenbrock_config.optimal, eval_budget, num_trials, 300.0)
  
  println("")
  println("╔══════════════════════════════════════════════════════════╗")
  println("║                    测试完成                              ║")
  println("╚══════════════════════════════════════════════════════════╝")
}
