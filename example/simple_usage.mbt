// 简单使用示例

import lib

fn main {
  println("=== Black-Box Optimization 使用示例 ===")
  println("")

  let dim = 3
  let bounds = lib.Bounds::{ min: -5.0, max: 5.0 }
  let rng = lib.RNG::{ seed: 123.0 }

  println("1. 运行随机搜索算法:")
  let rs_config = lib.RandomSearchConfig::{
    max_evaluations: 1000,
    dimension: dim,
    bounds: bounds
  }
  let rs_result = lib.random_search(rs_config, lib.sphere_function, rng)
  println("   最优值: \{rs_result.best_value}")
  println("   评估次数: \{rs_result.evaluations}")
  println("")

  println("2. 运行爬山法算法:")
  let hc_config = lib.HillClimbingConfig::{
    max_evaluations: 1000,
    dimension: dim,
    bounds: bounds,
    step_size: 0.5,
    restarts: 10
  }
  let hc_result = lib.hill_climbing(hc_config, lib.sphere_function, rng)
  println("   最优值: \{hc_result.best_value}")
  println("   评估次数: \{hc_result.evaluations}")
  println("")

  println("3. 运行模拟退火算法:")
  let sa_config = lib.SimulatedAnnealingConfig::{
    max_evaluations: 1000,
    dimension: dim,
    bounds: bounds,
    initial_temperature: 100.0,
    cooling_rate: 0.995,
    step_size: 0.5
  }
  let sa_result = lib.simulated_annealing(sa_config, lib.sphere_function, rng)
  println("   最优值: \{sa_result.best_value}")
  println("   评估次数: \{sa_result.evaluations}")
  println("")

  println("4. 比较不同评估预算下的性能:")
  let results = lib.run_benchmark_at_evals(
    lib.sphere_function,
    "Sphere Function",
    dim,
    bounds,
    0.0,
    5000,
    rng
  )
  lib.print_benchmark_table(results, "Sphere Function", 5000)
}