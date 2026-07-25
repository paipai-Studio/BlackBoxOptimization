// Black-Box Optimization Framework - 核心库

// ===== 数学常量 ===== 

pub const PI: Double = 3.14159265358979323846
pub const TWO_PI: Double = 6.28318530717958647692

// ===== 工具函数 ===== 

pub struct RNG {
  pub seed: Double
}

pub fn RNG::next(self: RNG) -> (RNG, Double) {
  let new_seed = (self.seed * 9301.0 + 49297.0) % 233280.0
  ({ seed: new_seed }, new_seed / 233280.0)
}

pub fn reduce_to_pi(x: Double) -> Double {
  let k = (x / TWO_PI).round()
  x - k * TWO_PI
}

pub fn cos_approx(x: Double) -> Double {
  let rx = reduce_to_pi(x)
  let x2 = rx * rx
  let x4 = x2 * x2
  let x6 = x4 * x2
  let x8 = x4 * x4
  1.0 - x2 / 2.0 + x4 / 24.0 - x6 / 720.0 + x8 / 40320.0 - x6 * x2 / 3628800.0 + x8 * x2 / 479001600.0
}

pub fn exp_approx(x: Double) -> Double {
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

pub fn max(a: Double, b: Double) -> Double {
  if a > b { a } else { b }
}

pub fn min(a: Double, b: Double) -> Double {
  if a < b { a } else { b }
}

pub fn abs_val(x: Double) -> Double {
  if x < 0.0 { -x } else { x }
}

// ===== 核心类型定义 ===== 

pub struct Bounds {
  pub min: Double
  pub max: Double
}

pub struct OptimizationResult {
  pub best_solution: Array[Double]
  pub best_value: Double
  pub evaluations: Int
}

pub struct BenchmarkResult {
  pub algorithm_name: String
  pub best_value: Double
  pub evaluations: Int
  pub relative_error: Double
}

pub struct Stats {
  pub min_value: Double
  pub max_value: Double
  pub avg_value: Double
}

pub struct TestConfig {
  pub f: (Array[Double]) -> Double
  pub name: String
  pub bounds: Bounds
  pub optimal: Double
}

pub fn make_test_config(f: (Array[Double]) -> Double, name: String, bounds: Bounds, optimal: Double) -> TestConfig {
  TestConfig::{ f: f, name: name, bounds: bounds, optimal: optimal }
}

// ===== 辅助函数 ===== 

pub fn random_solution(dim: Int, bounds: Bounds, rng: RNG) -> (RNG, Array[Double]) {
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

pub fn evaluate(f: (Array[Double]) -> Double, solution: Array[Double]) -> Double {
  f(solution)
}


// ===== 构造函数 ===== 

pub fn make_bounds(min: Double, max: Double) -> Bounds {
  Bounds::{ min: min, max: max }
}

pub fn make_rng(seed: Double) -> RNG {
  RNG::{ seed: seed }
}

pub fn make_random_search_config(max_evaluations: Int, dimension: Int, bounds: Bounds) -> RandomSearchConfig {
  RandomSearchConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds }
}

pub fn make_hill_climbing_config(max_evaluations: Int, dimension: Int, bounds: Bounds, step_size: Double, restarts: Int) -> HillClimbingConfig {
  HillClimbingConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, step_size: step_size, restarts: restarts }
}

pub fn make_simulated_annealing_config(max_evaluations: Int, dimension: Int, bounds: Bounds, initial_temperature: Double, cooling_rate: Double, step_size: Double) -> SimulatedAnnealingConfig {
  SimulatedAnnealingConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, initial_temperature: initial_temperature, cooling_rate: cooling_rate, step_size: step_size }
}

pub fn make_pso_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, inertia_weight: Double, cognitive_coeff: Double, social_coeff: Double) -> PSOConfig {
  PSOConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, inertia_weight: inertia_weight, cognitive_coeff: cognitive_coeff, social_coeff: social_coeff }
}

pub fn make_de_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, mutation_factor: Double, crossover_probability: Double) -> DEConfig {
  DEConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, mutation_factor: mutation_factor, crossover_probability: crossover_probability }
}

pub fn make_ga_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, crossover_rate: Double, mutation_rate: Double) -> GAConfig {
  GAConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, crossover_rate: crossover_rate, mutation_rate: mutation_rate }
}

pub fn make_nm_config(max_evaluations: Int, dimension: Int, bounds: Bounds, reflection_coeff: Double, expansion_coeff: Double, contraction_coeff: Double, shrink_coeff: Double) -> NMConfig {
  NMConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, reflection_coeff: reflection_coeff, expansion_coeff: expansion_coeff, contraction_coeff: contraction_coeff, shrink_coeff: shrink_coeff }
}

pub fn make_firefly_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, initial_attractiveness: Double, light_absorption: Double, step_size: Double) -> FireflyConfig {
  FireflyConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, initial_attractiveness: initial_attractiveness, light_absorption: light_absorption, step_size: step_size }
}

// ===== 随机搜索算法 ===== 

pub struct RandomSearchConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
}

pub fn random_search(config: RandomSearchConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
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

// ===== 爬山法算法 ======

pub struct HillClimbingConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub step_size: Double
  pub restarts: Int
}

pub fn generate_neighbor(solution: Array[Double], step_size: Double, bounds: Bounds, rng: RNG) -> (RNG, Array[Double]) {
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

pub fn hill_climb(f: (Array[Double]) -> Double, config: HillClimbingConfig, start_solution: Array[Double], start_value: Double, remaining_evals: Int, rng: RNG) -> (Array[Double], Double, Int, RNG) {
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

pub fn hill_climbing(config: HillClimbingConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
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

// ===== 模拟退火算法 ======

pub struct SimulatedAnnealingConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub initial_temperature: Double
  pub cooling_rate: Double
  pub step_size: Double
}

pub fn acceptance_probability(current_value: Double, new_value: Double, temperature: Double) -> Double {
  if new_value < current_value {
    return 1.0
  }
  let delta = new_value - current_value
  exp_approx(-delta / temperature)
}

pub fn generate_neighbor_sa(solution: Array[Double], step_size: Double, bounds: Bounds, rng: RNG) -> (RNG, Array[Double]) {
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

pub fn simulated_annealing(config: SimulatedAnnealingConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
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

// ===== 粒子群优化算法 ===== 

pub struct PSOConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub population_size: Int
  pub inertia_weight: Double
  pub cognitive_coeff: Double
  pub social_coeff: Double
}

pub fn pso(config: PSOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let w = config.inertia_weight
  let c1 = config.cognitive_coeff
  let c2 = config.social_coeff

  let mut current_rng = rng

  let positions: Array[Array[Double]] = []
  let velocities: Array[Array[Double]] = []
  let pbest_pos: Array[Array[Double]] = []
  let pbest_val: Array[Double] = []

  for i in 0..<pop_size {
    let (rng1, pos) = random_solution(dim, bounds, current_rng)
    current_rng = rng1

    let vel: Array[Double] = []
    for j in 0..<dim {
      let (rng2, r) = RNG::next(current_rng)
      current_rng = rng2
      vel.push((r * 2.0 - 1.0) * (bounds.max - bounds.min) * 0.1)
    }

    positions.push(pos)
    velocities.push(vel)
    pbest_pos.push(pos)
    let val = evaluate(f, pos)
    pbest_val.push(val)
  }

  let mut gbest_val = pbest_val[0]
  let mut gbest_pos = pbest_pos[0]
  let mut evaluations = pop_size

  for i in 1..<pop_size {
    if pbest_val[i] < gbest_val {
      gbest_val = pbest_val[i]
      gbest_pos = pbest_pos[i]
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals { break }

      let new_vel: Array[Double] = []
      let new_pos: Array[Double] = []

      for j in 0..<dim {
        let (rng3, r1) = RNG::next(current_rng)
        current_rng = rng3
        let (rng4, r2) = RNG::next(current_rng)
        current_rng = rng4

        let vel = velocities[i][j]
        let pos = positions[i][j]
        let pb = pbest_pos[i][j]
        let gb = gbest_pos[j]

        let updated_v = w * vel + c1 * r1 * (pb - pos) + c2 * r2 * (gb - pos)
        let clamped_v = max(bounds.min - pos, min(bounds.max - pos, updated_v))
        new_vel.push(clamped_v)

        let updated_p = max(bounds.min, min(bounds.max, pos + clamped_v))
        new_pos.push(updated_p)
      }

      positions[i] = new_pos
      velocities[i] = new_vel

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < pbest_val[i] {
        pbest_val[i] = new_val
        pbest_pos[i] = new_pos

        if new_val < gbest_val {
          gbest_val = new_val
          gbest_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: gbest_pos,
    best_value: gbest_val,
    evaluations: evaluations
  }
}

// ===== 差分进化算法 ===== 

pub struct DEConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub population_size: Int
  pub mutation_factor: Double
  pub crossover_probability: Double
}

pub fn differential_evolution(config: DEConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let f_mut = config.mutation_factor
  let cr = config.crossover_probability

  let mut current_rng = rng
  let population: Array[Array[Double]] = []

  for i in 0..<pop_size {
    let (rng1, pos) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(pos)
  }

  let fitness: Array[Double] = []
  let mut evaluations = 0
  let mut best_val = 1.0e308
  let mut best_pos: Array[Double] = []

  for i in 0..<pop_size {
    if evaluations >= max_evals { break }
    let val = evaluate(f, population[i])
    fitness.push(val)
    evaluations = evaluations + 1
    if val < best_val {
      best_val = val
      best_pos = population[i]
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals { break }

      let (rng_a, r1) = RNG::next(current_rng)
      current_rng = rng_a
      let mut a_idx = (r1 * pop_size.to_double()).to_int()
      if a_idx == i { a_idx = (a_idx + 1) % pop_size }

      let (rng_b, r2) = RNG::next(current_rng)
      current_rng = rng_b
      let mut b_idx = (r2 * pop_size.to_double()).to_int()
      if b_idx == i || b_idx == a_idx { b_idx = (b_idx + 1) % pop_size; if b_idx == i || b_idx == a_idx { b_idx = (b_idx + 1) % pop_size } }

      let (rng_c, r3) = RNG::next(current_rng)
      current_rng = rng_c
      let mut c_idx = (r3 * pop_size.to_double()).to_int()
      if c_idx == i || c_idx == a_idx || c_idx == b_idx { c_idx = (c_idx + 1) % pop_size; if c_idx == i || c_idx == a_idx || c_idx == b_idx { c_idx = (c_idx + 1) % pop_size } }

      let a = population[a_idx]
      let b = population[b_idx]
      let c = population[c_idx]

      let target = population[i]
      let trial: Array[Double] = []
      let (rng_jr, jr) = RNG::next(current_rng)
      current_rng = rng_jr
      let jrand = (jr * dim.to_double()).to_int() % dim

      for j in 0..<dim {
        let (rng_j, rj) = RNG::next(current_rng)
        current_rng = rng_j
        let use_mutant = rj < cr || j == jrand
        if use_mutant {
          let mutant = a[j] + f_mut * (b[j] - c[j])
          trial.push(max(bounds.min, min(bounds.max, mutant)))
        } else {
          trial.push(target[j])
        }
      }

      let trial_val = evaluate(f, trial)
      evaluations = evaluations + 1

      if trial_val < fitness[i] {
        population[i] = trial
        fitness[i] = trial_val
        if trial_val < best_val {
          best_val = trial_val
          best_pos = trial
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== 遗传算法 ===== 

pub struct GAConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub population_size: Int
  pub crossover_rate: Double
  pub mutation_rate: Double
}

pub fn genetic_algorithm(config: GAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let crossover_rate = config.crossover_rate
  let mutation_rate = config.mutation_rate

  let mut current_rng = rng
  let mut population: Array[Array[Double]] = []
  let mut fitness: Array[Double] = []

  for i in 0..<pop_size {
    let (rng1, pos) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(pos)
  }

  let mut evaluations = 0
  let mut best_val = 1.0e308
  let mut best_pos: Array[Double] = []

  for i in 0..<pop_size {
    if evaluations >= max_evals { break }
    let val = evaluate(f, population[i])
    fitness.push(val)
    evaluations = evaluations + 1
    if val < best_val {
      best_val = val
      best_pos = population[i]
    }
  }

  while evaluations < max_evals {
    let new_population: Array[Array[Double]] = []

    for i in 0..<pop_size {
      let (rng_t1, r_t1) = RNG::next(current_rng)
      current_rng = rng_t1
      let idx1 = (r_t1 * pop_size.to_double()).to_int()
      let (rng_t2, r_t2) = RNG::next(current_rng)
      current_rng = rng_t2
      let idx2 = (r_t2 * pop_size.to_double()).to_int()

      let mut parent: Array[Double] = []
      if fitness[idx1] < fitness[idx2] {
        parent = population[idx1]
      } else {
        parent = population[idx2]
      }

      let (rng_c, r_c) = RNG::next(current_rng)
      current_rng = rng_c
      if r_c < crossover_rate {
        let (rng_p2, r_p2) = RNG::next(current_rng)
        current_rng = rng_p2
        let idx3 = (r_p2 * pop_size.to_double()).to_int()
        let parent2 = population[idx3]

        let (rng_cp, r_cp) = RNG::next(current_rng)
        current_rng = rng_cp
        let crossover_point = (r_cp * dim.to_double()).to_int() % dim

        let child: Array[Double] = []
        for j in 0..<dim {
          if j < crossover_point {
            child.push(parent[j])
          } else {
            child.push(parent2[j])
          }
        }
        parent = child
      }

      let (rng_m, r_m) = RNG::next(current_rng)
      current_rng = rng_m
      if r_m < mutation_rate {
        let (rng_md, r_md) = RNG::next(current_rng)
        current_rng = rng_md
        let mut_idx = (r_md * dim.to_double()).to_int() % dim
        let (rng_mv, r_mv) = RNG::next(current_rng)
        current_rng = rng_mv
        parent[mut_idx] = max(bounds.min, min(bounds.max, parent[mut_idx] + (r_mv * 2.0 - 1.0) * (bounds.max - bounds.min) * 0.1))
      }

      new_population.push(parent)
    }

    let new_fitness: Array[Double] = []
    for i in 0..<pop_size {
      if evaluations >= max_evals { break }
      let val = evaluate(f, new_population[i])
      new_fitness.push(val)
      evaluations = evaluations + 1
      if val < best_val {
        best_val = val
        best_pos = new_population[i]
      }
    }

    let mut best_old_idx = 0
    for i in 1..<pop_size {
      if fitness[i] < fitness[best_old_idx] {
        best_old_idx = i
      }
    }

    let mut worst_new_idx = 0
    for i in 1..<new_fitness.length() {
      if new_fitness[i] > new_fitness[worst_new_idx] {
        worst_new_idx = i
      }
    }

    if fitness[best_old_idx] < new_fitness[worst_new_idx] {
      new_population[worst_new_idx] = population[best_old_idx]
      new_fitness[worst_new_idx] = fitness[best_old_idx]
    }

    population = new_population
    fitness = new_fitness
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Nelder-Mead 单纯形法 ===== 

pub struct NMConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub reflection_coeff: Double
  pub expansion_coeff: Double
  pub contraction_coeff: Double
  pub shrink_coeff: Double
}

pub fn nelder_mead(config: NMConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let alpha = config.reflection_coeff
  let gamma = config.expansion_coeff
  let rho = config.contraction_coeff
  let sigma = config.shrink_coeff

  let mut current_rng = rng

  let simplex_size = dim + 1
  let simplex: Array[Array[Double]] = []
  let simplex_val: Array[Double] = []

  let (rng1, first) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  simplex.push(first)
  let first_val = evaluate(f, first)
  simplex_val.push(first_val)

  for i in 1..<simplex_size {
    let vertex: Array[Double] = []
    for j in 0..<dim {
      let (rng_v, r_v) = RNG::next(current_rng)
      current_rng = rng_v
      if j == i - 1 {
        vertex.push(max(bounds.min, min(bounds.max, first[j] + (r_v * 2.0 - 1.0) * (bounds.max - bounds.min) * 0.05)))
      } else {
        vertex.push(first[j])
      }
    }
    simplex.push(vertex)
    let val = evaluate(f, vertex)
    simplex_val.push(val)
  }

  let mut evaluations = simplex_size
  let mut best_val = simplex_val[0]
  let mut best_pos = simplex[0]

  for i in 1..<simplex_size {
    if simplex_val[i] < best_val {
      best_val = simplex_val[i]
      best_pos = simplex[i]
    }
  }

  while evaluations < max_evals {
    let mut best_idx = 0
    let mut worst_idx = 0
    for i in 1..<simplex_size {
      if simplex_val[i] < simplex_val[best_idx] { best_idx = i }
      if simplex_val[i] > simplex_val[worst_idx] { worst_idx = i }
    }

    if simplex_val[best_idx] < best_val {
      best_val = simplex_val[best_idx]
      best_pos = simplex[best_idx]
    }

    let centroid: Array[Double] = []
    for j in 0..<dim {
      let mut sum = 0.0
      for i in 0..<simplex_size {
        if i != worst_idx {
          sum = sum + simplex[i][j]
        }
      }
      centroid.push(sum / (simplex_size.to_double() - 1.0))
    }

    let reflected: Array[Double] = []
    for j in 0..<dim {
      let r = centroid[j] + alpha * (centroid[j] - simplex[worst_idx][j])
      reflected.push(max(bounds.min, min(bounds.max, r)))
    }
    let reflected_val = evaluate(f, reflected)
    evaluations = evaluations + 1

    if reflected_val < simplex_val[best_idx] {
      let expanded: Array[Double] = []
      for j in 0..<dim {
        let e = centroid[j] + gamma * (reflected[j] - centroid[j])
        expanded.push(max(bounds.min, min(bounds.max, e)))
      }
      let expanded_val = evaluate(f, expanded)
      evaluations = evaluations + 1

      if expanded_val < reflected_val {
        simplex[worst_idx] = expanded
        simplex_val[worst_idx] = expanded_val
      } else {
        simplex[worst_idx] = reflected
        simplex_val[worst_idx] = reflected_val
      }
    } else {
      let mut second_worst_idx = 0
      for i in 1..<simplex_size {
        if i != worst_idx {
          if second_worst_idx == worst_idx || simplex_val[i] > simplex_val[second_worst_idx] {
            second_worst_idx = i
          }
        }
      }

      if reflected_val < simplex_val[second_worst_idx] {
        simplex[worst_idx] = reflected
        simplex_val[worst_idx] = reflected_val
      } else {
        if reflected_val < simplex_val[worst_idx] {
          let contracted: Array[Double] = []
          for j in 0..<dim {
            let c = centroid[j] + rho * (reflected[j] - centroid[j])
            contracted.push(max(bounds.min, min(bounds.max, c)))
          }
          let contracted_val = evaluate(f, contracted)
          evaluations = evaluations + 1

          if contracted_val <= reflected_val {
            simplex[worst_idx] = contracted
            simplex_val[worst_idx] = contracted_val
          } else {
            for i in 0..<simplex_size {
              if i != best_idx {
                for j in 0..<dim {
                  simplex[i][j] = max(bounds.min, min(bounds.max, simplex[best_idx][j] + sigma * (simplex[i][j] - simplex[best_idx][j])))
                }
                simplex_val[i] = evaluate(f, simplex[i])
                evaluations = evaluations + 1
              }
            }
          }
        } else {
          let contracted: Array[Double] = []
          for j in 0..<dim {
            let c = centroid[j] + rho * (simplex[worst_idx][j] - centroid[j])
            contracted.push(max(bounds.min, min(bounds.max, c)))
          }
          let contracted_val = evaluate(f, contracted)
          evaluations = evaluations + 1

          if contracted_val < simplex_val[worst_idx] {
            simplex[worst_idx] = contracted
            simplex_val[worst_idx] = contracted_val
          } else {
            for i in 0..<simplex_size {
              if i != best_idx {
                for j in 0..<dim {
                  simplex[i][j] = max(bounds.min, min(bounds.max, simplex[best_idx][j] + sigma * (simplex[i][j] - simplex[best_idx][j])))
                }
                simplex_val[i] = evaluate(f, simplex[i])
                evaluations = evaluations + 1
              }
            }
          }
        }
      }
    }
  }

  for i in 0..<simplex_size {
    if simplex_val[i] < best_val {
      best_val = simplex_val[i]
      best_pos = simplex[i]
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== 萤火虫算法 ======

pub struct FireflyConfig {
  pub max_evaluations: Int
  pub dimension: Int
  pub bounds: Bounds
  pub population_size: Int
  pub initial_attractiveness: Double
  pub light_absorption: Double
  pub step_size: Double
}

pub fn firefly_algorithm(config: FireflyConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let beta0 = config.initial_attractiveness
  let gamma = config.light_absorption
  let alpha = config.step_size

  let mut current_rng = rng
  let population: Array[Array[Double]] = []
  let brightness: Array[Double] = []

  for i in 0..<pop_size {
    let (rng1, pos) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(pos)
  }

  let mut evaluations = 0
  let mut best_val = 1.0e308
  let mut best_pos: Array[Double] = []

  for i in 0..<pop_size {
    if evaluations >= max_evals { break }
    let val = evaluate(f, population[i])
    brightness.push(val)
    evaluations = evaluations + 1
    if val < best_val {
      best_val = val
      best_pos = population[i]
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals { break }

      for j in 0..<pop_size {
        if brightness[j] < brightness[i] {
          let mut dist_sq = 0.0
          for d in 0..<dim {
            let diff = population[j][d] - population[i][d]
            dist_sq = dist_sq + diff * diff
          }

          let beta = beta0 * exp_approx(-gamma * dist_sq)

          for d in 0..<dim {
            let (rng_f, r_f) = RNG::next(current_rng)
            current_rng = rng_f
            let movement = beta * (population[j][d] - population[i][d]) + alpha * (r_f * 2.0 - 1.0) * (bounds.max - bounds.min) * 0.1
            population[i][d] = max(bounds.min, min(bounds.max, population[i][d] + movement))
          }
        }
      }

      let new_val = evaluate(f, population[i])
      evaluations = evaluations + 1
      brightness[i] = new_val

      if new_val < best_val {
        best_val = new_val
        best_pos = population[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== 测试函数 ===== 

pub fn sphere_function(x: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<x.length() {
    sum = sum + x[i] * x[i]
  }
  sum
}

pub fn rastrigin_function(x: Array[Double]) -> Double {
  let n: Double = x.length().to_double()
  let mut sum = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum = sum + xi * xi - 10.0 * cos_approx(2.0 * PI * xi)
  }
  sum + 10.0 * n
}

pub fn rosenbrock_function(x: Array[Double]) -> Double {
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

pub fn calc_relative_error(best_value: Double, known_optimal: Double) -> Double {
  let diff = abs_val(best_value - known_optimal)
  if known_optimal == 0.0 {
    if best_value == 0.0 { 0.0 } else { diff }
  } else {
    diff / abs_val(known_optimal)
  }
}

pub fn array_sum(arr: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<arr.length() {
    sum = sum + arr[i]
  }
  sum
}

pub fn array_min(arr: Array[Double]) -> Double {
  let mut min_val = arr[0]
  for i in 1..<arr.length() {
    if arr[i] < min_val {
      min_val = arr[i]
    }
  }
  min_val
}

pub fn array_max(arr: Array[Double]) -> Double {
  let mut max_val = arr[0]
  for i in 1..<arr.length() {
    if arr[i] > max_val {
      max_val = arr[i]
    }
  }
  max_val
}

pub fn array_avg(arr: Array[Double]) -> Double {
  array_sum(arr) / arr.length().to_double()
}

pub fn run_single_benchmark(
  algorithm: String,
  f: (Array[Double]) -> Double,
  dim: Int,
  bounds: Bounds,
  known_optimal: Double,
  max_evals: Int,
  rng: RNG
) -> BenchmarkResult {
  if algorithm == "Random Search" {
    let config = make_random_search_config(max_evals, dim, bounds)
    let result = random_search(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "Hill Climbing" {
    let config = make_hill_climbing_config(max_evals, dim, bounds, 0.5, 10)
    let result = hill_climbing(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "Simulated Annealing" {
    let config = make_simulated_annealing_config(max_evals, dim, bounds, 100.0, 0.995, 0.5)
    let result = simulated_annealing(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "PSO" {
    let config = make_pso_config(max_evals, dim, bounds, 30, 0.7, 1.5, 1.5)
    let result = pso(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "DE" {
    let config = make_de_config(max_evals, dim, bounds, 30, 0.8, 0.9)
    let result = differential_evolution(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "GA" {
    let config = make_ga_config(max_evals, dim, bounds, 50, 0.8, 0.01)
    let result = genetic_algorithm(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "Nelder-Mead" {
    let config = make_nm_config(max_evals, dim, bounds, 1.0, 2.0, 0.5, 0.5)
    let result = nelder_mead(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  if algorithm == "Firefly" {
    let config = make_firefly_config(max_evals, dim, bounds, 30, 1.0, 1.0, 0.5)
    let result = firefly_algorithm(config, f, rng)
    let error = calc_relative_error(result.best_value, known_optimal)
    return BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
  }
  let config = make_random_search_config(max_evals, dim, bounds)
  let result = random_search(config, f, rng)
  let error = calc_relative_error(result.best_value, known_optimal)
  BenchmarkResult::{ algorithm_name: algorithm, best_value: result.best_value, evaluations: result.evaluations, relative_error: error }
}

pub fn run_benchmark_at_evals(
  f: (Array[Double]) -> Double,
  func_name: String,
  dim: Int,
  bounds: Bounds,
  known_optimal: Double,
  max_evals: Int,
  rng: RNG
) -> Array[BenchmarkResult] {
  let results: Array[BenchmarkResult] = []

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

  let pso_config = PSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    population_size: 30,
    inertia_weight: 0.7,
    cognitive_coeff: 1.5,
    social_coeff: 1.5
  }
  let pso_result = pso(pso_config, f, rng)
  let pso_error = calc_relative_error(pso_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "PSO",
    best_value: pso_result.best_value,
    evaluations: pso_result.evaluations,
    relative_error: pso_error
  })

  let de_config = DEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    population_size: 30,
    mutation_factor: 0.8,
    crossover_probability: 0.9
  }
  let de_result = differential_evolution(de_config, f, rng)
  let de_error = calc_relative_error(de_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "DE",
    best_value: de_result.best_value,
    evaluations: de_result.evaluations,
    relative_error: de_error
  })

  let ga_config = GAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    population_size: 50,
    crossover_rate: 0.8,
    mutation_rate: 0.01
  }
  let ga_result = genetic_algorithm(ga_config, f, rng)
  let ga_error = calc_relative_error(ga_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "GA",
    best_value: ga_result.best_value,
    evaluations: ga_result.evaluations,
    relative_error: ga_error
  })

  let nm_config = NMConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    reflection_coeff: 1.0,
    expansion_coeff: 2.0,
    contraction_coeff: 0.5,
    shrink_coeff: 0.5
  }
  let nm_result = nelder_mead(nm_config, f, rng)
  let nm_error = calc_relative_error(nm_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Nelder-Mead",
    best_value: nm_result.best_value,
    evaluations: nm_result.evaluations,
    relative_error: nm_error
  })

  let ff_config = FireflyConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: bounds,
    population_size: 30,
    initial_attractiveness: 1.0,
    light_absorption: 1.0,
    step_size: 0.5
  }
  let ff_result = firefly_algorithm(ff_config, f, rng)
  let ff_error = calc_relative_error(ff_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Firefly",
    best_value: ff_result.best_value,
    evaluations: ff_result.evaluations,
    relative_error: ff_error
  })

  results
}

pub fn print_benchmark_table(results: Array[BenchmarkResult], func_name: String, max_evals: Int) -> Unit {
  println("  \{func_name} (evals = \{max_evals})")
  println("  +--------------------+-------------------+------------+-------------------+")
  println("  | Algorithm          | Best Value        | Evals      | Relative Error    |")
  println("  +--------------------+-------------------+------------+-------------------+")
  for i in 0..<results.length() {
    let r = results[i]
    let err_str = "\{r.relative_error}"
    println("  | \{r.algorithm_name} | \{r.best_value} | \{r.evaluations} | \{err_str} |")
  }
  println("  +--------------------+-------------------+------------+-------------------+")
  println("")
}

pub fn run_multiple_trials(
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
  let pso_values: Array[Double] = []
  let de_values: Array[Double] = []
  let ga_values: Array[Double] = []
  let nm_values: Array[Double] = []
  let ff_values: Array[Double] = []

  for trial in 0..<num_trials {
    let seed = base_seed + trial.to_double() * 1000.0
    let rng = RNG::{ seed: seed }

    let rs_config = RandomSearchConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds
    }
    let rs_result = random_search(rs_config, f, rng)
    rs_values.push(rs_result.best_value)

    let hc_config = HillClimbingConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      step_size: 0.5,
      restarts: 10
    }
    let hc_result = hill_climbing(hc_config, f, rng)
    hc_values.push(hc_result.best_value)

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

    let pso_config = PSOConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      population_size: 30,
      inertia_weight: 0.7,
      cognitive_coeff: 1.5,
      social_coeff: 1.5
    }
    let pso_result = pso(pso_config, f, rng)
    pso_values.push(pso_result.best_value)

    let de_config = DEConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      population_size: 30,
      mutation_factor: 0.8,
      crossover_probability: 0.9
    }
    let de_result = differential_evolution(de_config, f, rng)
    de_values.push(de_result.best_value)

    let ga_config = GAConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      population_size: 50,
      crossover_rate: 0.8,
      mutation_rate: 0.01
    }
    let ga_result = genetic_algorithm(ga_config, f, rng)
    ga_values.push(ga_result.best_value)

    let nm_config = NMConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      reflection_coeff: 1.0,
      expansion_coeff: 2.0,
      contraction_coeff: 0.5,
      shrink_coeff: 0.5
    }
    let nm_result = nelder_mead(nm_config, f, rng)
    nm_values.push(nm_result.best_value)

    let ff_config = FireflyConfig::{
      max_evaluations: max_evals,
      dimension: dim,
      bounds: bounds,
      population_size: 30,
      initial_attractiveness: 1.0,
      light_absorption: 1.0,
      step_size: 0.5
    }
    let ff_result = firefly_algorithm(ff_config, f, rng)
    ff_values.push(ff_result.best_value)
  }

  println("  Algorithm          | Min Value           | Max Value           | Avg Value           | Avg Error")
  println("  +-------------------+---------------------+---------------------+---------------------+------------------+")

  let rs_min = array_min(rs_values)
  let rs_max = array_max(rs_values)
  let rs_avg = array_avg(rs_values)
  let rs_err = calc_relative_error(rs_avg, known_optimal)
  println("  | Random Search     | \{rs_min} | \{rs_max} | \{rs_avg} | \{rs_err} |")

  let hc_min = array_min(hc_values)
  let hc_max = array_max(hc_values)
  let hc_avg = array_avg(hc_values)
  let hc_err = calc_relative_error(hc_avg, known_optimal)
  println("  | Hill Climbing     | \{hc_min} | \{hc_max} | \{hc_avg} | \{hc_err} |")

  let sa_min = array_min(sa_values)
  let sa_max = array_max(sa_values)
  let sa_avg = array_avg(sa_values)
  let sa_err = calc_relative_error(sa_avg, known_optimal)
  println("  | Simulated Annealing | \{sa_min} | \{sa_max} | \{sa_avg} | \{sa_err} |")

  let pso_min = array_min(pso_values)
  let pso_max = array_max(pso_values)
  let pso_avg = array_avg(pso_values)
  let pso_err = calc_relative_error(pso_avg, known_optimal)
  println("  | PSO               | \{pso_min} | \{pso_max} | \{pso_avg} | \{pso_err} |")

  let de_min = array_min(de_values)
  let de_max = array_max(de_values)
  let de_avg = array_avg(de_values)
  let de_err = calc_relative_error(de_avg, known_optimal)
  println("  | DE                | \{de_min} | \{de_max} | \{de_avg} | \{de_err} |")

  let ga_min = array_min(ga_values)
  let ga_max = array_max(ga_values)
  let ga_avg = array_avg(ga_values)
  let ga_err = calc_relative_error(ga_avg, known_optimal)
  println("  | GA                | \{ga_min} | \{ga_max} | \{ga_avg} | \{ga_err} |")

  let nm_min = array_min(nm_values)
  let nm_max = array_max(nm_values)
  let nm_avg = array_avg(nm_values)
  let nm_err = calc_relative_error(nm_avg, known_optimal)
  println("  | Nelder-Mead       | \{nm_min} | \{nm_max} | \{nm_avg} | \{nm_err} |")

  let ff_min = array_min(ff_values)
  let ff_max = array_max(ff_values)
  let ff_avg = array_avg(ff_values)
  let ff_err = calc_relative_error(ff_avg, known_optimal)
  println("  | Firefly           | \{ff_min} | \{ff_max} | \{ff_avg} | \{ff_err} |")

  println("")
}

pub fn run_evals_comparison(
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

