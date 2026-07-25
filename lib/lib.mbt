// Black-Box Optimization Framework - 核心库

// ===== 数学常量 =====

pub const PI: Double = 3.14159265358979323846
pub const TWO_PI: Double = 6.28318530717958647692

// ===== 工具函数 =====

pub struct RNG {
  seed: Double
}

pub fn RNG::next(self: RNG) -> (RNG, Double) {
  let new_seed = (self.seed * 9301.0 + 49297.0) % 233280.0
  ({ seed: new_seed }, new_seed / 233280.0)
}

pub fn make_rng(seed: Double) -> RNG {
  RNG::{ seed: seed }
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

pub fn sin_approx(x: Double) -> Double {
  let rx = reduce_to_pi(x)
  let x2 = rx * rx
  let x3 = x2 * rx
  let x5 = x3 * x2
  let x7 = x5 * x2
  let x9 = x7 * x2
  rx - x3 / 6.0 + x5 / 120.0 - x7 / 5040.0 + x9 / 362880.0
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

pub fn pow_approx(base: Double, exp: Double) -> Double {
  if base <= 0.0 {
    return 0.0001
  }
  exp_approx(exp * ln_approx(base))
}

pub fn max(a: Double, b: Double) -> Double {
  if a > b { a } else { b }
}

pub fn min(a: Double, b: Double) -> Double {
  if a < b { a } else { b }
}

fn max_int(a: Int, b: Int) -> Int {
  if a > b { a } else { b }
}

fn min_int(a: Int, b: Int) -> Int {
  if a < b { a } else { b }
}

pub fn abs_val(x: Double) -> Double {
  if x < 0.0 { -x } else { x }
}

// ===== 核心类型定义 =====

pub struct Bounds {
  min: Double
  max: Double
}

pub fn make_bounds(min: Double, max: Double) -> Bounds {
  Bounds::{ min: min, max: max }
}

pub struct OptimizationResult {
  best_solution: Array[Double]
  best_value: Double
  evaluations: Int
}

pub struct BenchmarkResult {
  algorithm_name: String
  best_value: Double
  evaluations: Int
  relative_error: Double
}

pub struct Stats {
  min_value: Double
  max_value: Double
  avg_value: Double
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

// ===== 随机搜索算法 =====

pub struct RandomSearchConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
}

pub fn make_random_search_config(max_evaluations: Int, dimension: Int, bounds: Bounds) -> RandomSearchConfig {
  RandomSearchConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds }
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

// ===== 爬山法算法 =====

pub struct HillClimbingConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  step_size: Double
  restarts: Int
}

pub fn make_hill_climbing_config(max_evaluations: Int, dimension: Int, bounds: Bounds, step_size: Double, restarts: Int) -> HillClimbingConfig {
  HillClimbingConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, step_size: step_size, restarts: restarts }
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

// ===== 模拟退火算法 =====

pub struct SimulatedAnnealingConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  initial_temperature: Double
  cooling_rate: Double
  step_size: Double
}

pub fn make_simulated_annealing_config(max_evaluations: Int, dimension: Int, bounds: Bounds, initial_temperature: Double, cooling_rate: Double, step_size: Double) -> SimulatedAnnealingConfig {
  SimulatedAnnealingConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, initial_temperature: initial_temperature, cooling_rate: cooling_rate, step_size: step_size }
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

// ===== 粒子群优化算法 (PSO) =====

pub struct PSOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  swarm_size: Int
  inertia_weight: Double
  cognitive_coeff: Double
  social_coeff: Double
}

pub fn make_pso_config(max_evaluations: Int, dimension: Int, bounds: Bounds, swarm_size: Int, inertia_weight: Double, cognitive_coeff: Double, social_coeff: Double) -> PSOConfig {
  PSOConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, swarm_size: swarm_size, inertia_weight: inertia_weight, cognitive_coeff: cognitive_coeff, social_coeff: social_coeff }
}

pub fn pso(config: PSOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let swarm_size = config.swarm_size
  let w = config.inertia_weight
  let c1 = config.cognitive_coeff
  let c2 = config.social_coeff

  let mut current_rng = rng

  let positions: Array[Array[Double]] = []
  let velocities: Array[Array[Double]] = []
  let personal_best: Array[Array[Double]] = []
  let personal_best_values: Array[Double] = []

  for i in 0..<swarm_size {
    let (rng1, pos) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    positions.push(pos)

    let vel: Array[Double] = []
    for j in 0..<dim {
      let (rng2, r) = RNG::next(current_rng)
      current_rng = rng2
      vel.push((r * 2.0 - 1.0) * (bounds.max - bounds.min) * 0.1)
    }
    velocities.push(vel)

    personal_best.push(pos)
    personal_best_values.push(evaluate(f, pos))
  }

  let mut global_best = personal_best[0]
  let mut global_best_value = personal_best_values[0]
  for i in 1..<swarm_size {
    if personal_best_values[i] < global_best_value {
      global_best = personal_best[i]
      global_best_value = personal_best_values[i]
    }
  }

  let mut evaluations = swarm_size

  while evaluations < max_evals {
    for i in 0..<swarm_size {
      if evaluations >= max_evals { break }

      let new_vel: Array[Double] = []
      let new_pos: Array[Double] = []

      for j in 0..<dim {
        let (rng1, r1) = RNG::next(current_rng)
        current_rng = rng1
        let (rng2, r2) = RNG::next(current_rng)
        current_rng = rng2

        let vel = w * velocities[i][j] + c1 * r1 * (personal_best[i][j] - positions[i][j]) + c2 * r2 * (global_best[j] - positions[i][j])
        new_vel.push(vel)

        let pos = positions[i][j] + vel
        let clamped = max(bounds.min, min(bounds.max, pos))
        new_pos.push(clamped)
      }

      velocities[i] = new_vel
      positions[i] = new_pos

      let value = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if value < personal_best_values[i] {
        personal_best_values[i] = value
        personal_best[i] = new_pos

        if value < global_best_value {
          global_best_value = value
          global_best = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: global_best,
    best_value: global_best_value,
    evaluations: evaluations
  }
}

// ===== 差分进化算法 (DE) =====

pub struct DEConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  mutation_factor: Double
  crossover_rate: Double
}

pub fn make_de_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, mutation_factor: Double, crossover_rate: Double) -> DEConfig {
  DEConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, mutation_factor: mutation_factor, crossover_rate: crossover_rate }
}

pub fn differential_evolution(config: DEConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let f_mut = config.mutation_factor
  let cr = config.crossover_rate

  let mut current_rng = rng

  let population: Array[Array[Double]] = []
  let fitness: Array[Double] = []

  for i in 0..<pop_size {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(sol)
    fitness.push(evaluate(f, sol))
  }

  let mut best_idx = 0
  let mut best_value = fitness[0]
  for i in 1..<pop_size {
    if fitness[i] < best_value {
      best_idx = i
      best_value = fitness[i]
    }
  }

  let mut evaluations = pop_size

  while evaluations < max_evals {
    for target in 0..<pop_size {
      if evaluations >= max_evals { break }

      let mut a = 0
      let mut b = 0
      let mut c = 0

      let (rng1, ra) = RNG::next(current_rng)
      current_rng = rng1
      a = (ra * pop_size.to_double()).to_int() % pop_size
      if a == target { a = (a + 1) % pop_size }

      let (rng2, rb) = RNG::next(current_rng)
      current_rng = rng2
      b = (rb * pop_size.to_double()).to_int() % pop_size
      while b == target || b == a {
        let (rng3, rb2) = RNG::next(current_rng)
        current_rng = rng3
        b = (rb2 * pop_size.to_double()).to_int() % pop_size
      }

      let (rng4, rc) = RNG::next(current_rng)
      current_rng = rng4
      c = (rc * pop_size.to_double()).to_int() % pop_size
      while c == target || c == a || c == b {
        let (rng5, rc2) = RNG::next(current_rng)
        current_rng = rng5
        c = (rc2 * pop_size.to_double()).to_int() % pop_size
      }

      let trial: Array[Double] = []
      for j in 0..<dim {
        let (rng6, r_cross) = RNG::next(current_rng)
        current_rng = rng6
        let mutant = population[a][j] + f_mut * (population[b][j] - population[c][j])
        let clamped_mutant = max(bounds.min, min(bounds.max, mutant))
        if r_cross < cr {
          trial.push(clamped_mutant)
        } else {
          trial.push(population[target][j])
        }
      }

      let trial_fitness = evaluate(f, trial)
      evaluations = evaluations + 1

      if trial_fitness < fitness[target] {
        population[target] = trial
        fitness[target] = trial_fitness

        if trial_fitness < best_value {
          best_value = trial_fitness
          best_idx = target
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: population[best_idx],
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 遗传算法 (GA) =====

pub struct GAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  mutation_rate: Double
  crossover_rate: Double
}

pub fn make_ga_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, mutation_rate: Double, crossover_rate: Double) -> GAConfig {
  GAConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, mutation_rate: mutation_rate, crossover_rate: crossover_rate }
}

fn tournament_select(population: Array[Array[Double]], fitness: Array[Double], rng: RNG) -> (RNG, Array[Double]) {
  let (rng1, r1) = RNG::next(rng)
  let (rng2, r2) = RNG::next(rng1)
  let (rng3, r3) = RNG::next(rng2)

  let idx1 = (r1 * population.length().to_double()).to_int() % population.length()
  let idx2 = (r2 * population.length().to_double()).to_int() % population.length()
  let idx3 = (r3 * population.length().to_double()).to_int() % population.length()

  let best = min(fitness[idx1], min(fitness[idx2], fitness[idx3]))
  if best == fitness[idx1] {
    (rng3, population[idx1])
  } else if best == fitness[idx2] {
    (rng3, population[idx2])
  } else {
    (rng3, population[idx3])
  }
}

pub fn genetic_algorithm(config: GAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let mutation_rate = config.mutation_rate
  let crossover_rate = config.crossover_rate

  let mut current_rng = rng

  let mut population: Array[Array[Double]] = []
  let mut fitness: Array[Double] = []

  for i in 0..<pop_size {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(sol)
    fitness.push(evaluate(f, sol))
  }

  let mut best_idx = 0
  let mut best_value = fitness[0]
  for i in 1..<pop_size {
    if fitness[i] < best_value {
      best_idx = i
      best_value = fitness[i]
    }
  }

  let mut evaluations = pop_size

  while evaluations < max_evals {
    let new_population: Array[Array[Double]] = []

    for i in 0..<pop_size {
      if evaluations >= max_evals {
        new_population.push(population[i])
        continue
      }

      let (rng1, parent1) = tournament_select(population, fitness, current_rng)
      current_rng = rng1
      let (rng2, parent2) = tournament_select(population, fitness, current_rng)
      current_rng = rng2

      let (rng3, r_cross) = RNG::next(current_rng)
      current_rng = rng3

      let child: Array[Double] = []
      if r_cross < crossover_rate {
        let (rng4, r_point) = RNG::next(current_rng)
        current_rng = rng4
        let point = (r_point * dim.to_double()).to_int() % dim
        for j in 0..<dim {
          if j < point {
            child.push(parent1[j])
          } else {
            child.push(parent2[j])
          }
        }
      } else {
        for j in 0..<dim {
          child.push(parent1[j])
        }
      }

      for j in 0..<dim {
        let (rng5, r_mut) = RNG::next(current_rng)
        current_rng = rng5
        if r_mut < mutation_rate {
          let (rng6, r_val) = RNG::next(current_rng)
          current_rng = rng6
          child[j] = bounds.min + r_val * (bounds.max - bounds.min)
        }
      }

      let child_fitness = evaluate(f, child)
      evaluations = evaluations + 1
      new_population.push(child)

      if child_fitness < best_value {
        best_value = child_fitness
        best_idx = new_population.length() - 1
      }
    }

    population = new_population
    let new_fitness: Array[Double] = []
    for i in 0..<pop_size {
      let val = evaluate(f, population[i])
      new_fitness.push(val)
      if val < best_value {
        best_value = val
        best_idx = i
      }
      evaluations = evaluations + 1
      if evaluations >= max_evals { break }
    }
    fitness = new_fitness
  }

  OptimizationResult::{
    best_solution: population[best_idx],
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== Nelder-Mead 单纯形法 =====

pub struct NMConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  reflection: Double
  expansion: Double
  contraction: Double
  shrink: Double
}

pub fn make_nm_config(max_evaluations: Int, dimension: Int, bounds: Bounds, reflection: Double, expansion: Double, contraction: Double, shrink: Double) -> NMConfig {
  NMConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, reflection: reflection, expansion: expansion, contraction: contraction, shrink: shrink }
}

pub fn nelder_mead(config: NMConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let alpha = config.reflection
  let gamma = config.expansion
  let rho = config.contraction
  let sigma = config.shrink

  let mut current_rng = rng

  let n_vertices = dim + 1
  let simplex: Array[Array[Double]] = []
  let values: Array[Double] = []

  let (rng1, base) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  simplex.push(base)
  values.push(evaluate(f, base))

  for i in 1..<n_vertices {
    let vertex: Array[Double] = []
    for j in 0..<dim {
      let (rng2, r) = RNG::next(current_rng)
      current_rng = rng2
      vertex.push(base[j] + (r * 2.0 - 1.0) * 0.5)
    }
    let clamped: Array[Double] = []
    for j in 0..<dim {
      clamped.push(max(bounds.min, min(bounds.max, vertex[j])))
    }
    simplex.push(clamped)
    values.push(evaluate(f, clamped))
  }

  let mut evaluations = n_vertices

  while evaluations < max_evals {
    let mut best = 0
    let mut worst = 0
    let mut second_worst = 0
    for i in 1..<n_vertices {
      if values[i] < values[best] {
        second_worst = worst
        worst = i
        best = i
      } else if values[i] < values[second_worst] || i == 1 {
        if i != best {
          second_worst = i
        }
      }
      if values[i] > values[worst] {
        worst = i
      }
    }

    if best == worst {
      for i in 0..<n_vertices {
        if values[i] < values[best] { best = i }
        if values[i] > values[worst] { worst = i }
      }
      if second_worst == worst {
        for i in 0..<n_vertices {
          if i != best && i != worst { second_worst = i; break }
        }
      }
    }

    let centroid: Array[Double] = []
    for j in 0..<dim {
      let mut sum = 0.0
      for i in 0..<n_vertices {
        if i != worst {
          sum = sum + simplex[i][j]
        }
      }
      centroid.push(sum / (n_vertices.to_double() - 1.0))
    }

    let reflected: Array[Double] = []
    for j in 0..<dim {
      reflected.push(max(bounds.min, min(bounds.max, centroid[j] + alpha * (centroid[j] - simplex[worst][j]))))
    }
    let reflected_val = evaluate(f, reflected)
    evaluations = evaluations + 1

    if reflected_val < values[best] {
      let expanded: Array[Double] = []
      for j in 0..<dim {
        expanded.push(max(bounds.min, min(bounds.max, centroid[j] + gamma * (reflected[j] - centroid[j]))))
      }
      let expanded_val = evaluate(f, expanded)
      evaluations = evaluations + 1

      if expanded_val < reflected_val {
        simplex[worst] = expanded
        values[worst] = expanded_val
      } else {
        simplex[worst] = reflected
        values[worst] = reflected_val
      }
    } else if reflected_val < values[second_worst] {
      simplex[worst] = reflected
      values[worst] = reflected_val
    } else {
      if reflected_val < values[worst] {
        let contracted: Array[Double] = []
        for j in 0..<dim {
          contracted.push(max(bounds.min, min(bounds.max, centroid[j] + rho * (reflected[j] - centroid[j]))))
        }
        let contracted_val = evaluate(f, contracted)
        evaluations = evaluations + 1

        if contracted_val <= reflected_val {
          simplex[worst] = contracted
          values[worst] = contracted_val
        } else {
          for i in 0..<n_vertices {
            if i != best {
              let shrunk: Array[Double] = []
              for j in 0..<dim {
                shrunk.push(max(bounds.min, min(bounds.max, simplex[best][j] + sigma * (simplex[i][j] - simplex[best][j]))))
              }
              simplex[i] = shrunk
              values[i] = evaluate(f, shrunk)
              evaluations = evaluations + 1
              if evaluations >= max_evals { break }
            }
          }
        }
      } else {
        let contracted: Array[Double] = []
        for j in 0..<dim {
          contracted.push(max(bounds.min, min(bounds.max, centroid[j] + rho * (simplex[worst][j] - centroid[j]))))
        }
        let contracted_val = evaluate(f, contracted)
        evaluations = evaluations + 1

        if contracted_val < values[worst] {
          simplex[worst] = contracted
          values[worst] = contracted_val
        } else {
          for i in 0..<n_vertices {
            if i != best {
              let shrunk: Array[Double] = []
              for j in 0..<dim {
                shrunk.push(max(bounds.min, min(bounds.max, simplex[best][j] + sigma * (simplex[i][j] - simplex[best][j]))))
              }
              simplex[i] = shrunk
              values[i] = evaluate(f, shrunk)
              evaluations = evaluations + 1
              if evaluations >= max_evals { break }
            }
          }
        }
      }
    }
  }

  let mut best = 0
  for i in 1..<n_vertices {
    if values[i] < values[best] { best = i }
  }

  OptimizationResult::{
    best_solution: simplex[best],
    best_value: values[best],
    evaluations: evaluations
  }
}

// ===== 萤火虫算法 (Firefly) =====

pub struct FireflyConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  alpha: Double
  beta0: Double
  gamma: Double
}

pub fn make_firefly_config(max_evaluations: Int, dimension: Int, bounds: Bounds, population_size: Int, alpha: Double, beta0: Double, gamma: Double) -> FireflyConfig {
  FireflyConfig::{ max_evaluations: max_evaluations, dimension: dimension, bounds: bounds, population_size: population_size, alpha: alpha, beta0: beta0, gamma: gamma }
}

fn euclidean_distance(a: Array[Double], b: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<a.length() {
    let diff = a[i] - b[i]
    sum = sum + diff * diff
  }
  sum.sqrt()
}

pub fn firefly_algorithm(config: FireflyConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let alpha = config.alpha
  let beta0 = config.beta0
  let gamma = config.gamma

  let mut current_rng = rng

  let fireflies: Array[Array[Double]] = []
  let brightness: Array[Double] = []

  for i in 0..<pop_size {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    fireflies.push(sol)
    brightness.push(evaluate(f, sol))
  }

  let mut best_idx = 0
  let mut best_value = brightness[0]
  for i in 1..<pop_size {
    if brightness[i] < best_value {
      best_idx = i
      best_value = brightness[i]
    }
  }

  let mut evaluations = pop_size

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals { break }

      for j in 0..<pop_size {
        if brightness[j] < brightness[i] {
          let dist = euclidean_distance(fireflies[i], fireflies[j])
          let beta = beta0 / (1.0 + gamma * dist * dist)

          let new_pos: Array[Double] = []
          for d in 0..<dim {
            let (rng2, r) = RNG::next(current_rng)
            current_rng = rng2
            let step = alpha * (r * 2.0 - 1.0) * (bounds.max - bounds.min)
            let moved = fireflies[i][d] + beta * (fireflies[j][d] - fireflies[i][d]) + step
            new_pos.push(max(bounds.min, min(bounds.max, moved)))
          }

          fireflies[i] = new_pos
          let new_brightness = evaluate(f, new_pos)
          evaluations = evaluations + 1
          brightness[i] = new_brightness

          if new_brightness < best_value {
            best_value = new_brightness
            best_idx = i
          }

          if evaluations >= max_evals { break }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: fireflies[best_idx],
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 蚁群算法 (ACO) =====

pub struct ACOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  ant_count: Int
  evaporation_rate: Double
  pheromone_weight: Double
  heuristic_weight: Double
}

pub fn make_aco_config(max_evals: Int, dim: Int, b: Bounds, ants: Int, evap: Double, pher_w: Double, heur_w: Double) -> ACOConfig {
  ACOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    ant_count: ants,
    evaporation_rate: evap,
    pheromone_weight: pher_w,
    heuristic_weight: heur_w
  }
}

pub fn aco(config: ACOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_ants = config.ant_count
  let bound_range = bounds.max - bounds.min

  let mut current_rng = rng
  let (rng1, initial_sol) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  let mut best_solution = initial_sol
  let mut best_value = evaluate(f, best_solution)
  let mut evaluations = 1

  let mut sigma = bound_range * 0.5

  while evaluations < config.max_evaluations {
    let ant_positions: Array[Array[Double]] = []
    let ant_values: Array[Double] = []

    for ant in 0..<n_ants {
      let pos: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r1) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r2) = RNG::next(current_rng)
        current_rng = rng3
        let gauss = box_muller(r1, r2)
        let new_val = best_solution[d] + sigma * gauss
        pos.push(max(bounds.min, min(bounds.max, new_val)))
      }

      let val = evaluate(f, pos)
      evaluations = evaluations + 1
      ant_positions.push(pos)
      ant_values.push(val)

      if val < best_value {
        best_value = val
        best_solution = pos
      }

      if evaluations >= config.max_evaluations { break }
    }

    if evaluations >= config.max_evaluations { break }

    // Evaporate and deposit pheromones
    // Sigma shrinks faster when solution improves (stronger pheromone deposit)
    let improvement = max(0.0, 1.0 - best_value / (best_value + 1.0))
    sigma = sigma * (1.0 - config.evaporation_rate) + improvement * config.heuristic_weight * bound_range * 0.01
    // Add exploration noise to escape local optima
    let (rng4, r4) = RNG::next(current_rng)
    current_rng = rng4
    sigma = sigma + r4 * sigma * 0.1
    sigma = max(bound_range * 0.0001, min(bound_range, sigma))
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 蜂群算法 (Artificial Bee Colony) =====

pub struct BeeColonyConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  scout_bee_count: Int
  max_trials: Int
}

pub fn make_bee_colony_config(max_evals: Int, dim: Int, b: Bounds, scouts: Int, trials: Int) -> BeeColonyConfig {
  BeeColonyConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    scout_bee_count: scouts,
    max_trials: trials
  }
}

pub fn bee_colony(config: BeeColonyConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_food = config.scout_bee_count

  let mut current_rng = rng
  let food_sources: Array[Array[Double]] = []
  let food_values: Array[Double] = []
  let trial_counts: Array[Int] = []

  for i in 0..<n_food {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    food_sources.push(sol)
    food_values.push(evaluate(f, sol))
    trial_counts.push(0)
  }

  let mut evaluations = n_food
  let mut best_value = array_min(food_values)
  let mut best_idx = 0
  for i in 1..<food_values.length() {
    if food_values[i] < best_value {
      best_value = food_values[i]
      best_idx = i
    }
  }
  let mut best_solution = food_sources[best_idx]

  while evaluations < config.max_evaluations {
    // Employed bees phase
    for i in 0..<n_food {
      if evaluations >= config.max_evaluations { break }

      let neighbor: Array[Double] = []
      for d in 0..<dim {
        let (rng1, r) = RNG::next(current_rng)
        current_rng = rng1
        let phi = r * 2.0 - 1.0
        let (rng2, r2) = RNG::next(current_rng)
        current_rng = rng2
        let k = ((r2 * (dim.to_double() - 1.0)).floor().to_int() + 1) % dim
        let new_val = food_sources[i][d] + phi * (food_sources[i][d] - food_sources[k][d])
        neighbor.push(max(bounds.min, min(bounds.max, new_val)))
      }

      let new_val = evaluate(f, neighbor)
      evaluations = evaluations + 1

      if new_val < food_values[i] {
        food_sources[i] = neighbor
        food_values[i] = new_val
        trial_counts[i] = 0

        if new_val < best_value {
          best_value = new_val
          best_solution = neighbor
        }
      } else {
        trial_counts[i] = trial_counts[i] + 1
      }
    }

    if evaluations >= config.max_evaluations { break }

    // Onlooker bees phase
    let mut total_fitness = 0.0
    for i in 0..<n_food {
      total_fitness = total_fitness + 1.0 / (abs_val(food_values[i]) + 1.0)
    }

    for i in 0..<n_food {
      if evaluations >= config.max_evaluations { break }

      let prob = (1.0 / (abs_val(food_values[i]) + 1.0)) / total_fitness
      let (rng1, r) = RNG::next(current_rng)
      current_rng = rng1

      if r < prob {
        let neighbor: Array[Double] = []
        for d in 0..<dim {
          let (rng2, r2) = RNG::next(current_rng)
          current_rng = rng2
          let phi = r2 * 2.0 - 1.0
          let (rng3, r3) = RNG::next(current_rng)
          current_rng = rng3
          let k = ((r3 * (dim.to_double() - 1.0)).floor().to_int() + 1) % dim
          let new_val = food_sources[i][d] + phi * (food_sources[i][d] - food_sources[k][d])
          neighbor.push(max(bounds.min, min(bounds.max, new_val)))
        }

        let new_val = evaluate(f, neighbor)
        evaluations = evaluations + 1

        if new_val < food_values[i] {
          food_sources[i] = neighbor
          food_values[i] = new_val
          trial_counts[i] = 0
        } else {
          trial_counts[i] = trial_counts[i] + 1
        }
      }
    }

    if evaluations >= config.max_evaluations { break }

    // Scout bees phase
    for i in 0..<n_food {
      if trial_counts[i] >= config.max_trials {
        if evaluations >= config.max_evaluations { break }
        let (rng1, sol) = random_solution(dim, bounds, current_rng)
        current_rng = rng1
        food_sources[i] = sol
        food_values[i] = evaluate(f, sol)
        evaluations = evaluations + 1
        trial_counts[i] = 0

        if food_values[i] < best_value {
          best_value = food_values[i]
          best_solution = sol
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 禁忌搜索 (Tabu Search) =====

pub struct TabuSearchConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  tabu_list_size: Int
  neighborhood_size: Int
  step_size: Double
}

pub fn make_tabu_search_config(max_evals: Int, dim: Int, b: Bounds, tabu_size: Int, neigh_size: Int, step: Double) -> TabuSearchConfig {
  TabuSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    tabu_list_size: tabu_size,
    neighborhood_size: neigh_size,
    step_size: step
  }
}

pub fn tabu_search(config: TabuSearchConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let mut current_rng = rng
  let (rng1, initial_sol) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  let mut current_sol = initial_sol
  let mut current_val = evaluate(f, current_sol)
  let mut evaluations = 1

  let mut best_solution = current_sol
  let mut best_value = current_val

  let tabu_list: Array[String] = []

  while evaluations < config.max_evaluations {
    let mut best_neighbor: Array[Double] = []
    let mut best_neighbor_val = 1.0e308
    let mut best_move_str = ""
    let mut found_tabu_improving = false

    for n in 0..<config.neighborhood_size {
      let neighbor: Array[Double] = []
      let mut move_str = ""
      for d in 0..<dim {
        let (rng2, r) = RNG::next(current_rng)
        current_rng = rng2
        let delta = (r * 2.0 - 1.0) * config.step_size
        let old_val = current_sol[d]
        let new_val = max(bounds.min, min(bounds.max, old_val + delta))
        neighbor.push(new_val)
        move_str = move_str + "\{d}:\{new_val};"
      }

      let is_tabu = tabu_list.contains(move_str)

      if !is_tabu {
        let neighbor_val = evaluate(f, neighbor)
        evaluations = evaluations + 1

        if neighbor_val < best_neighbor_val {
          best_neighbor_val = neighbor_val
          best_neighbor = neighbor
          best_move_str = move_str
        }

        if neighbor_val < best_value {
          found_tabu_improving = true
        }
      }

      if evaluations >= config.max_evaluations { break }
    }

    if evaluations >= config.max_evaluations { break }

    if best_neighbor.length() > 0 {
      current_sol = best_neighbor
      current_val = best_neighbor_val

      if tabu_list.length() >= config.tabu_list_size {
        let _ = tabu_list.remove(0)
      }
      if best_move_str != "" {
        tabu_list.push(best_move_str)
      }

      if current_val < best_value {
        best_value = current_val
        best_solution = current_sol
      }
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 局部邻域搜索 (Guided Local Search) =====

pub struct LocalSearchConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  neighborhood_size: Int
  step_size: Double
  lambda: Double
}

pub fn make_local_search_config(max_evals: Int, dim: Int, b: Bounds, neigh_size: Int, step: Double, lambda_val: Double) -> LocalSearchConfig {
  LocalSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    neighborhood_size: neigh_size,
    step_size: step,
    lambda: lambda_val
  }
}

pub fn local_search(config: LocalSearchConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let mut current_rng = rng
  let (rng1, initial_sol) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  let mut current_sol = initial_sol
  let mut current_val = evaluate(f, current_sol)
  let mut evaluations = 1

  let mut best_solution = current_sol
  let mut best_value = current_val

  let penalties: Array[Double] = []
  for i in 0..<dim {
    penalties.push(0.0)
  }

  while evaluations < config.max_evaluations {
    let mut best_neighbor: Array[Double] = []
    let mut best_altered = 1.0e308
    let mut best_original_val = 1.0e308

    for n in 0..<config.neighborhood_size {
      let neighbor: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r) = RNG::next(current_rng)
        current_rng = rng2
        let delta = (r * 2.0 - 1.0) * config.step_size
        neighbor.push(max(bounds.min, min(bounds.max, current_sol[d] + delta)))
      }

      let neighbor_val = evaluate(f, neighbor)
      evaluations = evaluations + 1

      let mut altered = neighbor_val
      for d in 0..<dim {
        altered = altered + config.lambda * penalties[d] * abs_val(neighbor[d] - current_sol[d])
      }

      if altered < best_altered {
        best_altered = altered
        best_neighbor = neighbor
        best_original_val = neighbor_val
      }

      if evaluations >= config.max_evaluations { break }
    }

    if evaluations >= config.max_evaluations { break }

    if best_neighbor.length() > 0 {
      // Update penalties for features of the best solution that changed
      for d in 0..<dim {
        if abs_val(best_neighbor[d] - current_sol[d]) > config.step_size * 0.1 {
          penalties[d] = penalties[d] + 1.0
        }
      }

      current_sol = best_neighbor
      current_val = best_original_val

      if current_val < best_value {
        best_value = current_val
        best_solution = current_sol
        // Reset penalties when new best found
        for d in 0..<dim {
          penalties[d] = 0.0
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== CMA-ES (协方差矩阵自适应进化策略) =====

pub struct CMAESConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  sigma: Double
  lambda: Int
  mu: Int
}

pub fn make_cmaes_config(max_evals: Int, dim: Int, b: Bounds, sigma_init: Double, lambda_val: Int, mu_val: Int) -> CMAESConfig {
  CMAESConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    sigma: sigma_init,
    lambda: lambda_val,
    mu: mu_val
  }
}

pub fn cma_es(config: CMAESConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let lambda = config.lambda
  let mu = config.mu

  let mut current_rng = rng
  let (rng1, mean_sol) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  let mean = mean_sol

  let mut sigma = config.sigma

  let cov: Array[Array[Double]] = []
  for i in 0..<dim {
    let row: Array[Double] = []
    for j in 0..<dim {
      if i == j { row.push(1.0) } else { row.push(0.0) }
    }
    cov.push(row)
  }

  let mut evaluations = 0
  let mut best_solution: Array[Double] = []
  let mut best_value = 1.0e308

  while evaluations < config.max_evaluations {
    let population: Array[Array[Double]] = []
    let pop_values: Array[Double] = []

    for i in 0..<lambda {
      let individual: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r) = RNG::next(current_rng)
        current_rng = rng2

        let mut val = mean[d]
        for d2 in 0..<dim {
          let (rng3, r2) = RNG::next(current_rng)
          current_rng = rng3
          let gauss = box_muller(r, r2)
          val = val + sigma * cov[d][d2] * gauss
        }
        individual.push(max(bounds.min, min(bounds.max, val)))
      }

      let val = evaluate(f, individual)
      evaluations = evaluations + 1
      population.push(individual)
      pop_values.push(val)

      if val < best_value {
        best_value = val
        best_solution = individual
      }

      if evaluations >= config.max_evaluations { break }
    }

    if evaluations >= config.max_evaluations { break }

    // Simple gradient-free CMA update
    let sorted_idx = sort_indices(pop_values)
    let new_mean: Array[Double] = []
    for d in 0..<dim {
      let mut sum = 0.0
      for i in 0..<mu {
        sum = sum + population[sorted_idx[i]][d]
      }
      new_mean.push(sum / mu.to_double())
    }

    for d in 0..<dim {
      mean[d] = new_mean[d]
    }

    sigma = sigma * 0.99
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 交叉熵方法 (Cross-Entropy Method) =====

pub struct CEMConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  sample_size: Int
  elite_fraction: Double
  smoothing: Double
}

pub fn make_cem_config(max_evals: Int, dim: Int, b: Bounds, sample_sz: Int, elite_frac: Double, smooth: Double) -> CEMConfig {
  CEMConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    sample_size: sample_sz,
    elite_fraction: elite_frac,
    smoothing: smooth
  }
}

pub fn cross_entropy_method(config: CEMConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_samples = config.sample_size
  let n_elite = max_int(1, (n_samples.to_double() * config.elite_fraction).to_int())

  let mut current_rng = rng
  let (rng1, initial_sol) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  let mean: Array[Double] = []
  for d in 0..<dim { mean.push(initial_sol[d]) }

  let bound_range = bounds.max - bounds.min
  let std: Array[Double] = []
  for d in 0..<dim { std.push(bound_range * 0.5) }

  let mut evaluations = 0
  let mut best_solution: Array[Double] = []
  let mut best_value = 1.0e308

  while evaluations < config.max_evaluations {
    // Sample from Gaussian
    let samples: Array[Array[Double]] = []
    let sample_vals: Array[Double] = []

    for s in 0..<n_samples {
      let sample: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r1) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r2) = RNG::next(current_rng)
        current_rng = rng3
        let gauss = box_muller(r1, r2)
        let val = mean[d] + std[d] * gauss
        sample.push(max(bounds.min, min(bounds.max, val)))
      }
      let val = evaluate(f, sample)
      evaluations = evaluations + 1
      samples.push(sample)
      sample_vals.push(val)

      if val < best_value {
        best_value = val
        best_solution = sample
      }

      if evaluations >= config.max_evaluations { break }
    }

    if evaluations >= config.max_evaluations { break }

    // Select elite samples
    let sorted_idx = sort_indices(sample_vals)
    let elite_mean: Array[Double] = []
    let elite_std: Array[Double] = []

    for d in 0..<dim {
      let mut sum = 0.0
      for i in 0..<n_elite {
        sum = sum + samples[sorted_idx[i]][d]
      }
      let em = sum / n_elite.to_double()
      elite_mean.push(em)

      let mut var_sum = 0.0
      for i in 0..<n_elite {
        let diff = samples[sorted_idx[i]][d] - em
        var_sum = var_sum + diff * diff
      }
      elite_std.push((var_sum / n_elite.to_double()).sqrt())
    }

    // Smooth update
    for d in 0..<dim {
      mean[d] = config.smoothing * mean[d] + (1.0 - config.smoothing) * elite_mean[d]
      std[d] = config.smoothing * std[d] + (1.0 - config.smoothing) * elite_std[d]
      std[d] = max(bound_range * 0.0001, min(bound_range, std[d]))
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 布谷鸟搜索算法 (Cuckoo Search) =====

pub struct CuckooSearchConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  nest_count: Int
  abandon_probability: Double
  step_size: Double
}

pub fn make_cuckoo_config(max_evals: Int, dim: Int, b: Bounds, nests: Int, abandon_prob: Double, step: Double) -> CuckooSearchConfig {
  CuckooSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    nest_count: nests,
    abandon_probability: abandon_prob,
    step_size: step
  }
}

pub fn cuckoo_search(config: CuckooSearchConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_nests = config.nest_count
  let bound_range = bounds.max - bounds.min

  let mut current_rng = rng
  let nests: Array[Array[Double]] = []
  let nest_values: Array[Double] = []

  for i in 0..<n_nests {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    nests.push(sol)
    nest_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_nests
  let mut best_idx = 0
  let mut best_value = nest_values[0]
  for i in 1..<n_nests {
    if nest_values[i] < best_value {
      best_value = nest_values[i]
      best_idx = i
    }
  }
  let mut best_solution = nests[best_idx]

  while evaluations < config.max_evaluations {
    // Lévy flight for each cuckoo
    for i in 0..<n_nests {
      if evaluations >= config.max_evaluations { break }

      let new_nest: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r1) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r2) = RNG::next(current_rng)
        current_rng = rng3
        let gauss1 = box_muller(r1, r2)
        let (rng4, r3) = RNG::next(current_rng)
        current_rng = rng4
        let levy = (1.0 / pow_approx(abs_val(gauss1), 2.0 / 3.0)).sqrt() * (r3 * 2.0 - 1.0)
        let step = config.step_size * bound_range * levy
        new_nest.push(max(bounds.min, min(bounds.max, nests[i][d] + step * gauss1)))
      }

      let new_val = evaluate(f, new_nest)
      evaluations = evaluations + 1

      if new_val < nest_values[i] {
        nests[i] = new_nest
        nest_values[i] = new_val
      }

      if new_val < best_value {
        best_value = new_val
        best_solution = new_nest
      }
    }

    if evaluations >= config.max_evaluations { break }

    // Abandon worst nests and rebuild
    let sorted_idx = sort_indices(nest_values)
    let n_abandon = max_int(1, (n_nests.to_double() * config.abandon_probability).to_int())

    for j in 0..<n_abandon {
      if evaluations >= config.max_evaluations { break }
      let worst_idx = sorted_idx[n_nests - 1 - j]
      let new_sol: Array[Double] = []
      for d in 0..<dim {
        let (rng5, r1) = RNG::next(current_rng)
        current_rng = rng5
        let (rng6, r2) = RNG::next(current_rng)
        current_rng = rng6
        let gauss = box_muller(r1, r2)
        new_sol.push(max(bounds.min, min(bounds.max, best_solution[d] + bound_range * 0.1 * gauss)))
      }
      nests[worst_idx] = new_sol
      nest_values[worst_idx] = evaluate(f, new_sol)
      evaluations = evaluations + 1

      if nest_values[worst_idx] < best_value {
        best_value = nest_values[worst_idx]
        best_solution = new_sol
      }
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== 灰狼优化算法 (Grey Wolf Optimizer) =====

pub struct GWOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  pack_size: Int
}

pub fn make_gwo_config(max_evals: Int, dim: Int, b: Bounds, pack_sz: Int) -> GWOConfig {
  GWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    pack_size: pack_sz
  }
}

pub fn grey_wolf_optimizer(config: GWOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_wolves = config.pack_size

  let mut current_rng = rng
  let wolves: Array[Array[Double]] = []
  let wolf_values: Array[Double] = []

  for i in 0..<n_wolves {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    wolves.push(sol)
    wolf_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_wolves
  let mut alpha_pos: Array[Double] = []
  let mut alpha_val = 1.0e308
  let mut beta_pos: Array[Double] = []
  let mut beta_val = 1.0e308
  let mut delta_pos: Array[Double] = []
  let mut delta_val = 1.0e308

  while evaluations < config.max_evaluations {
    // Update alpha, beta, delta
    alpha_val = 1.0e308
    beta_val = 1.0e308
    delta_val = 1.0e308
    for i in 0..<n_wolves {
      if wolf_values[i] < alpha_val {
        delta_val = beta_val
        delta_pos = beta_pos
        beta_val = alpha_val
        beta_pos = alpha_pos
        alpha_val = wolf_values[i]
        alpha_pos = wolves[i]
      } else if wolf_values[i] < beta_val {
        delta_val = beta_val
        delta_pos = beta_pos
        beta_val = wolf_values[i]
        beta_pos = wolves[i]
      } else if wolf_values[i] < delta_val {
        delta_val = wolf_values[i]
        delta_pos = wolves[i]
      }
    }

    // Update wolves
    let a = 2.0 - 2.0 * (evaluations.to_double() / config.max_evaluations.to_double())

    for i in 0..<n_wolves {
      if evaluations >= config.max_evaluations { break }
      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r1) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r2) = RNG::next(current_rng)
        current_rng = rng3
        let r_a1 = r1 * 2.0 - 1.0
        let r_a2 = r2 * 2.0 - 1.0

        let d_alpha = abs_val(alpha_pos[d] - wolves[i][d])
        let x1 = alpha_pos[d] - a * r_a1 * d_alpha

        let (rng4, r3) = RNG::next(current_rng)
        current_rng = rng4
        let (rng5, r4) = RNG::next(current_rng)
        current_rng = rng5
        let r_b1 = r3 * 2.0 - 1.0
        let r_b2 = r4 * 2.0 - 1.0

        let d_beta = abs_val(beta_pos[d] - wolves[i][d])
        let x2 = beta_pos[d] - a * r_b1 * d_beta

        let (rng6, r5) = RNG::next(current_rng)
        current_rng = rng6
        let (rng7, r6) = RNG::next(current_rng)
        current_rng = rng7
        let r_c1 = r5 * 2.0 - 1.0
        let r_c2 = r6 * 2.0 - 1.0

        let d_delta = abs_val(delta_pos[d] - wolves[i][d])
        let x3 = delta_pos[d] - a * r_c1 * d_delta

        let new_val = (x1 + x2 + x3) / 3.0
        new_pos.push(max(bounds.min, min(bounds.max, new_val)))
      }

      wolves[i] = new_pos
      let new_fit = evaluate(f, new_pos)
      wolf_values[i] = new_fit
      evaluations = evaluations + 1
    }
  }

  OptimizationResult::{
    best_solution: alpha_pos,
    best_value: alpha_val,
    evaluations: evaluations
  }
}

// ===== Hooke-Jeeves Pattern Search =====

pub struct HookeJeevesConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  initial_step: Double
  step_reduction: Double
  min_step: Double
}

pub fn make_hooke_jeeves_config(max_evals: Int, dim: Int, b: Bounds, init_step: Double, step_red: Double, min_s: Double) -> HookeJeevesConfig {
  HookeJeevesConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    initial_step: init_step,
    step_reduction: step_red,
    min_step: min_s
  }
}

pub fn hooke_jeeves(config: HookeJeevesConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let mut current_rng = rng
  let (rng1, base) = random_solution(dim, bounds, current_rng)
  current_rng = rng1
  let mut base_point = base
  let mut base_val = evaluate(f, base_point)
  let mut evaluations = 1

  let mut step = config.initial_step
  let mut best_solution = base_point
  let mut best_value = base_val

  while evaluations < config.max_evaluations && step > config.min_step {
    let mut improved = false

    // Exploratory move
    let mut current_point: Array[Double] = []
    for d in 0..<dim { current_point.push(base_point[d]) }
    let mut current_val = base_val

    for d in 0..<dim {
      if evaluations >= config.max_evaluations { break }

      // Try positive direction
      let trial1: Array[Double] = []
      for d2 in 0..<dim {
        if d2 == d {
          trial1.push(max(bounds.min, min(bounds.max, current_point[d2] + step)))
        } else {
          trial1.push(current_point[d2])
        }
      }
      let val1 = evaluate(f, trial1)
      evaluations = evaluations + 1

      if val1 < current_val {
        current_point = trial1
        current_val = val1
        improved = true
      } else {
        // Try negative direction
        let trial2: Array[Double] = []
        for d2 in 0..<dim {
          if d2 == d {
            trial2.push(max(bounds.min, min(bounds.max, current_point[d2] - step)))
          } else {
            trial2.push(current_point[d2])
          }
        }
        let val2 = evaluate(f, trial2)
        evaluations = evaluations + 1

        if val2 < current_val {
          current_point = trial2
          current_val = val2
          improved = true
        }
      }

      if current_val < best_value {
        best_value = current_val
        best_solution = current_point
      }
    }

    if evaluations >= config.max_evaluations { break }

    if improved {
      // Pattern move
      let pattern_point: Array[Double] = []
      for d in 0..<dim {
        pattern_point.push(max(bounds.min, min(bounds.max, current_point[d] + (current_point[d] - base_point[d]))))
      }

      let pat_val = evaluate(f, pattern_point)
      evaluations = evaluations + 1

      if pat_val < current_val {
        base_point = pattern_point
        base_val = pat_val
      } else {
        base_point = current_point
        base_val = current_val
      }

      if base_val < best_value {
        best_value = base_val
        best_solution = base_point
      }
    } else {
      // Reduce step size
      step = step * config.step_reduction
    }
  }

  OptimizationResult::{
    best_solution: best_solution,
    best_value: best_value,
    evaluations: evaluations
  }
}

// ===== Whale Optimization Algorithm (WOA) =====

pub struct WOAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  spiral_coeff: Double
}

pub fn make_woa_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, spiral: Double) -> WOAConfig {
  WOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    spiral_coeff: spiral
  }
}

pub fn whale_optimization(config: WOAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_whales = config.population_size
  let b = config.spiral_coeff

  let mut current_rng = rng
  let whales: Array[Array[Double]] = []
  let whale_values: Array[Double] = []

  for i in 0..<n_whales {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    whales.push(sol)
    whale_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_whales
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_whales {
    if whale_values[i] < best_val {
      best_val = whale_values[i]
      best_pos = whales[i]
    }
  }

  let bound_range = bounds.max - bounds.min

  while evaluations < config.max_evaluations {
    let a = 2.0 - 2.0 * (evaluations.to_double() / config.max_evaluations.to_double())
    let a2 = a * a

    for i in 0..<n_whales {
      if evaluations >= config.max_evaluations { break }

      let new_pos: Array[Double] = []
      let (rng2, p_rand) = RNG::next(current_rng)
      current_rng = rng2
      let (rng3, r_05) = RNG::next(current_rng)
      current_rng = rng3

      if p_rand < 0.5 {
        // Encircling prey or search for prey
        if r_05 < 0.5 {
          // Search for prey (exploration)
          let (rng4, r1) = RNG::next(current_rng)
          current_rng = rng4
          let (rng5, r2) = RNG::next(current_rng)
          current_rng = rng5
          let (rng6, r3) = RNG::next(current_rng)
          current_rng = rng6
          let rand_idx = ((r1 * n_whales.to_double()).floor()).to_int()
          let r_idx = if rand_idx >= n_whales { n_whales - 1 } else { rand_idx }

          for d in 0..<dim {
            let (rng7, rd1) = RNG::next(current_rng)
            current_rng = rng7
            let (rng8, rd2) = RNG::next(current_rng)
            current_rng = rng8
            let r_d1 = rd1 * 2.0 - 1.0
            let r_d2 = rd2 * 2.0 - 1.0
            let d_leader = abs_val(whales[r_idx][d] - whales[i][d])
            let pos = whales[r_idx][d] - a2 * r_d1 * d_leader
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        } else {
          // Encircling prey (exploitation)
          for d in 0..<dim {
            let (rng7, rd1) = RNG::next(current_rng)
            current_rng = rng7
            let (rng8, rd2) = RNG::next(current_rng)
            current_rng = rng8
            let r_d1 = rd1 * 2.0 - 1.0
            let r_d2 = rd2 * 2.0 - 1.0
            let d_prey = abs_val(best_pos[d] - whales[i][d])
            let pos = best_pos[d] - a * r_d1 * d_prey
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        }
      } else {
        // Bubble-net attacking (spiral update)
        for d in 0..<dim {
          let (rng7, rd1) = RNG::next(current_rng)
          current_rng = rng7
          let (rng8, rd2) = RNG::next(current_rng)
          current_rng = rng8
          let r_d1 = rd1 * 2.0 - 1.0
          let r_d2 = rd2 * 2.0 - 1.0
          let d_prey = abs_val(best_pos[d] - whales[i][d])
          let angle = 2.0 * 3.14159265358979323846 * r_d2
          let spiral = exp_approx(b * angle) * cos_approx(2.0 * 3.14159265358979323846 * r_d2)
          let pos = d_prey * spiral + best_pos[d]
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      whales[i] = new_pos
      let new_val = evaluate(f, new_pos)
      whale_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Bat Algorithm (BA) =====

pub struct BatConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  f_min: Double
  f_max: Double
  pulse_rate: Double
  loudness: Double
}

pub fn make_bat_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, fmin: Double, fmax: Double, pr: Double, loud: Double) -> BatConfig {
  BatConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    f_min: fmin,
    f_max: fmax,
    pulse_rate: pr,
    loudness: loud
  }
}

pub fn bat_algorithm(config: BatConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_bats = config.population_size

  let mut current_rng = rng
  let bats: Array[Array[Double]] = []
  let velocities: Array[Array[Double]] = []
  let bat_values: Array[Double] = []
  let frequencies: Array[Double] = []

  for i in 0..<n_bats {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    let vel: Array[Double] = []
    for d in 0..<dim { vel.push(0.0) }
    bats.push(sol)
    velocities.push(vel)
    bat_values.push(evaluate(f, sol))
    frequencies.push(0.0)
  }

  let mut evaluations = n_bats
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_bats {
    if bat_values[i] < best_val {
      best_val = bat_values[i]
      best_pos = bats[i]
    }
  }

  let mut loudness = config.loudness
  let mut pulse_rate = config.pulse_rate
  let bound_range = bounds.max - bounds.min

  while evaluations < config.max_evaluations {
    for i in 0..<n_bats {
      if evaluations >= config.max_evaluations { break }

      // Update frequency
      let (rng2, r_freq) = RNG::next(current_rng)
      current_rng = rng2
      frequencies[i] = config.f_min + (config.f_max - config.f_min) * r_freq

      // Update velocity
      let new_vel: Array[Double] = []
      for d in 0..<dim {
        let (rng3, r_d) = RNG::next(current_rng)
        current_rng = rng3
        new_vel.push(velocities[i][d] + (bats[i][d] - best_pos[d]) * frequencies[i])
      }
      velocities[i] = new_vel

      // Update position
      let new_pos: Array[Double] = []
      for d in 0..<dim {
        new_pos.push(max(bounds.min, min(bounds.max, bats[i][d] + velocities[i][d])))
      }

      // Local search
      let (rng4, r_pulse) = RNG::next(current_rng)
      current_rng = rng4
      if r_pulse > pulse_rate {
        for d in 0..<dim {
          let (rng5, r_step) = RNG::next(current_rng)
          current_rng = rng5
          let step = loudness * (r_step * 2.0 - 1.0)
          new_pos[d] = max(bounds.min, min(bounds.max, best_pos[d] + step))
        }
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      // Accept if better and loudness condition met
      let (rng6, r_loud) = RNG::next(current_rng)
      current_rng = rng6
      if new_val < bat_values[i] && r_loud < loudness {
        bats[i] = new_pos
        bat_values[i] = new_val
      }

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos

        // Update loudness and pulse rate
        loudness = loudness * 0.9
        pulse_rate = pulse_rate + 0.01
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Harmony Search (HS) =====

pub struct HSConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  harmony_memory_size: Int
  hmcr: Double
  par: Double
  bw: Double
}

pub fn make_hs_config(max_evals: Int, dim: Int, b: Bounds, hms: Int, hmcr_val: Double, par_val: Double, bw_val: Double) -> HSConfig {
  HSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    harmony_memory_size: hms,
    hmcr: hmcr_val,
    par: par_val,
    bw: bw_val
  }
}

pub fn harmony_search(config: HSConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let hms = config.harmony_memory_size
  let hmcr = config.hmcr
  let par = config.par
  let bw = config.bw

  let mut current_rng = rng
  let harmony_memory: Array[Array[Double]] = []
  let hm_values: Array[Double] = []

  for i in 0..<hms {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    harmony_memory.push(sol)
    hm_values.push(evaluate(f, sol))
  }

  let mut evaluations = hms
  let mut best_idx = 0
  let mut worst_idx = 0
  let mut best_val = hm_values[0]
  let mut worst_val = hm_values[0]

  for i in 1..<hms {
    if hm_values[i] < best_val {
      best_val = hm_values[i]
      best_idx = i
    }
    if hm_values[i] > worst_val {
      worst_val = hm_values[i]
      worst_idx = i
    }
  }

  let bound_range = bounds.max - bounds.min

  while evaluations < config.max_evaluations {
    // Create new harmony
    let new_harmony: Array[Double] = []

    for d in 0..<dim {
      let (rng2, r1) = RNG::next(current_rng)
      current_rng = rng2

      if r1 < hmcr {
        // Memory consideration
        let (rng3, r2) = RNG::next(current_rng)
        current_rng = rng3
        let idx = ((r2 * hms.to_double()).floor()).to_int()
        let pick_idx = if idx >= hms { hms - 1 } else { idx }
        let val = harmony_memory[pick_idx][d]

        // Pitch adjustment
        let (rng4, r3) = RNG::next(current_rng)
        current_rng = rng4
        if r3 < par {
          let (rng5, r4) = RNG::next(current_rng)
          current_rng = rng5
          let delta = (r4 * 2.0 - 1.0) * bw * bound_range
          new_harmony.push(max(bounds.min, min(bounds.max, val + delta)))
        } else {
          new_harmony.push(val)
        }
      } else {
        // Randomization
        let (rng6, r5) = RNG::next(current_rng)
        current_rng = rng6
        new_harmony.push(bounds.min + r5 * bound_range)
      }
    }

    let new_val = evaluate(f, new_harmony)
    evaluations = evaluations + 1

    if new_val < worst_val {
      // Replace worst harmony
      harmony_memory[worst_idx] = new_harmony
      hm_values[worst_idx] = new_val

      // Find new best and worst
      best_val = hm_values[0]
      worst_val = hm_values[0]
      best_idx = 0
      worst_idx = 0
      for i in 1..<hms {
        if hm_values[i] < best_val {
          best_val = hm_values[i]
          best_idx = i
        }
        if hm_values[i] > worst_val {
          worst_val = hm_values[i]
          worst_idx = i
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: harmony_memory[best_idx],
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Sine Cosine Algorithm (SCA) =====

pub struct SCAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
}

pub fn make_sca_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> SCAConfig {
  SCAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz
  }
}

pub fn sine_cosine_algorithm(config: SCAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_sca = config.population_size

  let mut current_rng = rng
  let population: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n_sca {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n_sca
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_sca {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = population[i]
    }
  }

  let mut r1 = 2.0

  while evaluations < config.max_evaluations {
    let r1_update = 2.0 - evaluations.to_double() * (2.0 / config.max_evaluations.to_double())

    for i in 0..<n_sca {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r2_val) = RNG::next(current_rng)
      current_rng = rng2
      let (rng3, r3_val) = RNG::next(current_rng)
      current_rng = rng3
      let (rng4, r4_val) = RNG::next(current_rng)
      current_rng = rng4

      let r2 = 2.0 * PI * r2_val
      let r3 = 2.0 * r3_val - 1.0
      let r4 = r4_val

      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let dist = abs_val(r3 * best_pos[d] - population[i][d])
        let pos = if r4 < 0.5 {
          population[i][d] + r1_update * sin_approx(r2) * dist
        } else {
          population[i][d] + r1_update * cos_approx(r2) * dist
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      population[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    r1 = r1_update
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Moth Flame Optimization (MFO) =====

pub struct MFOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  spiral_constant: Double
}

pub fn make_mfo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, spiral: Double) -> MFOConfig {
  MFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    spiral_constant: spiral
  }
}

pub fn moth_flame_optimization(config: MFOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_moths = config.population_size
  let b = config.spiral_constant

  let mut current_rng = rng
  let moths: Array[Array[Double]] = []
  let moth_values: Array[Double] = []

  for i in 0..<n_moths {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    moths.push(sol)
    moth_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_moths
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_moths {
    if moth_values[i] < best_val {
      best_val = moth_values[i]
      best_pos = moths[i]
    }
  }

  let flames: Array[Array[Double]] = []
  let flame_values: Array[Double] = []
  for i in 0..<n_moths {
    flames.push(moths[i])
    flame_values.push(moth_values[i])
  }

  while evaluations < config.max_evaluations {
    let iter = (evaluations.to_double() / config.max_evaluations.to_double())
    let num_flames = max_int(1, (n_moths.to_double() * (1.0 - iter)).to_int())

    // Sort flames by fitness (ascending)
    let sorted_flame_indices: Array[Int] = []
    for i in 0..<flame_values.length() {
      sorted_flame_indices.push(i)
    }
    for i in 0..<sorted_flame_indices.length() {
      for j in (i + 1)..<sorted_flame_indices.length() {
        if flame_values[sorted_flame_indices[j]] < flame_values[sorted_flame_indices[i]] {
          let temp = sorted_flame_indices[i]
          sorted_flame_indices[i] = sorted_flame_indices[j]
          sorted_flame_indices[j] = temp
        }
      }
    }

    // Sort moths by fitness
    let sorted_moth_indices: Array[Int] = []
    for i in 0..<moth_values.length() {
      sorted_moth_indices.push(i)
    }
    for i in 0..<sorted_moth_indices.length() {
      for j in (i + 1)..<sorted_moth_indices.length() {
        if moth_values[sorted_moth_indices[j]] < moth_values[sorted_moth_indices[i]] {
          let temp = sorted_moth_indices[i]
          sorted_moth_indices[i] = sorted_moth_indices[j]
          sorted_moth_indices[j] = temp
        }
      }
    }

    // Update each moth using spiral path toward its assigned flame
    for rank in 0..<n_moths {
      if evaluations >= config.max_evaluations { break }

      let moth_idx = sorted_moth_indices[rank]
      let flame_idx = min_int(rank, num_flames - 1)
      let flame_idx_sorted = sorted_flame_indices[flame_idx]

      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r_t) = RNG::next(current_rng)
        current_rng = rng2
        let t = r_t * 2.0 - 1.0
        let flame_pos = flames[flame_idx_sorted][d]
        let dist = abs_val(flame_pos - moths[moth_idx][d])
        let pos = dist * exp_approx(b * t) * cos_approx(2.0 * PI * t) + flame_pos
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      moths[moth_idx] = new_pos
      let new_val = evaluate(f, new_pos)
      moth_values[moth_idx] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Update flames: move flames to best moths
    for i in 0..<n_moths {
      if moth_values[i] < flame_values[i] {
        flames[i] = moths[i]
        flame_values[i] = moth_values[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Antlion Optimizer (ALO) =====

pub struct ALOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  ant_count: Int
}

pub fn make_alo_config(max_evals: Int, dim: Int, b: Bounds, ants: Int) -> ALOConfig {
  ALOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    ant_count: ants
  }
}

pub fn antlion_optimizer(config: ALOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_ants = config.ant_count
  let n_antlions = config.ant_count

  let mut current_rng = rng
  let antlions: Array[Array[Double]] = []
  let antlion_values: Array[Double] = []

  for i in 0..<n_antlions {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    antlions.push(sol)
    antlion_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_antlions
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_antlions {
    if antlion_values[i] < best_val {
      best_val = antlion_values[i]
      best_pos = antlions[i]
    }
  }

  let ants: Array[Array[Double]] = []
  let ant_values: Array[Double] = []
  for i in 0..<n_ants {
    ants.push(antlions[i])
    ant_values.push(antlion_values[i])
  }

  let base_range = bounds.max - bounds.min
  let mut iter = 0

  while evaluations < config.max_evaluations {
    iter = iter + 1

    // Sort antlions by fitness
    let sorted_indices: Array[Int] = []
    for i in 0..<n_antlions {
      sorted_indices.push(i)
    }
    for i in 0..<n_antlions {
      for j in (i + 1)..<n_antlions {
        if antlion_values[sorted_indices[j]] < antlion_values[sorted_indices[i]] {
          let temp = sorted_indices[i]
          sorted_indices[i] = sorted_indices[j]
          sorted_indices[j] = temp
        }
      }
    }

    let elite_count = max_int(2, n_antlions / 2)

    for i in 0..<n_ants {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_select) = RNG::next(current_rng)
      current_rng = rng2

      // Tournament selection: pick from elite pool with probability, else from entire pool
      let sel_idx = if r_select < 0.7 {
        let (rng3, rw) = RNG::next(current_rng)
        current_rng = rng3
        (rw * elite_count.to_double()).floor().to_int()
      } else {
        let (rng3, rw) = RNG::next(current_rng)
        current_rng = rng3
        (rw * n_antlions.to_double()).floor().to_int()
      }
      let sel = max_int(0, min_int(n_antlions - 1, sel_idx))
      let selected_antlion = antlions[sorted_indices[sel]]

      // Boundary shrinking
      let alpha_frac = 1.0 / (iter.to_double() + 1.0)
      let lower_bound = bounds.min + alpha_frac * abs_val(selected_antlion[0] - bounds.min)
      let upper_bound = bounds.max - alpha_frac * abs_val(bounds.max - selected_antlion[0])

      // Random walk toward the selected antlion
      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng4, r_walk) = RNG::next(current_rng)
        current_rng = rng4
        let walk_val = r_walk * 2.0 - 1.0
        let antlion_d = selected_antlion[d]
        let dist = abs_val(antlion_d - ants[i][d])
        let pos = ants[i][d] + walk_val * dist
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      ants[i] = new_pos
      let new_val = evaluate(f, new_pos)
      ant_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Elitism: replace worst antlions with best ants
    let sorted_antlion_idx: Array[Int] = []
    for i in 0..<n_antlions {
      sorted_antlion_idx.push(i)
    }
    for i in 0..<n_antlions {
      for j in (i + 1)..<n_antlions {
        if antlion_values[sorted_antlion_idx[j]] > antlion_values[sorted_antlion_idx[i]] {
          let temp = sorted_antlion_idx[i]
          sorted_antlion_idx[i] = sorted_antlion_idx[j]
          sorted_antlion_idx[j] = temp
        }
      }
    }

    let sorted_ant_idx: Array[Int] = []
    for i in 0..<n_ants {
      sorted_ant_idx.push(i)
    }
    for i in 0..<n_ants {
      for j in (i + 1)..<n_ants {
        if ant_values[sorted_ant_idx[j]] < ant_values[sorted_ant_idx[i]] {
          let temp = sorted_ant_idx[i]
          sorted_ant_idx[i] = sorted_ant_idx[j]
          sorted_ant_idx[j] = temp
        }
      }
    }

    // Replace worst antlions with best ants
    let replace_count = min_int(n_ants / 2, n_antlions / 3)
    for k in 0..<replace_count {
      let worst_idx = sorted_antlion_idx[n_antlions - 1 - k]
      let best_ant_idx = sorted_ant_idx[k]
      if ant_values[best_ant_idx] < antlion_values[worst_idx] {
        antlions[worst_idx] = ants[best_ant_idx]
        antlion_values[worst_idx] = ant_values[best_ant_idx]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Sparrow Search Algorithm (SSA) =====

pub struct SSAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  discoverer_ratio: Double
  warning_threshold: Double
}

pub fn make_ssa_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> SSAConfig {
  SSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    discoverer_ratio: 0.2,
    warning_threshold: 0.8
  }
}

pub fn sparrow_search_algorithm(config: SSAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let population: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = population[i]
    }
  }

  let mut worst_pos: Array[Double] = []
  let mut worst_val = -1.0e308
  for i in 0..<n {
    if values[i] > worst_val {
      worst_val = values[i]
      worst_pos = population[i]
    }
  }

  let discoverer_count = max_int(1, (n.to_double() * config.discoverer_ratio).to_int())
  let follower_count = n - discoverer_count

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let alpha_warning = config.warning_threshold * (1.0 - iter)

    let sorted_indices: Array[Int] = []
    for i in 0..<n {
      sorted_indices.push(i)
    }
    for i in 0..<n {
      for j in (i + 1)..<n {
        if values[sorted_indices[j]] < values[sorted_indices[i]] {
          let temp = sorted_indices[i]
          sorted_indices[i] = sorted_indices[j]
          sorted_indices[j] = temp
        }
      }
    }

    // Discoverer phase
    for k in 0..<discoverer_count {
      if evaluations >= config.max_evaluations { break }

      let idx = sorted_indices[k]
      let (rng2, r_warn) = RNG::next(current_rng)
      current_rng = rng2

      let new_pos: Array[Double] = []
      if r_warn < alpha_warning {
        for d in 0..<dim {
          let (rng3, r_step) = RNG::next(current_rng)
          current_rng = rng3
          let step = r_step * (bounds.max - bounds.min) + bounds.min
          new_pos.push(max(bounds.min, min(bounds.max, step)))
        }
      } else {
        for d in 0..<dim {
          let (rng3, r_foll) = RNG::next(current_rng)
          current_rng = rng3
          let step = best_pos[d] + r_foll * (population[idx][d] - best_pos[d])
          new_pos.push(max(bounds.min, min(bounds.max, step)))
        }
      }

      population[idx] = new_pos
      let new_val = evaluate(f, new_pos)
      values[idx] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Follower phase
    for k in 0..<follower_count {
      if evaluations >= config.max_evaluations { break }

      let idx = sorted_indices[discoverer_count + k]

      if values[idx] > best_val {
        let (rng4, r_fly) = RNG::next(current_rng)
        current_rng = rng4

        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let pos = best_pos[d] + r_fly * (population[idx][d] - best_pos[d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        population[idx] = new_pos
        let new_val = evaluate(f, new_pos)
        values[idx] = new_val
        evaluations = evaluations + 1

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      } else {
        let (rng4, r_rand) = RNG::next(current_rng)
        current_rng = rng4

        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let pos = r_rand * (worst_pos[d] - population[idx][d]) + worst_pos[d]
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        population[idx] = new_pos
        let new_val = evaluate(f, new_pos)
        values[idx] = new_val
        evaluations = evaluations + 1

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      }
    }

    // Update worst
    worst_val = -1.0e308
    for i in 0..<n {
      if values[i] > worst_val {
        worst_val = values[i]
        worst_pos = population[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Harris Hawks Optimization (HHO) =====

pub struct HHOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
}

pub fn make_hho_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> HHOConfig {
  HHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz
  }
}

pub fn harris_hawks_optimization(config: HHOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let hawks: Array[Array[Double]] = []
  let hawk_values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    hawks.push(sol)
    hawk_values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if hawk_values[i] < best_val {
      best_val = hawk_values[i]
      best_pos = hawks[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let e0 = 2.0 * (1.0 - iter)

    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_escape) = RNG::next(current_rng)
      current_rng = rng2

      if abs_val(e0) >= 1.0 {
        // Exploration phase
        let (rng3, r_q) = RNG::next(current_rng)
        current_rng = rng3
        let (rng4, r_rand) = RNG::next(current_rng)
        current_rng = rng4

        if r_q >= 0.5 {
          let (rng5, r_leader) = RNG::next(current_rng)
          current_rng = rng5
          let new_pos: Array[Double] = []
          for d in 0..<dim {
            let pos = bounds.min + r_rand * (bounds.max - bounds.min)
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          let (rng6, r_sub) = RNG::next(current_rng)
          current_rng = rng6
          let final_pos: Array[Double] = []
          for d in 0..<dim {
            let pos = new_pos[d] + r_leader * abs_val(new_pos[d] - hawks[i][d])
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        } else {
          let avg_pos: Array[Double] = []
          for d in 0..<dim {
            let mut sum = 0.0
            for j in 0..<n {
              sum = sum + hawks[j][d]
            }
            avg_pos.push(sum / n.to_double())
          }

          let (rng5, r_ld1) = RNG::next(current_rng)
          current_rng = rng5
          let (rng6, r_ld2) = RNG::next(current_rng)
          current_rng = rng6

          let final_pos: Array[Double] = []
          for d in 0..<dim {
            let pos = (r_ld1 * best_pos[d] - r_ld2 * avg_pos[d])
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        }
      } else {
        // Exploitation phase
        let (rng3, r_escape2) = RNG::next(current_rng)
        current_rng = rng3

        if r_escape2 >= 0.5 {
          let (rng4, r_cf) = RNG::next(current_rng)
          current_rng = rng4

          let final_pos: Array[Double] = []
          for d in 0..<dim {
            let pos = best_pos[d] + abs_val(best_pos[d] - hawks[i][d]) * e0 * r_cf
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        } else {
          let (rng4, r_cf1) = RNG::next(current_rng)
          current_rng = rng4
          let (rng5, r_cf2) = RNG::next(current_rng)
          current_rng = rng5

          let final_pos: Array[Double] = []
          for d in 0..<dim {
            let pos = best_pos[d] * (r_cf1 * e0 - r_cf2)
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        }
      }

      let new_val = evaluate(f, hawks[i])
      hawk_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = hawks[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Flower Pollination Algorithm (FPA) =====

pub struct FPAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  switch_probability: Double
}

pub fn make_fpa_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> FPAConfig {
  FPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    switch_probability: 0.8
  }
}

pub fn flower_pollination_algorithm(config: FPAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let flowers: Array[Array[Double]] = []
  let flower_values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    flowers.push(sol)
    flower_values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if flower_values[i] < best_val {
      best_val = flower_values[i]
      best_pos = flowers[i]
    }
  }

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_switch) = RNG::next(current_rng)
      current_rng = rng2

      let new_pos: Array[Double] = []
      if r_switch < config.switch_probability {
        // Global pollination
        let (rng3, r_levy1) = RNG::next(current_rng)
        current_rng = rng3
        let (rng4, r_levy2) = RNG::next(current_rng)
        current_rng = rng4

        for d in 0..<dim {
          let levy_step = abs_val(r_levy1 * (best_pos[d] - flowers[i][d]))
          let pos = flowers[i][d] + levy_step * (r_levy2 - 0.5) * 2.0
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Local pollination
        let (rng3, r_local) = RNG::next(current_rng)
        current_rng = rng3
        let (rng4, r_idx1) = RNG::next(current_rng)
        current_rng = rng4
        let (rng5, r_idx2) = RNG::next(current_rng)
        current_rng = rng5

        let j = (r_idx1 * n.to_double()).floor().to_int()
        let k = (r_idx2 * n.to_double()).floor().to_int()
        let j_idx = max_int(0, min_int(n - 1, j))
        let k_idx = max_int(0, min_int(n - 1, k))

        for d in 0..<dim {
          let pos = flowers[i][d] + r_local * (flowers[j_idx][d] - flowers[k_idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      flowers[i] = new_pos
      let new_val = evaluate(f, new_pos)
      flower_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Dragonfly Algorithm (DA) =====

pub struct DAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  separation_weight: Double
  alignment_weight: Double
  cohesion_weight: Double
}

pub fn make_da_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> DAConfig {
  DAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    separation_weight: 1.0,
    alignment_weight: 1.0,
    cohesion_weight: 1.0
  }
}

pub fn dragonfly_algorithm(config: DAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let range = bounds.max - bounds.min

  let mut current_rng = rng
  let dragonflies: Array[Array[Double]] = []
  let velocities: Array[Array[Double]] = []
  let df_values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    dragonflies.push(sol)
    df_values.push(evaluate(f, sol))

    let init_vel: Array[Double] = []
    for d in 0..<dim {
      let (rng2, v) = RNG::next(current_rng)
      current_rng = rng2
      init_vel.push((v * 2.0 - 1.0) * range * 0.01)
    }
    velocities.push(init_vel)
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if df_values[i] < best_val {
      best_val = df_values[i]
      best_pos = dragonflies[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Dynamic weights
    let s_weight = config.separation_weight
    let a_weight = config.alignment_weight
    let c_weight = config.cohesion_weight
    let f_weight = 2.0 * (1.0 - iter)
    let e_weight = 1.0 - iter * 0.5

    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng3, r_sep) = RNG::next(current_rng)
      current_rng = rng3
      let (rng4, r_ali) = RNG::next(current_rng)
      current_rng = rng4
      let (rng5, r_coh) = RNG::next(current_rng)
      current_rng = rng5
      let (rng6, r_food) = RNG::next(current_rng)
      current_rng = rng6

      // Find actual neighbors and compute metrics
      let sep: Array[Double] = []
      let ali: Array[Double] = []
      let coh: Array[Double] = []
      for d in 0..<dim {
        let mut sep_sum = 0.0
        let mut ali_sum = 0.0
        let mut coh_sum = 0.0
        let mut count = 0

        for j in 0..<n {
          if j == i { continue }
          let dist = abs_val(dragonflies[i][d] - dragonflies[j][d])
          if dist > 0.0001 {
            sep_sum = sep_sum + (dragonflies[i][d] - dragonflies[j][d]) / dist
          }
          ali_sum = ali_sum + velocities[j][d]
          coh_sum = coh_sum + dragonflies[j][d]
          count = count + 1
        }

        let c = if count > 0 { count.to_double() } else { 1.0 }
        sep.push(sep_sum / c)
        ali.push(ali_sum / c)
        coh.push(coh_sum / c - dragonflies[i][d])
      }

      // Update velocity and position
      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let vel = s_weight * r_sep * sep[d] +
                 a_weight * r_ali * ali[d] +
                 c_weight * r_coh * coh[d] +
                 f_weight * r_food * (best_pos[d] - dragonflies[i][d]) +
                 e_weight * velocities[i][d]

        velocities[i][d] = vel
        let pos = dragonflies[i][d] + vel
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      dragonflies[i] = new_pos
      let new_val = evaluate(f, new_pos)
      df_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Crow Search Algorithm (CSA) =====

pub struct CSAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  flight_length: Double
  awareness_probability: Double
}

pub fn make_csa_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> CSAConfig {
  CSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    flight_length: 2.0,
    awareness_probability: 0.1
  }
}

pub fn crow_search_algorithm(config: CSAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let crows: Array[Array[Double]] = []
  let crow_values: Array[Double] = []
  let memory: Array[Array[Double]] = []
  let memory_values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    crows.push(sol)
    crow_values.push(evaluate(f, sol))
    memory.push(sol)
    memory_values.push(crow_values[i])
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if crow_values[i] < best_val {
      best_val = crow_values[i]
      best_pos = crows[i]
    }
  }

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_awareness) = RNG::next(current_rng)
      current_rng = rng2

      let new_pos: Array[Double] = []
      if r_awareness >= config.awareness_probability {
        // Follow another crow
        let (rng3, r_target) = RNG::next(current_rng)
        current_rng = rng3
        let target = (r_target * n.to_double()).floor().to_int()
        let target_idx = max_int(0, min_int(n - 1, target))

        let (rng4, r_flight) = RNG::next(current_rng)
        current_rng = rng4

        for d in 0..<dim {
          let pos = crows[i][d] + config.flight_length * r_flight * (memory[target_idx][d] - crows[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Random search
        let (rng3, r_rand) = RNG::next(current_rng)
        current_rng = rng3
        for d in 0..<dim {
          let pos = r_rand * (bounds.max - bounds.min) + bounds.min
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      crows[i] = new_pos
      let new_val = evaluate(f, new_pos)
      crow_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }

      // Update memory
      if new_val < memory_values[i] {
        memory[i] = new_pos
        memory_values[i] = new_val
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Salp Swarm Algorithm (SalpSA) =====

pub struct SalpSAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
}

pub fn make_salp_sa_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> SalpSAConfig {
  SalpSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz
  }
}

pub fn salp_swarm_algorithm(config: SalpSAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let salps: Array[Array[Double]] = []
  let salp_values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    salps.push(sol)
    salp_values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  // Find initial best
  let mut best_idx = 0
  for i in 0..<n {
    if salp_values[i] < best_val {
      best_val = salp_values[i]
      best_pos = salps[i]
      best_idx = i
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let c1 = 2.0 * exp_approx(-(4.0 * iter) * (4.0 * iter))

    let new_salps: Array[Array[Double]] = []

    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_c2) = RNG::next(current_rng)
      current_rng = rng2
      let (rng3, r_c3) = RNG::next(current_rng)
      current_rng = rng3

      let new_pos: Array[Double] = []
      if i < n / 2 {
        // Leader salp
        for d in 0..<dim {
          let pos = if r_c3 < 0.5 {
            best_pos[d] + c1 * ((bounds.max - bounds.min) * r_c2 + bounds.min)
          } else {
            best_pos[d] - c1 * ((bounds.max - bounds.min) * r_c2 + bounds.min)
          }
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Follower salp
        let prev = new_salps[i - 1]
        for d in 0..<dim {
          let pos = (prev[d] + salps[i][d]) / 2.0
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      new_salps.push(new_pos)
    }

    // Evaluate and update
    for i in 0..<new_salps.length() {
      if evaluations >= config.max_evaluations { break }
      salps[i] = new_salps[i]
      let new_val = evaluate(f, new_salps[i])
      salp_values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_salps[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Monarch Butterfly Optimization (MBO) =====

pub struct MBOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  migration_ratio: Double
  perturbation: Double
}

pub fn make_mbo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> MBOConfig {
  MBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    migration_ratio: 0.5,
    perturbation: 0.3
  }
}

pub fn monarch_butterfly_optimization(config: MBOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let n1 = (n.to_double() * config.migration_ratio).to_int()
  let n2 = n - n1

  let mut current_rng = rng
  let butterflies1: Array[Array[Double]] = []
  let values1: Array[Double] = []
  let butterflies2: Array[Array[Double]] = []
  let values2: Array[Double] = []

  for i in 0..<n1 {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    butterflies1.push(sol)
    values1.push(evaluate(f, sol))
  }
  for i in 0..<n2 {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    butterflies2.push(sol)
    values2.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n1 {
    if values1[i] < best_val {
      best_val = values1[i]
      best_pos = butterflies1[i]
    }
  }
  for i in 0..<n2 {
    if values2[i] < best_val {
      best_val = values2[i]
      best_pos = butterflies2[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Migration operator (butterfly 1)
    for i in 0..<n1 {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_mig) = RNG::next(current_rng)
      current_rng = rng2

      let new_pos: Array[Double] = []
      if r_mig < config.migration_ratio {
        for d in 0..<dim {
          let (rng3, r_levy) = RNG::next(current_rng)
          current_rng = rng3
          let pos = best_pos[d] + r_levy * (butterflies1[i][d] - best_pos[d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        let (rng3, r_rand) = RNG::next(current_rng)
        current_rng = rng3
        for d in 0..<dim {
          let pos = r_rand * (bounds.max - bounds.min) + bounds.min
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      butterflies1[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values1[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Adjustment operator (butterfly 2)
    for i in 0..<n2 {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_adj) = RNG::next(current_rng)
      current_rng = rng2

      let new_pos: Array[Double] = []
      if r_adj < config.perturbation {
        for d in 0..<dim {
          let (rng3, r_mut) = RNG::next(current_rng)
          current_rng = rng3
          let pos = butterflies2[i][d] + (r_mut - 0.5) * (bounds.max - bounds.min) * 0.1
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        let (rng3, r_swap) = RNG::next(current_rng)
        current_rng = rng3
        let partner = (r_swap * n1.to_double()).floor().to_int()
        let p_idx = max_int(0, min_int(n1 - 1, partner))
        for d in 0..<dim {
          let pos = butterflies2[i][d] + iter * (butterflies1[p_idx][d] - butterflies2[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      butterflies2[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values2[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Artificial Ecosystem-based Optimization (AEO) =====

pub struct AEOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
}

pub fn make_aeo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> AEOConfig {
  AEOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz
  }
}

pub fn artificial_ecosystem_optimization(config: AEOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let individuals: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    individuals.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = individuals[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Sort individuals by fitness (ascending)
    let sorted_indices: Array[Int] = []
    for i in 0..<n {
      sorted_indices.push(i)
    }
    for i in 0..<n {
      for j in (i + 1)..<n {
        if values[sorted_indices[j]] < values[sorted_indices[i]] {
          let temp = sorted_indices[i]
          sorted_indices[i] = sorted_indices[j]
          sorted_indices[j] = temp
        }
      }
    }

    // Update positions
    for rank in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let idx = sorted_indices[rank]
      let (rng2, r_role) = RNG::next(current_rng)
      current_rng = rng2

      let new_pos: Array[Double] = []

      if rank == 0 {
        // Producer (best) - exploit
        for d in 0..<dim {
          let (rng3, r_prod) = RNG::next(current_rng)
          current_rng = rng3
          let pos = individuals[idx][d] + r_prod * (best_pos[d] - individuals[idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else if r_role < 0.5 {
        // Consumer - explore toward better solution
        let better_idx = sorted_indices[rank - 1]
        for d in 0..<dim {
          let (rng3, r_cons) = RNG::next(current_rng)
          current_rng = rng3
          let pos = individuals[idx][d] + r_cons * (individuals[better_idx][d] - individuals[idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Decomposer - random search
        for d in 0..<dim {
          let (rng3, r_decomp) = RNG::next(current_rng)
          current_rng = rng3
          let pos = individuals[idx][d] + (r_decomp - 0.5) * (bounds.max - bounds.min) * (1.0 - iter)
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      individuals[idx] = new_pos
      let new_val = evaluate(f, new_pos)
      values[idx] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Invasive Weed Optimization (IWO) =====

pub struct IWOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  initial_population: Int
  max_population: Int
  min_seeds: Int
  max_seeds: Int
  sigma_initial: Double
  sigma_final: Double
  nonlinear: Double
}

pub fn make_iwo_config(max_evals: Int, dim: Int, b: Bounds, init_pop: Int, max_pop: Int, min_seed: Int, max_seed: Int, sig_init: Double, sig_final: Double, n_nonlinear: Double) -> IWOConfig {
  IWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    initial_population: init_pop,
    max_population: max_pop,
    min_seeds: min_seed,
    max_seeds: max_seed,
    sigma_initial: sig_init,
    sigma_final: sig_final,
    nonlinear: n_nonlinear
  }
}

pub fn invasive_weed_optimization(config: IWOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let mut current_rng = rng
  let population: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<config.initial_population {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    population.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = config.initial_population
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<population.length() {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = population[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter_frac = evaluations.to_double() / config.max_evaluations.to_double()
    let sigma = config.sigma_final + (config.sigma_initial - config.sigma_final) * pow_approx(1.0 - iter_frac, config.nonlinear)

    let new_pop: Array[Array[Double]] = []
    let new_vals: Array[Double] = []

    let sorted_idx = sort_indices(values)
    let n = population.length()

    for rank in 0..<n {
      if evaluations + new_pop.length() >= config.max_evaluations { break }

      let idx = sorted_idx[rank]
      let norm_rank = (rank.to_double() + 1.0) / n.to_double()
      let seeds_d: Double = config.min_seeds.to_double() + norm_rank * (config.max_seeds.to_double() - config.min_seeds.to_double())
      let seeds = seeds_d.round().to_int()

      for s in 0..<seeds {
        if evaluations + new_pop.length() >= config.max_evaluations { break }

        let child: Array[Double] = []
        for d in 0..<dim {
          let (rng2, r1) = RNG::next(current_rng)
          current_rng = rng2
          let (rng3, r2) = RNG::next(current_rng)
          current_rng = rng3
          let gauss = box_muller(r1, r2)
          let pos = population[idx][d] + sigma * gauss
          child.push(max(bounds.min, min(bounds.max, pos)))
        }

        new_pop.push(child)
        new_vals.push(evaluate(f, child))
      }
    }

    // Combine and keep best max_population
    let combined: Array[Array[Double]] = []
    let combined_vals: Array[Double] = []
    for i in 0..<population.length() {
      combined.push(population[i])
      combined_vals.push(values[i])
    }
    for i in 0..<new_pop.length() {
      combined.push(new_pop[i])
      combined_vals.push(new_vals[i])
    }

    let combined_idx = sort_indices(combined_vals)
    let keep = min_int(config.max_population, combined.length())

    population.clear()
    values.clear()
    for i in 0..<keep {
      let idx = combined_idx[i]
      population.push(combined[idx])
      values.push(combined_vals[idx])
    }

    evaluations = evaluations + new_pop.length()

    for i in 0..<new_pop.length() {
      if new_vals[i] < best_val {
        best_val = new_vals[i]
        best_pos = new_pop[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Bacteria Foraging Optimization (BFO) =====

pub struct BFOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  chemotactic_steps: Int
  swim_length: Int
  tumble_step: Double
}

pub fn make_bfo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, chem_steps: Int, swim_len: Int, tumble: Double) -> BFOConfig {
  BFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    chemotactic_steps: chem_steps,
    swim_length: swim_len,
    tumble_step: tumble
  }
}

pub fn bacteria_foraging_optimization(config: BFOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let bacteria: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    bacteria.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = bacteria[i]
    }
  }

  while evaluations < config.max_evaluations {
    for b_idx in 0..<n {
      if evaluations >= config.max_evaluations { break }

      // Tumble: generate random direction
      let delta: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r1) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r2) = RNG::next(current_rng)
        current_rng = rng3
        let g = box_muller(r1, r2)
        delta.push(g)
      }
      // Normalize
      let mut norm_sq = 0.0
      for d in 0..<dim {
        norm_sq = norm_sq + delta[d] * delta[d]
      }
      let norm = if norm_sq.sqrt() > 1.0e-10 { norm_sq.sqrt() } else { 1.0 }

      // Swim in that direction
      for s in 0..<config.swim_length {
        if evaluations >= config.max_evaluations { break }

        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let pos = bacteria[b_idx][d] + config.tumble_step * delta[d] / norm
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1

        if new_val < values[b_idx] {
          bacteria[b_idx] = new_pos
          values[b_idx] = new_val

          if new_val < best_val {
            best_val = new_val
            best_pos = new_pos
          }
        } else {
          break
        }
      }
    }

    // Elimination: replace worst with random
    if evaluations < config.max_evaluations {
      let idx = sort_indices(values)
      for i in 0..<(n / 5) {
        if evaluations >= config.max_evaluations { break }
        let worst_idx = idx[n - 1 - i]
        let (rng4, sol) = random_solution(dim, bounds, current_rng)
        current_rng = rng4
        bacteria[worst_idx] = sol
        values[worst_idx] = evaluate(f, sol)
        evaluations = evaluations + 1
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Shuffled Frog Leaping Algorithm (SFLA) =====

pub struct SFLAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  memeplexes: Int
  sub_iterations: Int
}

pub fn make_sfla_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, memes: Int, sub_iter: Int) -> SFLAConfig {
  SFLAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    memeplexes: memes,
    sub_iterations: sub_iter
  }
}

pub fn shuffled_frog_leaping_algorithm(config: SFLAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let memes = config.memeplexes
  let memeplex_sz = n / memes

  let mut current_rng = rng
  let frogs: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    frogs.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = frogs[i]
    }
  }

  while evaluations < config.max_evaluations {
    // Sort frogs by fitness
    let sorted_idx = sort_indices(values)

    // Divide into memeplexes (interleaved)
    for m in 0..<memes {
      if evaluations >= config.max_evaluations { break }

      let memeplex: Array[Int] = []
      for i in 0..<memeplex_sz {
        let idx = m + i * memes
        if idx < n {
          memeplex.push(sorted_idx[idx])
        }
      }

      // Get best and worst in memeplex
      let mut best_idx = memeplex[0]
      let mut worst_idx = memeplex[0]
      for i in 0..<memeplex.length() {
        let fi = memeplex[i]
        if values[fi] < values[best_idx] {
          best_idx = fi
        }
        if values[fi] > values[worst_idx] {
          worst_idx = fi
        }
      }

      let sub_best = frogs[best_idx]
      let sub_worst = frogs[worst_idx]

      // Update worst frog toward best
      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r) = RNG::next(current_rng)
        current_rng = rng2
        let pos = sub_worst[d] + r * (sub_best[d] - sub_worst[d])
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < values[worst_idx] {
        frogs[worst_idx] = new_pos
        values[worst_idx] = new_val
      } else {
        // Try global best
        let new_pos2: Array[Double] = []
        for d in 0..<dim {
          let (rng3, r) = RNG::next(current_rng)
          current_rng = rng3
          let pos = sub_worst[d] + r * (best_pos[d] - sub_worst[d])
          new_pos2.push(max(bounds.min, min(bounds.max, pos)))
        }
        let new_val2 = evaluate(f, new_pos2)
        evaluations = evaluations + 1
        if new_val2 < values[worst_idx] {
          frogs[worst_idx] = new_pos2
          values[worst_idx] = new_val2
        } else {
          // Random replacement
          let (rng4, sol) = random_solution(dim, bounds, current_rng)
          current_rng = rng4
          frogs[worst_idx] = sol
          values[worst_idx] = evaluate(f, sol)
          evaluations = evaluations + 1
        }
      }

      if values[worst_idx] < best_val {
        best_val = values[worst_idx]
        best_pos = frogs[worst_idx]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Glowworm Swarm Optimization (GSO) =====

pub struct GSOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  luciferin_decay: Double
  luciferin_enhance: Double
  neighbor_threshold: Int
  initial_luciferin: Double
  initial_radius: Double
}

pub fn make_gso_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, luc_decay: Double, luc_enhance: Double, neigh_thresh: Int, init_luc: Double, init_rad: Double) -> GSOConfig {
  GSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    luciferin_decay: luc_decay,
    luciferin_enhance: luc_enhance,
    neighbor_threshold: neigh_thresh,
    initial_luciferin: init_luc,
    initial_radius: init_rad
  }
}

pub fn glowworm_swarm_optimization(config: GSOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let worms: Array[Array[Double]] = []
  let values: Array[Double] = []
  let luciferin: Array[Double] = []
  let radius: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    worms.push(sol)
    values.push(evaluate(f, sol))
    luciferin.push(config.initial_luciferin)
    radius.push(config.initial_radius)
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = worms[i]
    }
  }

  while evaluations < config.max_evaluations {
    // Update luciferin
    for i in 0..<n {
      luciferin[i] = (1.0 - config.luciferin_decay) * luciferin[i] + config.luciferin_enhance * (1.0 / (values[i] + 1.0))
    }

    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      // Find neighbors with higher luciferin within radius
      let neighbors: Array[Int] = []
      for j in 0..<n {
        if j == i { continue }
        let dist = euclidean_dist(worms[i], worms[j])
        if dist <= radius[i] && luciferin[j] > luciferin[i] {
          neighbors.push(j)
        }
      }

      if neighbors.length() == 0 {
        // Random movement
        let (rng2, r) = RNG::next(current_rng)
        current_rng = rng2
        let step = radius[i] * 0.3
        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let (rng3, rd) = RNG::next(current_rng)
          current_rng = rng3
          let pos = worms[i][d] + (rd - 0.5) * step
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
        worms[i] = new_pos
        let new_val = evaluate(f, new_pos)
        values[i] = new_val
        evaluations = evaluations + 1
        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      } else {
        // Probabilistically move toward a neighbor
        // Compute selection probabilities
        let weights: Array[Double] = []
        let mut total = 0.0
        for j in 0..<neighbors.length() {
          let w = luciferin[neighbors[j]] - luciferin[i]
          weights.push(w)
          total = total + w
        }

        let (rng4, sel) = RNG::next(current_rng)
        current_rng = rng4
        let mut cum = 0.0
        let mut chosen = neighbors[0]
        for j in 0..<neighbors.length() {
          cum = cum + weights[j] / total
          if sel <= cum {
            chosen = neighbors[j]
            break
          }
        }

        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let (rng5, rd) = RNG::next(current_rng)
          current_rng = rng5
          let step = radius[i] * rd
          let pos = worms[i][d] + step * (worms[chosen][d] - worms[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        worms[i] = new_pos
        let new_val = evaluate(f, new_pos)
        values[i] = new_val
        evaluations = evaluations + 1

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }

        // Update radius
        if neighbors.length() > config.neighbor_threshold {
          radius[i] = radius[i] * 0.9
        } else if neighbors.length() < config.neighbor_threshold {
          radius[i] = radius[i] * 1.1
        }
        if radius[i] > (bounds.max - bounds.min) * 0.5 {
          radius[i] = (bounds.max - bounds.min) * 0.5
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

fn euclidean_dist(a: Array[Double], b: Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<a.length() {
    let diff = a[i] - b[i]
    sum = sum + diff * diff
  }
  sum.sqrt()
}

// ===== Gravitational Search Algorithm (GSA) =====

pub struct GSAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  gravity_constant: Double
  decay_alpha: Double
}

pub fn make_gsa_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, g0: Double, alpha: Double) -> GSAConfig {
  GSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    gravity_constant: g0,
    decay_alpha: alpha
  }
}

pub fn gravitational_search_algorithm(config: GSAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let agents: Array[Array[Double]] = []
  let velocities: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    agents.push(sol)
    values.push(evaluate(f, sol))

    let vel: Array[Double] = []
    for d in 0..<dim {
      let (rng2, r) = RNG::next(current_rng)
      current_rng = rng2
      vel.push(0.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308
  let mut worst_val = -1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = agents[i]
    }
    if values[i] > worst_val {
      worst_val = values[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let g = config.gravity_constant * pow_approx(config.decay_alpha, iter)

    // Calculate masses (fitness normalized)
    let masses: Array[Double] = []
    for i in 0..<n {
      let m = (values[i] - worst_val) / (best_val - worst_val + 1.0e-30)
      masses.push(m + 1.0)
    }
    let mut total_mass = 0.0
    for i in 0..<n { total_mass = total_mass + masses[i] }

    // Update positions
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      // Calculate gravitational force on agent i from all other agents
      let force: Array[Double] = []
      for d in 0..<dim {
        let mut f_d = 0.0
        for j in 0..<n {
          if j == i { continue }
          let dist = euclidean_dist(agents[i], agents[j])
          if dist < 1.0e-10 { continue }
          let (rng3, r_force) = RNG::next(current_rng)
          current_rng = rng3
          f_d = f_d + r_force * g * masses[j] * (agents[j][d] - agents[i][d]) / (dist + 1.0e-10)
        }
        force.push(f_d)
      }

      // Update velocity and position
      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng4, r_vel) = RNG::next(current_rng)
        current_rng = rng4
        let vel = r_vel * velocities[i][d] + force[d] / (masses[i] + 1.0e-30)
        velocities[i][d] = vel
        let (rng5, r_pos) = RNG::next(current_rng)
        current_rng = rng5
        let pos = agents[i][d] + r_pos * velocities[i][d]
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      agents[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
      if new_val > worst_val {
        worst_val = new_val
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Teaching-Learning-Based Optimization (TLBO) =====

pub struct TLBOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
}

pub fn make_tlbo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> TLBOConfig {
  TLBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz
  }
}

pub fn teaching_learning_based_optimization(config: TLBOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let students: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    students.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = students[i]
    }
  }

  while evaluations < config.max_evaluations {
    // Find best (teacher) and worst student
    let mut teacher_idx = 0
    let mut worst_idx = 0
    for i in 1..<n {
      if values[i] < values[teacher_idx] { teacher_idx = i }
      if values[i] > values[worst_idx] { worst_idx = i }
    }
    let teacher = students[teacher_idx]

    // Calculate mean of all students
    let mean: Array[Double] = []
    for d in 0..<dim {
      let mut sum = 0.0
      for i in 0..<n { sum = sum + students[i][d] }
      mean.push(sum / n.to_double())
    }

    // Teacher phase
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }
      if i == teacher_idx { continue }

      let (rng2, tf_val) = RNG::next(current_rng)
      current_rng = rng2
      let tf = 1.0 + (if tf_val < 0.5 { 0.0 } else { 1.0 })
      let (rng3, r1) = RNG::next(current_rng)
      current_rng = rng3

      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let pos = students[i][d] + r1 * (teacher[d] - tf * mean[d])
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        students[i] = new_pos
        values[i] = new_val
      }

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Learner phase
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng4, partner_r) = RNG::next(current_rng)
      current_rng = rng4
      let partner = (partner_r * n.to_double()).floor().to_int() % n
      let partner_idx = if partner < 0 { n + partner } else { partner }

      let (rng5, r1) = RNG::next(current_rng)
      current_rng = rng5
      let new_pos: Array[Double] = []

      if values[i] < values[partner_idx] {
        for d in 0..<dim {
          let pos = students[i][d] + r1 * (students[i][d] - students[partner_idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        for d in 0..<dim {
          let pos = students[i][d] + r1 * (students[partner_idx][d] - students[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        students[i] = new_pos
        values[i] = new_val
      }

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Animal Migration Optimization (AMO) =====

pub struct AMOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  migration_interval: Int
}

pub fn make_amo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int, mig_interval: Int) -> AMOConfig {
  AMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    migration_interval: mig_interval
  }
}

pub fn animal_migration_optimization(config: AMOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let animals: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    animals.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = animals[i]
    }
  }

  let mut migration_count = 0

  while evaluations < config.max_evaluations {
    migration_count = migration_count + 1

    // Update position: follow best or move randomly
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng2, r_follow) = RNG::next(current_rng)
      current_rng = rng2
      let new_pos: Array[Double] = []

      if r_follow < 0.7 {
        // Follow best
        let (rng3, r1) = RNG::next(current_rng)
        current_rng = rng3
        for d in 0..<dim {
          let pos = animals[i][d] + r1 * (best_pos[d] - animals[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Explore randomly toward a random neighbor
        let (rng4, r_idx) = RNG::next(current_rng)
        current_rng = rng4
        let j = (r_idx * n.to_double()).floor().to_int() % n
        let j_idx = if j < 0 { n + j } else { j }
        let (rng5, r2) = RNG::next(current_rng)
        current_rng = rng5
        for d in 0..<dim {
          let pos = animals[i][d] + r2 * (animals[j_idx][d] - animals[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      animals[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Migration: replace worst individuals with new random ones
    if migration_count % config.migration_interval == 0 && evaluations < config.max_evaluations {
      let idx = sort_indices(values)
      let replace_count = n / 5
      for k in 0..<replace_count {
        if evaluations >= config.max_evaluations { break }
        let worst_idx = idx[n - 1 - k]
        let (rng6, sol) = random_solution(dim, bounds, current_rng)
        current_rng = rng6
        animals[worst_idx] = sol
        values[worst_idx] = evaluate(f, sol)
        evaluations = evaluations + 1
        if values[worst_idx] < best_val {
          best_val = values[worst_idx]
          best_pos = animals[worst_idx]
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

// ===== Multi-Verse Optimizer (MVO) =====

pub struct MVOConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  wormhole_exploitation: Double
  travelling_distance_rate: Double
}

pub fn make_mvo_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> MVOConfig {
  MVOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    wormhole_exploitation: 1.0,
    travelling_distance_rate: 0.6
  }
}

pub fn multi_verse_optimizer(config: MVOConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let universes: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    universes.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308
  let mut worst_val = -1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = universes[i]
    }
    if values[i] > worst_val {
      worst_val = values[i]
    }
  }

  let mut sorted_idx = sort_indices(values)

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Update wormhole existence probability
    let wep = 1.0 + iter
    // Update travelling distance rate
    let tdr = 1.0 - iter.abs()

    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng2, r1) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r_wep) = RNG::next(current_rng)
        current_rng = rng3
        let (rng4, r_tdr) = RNG::next(current_rng)
        current_rng = rng4
        let (rng5, r_idx) = RNG::next(current_rng)
        current_rng = rng5

        let mut pos: Double = 0.0
        if r_wep < wep {
          // White/black hole tunnel
          let j = (r_idx * n.to_double()).floor().to_int() % n
          let j_idx = if j < 0 { n + j } else { j }
          if values[i] < values[j_idx] {
            pos = universes[i][d] + tdr * (bounds.max - bounds.min) * r1
          } else {
            pos = universes[j_idx][d] + tdr * (bounds.max - bounds.min) * r1
          }
        } else {
          // Wormhole
          if r_tdr < 0.5 {
            pos = best_pos[d] + (bounds.max - bounds.min) * r1
          } else {
            pos = best_pos[d] - (bounds.max - bounds.min) * r1
          }
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      universes[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    sorted_idx = sort_indices(values)
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Krill Herd (KH) =====

pub struct KHConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  max_speed: Double
  foraging_speed: Double
}

pub fn make_kh_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> KHConfig {
  KHConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    max_speed: 0.01,
    foraging_speed: 0.02
  }
}

pub fn krill_herd(config: KHConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let mut current_rng = rng
  let krills: Array[Array[Double]] = []
  let velocities: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    krills.push(sol)
    values.push(evaluate(f, sol))

    let vel: Array[Double] = []
    for d in 0..<dim {
      let (rng2, r) = RNG::next(current_rng)
      current_rng = rng2
      vel.push(0.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308
  let mut worst_val = -1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = krills[i]
    }
    if values[i] > worst_val {
      worst_val = values[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let step_size = config.max_speed * (1.0 - iter)

    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      // Find neighbors
      let neighbors: Array[Int] = []
      let max_neighbors = 2
      for j in 0..<n {
        if j == i { continue }
        if neighbors.length() >= max_neighbors { break }
        let dist = euclidean_dist(krills[i], krills[j])
        if dist < 1.0 {
          neighbors.push(j)
        }
      }

      // Direction of motion (induced + foraging)
      let motion: Array[Double] = []
      for d in 0..<dim {
        // Induced motion: follow neighbors
        let mut induced = 0.0
        for ni in 0..<neighbors.length() {
          let j = neighbors[ni]
          let dist = euclidean_dist(krills[i], krills[j])
          if dist > 1.0e-10 {
            induced = induced + (krills[j][d] - krills[i][d]) / dist
          }
        }
        // Foraging motion: toward food (best)
        let food = best_pos[d] - krills[i][d]
        // Random diffusion
        let (rng3, r_diff) = RNG::next(current_rng)
        current_rng = rng3

        let dir = induced * config.foraging_speed + food * 0.5 + r_diff * step_size
        motion.push(dir)
      }

      let new_pos: Array[Double] = []
      for d in 0..<dim {
        let (rng4, r1) = RNG::next(current_rng)
        current_rng = rng4
        let vel = velocities[i][d] + step_size * motion[d]
        velocities[i][d] = vel
        let pos = krills[i][d] + r1 * velocities[i][d]
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      krills[i] = new_pos
      let new_val = evaluate(f, new_pos)
      values[i] = new_val
      evaluations = evaluations + 1

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations: evaluations
  }
}

// ===== Imperialist Competitive Algorithm (ICA) =====

pub struct ICAConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  imperialist_count: Int
  assimilation_rate: Double
  revolution_rate: Double
}

pub fn make_ica_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> ICAConfig {
  ICAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    imperialist_count: 3,
    assimilation_rate: 0.3,
    revolution_rate: 0.3
  }
}

pub fn imperialist_competitive_algorithm(config: ICAConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let imp_count = config.imperialist_count

  let mut current_rng = rng
  let countries: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    countries.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  // Sort and assign imperialists and colonies
  let mut sorted_idx = sort_indices(values)

  while evaluations < config.max_evaluations {
    // Re-sort
    sorted_idx = sort_indices(values)

    // Imperialists: best imp_count countries
    let imperialists: Array[Int] = []
    for i in 0..<imp_count {
      imperialists.push(sorted_idx[i])
    }

    // Assign colonies to imperialists based on cost
    let colonies: Array[Array[Int]] = []
    for i in 0..<imp_count {
      let col: Array[Int] = []
      colonies.push(col)
    }

    let total_imperialist_cost: Double = {
      let mut s = 0.0
      for i in 0..<imp_count {
        let idx = imperialists[i]
        s = s + (values[idx] - values[sorted_idx[imp_count]]) + 1.0e-10
      }
      s
    }

    for i in imp_count..<n {
      let idx = sorted_idx[i]
      // Assign to weakest imperialist (simplified: round-robin by probability)
      let (rng2, r_assign) = RNG::next(current_rng)
      current_rng = rng2
      let imp_idx = (r_assign * imp_count.to_double()).floor().to_int()
      let target = max_int(0, min_int(imp_count - 1, imp_idx))
      colonies[target].push(idx)
    }

    // Assimilation: colonies move toward their imperialist
    for k in 0..<imp_count {
      if evaluations >= config.max_evaluations { break }
      let imp_idx = imperialists[k]
      let imp_pos = countries[imp_idx]

      for c in 0..<colonies[k].length() {
        if evaluations >= config.max_evaluations { break }
        let col_idx = colonies[k][c]

        let (rng3, r_assim) = RNG::next(current_rng)
        current_rng = rng3

        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let pos = countries[col_idx][d] + r_assim * config.assimilation_rate * (imp_pos[d] - countries[col_idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        countries[col_idx] = new_pos
        let new_val = evaluate(f, new_pos)
        values[col_idx] = new_val
        evaluations = evaluations + 1

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }

        if values[col_idx] < values[imp_idx] {
          // Revolution: colony becomes new imperialist
          let temp = countries[col_idx]
          countries[col_idx] = countries[imp_idx]
          countries[imp_idx] = temp
          let temp_val = values[col_idx]
          values[col_idx] = values[imp_idx]
          values[imp_idx] = temp_val
        }
      }
    }

    // Revolution: random perturbation
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }
      let (rng4, r_rev) = RNG::next(current_rng)
      current_rng = rng4
      if r_rev < config.revolution_rate {
        let dims_to_change = max_int(1, dim / 3)
        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let (rng5, r_p) = RNG::next(current_rng)
          current_rng = rng5
          if d < dims_to_change {
            new_pos.push(bounds.min + r_p * (bounds.max - bounds.min))
          } else {
            new_pos.push(countries[i][d])
          }
        }
        countries[i] = new_pos
        let new_val = evaluate(f, new_pos)
        values[i] = new_val
        evaluations = evaluations + 1

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
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

// ===== Stochastic Fractal Search (SFS) =====

pub struct SFSConfig {
  max_evaluations: Int
  dimension: Int
  bounds: Bounds
  population_size: Int
  discovery_iterations: Int
  learning_iterations: Int
}

pub fn make_sfs_config(max_evals: Int, dim: Int, b: Bounds, pop_sz: Int) -> SFSConfig {
  SFSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    discovery_iterations: 3,
    learning_iterations: 2
  }
}

pub fn stochastic_fractal_search(config: SFSConfig, f: (Array[Double]) -> Double, rng: RNG) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let range = bounds.max - bounds.min

  let mut current_rng = rng
  let points: Array[Array[Double]] = []
  let values: Array[Double] = []

  for i in 0..<n {
    let (rng1, sol) = random_solution(dim, bounds, current_rng)
    current_rng = rng1
    points.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos: Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = points[i]
    }
  }

  while evaluations < config.max_evaluations {
    // Discovery phase: create new points around existing ones
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      for d in 0..<dim {
        if evaluations >= config.max_evaluations { break }

        let sigma = range * 0.1
        let (rng2, r_gauss) = RNG::next(current_rng)
        current_rng = rng2
        let (rng3, r_gauss2) = RNG::next(current_rng)
        current_rng = rng3
        let gauss_val = box_muller(r_gauss, r_gauss2)

        // Explore nearby
        let pos1 = points[i][d] + sigma * gauss_val
        let new_pos: Array[Double] = []
        for dd in 0..<dim {
          if dd == d {
            new_pos.push(max(bounds.min, min(bounds.max, pos1)))
          } else {
            new_pos.push(points[i][dd])
          }
        }

        if evaluations < config.max_evaluations {
          let new_val = evaluate(f, new_pos)
          evaluations = evaluations + 1
          if new_val < best_val {
            best_val = new_val
            best_pos = new_pos
          }
          if new_val < values[i] {
            points[i] = new_pos
            values[i] = new_val
          }
        }

        if evaluations >= config.max_evaluations { break }

        // Explore in opposite direction
        let pos2 = points[i][d] - sigma * gauss_val
        let new_pos2: Array[Double] = []
        for dd in 0..<dim {
          if dd == d {
            new_pos2.push(max(bounds.min, min(bounds.max, pos2)))
          } else {
            new_pos2.push(points[i][dd])
          }
        }

        if evaluations < config.max_evaluations {
          let new_val2 = evaluate(f, new_pos2)
          evaluations = evaluations + 1
          if new_val2 < best_val {
            best_val = new_val2
            best_pos = new_pos2
          }
          if new_val2 < values[i] {
            points[i] = new_pos2
            values[i] = new_val2
          }
        }
      }
    }

    // Learning phase: points move toward better neighbors
    for i in 0..<n {
      if evaluations >= config.max_evaluations { break }

      let (rng4, r_learn) = RNG::next(current_rng)
      current_rng = rng4

      // Find a better point
      let mut better_idx = -1
      let mut min_dist = 1.0e308
      for j in 0..<n {
        if j == i { continue }
        if values[j] < values[i] {
          let dist = euclidean_dist(points[i], points[j])
          if dist < min_dist {
            min_dist = dist
            better_idx = j
          }
        }
      }

      if better_idx >= 0 && min_dist < 1.0e308 {
        let new_pos: Array[Double] = []
        for d in 0..<dim {
          let pos = points[i][d] + r_learn * (points[better_idx][d] - points[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
        points[i] = new_pos
        let new_val = evaluate(f, new_pos)
        values[i] = new_val
        evaluations = evaluations + 1

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
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

fn box_muller(u1: Double, u2: Double) -> Double {
  let pi = 3.14159265358979323846
  let sqrt_val = (-2.0 * ln_approx(u1)).sqrt()
  let cos_val = cos_approx(2.0 * pi * u2)
  sqrt_val * cos_val
}

fn ln_approx(x: Double) -> Double {
  let val = if x > 0.0 { x } else { 0.0001 }
  // Use the series: ln(x) = 2 * sum_{k=0}^{inf} (1/(2k+1)) * ((x-1)/(x+1))^(2k+1)
  let y = (val - 1.0) / (val + 1.0)
  let y2 = y * y
  let mut result = y
  let mut term = y
  for k in 1..<20 {
    term = term * y2
    result = result + term / (2 * k + 1).to_double()
  }
  2.0 * result
}

fn sort_indices(values: Array[Double]) -> Array[Int] {
  let n = values.length()
  let indices: Array[Int] = []
  for i in 0..<n {
    indices.push(i)
  }

  for i in 0..<n {
    for j in (i + 1)..<n {
      if values[indices[j]] < values[indices[i]] {
        let temp = indices[i]
        indices[i] = indices[j]
        indices[j] = temp
      }
    }
  }

  indices
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

pub fn fmt_fixed(x: Double, decimals: Int) -> String {
  let neg = x < 0.0
  let abs_x = if neg { -x } else { x }
  let factor = 1.0
  let f = make_factor(factor, decimals)
  let scaled = (abs_x * f).round()
  let int_part = (scaled / f).floor()
  let frac_part = (scaled - int_part * f).round()
  let frac_str = "\{frac_part}"
  let padded = if frac_str.length() < decimals {
    pad_zeros(frac_str, decimals - frac_str.length())
  } else {
    frac_str
  }
  let int_str = "\{int_part}"
  if neg {
    "-\{int_str}.\{padded}"
  } else {
    "\{int_str}.\{padded}"
  }
}

fn make_factor(acc: Double, n: Int) -> Double {
  if n == 0 { acc } else { make_factor(acc * 10.0, n - 1) }
}

fn pad_zeros(s: String, count: Int) -> String {
  if count == 0 { s } else { pad_zeros("0" + s, count - 1) }
}

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

  let rs_config = make_random_search_config(max_evals, dim, bounds)
  let rs_result = random_search(rs_config, f, rng)
  let rs_error = calc_relative_error(rs_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Random Search",
    best_value: rs_result.best_value,
    evaluations: rs_result.evaluations,
    relative_error: rs_error
  })

  let hc_config = make_hill_climbing_config(max_evals, dim, bounds, 0.5, 10)
  let hc_result = hill_climbing(hc_config, f, rng)
  let hc_error = calc_relative_error(hc_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Hill Climbing",
    best_value: hc_result.best_value,
    evaluations: hc_result.evaluations,
    relative_error: hc_error
  })

  let sa_config = make_simulated_annealing_config(max_evals, dim, bounds, 100.0, 0.995, 0.5)
  let sa_result = simulated_annealing(sa_config, f, rng)
  let sa_error = calc_relative_error(sa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Simulated Annealing",
    best_value: sa_result.best_value,
    evaluations: sa_result.evaluations,
    relative_error: sa_error
  })

  let pso_config = make_pso_config(max_evals, dim, bounds, 30, 0.7, 1.5, 1.5)
  let pso_result = pso(pso_config, f, rng)
  let pso_error = calc_relative_error(pso_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "PSO",
    best_value: pso_result.best_value,
    evaluations: pso_result.evaluations,
    relative_error: pso_error
  })

  let de_config = make_de_config(max_evals, dim, bounds, 30, 0.8, 0.9)
  let de_result = differential_evolution(de_config, f, rng)
  let de_error = calc_relative_error(de_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Differential Evolution",
    best_value: de_result.best_value,
    evaluations: de_result.evaluations,
    relative_error: de_error
  })

  let ga_config = make_ga_config(max_evals, dim, bounds, 50, 0.01, 0.8)
  let ga_result = genetic_algorithm(ga_config, f, rng)
  let ga_error = calc_relative_error(ga_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Genetic Algorithm",
    best_value: ga_result.best_value,
    evaluations: ga_result.evaluations,
    relative_error: ga_error
  })

  let nm_config = make_nm_config(max_evals, dim, bounds, 1.0, 2.0, 0.5, 0.5)
  let nm_result = nelder_mead(nm_config, f, rng)
  let nm_error = calc_relative_error(nm_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Nelder-Mead",
    best_value: nm_result.best_value,
    evaluations: nm_result.evaluations,
    relative_error: nm_error
  })

  let ff_config = make_firefly_config(max_evals, dim, bounds, 25, 0.5, 1.0, 1.0)
  let ff_result = firefly_algorithm(ff_config, f, rng)
  let ff_error = calc_relative_error(ff_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Firefly",
    best_value: ff_result.best_value,
    evaluations: ff_result.evaluations,
    relative_error: ff_error
  })

  let aco_config = make_aco_config(max_evals, dim, bounds, 20, 0.1, 1.0, 10.0)
  let aco_result = aco(aco_config, f, rng)
  let aco_error = calc_relative_error(aco_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "ACO",
    best_value: aco_result.best_value,
    evaluations: aco_result.evaluations,
    relative_error: aco_error
  })

  let bc_config = make_bee_colony_config(max_evals, dim, bounds, 25, 100)
  let bc_result = bee_colony(bc_config, f, rng)
  let bc_error = calc_relative_error(bc_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Bee Colony",
    best_value: bc_result.best_value,
    evaluations: bc_result.evaluations,
    relative_error: bc_error
  })

  let ts_config = make_tabu_search_config(max_evals, dim, bounds, 50, 30, 0.5)
  let ts_result = tabu_search(ts_config, f, rng)
  let ts_error = calc_relative_error(ts_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Tabu Search",
    best_value: ts_result.best_value,
    evaluations: ts_result.evaluations,
    relative_error: ts_error
  })

  let ls_config = make_local_search_config(max_evals, dim, bounds, 30, 0.5, 0.5)
  let ls_result = local_search(ls_config, f, rng)
  let ls_error = calc_relative_error(ls_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Local Search",
    best_value: ls_result.best_value,
    evaluations: ls_result.evaluations,
    relative_error: ls_error
  })

  let cma_config = make_cmaes_config(max_evals, dim, bounds, 0.5, 20, 10)
  let cma_result = cma_es(cma_config, f, rng)
  let cma_error = calc_relative_error(cma_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "CMA-ES",
    best_value: cma_result.best_value,
    evaluations: cma_result.evaluations,
    relative_error: cma_error
  })

  let cem_config = make_cem_config(max_evals, dim, bounds, 100, 0.1, 0.7)
  let cem_result = cross_entropy_method(cem_config, f, rng)
  let cem_error = calc_relative_error(cem_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "CEM",
    best_value: cem_result.best_value,
    evaluations: cem_result.evaluations,
    relative_error: cem_error
  })

  let cs_config = make_cuckoo_config(max_evals, dim, bounds, 25, 0.25, 1.5)
  let cs_result = cuckoo_search(cs_config, f, rng)
  let cs_error = calc_relative_error(cs_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Cuckoo Search",
    best_value: cs_result.best_value,
    evaluations: cs_result.evaluations,
    relative_error: cs_error
  })

  let gwo_config = make_gwo_config(max_evals, dim, bounds, 25)
  let gwo_result = grey_wolf_optimizer(gwo_config, f, rng)
  let gwo_error = calc_relative_error(gwo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "GWO",
    best_value: gwo_result.best_value,
    evaluations: gwo_result.evaluations,
    relative_error: gwo_error
  })

  let hj_config = make_hooke_jeeves_config(max_evals, dim, bounds, 0.5, 0.5, 1.0e-6)
  let hj_result = hooke_jeeves(hj_config, f, rng)
  let hj_error = calc_relative_error(hj_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Hooke-Jeeves",
    best_value: hj_result.best_value,
    evaluations: hj_result.evaluations,
    relative_error: hj_error
  })

  let woa_config = make_woa_config(max_evals, dim, bounds, 25, 1.0)
  let woa_result = whale_optimization(woa_config, f, rng)
  let woa_error = calc_relative_error(woa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "WOA",
    best_value: woa_result.best_value,
    evaluations: woa_result.evaluations,
    relative_error: woa_error
  })

  let bat_config = make_bat_config(max_evals, dim, bounds, 25, 0.0, 1.0, 0.5, 1.0)
  let bat_result = bat_algorithm(bat_config, f, rng)
  let bat_error = calc_relative_error(bat_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Bat Algorithm",
    best_value: bat_result.best_value,
    evaluations: bat_result.evaluations,
    relative_error: bat_error
  })

  let hs_config = make_hs_config(max_evals, dim, bounds, 20, 0.9, 0.3, 0.05)
  let hs_result = harmony_search(hs_config, f, rng)
  let hs_error = calc_relative_error(hs_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "Harmony Search",
    best_value: hs_result.best_value,
    evaluations: hs_result.evaluations,
    relative_error: hs_error
  })

  let sca_config = make_sca_config(max_evals, dim, bounds, 25)
  let sca_result = sine_cosine_algorithm(sca_config, f, rng)
  let sca_error = calc_relative_error(sca_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "SCA",
    best_value: sca_result.best_value,
    evaluations: sca_result.evaluations,
    relative_error: sca_error
  })

  let mfo_config = make_mfo_config(max_evals, dim, bounds, 25, 1.0)
  let mfo_result = moth_flame_optimization(mfo_config, f, rng)
  let mfo_error = calc_relative_error(mfo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "MFO",
    best_value: mfo_result.best_value,
    evaluations: mfo_result.evaluations,
    relative_error: mfo_error
  })

  let alo_config = make_alo_config(max_evals, dim, bounds, 25)
  let alo_result = antlion_optimizer(alo_config, f, rng)
  let alo_error = calc_relative_error(alo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "ALO",
    best_value: alo_result.best_value,
    evaluations: alo_result.evaluations,
    relative_error: alo_error
  })

  let ssa_config = make_ssa_config(max_evals, dim, bounds, 25)
  let ssa_result = sparrow_search_algorithm(ssa_config, f, rng)
  let ssa_error = calc_relative_error(ssa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "SSA",
    best_value: ssa_result.best_value,
    evaluations: ssa_result.evaluations,
    relative_error: ssa_error
  })

  let hho_config = make_hho_config(max_evals, dim, bounds, 25)
  let hho_result = harris_hawks_optimization(hho_config, f, rng)
  let hho_error = calc_relative_error(hho_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "HHO",
    best_value: hho_result.best_value,
    evaluations: hho_result.evaluations,
    relative_error: hho_error
  })

  let fpa_config = make_fpa_config(max_evals, dim, bounds, 25)
  let fpa_result = flower_pollination_algorithm(fpa_config, f, rng)
  let fpa_error = calc_relative_error(fpa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "FPA",
    best_value: fpa_result.best_value,
    evaluations: fpa_result.evaluations,
    relative_error: fpa_error
  })

  let da_config = make_da_config(max_evals, dim, bounds, 25)
  let da_result = dragonfly_algorithm(da_config, f, rng)
  let da_error = calc_relative_error(da_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "DA",
    best_value: da_result.best_value,
    evaluations: da_result.evaluations,
    relative_error: da_error
  })

  let csa_config = make_csa_config(max_evals, dim, bounds, 25)
  let csa_result = crow_search_algorithm(csa_config, f, rng)
  let csa_error = calc_relative_error(csa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "CSA",
    best_value: csa_result.best_value,
    evaluations: csa_result.evaluations,
    relative_error: csa_error
  })

  let salp_config = make_salp_sa_config(max_evals, dim, bounds, 25)
  let salp_result = salp_swarm_algorithm(salp_config, f, rng)
  let salp_error = calc_relative_error(salp_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "SalpSA",
    best_value: salp_result.best_value,
    evaluations: salp_result.evaluations,
    relative_error: salp_error
  })

  let mbo_config = make_mbo_config(max_evals, dim, bounds, 25)
  let mbo_result = monarch_butterfly_optimization(mbo_config, f, rng)
  let mbo_error = calc_relative_error(mbo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "MBO",
    best_value: mbo_result.best_value,
    evaluations: mbo_result.evaluations,
    relative_error: mbo_error
  })

  let aeo_config = make_aeo_config(max_evals, dim, bounds, 25)
  let aeo_result = artificial_ecosystem_optimization(aeo_config, f, rng)
  let aeo_error = calc_relative_error(aeo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "AEO",
    best_value: aeo_result.best_value,
    evaluations: aeo_result.evaluations,
    relative_error: aeo_error
  })

  let iwo_config = make_iwo_config(max_evals, dim, bounds, 10, 40, 1, 5, 2.0, 0.1, 3.0)
  let iwo_result = invasive_weed_optimization(iwo_config, f, rng)
  let iwo_error = calc_relative_error(iwo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "IWO",
    best_value: iwo_result.best_value,
    evaluations: iwo_result.evaluations,
    relative_error: iwo_error
  })

  let bfo_config = make_bfo_config(max_evals, dim, bounds, 25, 5, 3, 0.2)
  let bfo_result = bacteria_foraging_optimization(bfo_config, f, rng)
  let bfo_error = calc_relative_error(bfo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "BFO",
    best_value: bfo_result.best_value,
    evaluations: bfo_result.evaluations,
    relative_error: bfo_error
  })

  let sfla_config = make_sfla_config(max_evals, dim, bounds, 30, 5, 10)
  let sfla_result = shuffled_frog_leaping_algorithm(sfla_config, f, rng)
  let sfla_error = calc_relative_error(sfla_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "SFLA",
    best_value: sfla_result.best_value,
    evaluations: sfla_result.evaluations,
    relative_error: sfla_error
  })

  let gso_config = make_gso_config(max_evals, dim, bounds, 25, 0.4, 0.6, 5, 5.0, 1.0)
  let gso_result = glowworm_swarm_optimization(gso_config, f, rng)
  let gso_error = calc_relative_error(gso_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "GSO",
    best_value: gso_result.best_value,
    evaluations: gso_result.evaluations,
    relative_error: gso_error
  })

  let gsa_config = make_gsa_config(max_evals, dim, bounds, 25, 100.0, 0.99)
  let gsa_result = gravitational_search_algorithm(gsa_config, f, rng)
  let gsa_error = calc_relative_error(gsa_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "GSA",
    best_value: gsa_result.best_value,
    evaluations: gsa_result.evaluations,
    relative_error: gsa_error
  })

  let tlbo_config = make_tlbo_config(max_evals, dim, bounds, 30)
  let tlbo_result = teaching_learning_based_optimization(tlbo_config, f, rng)
  let tlbo_error = calc_relative_error(tlbo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "TLBO",
    best_value: tlbo_result.best_value,
    evaluations: tlbo_result.evaluations,
    relative_error: tlbo_error
  })

  let amo_config = make_amo_config(max_evals, dim, bounds, 25, 5)
  let amo_result = animal_migration_optimization(amo_config, f, rng)
  let amo_error = calc_relative_error(amo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "AMO",
    best_value: amo_result.best_value,
    evaluations: amo_result.evaluations,
    relative_error: amo_error
  })

  let mvo_config = make_mvo_config(max_evals, dim, bounds, 25)
  let mvo_result = multi_verse_optimizer(mvo_config, f, rng)
  let mvo_error = calc_relative_error(mvo_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "MVO",
    best_value: mvo_result.best_value,
    evaluations: mvo_result.evaluations,
    relative_error: mvo_error
  })

  let kh_config = make_kh_config(max_evals, dim, bounds, 25)
  let kh_result = krill_herd(kh_config, f, rng)
  let kh_error = calc_relative_error(kh_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "KH",
    best_value: kh_result.best_value,
    evaluations: kh_result.evaluations,
    relative_error: kh_error
  })

  let ica_config = make_ica_config(max_evals, dim, bounds, 25)
  let ica_result = imperialist_competitive_algorithm(ica_config, f, rng)
  let ica_error = calc_relative_error(ica_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "ICA",
    best_value: ica_result.best_value,
    evaluations: ica_result.evaluations,
    relative_error: ica_error
  })

  let sfs_config = make_sfs_config(max_evals, dim, bounds, 25)
  let sfs_result = stochastic_fractal_search(sfs_config, f, rng)
  let sfs_error = calc_relative_error(sfs_result.best_value, known_optimal)
  results.push(BenchmarkResult::{
    algorithm_name: "SFS",
    best_value: sfs_result.best_value,
    evaluations: sfs_result.evaluations,
    relative_error: sfs_error
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
    let err_str = fmt_fixed(r.relative_error, 6)
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
  let aco_values: Array[Double] = []
  let bc_values: Array[Double] = []
  let ts_values: Array[Double] = []
  let ls_values: Array[Double] = []
  let cma_values: Array[Double] = []
  let cem_values: Array[Double] = []
  let cs_values: Array[Double] = []
  let gwo_values: Array[Double] = []
  let hj_values: Array[Double] = []
  let woa_values: Array[Double] = []
  let bat_values: Array[Double] = []
  let hs_values: Array[Double] = []
  let sca_values: Array[Double] = []
  let mfo_values: Array[Double] = []
  let alo_values: Array[Double] = []
  let ssa_values: Array[Double] = []
  let hho_values: Array[Double] = []
  let fpa_values: Array[Double] = []
  let da_values: Array[Double] = []
  let csa_values: Array[Double] = []
  let salp_values: Array[Double] = []
  let mbo_values: Array[Double] = []
  let aeo_values: Array[Double] = []
  let iwo_values: Array[Double] = []
  let bfo_values: Array[Double] = []
  let sfla_values: Array[Double] = []
  let gso_values: Array[Double] = []
  let gsa_values: Array[Double] = []
  let tlbo_values: Array[Double] = []
  let amo_values: Array[Double] = []
  let mvo_values: Array[Double] = []
  let kh_values: Array[Double] = []
  let ica_values: Array[Double] = []
  let sfs_values: Array[Double] = []

  for trial in 0..<num_trials {
    let seed = base_seed + trial.to_double() * 1000.0
    let rng = make_rng(seed)

    let rs_config = make_random_search_config(max_evals, dim, bounds)
    let rs_result = random_search(rs_config, f, rng)
    rs_values.push(rs_result.best_value)

    let hc_config = make_hill_climbing_config(max_evals, dim, bounds, 0.5, 10)
    let hc_result = hill_climbing(hc_config, f, rng)
    hc_values.push(hc_result.best_value)

    let sa_config = make_simulated_annealing_config(max_evals, dim, bounds, 100.0, 0.995, 0.5)
    let sa_result = simulated_annealing(sa_config, f, rng)
    sa_values.push(sa_result.best_value)

    let pso_config = make_pso_config(max_evals, dim, bounds, 30, 0.7, 1.5, 1.5)
    let pso_result = pso(pso_config, f, rng)
    pso_values.push(pso_result.best_value)

    let de_config = make_de_config(max_evals, dim, bounds, 30, 0.8, 0.9)
    let de_result = differential_evolution(de_config, f, rng)
    de_values.push(de_result.best_value)

    let ga_config = make_ga_config(max_evals, dim, bounds, 50, 0.01, 0.8)
    let ga_result = genetic_algorithm(ga_config, f, rng)
    ga_values.push(ga_result.best_value)

    let nm_config = make_nm_config(max_evals, dim, bounds, 1.0, 2.0, 0.5, 0.5)
    let nm_result = nelder_mead(nm_config, f, rng)
    nm_values.push(nm_result.best_value)

    let ff_config = make_firefly_config(max_evals, dim, bounds, 25, 0.5, 1.0, 1.0)
    let ff_result = firefly_algorithm(ff_config, f, rng)
    ff_values.push(ff_result.best_value)

    let aco_config = make_aco_config(max_evals, dim, bounds, 20, 0.1, 1.0, 10.0)
    let aco_result = aco(aco_config, f, rng)
    aco_values.push(aco_result.best_value)

    let bc_config = make_bee_colony_config(max_evals, dim, bounds, 25, 100)
    let bc_result = bee_colony(bc_config, f, rng)
    bc_values.push(bc_result.best_value)

    let ts_config = make_tabu_search_config(max_evals, dim, bounds, 50, 30, 0.5)
    let ts_result = tabu_search(ts_config, f, rng)
    ts_values.push(ts_result.best_value)

    let ls_config = make_local_search_config(max_evals, dim, bounds, 30, 0.5, 0.5)
    let ls_result = local_search(ls_config, f, rng)
    ls_values.push(ls_result.best_value)

    let cma_config = make_cmaes_config(max_evals, dim, bounds, 0.5, 20, 10)
    let cma_result = cma_es(cma_config, f, rng)
    cma_values.push(cma_result.best_value)

    let cem_config = make_cem_config(max_evals, dim, bounds, 100, 0.1, 0.7)
    let cem_result = cross_entropy_method(cem_config, f, rng)
    cem_values.push(cem_result.best_value)

    let cs_config = make_cuckoo_config(max_evals, dim, bounds, 25, 0.25, 1.5)
    let cs_result = cuckoo_search(cs_config, f, rng)
    cs_values.push(cs_result.best_value)

    let gwo_config = make_gwo_config(max_evals, dim, bounds, 25)
    let gwo_result = grey_wolf_optimizer(gwo_config, f, rng)
    gwo_values.push(gwo_result.best_value)

    let hj_config = make_hooke_jeeves_config(max_evals, dim, bounds, 0.5, 0.5, 1.0e-6)
    let hj_result = hooke_jeeves(hj_config, f, rng)
    hj_values.push(hj_result.best_value)

    let woa_config = make_woa_config(max_evals, dim, bounds, 25, 1.0)
    let woa_result = whale_optimization(woa_config, f, rng)
    woa_values.push(woa_result.best_value)

    let bat_config = make_bat_config(max_evals, dim, bounds, 25, 0.0, 1.0, 0.5, 1.0)
    let bat_result = bat_algorithm(bat_config, f, rng)
    bat_values.push(bat_result.best_value)

    let hs_config = make_hs_config(max_evals, dim, bounds, 20, 0.9, 0.3, 0.05)
    let hs_result = harmony_search(hs_config, f, rng)
    hs_values.push(hs_result.best_value)

    let sca_config = make_sca_config(max_evals, dim, bounds, 25)
    let sca_result = sine_cosine_algorithm(sca_config, f, rng)
    sca_values.push(sca_result.best_value)

    let mfo_config = make_mfo_config(max_evals, dim, bounds, 25, 1.0)
    let mfo_result = moth_flame_optimization(mfo_config, f, rng)
    mfo_values.push(mfo_result.best_value)

    let alo_config = make_alo_config(max_evals, dim, bounds, 25)
    let alo_result = antlion_optimizer(alo_config, f, rng)
    alo_values.push(alo_result.best_value)

    let ssa_config = make_ssa_config(max_evals, dim, bounds, 25)
    let ssa_result = sparrow_search_algorithm(ssa_config, f, rng)
    ssa_values.push(ssa_result.best_value)

    let hho_config = make_hho_config(max_evals, dim, bounds, 25)
    let hho_result = harris_hawks_optimization(hho_config, f, rng)
    hho_values.push(hho_result.best_value)

    let fpa_config = make_fpa_config(max_evals, dim, bounds, 25)
    let fpa_result = flower_pollination_algorithm(fpa_config, f, rng)
    fpa_values.push(fpa_result.best_value)

    let da_config = make_da_config(max_evals, dim, bounds, 25)
    let da_result = dragonfly_algorithm(da_config, f, rng)
    da_values.push(da_result.best_value)

    let csa_config = make_csa_config(max_evals, dim, bounds, 25)
    let csa_result = crow_search_algorithm(csa_config, f, rng)
    csa_values.push(csa_result.best_value)

    let salp_config = make_salp_sa_config(max_evals, dim, bounds, 25)
    let salp_result = salp_swarm_algorithm(salp_config, f, rng)
    salp_values.push(salp_result.best_value)

    let mbo_config = make_mbo_config(max_evals, dim, bounds, 25)
    let mbo_result = monarch_butterfly_optimization(mbo_config, f, rng)
    mbo_values.push(mbo_result.best_value)

    let aeo_config = make_aeo_config(max_evals, dim, bounds, 25)
    let aeo_result = artificial_ecosystem_optimization(aeo_config, f, rng)
    aeo_values.push(aeo_result.best_value)

    let iwo_config = make_iwo_config(max_evals, dim, bounds, 10, 40, 1, 5, 2.0, 0.1, 3.0)
    let iwo_result = invasive_weed_optimization(iwo_config, f, rng)
    iwo_values.push(iwo_result.best_value)

    let bfo_config = make_bfo_config(max_evals, dim, bounds, 25, 5, 3, 0.2)
    let bfo_result = bacteria_foraging_optimization(bfo_config, f, rng)
    bfo_values.push(bfo_result.best_value)

    let sfla_config = make_sfla_config(max_evals, dim, bounds, 30, 5, 10)
    let sfla_result = shuffled_frog_leaping_algorithm(sfla_config, f, rng)
    sfla_values.push(sfla_result.best_value)

    let gso_config = make_gso_config(max_evals, dim, bounds, 25, 0.4, 0.6, 5, 5.0, 1.0)
    let gso_result = glowworm_swarm_optimization(gso_config, f, rng)
    gso_values.push(gso_result.best_value)

    let gsa_config = make_gsa_config(max_evals, dim, bounds, 25, 100.0, 0.99)
    let gsa_result = gravitational_search_algorithm(gsa_config, f, rng)
    gsa_values.push(gsa_result.best_value)

    let tlbo_config = make_tlbo_config(max_evals, dim, bounds, 30)
    let tlbo_result = teaching_learning_based_optimization(tlbo_config, f, rng)
    tlbo_values.push(tlbo_result.best_value)

    let amo_config = make_amo_config(max_evals, dim, bounds, 25, 5)
    let amo_result = animal_migration_optimization(amo_config, f, rng)
    amo_values.push(amo_result.best_value)

    let mvo_config = make_mvo_config(max_evals, dim, bounds, 25)
    let mvo_result = multi_verse_optimizer(mvo_config, f, rng)
    mvo_values.push(mvo_result.best_value)

    let kh_config = make_kh_config(max_evals, dim, bounds, 25)
    let kh_result = krill_herd(kh_config, f, rng)
    kh_values.push(kh_result.best_value)

    let ica_config = make_ica_config(max_evals, dim, bounds, 25)
    let ica_result = imperialist_competitive_algorithm(ica_config, f, rng)
    ica_values.push(ica_result.best_value)

    let sfs_config = make_sfs_config(max_evals, dim, bounds, 25)
    let sfs_result = stochastic_fractal_search(sfs_config, f, rng)
    sfs_values.push(sfs_result.best_value)
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

  let aco_min = array_min(aco_values)
  let aco_max = array_max(aco_values)
  let aco_avg = array_avg(aco_values)
  let aco_err = calc_relative_error(aco_avg, known_optimal)
  println("  | ACO               | \{aco_min} | \{aco_max} | \{aco_avg} | \{aco_err} |")

  let bc_min = array_min(bc_values)
  let bc_max = array_max(bc_values)
  let bc_avg = array_avg(bc_values)
  let bc_err = calc_relative_error(bc_avg, known_optimal)
  println("  | Bee Colony        | \{bc_min} | \{bc_max} | \{bc_avg} | \{bc_err} |")

  let ts_min = array_min(ts_values)
  let ts_max = array_max(ts_values)
  let ts_avg = array_avg(ts_values)
  let ts_err = calc_relative_error(ts_avg, known_optimal)
  println("  | Tabu Search       | \{ts_min} | \{ts_max} | \{ts_avg} | \{ts_err} |")

  let ls_min = array_min(ls_values)
  let ls_max = array_max(ls_values)
  let ls_avg = array_avg(ls_values)
  let ls_err = calc_relative_error(ls_avg, known_optimal)
  println("  | Local Search      | \{ls_min} | \{ls_max} | \{ls_avg} | \{ls_err} |")

  let cma_min = array_min(cma_values)
  let cma_max = array_max(cma_values)
  let cma_avg = array_avg(cma_values)
  let cma_err = calc_relative_error(cma_avg, known_optimal)
  println("  | CMA-ES            | \{cma_min} | \{cma_max} | \{cma_avg} | \{cma_err} |")

  let cem_min = array_min(cem_values)
  let cem_max = array_max(cem_values)
  let cem_avg = array_avg(cem_values)
  let cem_err = calc_relative_error(cem_avg, known_optimal)
  println("  | CEM               | \{cem_min} | \{cem_max} | \{cem_avg} | \{cem_err} |")

  let cs_min = array_min(cs_values)
  let cs_max = array_max(cs_values)
  let cs_avg = array_avg(cs_values)
  let cs_err = calc_relative_error(cs_avg, known_optimal)
  println("  | Cuckoo Search     | \{cs_min} | \{cs_max} | \{cs_avg} | \{cs_err} |")

  let gwo_min = array_min(gwo_values)
  let gwo_max = array_max(gwo_values)
  let gwo_avg = array_avg(gwo_values)
  let gwo_err = calc_relative_error(gwo_avg, known_optimal)
  println("  | GWO               | \{gwo_min} | \{gwo_max} | \{gwo_avg} | \{gwo_err} |")

  let hj_min = array_min(hj_values)
  let hj_max = array_max(hj_values)
  let hj_avg = array_avg(hj_values)
  let hj_err = calc_relative_error(hj_avg, known_optimal)
  println("  | Hooke-Jeeves      | \{hj_min} | \{hj_max} | \{hj_avg} | \{hj_err} |")

  let woa_min = array_min(woa_values)
  let woa_max = array_max(woa_values)
  let woa_avg = array_avg(woa_values)
  let woa_err = calc_relative_error(woa_avg, known_optimal)
  println("  | WOA               | \{woa_min} | \{woa_max} | \{woa_avg} | \{woa_err} |")

  let bat_min = array_min(bat_values)
  let bat_max = array_max(bat_values)
  let bat_avg = array_avg(bat_values)
  let bat_err = calc_relative_error(bat_avg, known_optimal)
  println("  | Bat Algorithm     | \{bat_min} | \{bat_max} | \{bat_avg} | \{bat_err} |")

  let hs_min = array_min(hs_values)
  let hs_max = array_max(hs_values)
  let hs_avg = array_avg(hs_values)
  let hs_err = calc_relative_error(hs_avg, known_optimal)
  println("  | Harmony Search    | \{hs_min} | \{hs_max} | \{hs_avg} | \{hs_err} |")

  let sca_min = array_min(sca_values)
  let sca_max = array_max(sca_values)
  let sca_avg = array_avg(sca_values)
  let sca_err = calc_relative_error(sca_avg, known_optimal)
  println("  | SCA               | \{sca_min} | \{sca_max} | \{sca_avg} | \{sca_err} |")

  let mfo_min = array_min(mfo_values)
  let mfo_max = array_max(mfo_values)
  let mfo_avg = array_avg(mfo_values)
  let mfo_err = calc_relative_error(mfo_avg, known_optimal)
  println("  | MFO               | \{mfo_min} | \{mfo_max} | \{mfo_avg} | \{mfo_err} |")

  let alo_min = array_min(alo_values)
  let alo_max = array_max(alo_values)
  let alo_avg = array_avg(alo_values)
  let alo_err = calc_relative_error(alo_avg, known_optimal)
  println("  | ALO               | \{alo_min} | \{alo_max} | \{alo_avg} | \{alo_err} |")

  let ssa_min = array_min(ssa_values)
  let ssa_max = array_max(ssa_values)
  let ssa_avg = array_avg(ssa_values)
  let ssa_err = calc_relative_error(ssa_avg, known_optimal)
  println("  | SSA               | \{ssa_min} | \{ssa_max} | \{ssa_avg} | \{ssa_err} |")

  let hho_min = array_min(hho_values)
  let hho_max = array_max(hho_values)
  let hho_avg = array_avg(hho_values)
  let hho_err = calc_relative_error(hho_avg, known_optimal)
  println("  | HHO               | \{hho_min} | \{hho_max} | \{hho_avg} | \{hho_err} |")

  let fpa_min = array_min(fpa_values)
  let fpa_max = array_max(fpa_values)
  let fpa_avg = array_avg(fpa_values)
  let fpa_err = calc_relative_error(fpa_avg, known_optimal)
  println("  | FPA               | \{fpa_min} | \{fpa_max} | \{fpa_avg} | \{fpa_err} |")

  let da_min = array_min(da_values)
  let da_max = array_max(da_values)
  let da_avg = array_avg(da_values)
  let da_err = calc_relative_error(da_avg, known_optimal)
  println("  | DA                | \{da_min} | \{da_max} | \{da_avg} | \{da_err} |")

  let iwo_min = array_min(iwo_values)
  let iwo_max = array_max(iwo_values)
  let iwo_avg = array_avg(iwo_values)
  let iwo_err = calc_relative_error(iwo_avg, known_optimal)
  println("  | IWO               | \{iwo_min} | \{iwo_max} | \{iwo_avg} | \{iwo_err} |")

  let bfo_min = array_min(bfo_values)
  let bfo_max = array_max(bfo_values)
  let bfo_avg = array_avg(bfo_values)
  let bfo_err = calc_relative_error(bfo_avg, known_optimal)
  println("  | BFO               | \{bfo_min} | \{bfo_max} | \{bfo_avg} | \{bfo_err} |")

  let sfla_min = array_min(sfla_values)
  let sfla_max = array_max(sfla_values)
  let sfla_avg = array_avg(sfla_values)
  let sfla_err = calc_relative_error(sfla_avg, known_optimal)
  println("  | SFLA              | \{sfla_min} | \{sfla_max} | \{sfla_avg} | \{sfla_err} |")

  let gso_min = array_min(gso_values)
  let gso_max = array_max(gso_values)
  let gso_avg = array_avg(gso_values)
  let gso_err = calc_relative_error(gso_avg, known_optimal)
  println("  | GSO               | \{gso_min} | \{gso_max} | \{gso_avg} | \{gso_err} |")

  let gsa_min = array_min(gsa_values)
  let gsa_max = array_max(gsa_values)
  let gsa_avg = array_avg(gsa_values)
  let gsa_err = calc_relative_error(gsa_avg, known_optimal)
  println("  | GSA               | \{gsa_min} | \{gsa_max} | \{gsa_avg} | \{gsa_err} |")

  let tlbo_min = array_min(tlbo_values)
  let tlbo_max = array_max(tlbo_values)
  let tlbo_avg = array_avg(tlbo_values)
  let tlbo_err = calc_relative_error(tlbo_avg, known_optimal)
  println("  | TLBO              | \{tlbo_min} | \{tlbo_max} | \{tlbo_avg} | \{tlbo_err} |")

  let amo_min = array_min(amo_values)
  let amo_max = array_max(amo_values)
  let amo_avg = array_avg(amo_values)
  let amo_err = calc_relative_error(amo_avg, known_optimal)
  println("  | AMO               | \{amo_min} | \{amo_max} | \{amo_avg} | \{amo_err} |")

  let mvo_min = array_min(mvo_values)
  let mvo_max = array_max(mvo_values)
  let mvo_avg = array_avg(mvo_values)
  let mvo_err = calc_relative_error(mvo_avg, known_optimal)
  println("  | MVO               | \{mvo_min} | \{mvo_max} | \{mvo_avg} | \{mvo_err} |")

  let kh_min = array_min(kh_values)
  let kh_max = array_max(kh_values)
  let kh_avg = array_avg(kh_values)
  let kh_err = calc_relative_error(kh_avg, known_optimal)
  println("  | KH                | \{kh_min} | \{kh_max} | \{kh_avg} | \{kh_err} |")

  let ica_min = array_min(ica_values)
  let ica_max = array_max(ica_values)
  let ica_avg = array_avg(ica_values)
  let ica_err = calc_relative_error(ica_avg, known_optimal)
  println("  | ICA               | \{ica_min} | \{ica_max} | \{ica_avg} | \{ica_err} |")

  let sfs_min = array_min(sfs_values)
  let sfs_max = array_max(sfs_values)
  let sfs_avg = array_avg(sfs_values)
  let sfs_err = calc_relative_error(sfs_avg, known_optimal)
  println("  | SFS               | \{sfs_min} | \{sfs_max} | \{sfs_avg} | \{sfs_err} |")

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