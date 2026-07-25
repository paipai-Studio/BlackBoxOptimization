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