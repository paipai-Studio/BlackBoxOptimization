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
        let levy = (1.0 / (abs_val(gauss1).pow(2.0 / 3.0))).sqrt() * (r3 * 2.0 - 1.0)
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