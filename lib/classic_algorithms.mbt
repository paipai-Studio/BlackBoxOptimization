// ===== 随机搜索算法 =====

///|
pub struct RandomSearchConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
}

///|
pub fn make_random_search_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
) -> RandomSearchConfig {
  RandomSearchConfig::{ max_evaluations, dimension, bounds, }
}

///|
pub fn random_search(
  config : RandomSearchConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations

  let current_rng = rng
  let initial_solution = random_solution(dim, bounds, current_rng)

  let mut best_solution = initial_solution
  let mut best_value = evaluate(f, best_solution)
  let mut evaluations = 1

  for _i in 1..<max_evals {
    let candidate = random_solution(dim, bounds, current_rng)
    let value = evaluate(f, candidate)
    evaluations = evaluations + 1

    if value < best_value {
      best_value = value
      best_solution = candidate
    }
  }

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 爬山法算法 =====

///|
pub struct HillClimbingConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  step_size : Double
  restarts : Int32
}

///|
pub fn make_hill_climbing_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  step_size : Double,
  restarts : Int32,
) -> HillClimbingConfig {
  HillClimbingConfig::{
    max_evaluations,
    dimension,
    bounds,
    step_size,
    restarts,
  }
}

///|
fn generate_neighbor(
  solution : Array[Double],
  step_size : Double,
  bounds : Bounds,
  rng : RNG,
) -> Array[Double] {
  let neighbor : Array[Double] = []
  for i in 0..<solution.length() {
    let r = RNG::next(rng)
    let perturbation : Double = (r * 2.0 - 1.0) * step_size
    let new_val : Double = max(
      bounds.min,
      min(bounds.max, solution[i] + perturbation),
    )
    neighbor.push(new_val)
  }
  neighbor
}

///|
fn hill_climb(
  f : (Array[Double]) -> Double,
  config : HillClimbingConfig,
  start_solution : Array[Double],
  start_value : Double,
  remaining_evals : Int32,
  rng : RNG,
) -> (Array[Double], Double, Int32) {
  let mut current = start_solution
  let mut current_value = start_value
  let mut evaluations = 0

  while evaluations < remaining_evals {
    let neighbor = generate_neighbor(
      current,
      config.step_size,
      config.bounds,
      rng,
    )
    let neighbor_value = evaluate(f, neighbor)
    evaluations = evaluations + 1

    if neighbor_value < current_value {
      current = neighbor
      current_value = neighbor_value
    }
  }

  (current, current_value, evaluations)
}

///|
pub fn hill_climbing(
  config : HillClimbingConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations

  let mut best_solution : Array[Double] = []
  let mut best_value = 1.0e308
  let mut total_evaluations = 0
  let current_rng = rng

  for _restart in 0..<config.restarts {
    if total_evaluations >= max_evals {
      break
    }

    let solution = random_solution(dim, bounds, current_rng)
    let value = evaluate(f, solution)
    total_evaluations = total_evaluations + 1

    let remaining = max_evals - total_evaluations
    let (final_sol, final_val, evals) = hill_climb(
      f, config, solution, value, remaining, current_rng,
    )
    total_evaluations = total_evaluations + evals

    if final_val < best_value {
      best_value = final_val
      best_solution = final_sol
    }
  }

  OptimizationResult::{
    best_solution,
    best_value,
    evaluations: total_evaluations,
  }
}

// ===== 模拟退火算法 =====

///|
pub struct SimulatedAnnealingConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  initial_temperature : Double
  cooling_rate : Double
  step_size : Double
}

///|
pub fn make_simulated_annealing_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  initial_temperature : Double,
  cooling_rate : Double,
  step_size : Double,
) -> SimulatedAnnealingConfig {
  SimulatedAnnealingConfig::{
    max_evaluations,
    dimension,
    bounds,
    initial_temperature,
    cooling_rate,
    step_size,
  }
}

///|
fn acceptance_probability(
  current_value : Double,
  new_value : Double,
  temperature : Double,
) -> Double {
  if new_value < current_value {
    return 1.0
  }
  let delta = new_value - current_value
  exp_approx(-delta / temperature)
}

///|
fn generate_neighbor_sa(
  solution : Array[Double],
  step_size : Double,
  bounds : Bounds,
  rng : RNG,
) -> Array[Double] {
  let neighbor : Array[Double] = []
  for i in 0..<solution.length() {
    let r = RNG::next(rng)
    let perturbation : Double = (r * 2.0 - 1.0) * step_size
    let new_val : Double = max(
      bounds.min,
      min(bounds.max, solution[i] + perturbation),
    )
    neighbor.push(new_val)
  }
  neighbor
}

///|
pub fn simulated_annealing(
  config : SimulatedAnnealingConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let initial_temp = config.initial_temperature
  let cooling_rate = config.cooling_rate
  let step_size = config.step_size

  let current_rng = rng
  let initial_solution = random_solution(dim, bounds, current_rng)

  let mut current_solution = initial_solution
  let mut current_value = evaluate(f, current_solution)
  let mut best_solution = current_solution
  let mut best_value = current_value
  let mut temperature = initial_temp
  let mut evaluations = 1

  while evaluations < max_evals {
    let neighbor = generate_neighbor_sa(
      current_solution, step_size, bounds, current_rng,
    )
    let neighbor_value = evaluate(f, neighbor)
    evaluations = evaluations + 1

    let prob = acceptance_probability(
      current_value, neighbor_value, temperature,
    )
    let rand_val = RNG::next(current_rng)

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

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 粒子群优化算法 (PSO) =====

///|
pub struct PSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  swarm_size : Int32
  inertia_weight : Double
  cognitive_coeff : Double
  social_coeff : Double
}

///|
pub fn make_pso_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  swarm_size : Int32,
  inertia_weight : Double,
  cognitive_coeff : Double,
  social_coeff : Double,
) -> PSOConfig {
  PSOConfig::{
    max_evaluations,
    dimension,
    bounds,
    swarm_size,
    inertia_weight,
    cognitive_coeff,
    social_coeff,
  }
}

///|
pub fn pso(
  config : PSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let swarm_size = config.swarm_size
  let w = config.inertia_weight
  let c1 = config.cognitive_coeff
  let c2 = config.social_coeff

  let current_rng = rng

  let positions : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let personal_best : Array[Array[Double]] = []
  let personal_best_values : Array[Double] = []

  for _i in 0..<swarm_size {
    let pos = random_solution(dim, bounds, current_rng)
    positions.push(pos)

    let vel : Array[Double] = []
    for _j in 0..<dim {
      let r = RNG::next(current_rng)
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
      if evaluations >= max_evals {
        break
      }

      let new_vel : Array[Double] = []
      let new_pos : Array[Double] = []

      for j in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)

        let vel = w * velocities[i][j] +
          c1 * r1 * (personal_best[i][j] - positions[i][j]) +
          c2 * r2 * (global_best[j] - positions[i][j])
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
    evaluations,
  }
}

// ===== 差分进化算法 (DE) =====

///|
pub struct DEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  mutation_factor : Double
  crossover_rate : Double
}

///|
pub fn make_de_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  population_size : Int32,
  mutation_factor : Double,
  crossover_rate : Double,
) -> DEConfig {
  DEConfig::{
    max_evaluations,
    dimension,
    bounds,
    population_size,
    mutation_factor,
    crossover_rate,
  }
}

///|
pub fn differential_evolution(
  config : DEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let f_mut = config.mutation_factor
  let cr = config.crossover_rate

  let current_rng = rng

  let population : Array[Array[Double]] = []
  let fitness : Array[Double] = []

  for _i in 0..<pop_size {
    let sol = random_solution(dim, bounds, current_rng)
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
      if evaluations >= max_evals {
        break
      }

      let mut a = 0
      let mut b = 0
      let mut c = 0

      let ra = RNG::next(current_rng)
      a = (ra * pop_size.to_double()).to_int() % pop_size
      if a == target {
        a = (a + 1) % pop_size
      }

      let rb = RNG::next(current_rng)
      b = (rb * pop_size.to_double()).to_int() % pop_size
      while b == target || b == a {
        let rb2 = RNG::next(current_rng)
        b = (rb2 * pop_size.to_double()).to_int() % pop_size
      }

      let rc = RNG::next(current_rng)
      c = (rc * pop_size.to_double()).to_int() % pop_size
      while c == target || c == a || c == b {
        let rc2 = RNG::next(current_rng)
        c = (rc2 * pop_size.to_double()).to_int() % pop_size
      }

      let trial : Array[Double] = []
      for j in 0..<dim {
        let r_cross = RNG::next(current_rng)
        let mutant = population[a][j] +
          f_mut * (population[b][j] - population[c][j])
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
    best_value,
    evaluations,
  }
}

// ===== 遗传算法 (GA) =====

///|
pub struct GAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  mutation_rate : Double
  crossover_rate : Double
}

///|
pub fn make_ga_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  population_size : Int32,
  mutation_rate : Double,
  crossover_rate : Double,
) -> GAConfig {
  GAConfig::{
    max_evaluations,
    dimension,
    bounds,
    population_size,
    mutation_rate,
    crossover_rate,
  }
}

///|
fn tournament_select(
  population : Array[Array[Double]],
  fitness : Array[Double],
  rng : RNG,
) -> Array[Double] {
  let r1 = RNG::next(rng)
  let r2 = RNG::next(rng)
  let r3 = RNG::next(rng)

  let idx1 = (r1 * population.length().to_double()).to_int() %
    population.length()
  let idx2 = (r2 * population.length().to_double()).to_int() %
    population.length()
  let idx3 = (r3 * population.length().to_double()).to_int() %
    population.length()

  let best = min(fitness[idx1], min(fitness[idx2], fitness[idx3]))
  if best == fitness[idx1] {
    population[idx1]
  } else if best == fitness[idx2] {
    population[idx2]
  } else {
    population[idx3]
  }
}

///|
pub fn genetic_algorithm(
  config : GAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let mutation_rate = config.mutation_rate
  let crossover_rate = config.crossover_rate

  let current_rng = rng

  let mut population : Array[Array[Double]] = []
  let mut fitness : Array[Double] = []

  for _i in 0..<pop_size {
    let sol = random_solution(dim, bounds, current_rng)
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
    let new_population : Array[Array[Double]] = []
    let new_fitness : Array[Double] = []

    for i in 0..<pop_size {
      if evaluations >= max_evals {
        new_population.push(population[i])
        new_fitness.push(fitness[i])
        continue
      }

      let parent1 = tournament_select(population, fitness, current_rng)
      let parent2 = tournament_select(population, fitness, current_rng)

      let r_cross = RNG::next(current_rng)

      let child : Array[Double] = []
      if r_cross < crossover_rate {
        let r_point = RNG::next(current_rng)
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
        let r_mut = RNG::next(current_rng)
        if r_mut < mutation_rate {
          let r_val = RNG::next(current_rng)
          child[j] = bounds.min + r_val * (bounds.max - bounds.min)
        }
      }

      let child_fitness = evaluate(f, child)
      evaluations = evaluations + 1
      new_population.push(child)
      new_fitness.push(child_fitness)

      if child_fitness < best_value {
        best_value = child_fitness
        best_idx = new_population.length() - 1
      }
    }

    population = new_population
    fitness = new_fitness
  }

  OptimizationResult::{
    best_solution: population[best_idx],
    best_value,
    evaluations,
  }
}

// ===== Nelder-Mead 单纯形法 =====

///|
pub struct NMConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  reflection : Double
  expansion : Double
  contraction : Double
  shrink : Double
}

///|
pub fn make_nm_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  reflection : Double,
  expansion : Double,
  contraction : Double,
  shrink : Double,
) -> NMConfig {
  NMConfig::{
    max_evaluations,
    dimension,
    bounds,
    reflection,
    expansion,
    contraction,
    shrink,
  }
}

///|
pub fn nelder_mead(
  config : NMConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let alpha = config.reflection
  let gamma = config.expansion
  let rho = config.contraction
  let sigma = config.shrink

  let current_rng = rng

  let n_vertices = dim + 1
  let simplex : Array[Array[Double]] = []
  let values : Array[Double] = []

  let base = random_solution(dim, bounds, current_rng)
  simplex.push(base)
  values.push(evaluate(f, base))

  for _i in 1..<n_vertices {
    let vertex : Array[Double] = []
    for j in 0..<dim {
      let r = RNG::next(current_rng)
      vertex.push(base[j] + (r * 2.0 - 1.0) * 0.5)
    }
    let clamped : Array[Double] = []
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
    for i in 1..<n_vertices {
      if values[i] < values[best] {
        best = i
      }
      if values[i] > values[worst] {
        worst = i
      }
    }
    let mut second_worst = if worst == 0 { 1 } else { 0 }
    for i in 0..<n_vertices {
      if i != worst && values[i] > values[second_worst] {
        second_worst = i
      }
    }

    let centroid : Array[Double] = []
    for j in 0..<dim {
      let mut sum = 0.0
      for i in 0..<n_vertices {
        if i != worst {
          sum = sum + simplex[i][j]
        }
      }
      centroid.push(sum / (n_vertices.to_double() - 1.0))
    }

    let reflected : Array[Double] = []
    for j in 0..<dim {
      reflected.push(
        max(
          bounds.min,
          min(
            bounds.max,
            centroid[j] + alpha * (centroid[j] - simplex[worst][j]),
          ),
        ),
      )
    }
    let reflected_val = evaluate(f, reflected)
    evaluations = evaluations + 1

    if reflected_val < values[best] {
      let expanded : Array[Double] = []
      for j in 0..<dim {
        expanded.push(
          max(
            bounds.min,
            min(bounds.max, centroid[j] + gamma * (reflected[j] - centroid[j])),
          ),
        )
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
    } else if reflected_val < values[worst] {
      let contracted : Array[Double] = []
      for j in 0..<dim {
        contracted.push(
          max(
            bounds.min,
            min(bounds.max, centroid[j] + rho * (reflected[j] - centroid[j])),
          ),
        )
      }
      let contracted_val = evaluate(f, contracted)
      evaluations = evaluations + 1

      if contracted_val <= reflected_val {
        simplex[worst] = contracted
        values[worst] = contracted_val
      } else {
        for i in 0..<n_vertices {
          if i != best {
            let shrunk : Array[Double] = []
            for j in 0..<dim {
              shrunk.push(
                max(
                  bounds.min,
                  min(
                    bounds.max,
                    simplex[best][j] +
                    sigma * (simplex[i][j] - simplex[best][j]),
                  ),
                ),
              )
            }
            simplex[i] = shrunk
            values[i] = evaluate(f, shrunk)
            evaluations = evaluations + 1
            if evaluations >= max_evals {
              break
            }
          }
        }
      }
    } else {
      let contracted : Array[Double] = []
      for j in 0..<dim {
        contracted.push(
          max(
            bounds.min,
            min(
              bounds.max,
              centroid[j] + rho * (simplex[worst][j] - centroid[j]),
            ),
          ),
        )
      }
      let contracted_val = evaluate(f, contracted)
      evaluations = evaluations + 1

      if contracted_val < values[worst] {
        simplex[worst] = contracted
        values[worst] = contracted_val
      } else {
        for i in 0..<n_vertices {
          if i != best {
            let shrunk : Array[Double] = []
            for j in 0..<dim {
              shrunk.push(
                max(
                  bounds.min,
                  min(
                    bounds.max,
                    simplex[best][j] +
                    sigma * (simplex[i][j] - simplex[best][j]),
                  ),
                ),
              )
            }
            simplex[i] = shrunk
            values[i] = evaluate(f, shrunk)
            evaluations = evaluations + 1
            if evaluations >= max_evals {
              break
            }
          }
        }
      }
    }
  }

  let mut best = 0
  for i in 1..<n_vertices {
    if values[i] < values[best] {
      best = i
    }
  }

  OptimizationResult::{
    best_solution: simplex[best],
    best_value: values[best],
    evaluations,
  }
}

// ===== Evolution Strategy (mu+lambda)-ES =====

///|
pub struct ESConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  mu : Int32
  lambda : Int32
  sigma_init : Double
}

///|
pub fn make_es_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  mu_val : Int32,
  lambda_val : Int32,
  sigma_val : Double,
) -> ESConfig {
  ESConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    mu: mu_val,
    lambda: lambda_val,
    sigma_init: sigma_val,
  }
}

///|
/// (mu + lambda) Evolution Strategy: mu parents generate lambda offspring by
/// Gaussian mutation, then the best mu individuals survive. The mutation
/// strength sigma shrinks over the run so the search refines late.
pub fn evolution_strategy(
  config : ESConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let mu = config.mu
  let lambda = config.lambda
  let current_rng = rng
  let range = bounds.max - bounds.min

  let parents : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<mu {
    let sol = random_solution(dim, bounds, current_rng)
    parents.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = mu
  let sorted0 = sort_indices(values)
  let mut best_pos : Array[Double] = parents[sorted0[0]]
  let mut best_val = values[sorted0[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let step = config.sigma_init * range * (1.0 - 0.95 * iter)

    // sigma-scaled Gaussian mutation of randomly chosen parents
    let offspring : Array[Array[Double]] = []
    let child_values : Array[Double] = []
    for _l in 0..<lambda {
      if evaluations >= config.max_evaluations {
        break
      }
      let pi = (RNG::next(current_rng) * mu.to_double()).to_int() % mu
      let child : Array[Double] = []
      for d in 0..<dim {
        let g1 = RNG::next(current_rng)
        let g2 = RNG::next(current_rng)
        let gauss = box_muller(g1, g2)
        let pos = parents[pi][d] + gauss * step
        child.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, child)
      evaluations = evaluations + 1
      offspring.push(child)
      child_values.push(cv)
      if cv < best_val {
        best_val = cv
        best_pos = child
      }
    }

    // (mu + lambda) selection: keep the best mu from parents plus offspring
    let pool : Array[Array[Double]] = []
    let pool_values : Array[Double] = []
    for i in 0..<parents.length() {
      pool.push(parents[i])
      pool_values.push(values[i])
    }
    for i in 0..<offspring.length() {
      pool.push(offspring[i])
      pool_values.push(child_values[i])
    }
    let order = sort_indices(pool_values)
    let keep = min_int(mu, pool.length())
    for i in 0..<keep {
      parents[i] = pool[order[i]]
      values[i] = pool_values[order[i]]
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Population-Based Incremental Learning (PBIL) =====

///|
pub struct PBILConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  sample_size : Int32
  elite_size : Int32
  learning_rate : Double
}

///|
pub fn make_pbil_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  sample_sz : Int32,
  elite_sz : Int32,
  lr : Double,
) -> PBILConfig {
  PBILConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    sample_size: sample_sz,
    elite_size: elite_sz,
    learning_rate: lr,
  }
}

///|
/// Continuous PBIL (estimation-of-distribution): keeps a per-dimension mean
/// and standard deviation, samples a population from that distribution, then
/// moves the parameters toward the elite sample's mean/std. The distribution
/// narrows as the population converges, with a small floor to keep exploring.
pub fn population_based_incremental_learning(
  config : PBILConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.sample_size
  let m = max_int(1, min_int(config.elite_size, n))
  let lr = config.learning_rate
  let current_rng = rng
  let range = bounds.max - bounds.min

  let mean : Array[Double] = []
  let sigma : Array[Double] = []
  for _d in 0..<dim {
    mean.push((bounds.min + bounds.max) * 0.5)
    sigma.push(range * 0.25)
  }

  let mut evaluations = 0
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  while evaluations < config.max_evaluations {
    // Sample a population from the current distribution
    let pop : Array[Array[Double]] = []
    let pop_values : Array[Double] = []
    for _i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let sol : Array[Double] = []
      for d in 0..<dim {
        let g1 = RNG::next(current_rng)
        let g2 = RNG::next(current_rng)
        let gauss = box_muller(g1, g2)
        let pos = mean[d] + sigma[d] * gauss
        sol.push(max(bounds.min, min(bounds.max, pos)))
      }
      let v = evaluate(f, sol)
      evaluations = evaluations + 1
      pop.push(sol)
      pop_values.push(v)
      if v < best_val {
        best_val = v
        best_pos = sol
      }
    }

    if pop.length() == 0 {
      break
    }

    // Move mean/sigma toward the elite sample statistics
    let order = sort_indices(pop_values)
    let keep = min_int(m, pop.length())
    let elite_mean : Array[Double] = []
    let elite_var : Array[Double] = []
    for _d in 0..<dim {
      elite_mean.push(0.0)
      elite_var.push(0.0)
    }
    for i in 0..<keep {
      let idx = order[i]
      for d in 0..<dim {
        elite_mean[d] = elite_mean[d] + pop[idx][d]
      }
    }
    for d in 0..<dim {
      elite_mean[d] = elite_mean[d] / keep.to_double()
    }
    for i in 0..<keep {
      let idx = order[i]
      for d in 0..<dim {
        let diff = pop[idx][d] - elite_mean[d]
        elite_var[d] = elite_var[d] + diff * diff
      }
    }
    for d in 0..<dim {
      let sd = (elite_var[d] / keep.to_double()).sqrt()
      mean[d] = (1.0 - lr) * mean[d] + lr * elite_mean[d]
      sigma[d] = (1.0 - lr) * sigma[d] + lr * sd
      if sigma[d] < range * 0.0005 {
        sigma[d] = range * 0.0005
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}
