// ===== 萤火虫算法 (Firefly) =====

///|
pub struct FireflyConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  alpha : Double
  beta0 : Double
  gamma : Double
}

///|
pub fn make_firefly_config(
  max_evaluations : Int32,
  dimension : Int32,
  bounds : Bounds,
  population_size : Int32,
  alpha : Double,
  beta0 : Double,
  gamma : Double,
) -> FireflyConfig {
  FireflyConfig::{
    max_evaluations,
    dimension,
    bounds,
    population_size,
    alpha,
    beta0,
    gamma,
  }
}

///|
fn euclidean_distance(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<a.length() {
    let diff = a[i] - b[i]
    sum = sum + diff * diff
  }
  sum.sqrt()
}

///|
pub fn firefly_algorithm(
  config : FireflyConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let pop_size = config.population_size
  let alpha = config.alpha
  let beta0 = config.beta0
  let gamma = config.gamma

  let current_rng = rng

  let fireflies : Array[Array[Double]] = []
  let brightness : Array[Double] = []

  for _i in 0..<pop_size {
    let sol = random_solution(dim, bounds, current_rng)
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
    let iter = evaluations.to_double() / max_evals.to_double()
    let alpha_t = alpha * (1.0 - iter)
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }

      for j in 0..<pop_size {
        if brightness[j] < brightness[i] {
          let dist = euclidean_distance(fireflies[i], fireflies[j])
          let beta = beta0 / (1.0 + gamma * dist * dist)

          let new_pos : Array[Double] = []
          for d in 0..<dim {
            let r = RNG::next(current_rng)
            let step = alpha_t * (r * 2.0 - 1.0) * (bounds.max - bounds.min)
            let moved = fireflies[i][d] +
              beta * (fireflies[j][d] - fireflies[i][d]) +
              step
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

          if evaluations >= max_evals {
            break
          }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: fireflies[best_idx],
    best_value,
    evaluations,
  }
}

// ===== 蚁群算法 (ACO) =====

///|
pub struct ACOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  ant_count : Int32
  evaporation_rate : Double
  pheromone_weight : Double
  heuristic_weight : Double
}

///|
pub fn make_aco_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  ants : Int32,
  evap : Double,
  pher_w : Double,
  heur_w : Double,
) -> ACOConfig {
  ACOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    ant_count: ants,
    evaporation_rate: evap,
    pheromone_weight: pher_w,
    heuristic_weight: heur_w,
  }
}

///|
pub fn aco(
  config : ACOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_ants = config.ant_count
  let bound_range = bounds.max - bounds.min

  let current_rng = rng
  let initial_sol = random_solution(dim, bounds, current_rng)
  let mut best_solution = initial_sol
  let mut best_value = evaluate(f, best_solution)
  let mut evaluations = 1

  let mut sigma = bound_range * 0.5

  while evaluations < config.max_evaluations {
    let ant_positions : Array[Array[Double]] = []
    let ant_values : Array[Double] = []

    for _ant in 0..<n_ants {
      let pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
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

      if evaluations >= config.max_evaluations {
        break
      }
    }

    if evaluations >= config.max_evaluations {
      break
    }

    // Pheromone evaporation drives sigma down over time.
    // improvement is larger when the best value is closer to optimal, so
    // better progress shrinks the sampling radius faster.
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let improvement = 1.0 / (1.0 + abs_val(best_value))
    let decay = (0.98 - 0.01 * improvement) * (1.0 - 0.3 * iter)
    sigma = sigma * decay
    sigma = max(bound_range * 0.0001, sigma)
  }

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 蜂群算法 (Artificial Bee Colony) =====

///|
pub struct BeeColonyConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  scout_bee_count : Int32
  max_trials : Int32
}

///|
pub fn make_bee_colony_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  scouts : Int32,
  trials : Int32,
) -> BeeColonyConfig {
  BeeColonyConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    scout_bee_count: scouts,
    max_trials: trials,
  }
}

///|
pub fn bee_colony(
  config : BeeColonyConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_food = config.scout_bee_count

  let current_rng = rng
  let food_sources : Array[Array[Double]] = []
  let food_values : Array[Double] = []
  let trial_counts : Array[Int32] = []

  for _i in 0..<n_food {
    let sol = random_solution(dim, bounds, current_rng)
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
      if evaluations >= config.max_evaluations {
        break
      }

      let neighbor : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        let phi = r * 2.0 - 1.0
        let r2 = RNG::next(current_rng)
        let mut k = (r2 * n_food.to_double()).to_int() % n_food
        if k == i {
          k = (k + 1) % n_food
        }
        let new_val = food_sources[i][d] +
          phi * (food_sources[i][d] - food_sources[k][d])
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

    if evaluations >= config.max_evaluations {
      break
    }

    // Onlooker bees phase
    let mut total_fitness = 0.0
    for i in 0..<n_food {
      total_fitness = total_fitness + 1.0 / (abs_val(food_values[i]) + 1.0)
    }

    for i in 0..<n_food {
      if evaluations >= config.max_evaluations {
        break
      }

      let prob = 1.0 / (abs_val(food_values[i]) + 1.0) / total_fitness
      let r = RNG::next(current_rng)

      if r < prob {
        let neighbor : Array[Double] = []
        for d in 0..<dim {
          let r2 = RNG::next(current_rng)
          let phi = r2 * 2.0 - 1.0
          let r3 = RNG::next(current_rng)
          let mut k = (r3 * n_food.to_double()).to_int() % n_food
          if k == i {
            k = (k + 1) % n_food
          }
          let new_val = food_sources[i][d] +
            phi * (food_sources[i][d] - food_sources[k][d])
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

    if evaluations >= config.max_evaluations {
      break
    }

    // Scout bees phase
    for i in 0..<n_food {
      if trial_counts[i] >= config.max_trials {
        if evaluations >= config.max_evaluations {
          break
        }
        let sol = random_solution(dim, bounds, current_rng)
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

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 禁忌搜索 (Tabu Search) =====

///|
pub struct TabuSearchConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  tabu_list_size : Int32
  neighborhood_size : Int32
  step_size : Double
}

///|
pub fn make_tabu_search_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  tabu_size : Int32,
  neigh_size : Int32,
  step : Double,
) -> TabuSearchConfig {
  TabuSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    tabu_list_size: tabu_size,
    neighborhood_size: neigh_size,
    step_size: step,
  }
}

///|
pub fn tabu_search(
  config : TabuSearchConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let current_rng = rng
  let initial_sol = random_solution(dim, bounds, current_rng)
  let mut current_sol = initial_sol
  let mut current_val = evaluate(f, current_sol)
  let mut evaluations = 1

  let mut best_solution = current_sol
  let mut best_value = current_val

  let tabu_list : Array[String] = []

  while evaluations < config.max_evaluations {
    let mut best_neighbor : Array[Double] = []
    let mut best_neighbor_val = 1.0e308
    let mut best_move_str = ""

    for _n in 0..<config.neighborhood_size {
      let neighbor : Array[Double] = []
      let mut move_str = ""
      for d in 0..<dim {
        let r = RNG::next(current_rng)
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
      }

      if evaluations >= config.max_evaluations {
        break
      }
    }

    if evaluations >= config.max_evaluations {
      break
    }

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

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 局部邻域搜索 (Guided Local Search) =====

///|
pub struct LocalSearchConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  neighborhood_size : Int32
  step_size : Double
  lambda : Double
}

///|
pub fn make_local_search_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  neigh_size : Int32,
  step : Double,
  lambda_val : Double,
) -> LocalSearchConfig {
  LocalSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    neighborhood_size: neigh_size,
    step_size: step,
    lambda: lambda_val,
  }
}

///|
pub fn local_search(
  config : LocalSearchConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let current_rng = rng
  let initial_sol = random_solution(dim, bounds, current_rng)
  let mut best_solution = initial_sol
  let mut best_value = evaluate(f, initial_sol)
  let mut evaluations = 1

  let penalties : Array[Double] = []
  for _i in 0..<dim {
    penalties.push(0.0)
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Neighborhood radius shrinks over time so the search can refine
    let step = config.step_size * (1.0 - 0.99 * iter)
    let mut best_neighbor : Array[Double] = []
    let mut best_altered = 1.0e308
    let mut best_original_val = 1.0e308

    // Sample the neighborhood around the incumbent best solution
    for _n in 0..<config.neighborhood_size {
      let neighbor : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        let delta = (r * 2.0 - 1.0) * step
        neighbor.push(
          max(bounds.min, min(bounds.max, best_solution[d] + delta)),
        )
      }

      let neighbor_val = evaluate(f, neighbor)
      evaluations = evaluations + 1

      let mut altered = neighbor_val
      for d in 0..<dim {
        altered = altered +
          config.lambda * penalties[d] * abs_val(neighbor[d] - best_solution[d])
      }

      if altered < best_altered {
        best_altered = altered
        best_neighbor = neighbor
        best_original_val = neighbor_val
      }

      if evaluations >= config.max_evaluations {
        break
      }
    }

    if evaluations >= config.max_evaluations {
      break
    }

    if best_neighbor.length() > 0 {
      // Penalize dimensions used by this move so the search explores elsewhere
      for d in 0..<dim {
        if abs_val(best_neighbor[d] - best_solution[d]) > step * 0.1 {
          penalties[d] = penalties[d] + 1.0
        }
      }

      // Greedy acceptance: only replace the incumbent when strictly better
      if best_original_val < best_value {
        best_value = best_original_val
        best_solution = best_neighbor
        for d in 0..<dim {
          penalties[d] = 0.0
        }
      }
    }
  }

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== CMA-ES (协方差矩阵自适应进化策略) =====

///|
pub struct CMAESConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  sigma : Double
  lambda : Int32
  mu : Int32
}

///|
pub fn make_cmaes_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  sigma_init : Double,
  lambda_val : Int32,
  mu_val : Int32,
) -> CMAESConfig {
  CMAESConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    sigma: sigma_init,
    lambda: lambda_val,
    mu: mu_val,
  }
}

///|
pub fn cma_es(
  config : CMAESConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let lambda = config.lambda
  let mu = config.mu

  let current_rng = rng
  let mean_sol = random_solution(dim, bounds, current_rng)
  let mean = mean_sol

  let mut sigma = config.sigma

  let cov : Array[Array[Double]] = []
  for i in 0..<dim {
    let row : Array[Double] = []
    for j in 0..<dim {
      if i == j {
        row.push(1.0)
      } else {
        row.push(0.0)
      }
    }
    cov.push(row)
  }

  let mut evaluations = 0
  let mut best_solution : Array[Double] = []
  let mut best_value = 1.0e308

  while evaluations < config.max_evaluations {
    let population : Array[Array[Double]] = []
    let pop_values : Array[Double] = []

    for _i in 0..<lambda {
      let individual : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let gauss = box_muller(r1, r2)
        let val = mean[d] + sigma * cov[d][d].sqrt() * gauss
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

      if evaluations >= config.max_evaluations {
        break
      }
    }

    if evaluations >= config.max_evaluations {
      break
    }

    // Simple gradient-free CMA update
    let sorted_idx = sort_indices(pop_values)
    let new_mean : Array[Double] = []
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

    // Adaptive diagonal covariance from elite spread
    for d in 0..<dim {
      let mut var_sum = 0.0
      for i in 0..<mu {
        let diff = population[sorted_idx[i]][d] - new_mean[d]
        var_sum = var_sum + diff * diff
      }
      let elite_var = var_sum / mu.to_double()
      let target = elite_var / (sigma * sigma)
      if target > 1.0e-12 {
        cov[d][d] = cov[d][d] * 0.9 + 0.1 * target
      }
    }

    sigma = sigma * 0.99
  }

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 交叉熵方法 (Cross-Entropy Method) =====

///|
pub struct CEMConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  sample_size : Int32
  elite_fraction : Double
  smoothing : Double
}

///|
pub fn make_cem_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  sample_sz : Int32,
  elite_frac : Double,
  smooth : Double,
) -> CEMConfig {
  CEMConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    sample_size: sample_sz,
    elite_fraction: elite_frac,
    smoothing: smooth,
  }
}

///|
pub fn cross_entropy_method(
  config : CEMConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_samples = config.sample_size
  let n_elite = max_int(
    1,
    (n_samples.to_double() * config.elite_fraction).to_int(),
  )

  let current_rng = rng
  let initial_sol = random_solution(dim, bounds, current_rng)
  let mean : Array[Double] = []
  for d in 0..<dim {
    mean.push(initial_sol[d])
  }

  let bound_range = bounds.max - bounds.min
  let std : Array[Double] = []
  for _d in 0..<dim {
    std.push(bound_range * 0.5)
  }

  let mut evaluations = 0
  let mut best_solution : Array[Double] = []
  let mut best_value = 1.0e308

  while evaluations < config.max_evaluations {
    // Sample from Gaussian
    let samples : Array[Array[Double]] = []
    let sample_vals : Array[Double] = []

    for _s in 0..<n_samples {
      let sample : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
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

      if evaluations >= config.max_evaluations {
        break
      }
    }

    if evaluations >= config.max_evaluations {
      break
    }

    // Select elite samples
    let sorted_idx = sort_indices(sample_vals)
    let elite_mean : Array[Double] = []
    let elite_std : Array[Double] = []

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
      mean[d] = config.smoothing * mean[d] +
        (1.0 - config.smoothing) * elite_mean[d]
      std[d] = config.smoothing * std[d] +
        (1.0 - config.smoothing) * elite_std[d]
      std[d] = max(bound_range * 0.0001, min(bound_range, std[d]))
    }
  }

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 布谷鸟搜索算法 (Cuckoo Search) =====

///|
pub struct CuckooSearchConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  nest_count : Int32
  abandon_probability : Double
  step_size : Double
}

///|
pub fn make_cuckoo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  nests : Int32,
  abandon_prob : Double,
  step : Double,
) -> CuckooSearchConfig {
  CuckooSearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    nest_count: nests,
    abandon_probability: abandon_prob,
    step_size: step,
  }
}

///|
pub fn cuckoo_search(
  config : CuckooSearchConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_nests = config.nest_count
  let bound_range = bounds.max - bounds.min

  let current_rng = rng
  let nests : Array[Array[Double]] = []
  let nest_values : Array[Double] = []

  for _i in 0..<n_nests {
    let sol = random_solution(dim, bounds, current_rng)
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
      if evaluations >= config.max_evaluations {
        break
      }

      let new_nest : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let gauss1 = box_muller(r1, r2)
        let r3 = RNG::next(current_rng)
        let levy = (1.0 / pow_approx(abs_val(gauss1), 2.0 / 3.0)).sqrt() *
          (r3 * 2.0 - 1.0)
        let step = config.step_size * bound_range * levy
        new_nest.push(
          max(bounds.min, min(bounds.max, nests[i][d] + step * gauss1)),
        )
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

    if evaluations >= config.max_evaluations {
      break
    }

    // Abandon worst nests and rebuild
    let sorted_idx = sort_indices(nest_values)
    let n_abandon = max_int(
      1,
      (n_nests.to_double() * config.abandon_probability).to_int(),
    )

    for j in 0..<n_abandon {
      if evaluations >= config.max_evaluations {
        break
      }
      let worst_idx = sorted_idx[n_nests - 1 - j]
      let new_sol : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let gauss = box_muller(r1, r2)
        new_sol.push(
          max(
            bounds.min,
            min(bounds.max, best_solution[d] + bound_range * 0.1 * gauss),
          ),
        )
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

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== 灰狼优化算法 (Grey Wolf Optimizer) =====

///|
pub struct GWOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  pack_size : Int32
}

///|
pub fn make_gwo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pack_sz : Int32,
) -> GWOConfig {
  GWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    pack_size: pack_sz,
  }
}

///|
pub fn grey_wolf_optimizer(
  config : GWOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_wolves = config.pack_size

  let current_rng = rng
  let wolves : Array[Array[Double]] = []
  let wolf_values : Array[Double] = []

  for _i in 0..<n_wolves {
    let sol = random_solution(dim, bounds, current_rng)
    wolves.push(sol)
    wolf_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_wolves
  let mut alpha_pos : Array[Double] = []
  let mut alpha_val = 1.0e308
  let mut beta_pos : Array[Double] = []
  let mut beta_val = 1.0e308
  let mut delta_pos : Array[Double] = []
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
    let a = 2.0 -
      2.0 * (evaluations.to_double() / config.max_evaluations.to_double())

    for i in 0..<n_wolves {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r_a1 = r1 * 2.0 - 1.0
        let _r_a2 = r2 * 2.0 - 1.0

        let d_alpha = abs_val(alpha_pos[d] - wolves[i][d])
        let x1 = alpha_pos[d] - a * r_a1 * d_alpha

        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        let r_b1 = r3 * 2.0 - 1.0
        let _r_b2 = r4 * 2.0 - 1.0

        let d_beta = abs_val(beta_pos[d] - wolves[i][d])
        let x2 = beta_pos[d] - a * r_b1 * d_beta

        let r5 = RNG::next(current_rng)
        let r6 = RNG::next(current_rng)
        let r_c1 = r5 * 2.0 - 1.0
        let _r_c2 = r6 * 2.0 - 1.0

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
    evaluations,
  }
}

// ===== Hooke-Jeeves Pattern Search =====

///|
pub struct HookeJeevesConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  initial_step : Double
  step_reduction : Double
  min_step : Double
}

///|
pub fn make_hooke_jeeves_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  init_step : Double,
  step_red : Double,
  min_s : Double,
) -> HookeJeevesConfig {
  HookeJeevesConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    initial_step: init_step,
    step_reduction: step_red,
    min_step: min_s,
  }
}

///|
pub fn hooke_jeeves(
  config : HookeJeevesConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let current_rng = rng
  let base = random_solution(dim, bounds, current_rng)
  let mut base_point = base
  let mut base_val = evaluate(f, base_point)
  let mut evaluations = 1

  let mut step = config.initial_step
  let mut best_solution = base_point
  let mut best_value = base_val

  while evaluations < config.max_evaluations && step > config.min_step {
    let mut improved = false

    // Exploratory move
    let mut current_point : Array[Double] = []
    for d in 0..<dim {
      current_point.push(base_point[d])
    }
    let mut current_val = base_val

    for d in 0..<dim {
      if evaluations >= config.max_evaluations {
        break
      }

      // Try positive direction
      let trial1 : Array[Double] = []
      for d2 in 0..<dim {
        if d2 == d {
          trial1.push(
            max(bounds.min, min(bounds.max, current_point[d2] + step)),
          )
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
        let trial2 : Array[Double] = []
        for d2 in 0..<dim {
          if d2 == d {
            trial2.push(
              max(bounds.min, min(bounds.max, current_point[d2] - step)),
            )
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

    if evaluations >= config.max_evaluations {
      break
    }

    if improved {
      // Pattern move
      let pattern_point : Array[Double] = []
      for d in 0..<dim {
        pattern_point.push(
          max(
            bounds.min,
            min(
              bounds.max,
              current_point[d] + (current_point[d] - base_point[d]),
            ),
          ),
        )
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

  OptimizationResult::{ best_solution, best_value, evaluations, }
}

// ===== Whale Optimization Algorithm (WOA) =====

///|
pub struct WOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  spiral_coeff : Double
}

///|
pub fn make_woa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  spiral : Double,
) -> WOAConfig {
  WOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    spiral_coeff: spiral,
  }
}

///|
pub fn whale_optimization(
  config : WOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_whales = config.population_size
  let b = config.spiral_coeff

  let current_rng = rng
  let whales : Array[Array[Double]] = []
  let whale_values : Array[Double] = []

  for _i in 0..<n_whales {
    let sol = random_solution(dim, bounds, current_rng)
    whales.push(sol)
    whale_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_whales
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_whales {
    if whale_values[i] < best_val {
      best_val = whale_values[i]
      best_pos = whales[i]
    }
  }

  let _bound_range = bounds.max - bounds.min

  while evaluations < config.max_evaluations {
    let a = 2.0 -
      2.0 * (evaluations.to_double() / config.max_evaluations.to_double())
    let a2 = a * a

    for i in 0..<n_whales {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      let p_rand = RNG::next(current_rng)
      let r_05 = RNG::next(current_rng)

      if p_rand < 0.5 {
        // Encircling prey or search for prey
        if r_05 < 0.5 {
          // Search for prey (exploration)
          let r1 = RNG::next(current_rng)
          let _r2 = RNG::next(current_rng)
          let _r3 = RNG::next(current_rng)
          let rand_idx = (r1 * n_whales.to_double()).floor().to_int()
          let r_idx = if rand_idx >= n_whales { n_whales - 1 } else { rand_idx }

          for d in 0..<dim {
            let rd1 = RNG::next(current_rng)
            let rd2 = RNG::next(current_rng)
            let r_d1 = rd1 * 2.0 - 1.0
            let _r_d2 = rd2 * 2.0 - 1.0
            let d_leader = abs_val(whales[r_idx][d] - whales[i][d])
            let pos = whales[r_idx][d] - a2 * r_d1 * d_leader
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        } else {
          // Encircling prey (exploitation)
          for d in 0..<dim {
            let rd1 = RNG::next(current_rng)
            let rd2 = RNG::next(current_rng)
            let r_d1 = rd1 * 2.0 - 1.0
            let _r_d2 = rd2 * 2.0 - 1.0
            let d_prey = abs_val(best_pos[d] - whales[i][d])
            let pos = best_pos[d] - a * r_d1 * d_prey
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        }
      } else {
        // Bubble-net attacking (spiral update)
        for d in 0..<dim {
          let rd1 = RNG::next(current_rng)
          let rd2 = RNG::next(current_rng)
          let _r_d1 = rd1 * 2.0 - 1.0
          let r_d2 = rd2 * 2.0 - 1.0
          let d_prey = abs_val(best_pos[d] - whales[i][d])
          let angle = 2.0 * 3.14159265358979323846 * r_d2
          let spiral = exp_approx(b * angle) *
            cos_approx(2.0 * 3.14159265358979323846 * r_d2)
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
    evaluations,
  }
}

// ===== Bat Algorithm (BA) =====

///|
pub struct BatConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  f_min : Double
  f_max : Double
  pulse_rate : Double
  loudness : Double
}

///|
pub fn make_bat_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  fmin : Double,
  fmax : Double,
  pr : Double,
  loud : Double,
) -> BatConfig {
  BatConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    f_min: fmin,
    f_max: fmax,
    pulse_rate: pr,
    loudness: loud,
  }
}

///|
pub fn bat_algorithm(
  config : BatConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_bats = config.population_size

  let current_rng = rng
  let bats : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let bat_values : Array[Double] = []
  let frequencies : Array[Double] = []

  for _i in 0..<n_bats {
    let sol = random_solution(dim, bounds, current_rng)
    let vel : Array[Double] = []
    for _d in 0..<dim {
      vel.push(0.0)
    }
    bats.push(sol)
    velocities.push(vel)
    bat_values.push(evaluate(f, sol))
    frequencies.push(0.0)
  }

  let mut evaluations = n_bats
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_bats {
    if bat_values[i] < best_val {
      best_val = bat_values[i]
      best_pos = bats[i]
    }
  }

  let mut loudness = config.loudness
  let mut pulse_rate = config.pulse_rate
  let _bound_range = bounds.max - bounds.min

  while evaluations < config.max_evaluations {
    for i in 0..<n_bats {
      if evaluations >= config.max_evaluations {
        break
      }

      // Update frequency
      let r_freq = RNG::next(current_rng)
      frequencies[i] = config.f_min + (config.f_max - config.f_min) * r_freq

      // Update velocity
      let new_vel : Array[Double] = []
      for d in 0..<dim {
        let _r_d = RNG::next(current_rng)
        new_vel.push(
          velocities[i][d] + (bats[i][d] - best_pos[d]) * frequencies[i],
        )
      }
      velocities[i] = new_vel

      // Update position
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        new_pos.push(
          max(bounds.min, min(bounds.max, bats[i][d] + velocities[i][d])),
        )
      }

      // Local search
      let r_pulse = RNG::next(current_rng)
      if r_pulse > pulse_rate {
        for d in 0..<dim {
          let r_step = RNG::next(current_rng)
          let step = loudness * (r_step * 2.0 - 1.0)
          new_pos[d] = max(bounds.min, min(bounds.max, best_pos[d] + step))
        }
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      // Accept if better and loudness condition met
      let r_loud = RNG::next(current_rng)
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
    evaluations,
  }
}

// ===== Harmony Search (HS) =====

///|
pub struct HSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  harmony_memory_size : Int32
  hmcr : Double
  par : Double
  bw : Double
}

///|
pub fn make_hs_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  hms : Int32,
  hmcr_val : Double,
  par_val : Double,
  bw_val : Double,
) -> HSConfig {
  HSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    harmony_memory_size: hms,
    hmcr: hmcr_val,
    par: par_val,
    bw: bw_val,
  }
}

///|
pub fn harmony_search(
  config : HSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let hms = config.harmony_memory_size
  let hmcr = config.hmcr
  let par = config.par
  let bw = config.bw

  let current_rng = rng
  let harmony_memory : Array[Array[Double]] = []
  let hm_values : Array[Double] = []

  for _i in 0..<hms {
    let sol = random_solution(dim, bounds, current_rng)
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
    let new_harmony : Array[Double] = []

    for d in 0..<dim {
      let r1 = RNG::next(current_rng)

      if r1 < hmcr {
        // Memory consideration
        let r2 = RNG::next(current_rng)
        let idx = (r2 * hms.to_double()).floor().to_int()
        let pick_idx = if idx >= hms { hms - 1 } else { idx }
        let val = harmony_memory[pick_idx][d]

        // Pitch adjustment
        let r3 = RNG::next(current_rng)
        if r3 < par {
          let r4 = RNG::next(current_rng)
          let delta = (r4 * 2.0 - 1.0) * bw * bound_range
          new_harmony.push(max(bounds.min, min(bounds.max, val + delta)))
        } else {
          new_harmony.push(val)
        }
      } else {
        // Randomization
        let r5 = RNG::next(current_rng)
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
    evaluations,
  }
}

// ===== Sine Cosine Algorithm (SCA) =====

///|
pub struct SCAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sca_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SCAConfig {
  SCAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn sine_cosine_algorithm(
  config : SCAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_sca = config.population_size

  let current_rng = rng
  let population : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n_sca {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n_sca
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_sca {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = population[i]
    }
  }

  let mut _r1 = 2.0

  while evaluations < config.max_evaluations {
    let r1_update = 2.0 -
      evaluations.to_double() * (2.0 / config.max_evaluations.to_double())

    for i in 0..<n_sca {
      if evaluations >= config.max_evaluations {
        break
      }

      let r2_val = RNG::next(current_rng)
      let r3_val = RNG::next(current_rng)
      let r4_val = RNG::next(current_rng)

      let r2 = 2.0 * PI * r2_val
      let r3 = 2.0 * r3_val - 1.0
      let r4 = r4_val

      let new_pos : Array[Double] = []
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

    _r1 = r1_update
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Moth Flame Optimization (MFO) =====

///|
pub struct MFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  spiral_constant : Double
}

///|
pub fn make_mfo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  spiral : Double,
) -> MFOConfig {
  MFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    spiral_constant: spiral,
  }
}

///|
pub fn moth_flame_optimization(
  config : MFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_moths = config.population_size
  let b = config.spiral_constant

  let current_rng = rng
  let moths : Array[Array[Double]] = []
  let moth_values : Array[Double] = []

  for _i in 0..<n_moths {
    let sol = random_solution(dim, bounds, current_rng)
    moths.push(sol)
    moth_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_moths
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_moths {
    if moth_values[i] < best_val {
      best_val = moth_values[i]
      best_pos = moths[i]
    }
  }

  let flames : Array[Array[Double]] = []
  let flame_values : Array[Double] = []
  for i in 0..<n_moths {
    flames.push(moths[i])
    flame_values.push(moth_values[i])
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let num_flames = max_int(1, (n_moths.to_double() * (1.0 - iter)).to_int())

    // Sort flames by fitness (ascending)
    let sorted_flame_indices : Array[Int32] = []
    for i in 0..<flame_values.length() {
      sorted_flame_indices.push(i)
    }
    for i in 0..<sorted_flame_indices.length() {
      for j in (i + 1)..<sorted_flame_indices.length() {
        if flame_values[sorted_flame_indices[j]] <
          flame_values[sorted_flame_indices[i]] {
          let temp = sorted_flame_indices[i]
          sorted_flame_indices[i] = sorted_flame_indices[j]
          sorted_flame_indices[j] = temp
        }
      }
    }

    // Sort moths by fitness
    let sorted_moth_indices : Array[Int32] = []
    for i in 0..<moth_values.length() {
      sorted_moth_indices.push(i)
    }
    for i in 0..<sorted_moth_indices.length() {
      for j in (i + 1)..<sorted_moth_indices.length() {
        if moth_values[sorted_moth_indices[j]] <
          moth_values[sorted_moth_indices[i]] {
          let temp = sorted_moth_indices[i]
          sorted_moth_indices[i] = sorted_moth_indices[j]
          sorted_moth_indices[j] = temp
        }
      }
    }

    // Update each moth using spiral path toward its assigned flame
    for rank in 0..<n_moths {
      if evaluations >= config.max_evaluations {
        break
      }

      let moth_idx = sorted_moth_indices[rank]
      let flame_idx = min_int(rank, num_flames - 1)
      let flame_idx_sorted = sorted_flame_indices[flame_idx]

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r_t = RNG::next(current_rng)
        let t = r_t * 2.0 - 1.0
        let flame_pos = flames[flame_idx_sorted][d]
        let dist = abs_val(flame_pos - moths[moth_idx][d])
        let pos = dist * exp_approx(b * t) * cos_approx(2.0 * PI * t) +
          flame_pos
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
    evaluations,
  }
}

// ===== Antlion Optimizer (ALO) =====

///|
pub struct ALOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  ant_count : Int32
}

///|
pub fn make_alo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  ants : Int32,
) -> ALOConfig {
  ALOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    ant_count: ants,
  }
}

///|
pub fn antlion_optimizer(
  config : ALOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n_ants = config.ant_count
  let n_antlions = config.ant_count

  let current_rng = rng
  let antlions : Array[Array[Double]] = []
  let antlion_values : Array[Double] = []

  for _i in 0..<n_antlions {
    let sol = random_solution(dim, bounds, current_rng)
    antlions.push(sol)
    antlion_values.push(evaluate(f, sol))
  }

  let mut evaluations = n_antlions
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n_antlions {
    if antlion_values[i] < best_val {
      best_val = antlion_values[i]
      best_pos = antlions[i]
    }
  }

  let ants : Array[Array[Double]] = []
  let ant_values : Array[Double] = []
  for i in 0..<n_ants {
    ants.push(antlions[i])
    ant_values.push(antlion_values[i])
  }

  while evaluations < config.max_evaluations {
    // Sort antlions by fitness
    let sorted_indices : Array[Int32] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let r_select = RNG::next(current_rng)

      // Tournament selection: pick from elite pool with probability, else from entire pool
      let sel_idx = if r_select < 0.7 {
        let rw = RNG::next(current_rng)
        (rw * elite_count.to_double()).floor().to_int()
      } else {
        let rw = RNG::next(current_rng)
        (rw * n_antlions.to_double()).floor().to_int()
      }
      let sel = max_int(0, min_int(n_antlions - 1, sel_idx))
      let selected_antlion = antlions[sorted_indices[sel]]

      // Random walk toward the selected antlion (step scales with distance,
      // so ants converge onto their antlion)
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r_walk = RNG::next(current_rng)
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
    let sorted_antlion_idx : Array[Int32] = []
    for i in 0..<n_antlions {
      sorted_antlion_idx.push(i)
    }
    for i in 0..<n_antlions {
      for j in (i + 1)..<n_antlions {
        if antlion_values[sorted_antlion_idx[j]] >
          antlion_values[sorted_antlion_idx[i]] {
          let temp = sorted_antlion_idx[i]
          sorted_antlion_idx[i] = sorted_antlion_idx[j]
          sorted_antlion_idx[j] = temp
        }
      }
    }

    let sorted_ant_idx : Array[Int32] = []
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
    evaluations,
  }
}

// ===== Sparrow Search Algorithm (SSA) =====

///|
pub struct SSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  discoverer_ratio : Double
  warning_threshold : Double
}

///|
pub fn make_ssa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SSAConfig {
  SSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    discoverer_ratio: 0.2,
    warning_threshold: 0.8,
  }
}

///|
pub fn sparrow_search_algorithm(
  config : SSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let population : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = population[i]
    }
  }

  let mut worst_pos : Array[Double] = []
  let mut worst_val = -1.0e308
  for i in 0..<n {
    if values[i] > worst_val {
      worst_val = values[i]
      worst_pos = population[i]
    }
  }

  let discoverer_count = max_int(
    1,
    (n.to_double() * config.discoverer_ratio).to_int(),
  )
  let follower_count = n - discoverer_count

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let alpha_warning = config.warning_threshold * (1.0 - iter)

    let sorted_indices : Array[Int32] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let idx = sorted_indices[k]
      let r_warn = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      if r_warn < alpha_warning {
        for _d in 0..<dim {
          let r_step = RNG::next(current_rng)
          let step = r_step * (bounds.max - bounds.min) + bounds.min
          new_pos.push(max(bounds.min, min(bounds.max, step)))
        }
      } else {
        for d in 0..<dim {
          let r_foll = RNG::next(current_rng)
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
      if evaluations >= config.max_evaluations {
        break
      }

      let idx = sorted_indices[discoverer_count + k]

      if values[idx] > best_val {
        let r_fly = RNG::next(current_rng)

        let new_pos : Array[Double] = []
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
        let r_rand = RNG::next(current_rng)

        let new_pos : Array[Double] = []
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
    evaluations,
  }
}
