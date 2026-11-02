// ===== Moth Search Algorithm (MSA) =====

///|
pub struct MSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_msa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MSAConfig {
  MSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn moth_search_algorithm(
  config : MSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let moths : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    moths.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = moths[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Follow best moth via Lévy flight
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let g1 = box_muller(r1, r2)

        let pos = if r3 < 0.5 {
          // Lévy flight towards best
          moths[i][d] + 2.0 * g1 * (best_pos[d] - moths[i][d])
        } else {
          // Spiral search (convergent factor decays quadratically so larger
          // budgets also leave the exploration regime in time)
          moths[i][d] +
          (2.0 * r3 - 1.0) *
          (bounds.max - bounds.min) *
          ((1.0 - iter) * (1.0 - iter))
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      // Greedy acceptance: discarding worse moves keeps the swarm from being
      // randomized by the wide spiral steps (it stalled in high dimensions).
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        moths[i] = new_pos
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
    evaluations,
  }
}

// ===== Pathfinder Algorithm (PFA) =====

///|
pub struct PFAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pfa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PFAConfig {
  PFAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn pathfinder_algorithm(
  config : PFAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let animals : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    animals.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let init_sorted = sort_indices(values)
  let mut global_best_pos : Array[Double] = animals[init_sorted[0]]
  let mut global_best_val : Double = values[init_sorted[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let sorted_idx = sort_indices(values)
    let mut best_pos : Array[Double] = animals[sorted_idx[0]]
    let mut best_val = values[sorted_idx[0]]
    let worst_pos : Array[Double] = animals[sorted_idx[n - 1]]

    // Update pathfinder (best individual) via random perturbation
    let rp = RNG::next(current_rng)
    let pf_new : Array[Double] = []
    for d in 0..<dim {
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let pos = if rp < 0.5 {
        best_pos[d] +
        2.0 * (r2 - 0.5) * (bounds.max - bounds.min) * (1.0 - iter)
      } else {
        best_pos[d] + 2.0 * (r3 - 0.5) * (best_pos[d] - worst_pos[d])
      }
      pf_new.push(max(bounds.min, min(bounds.max, pos)))
    }

    let pf_val = evaluate(f, pf_new)
    evaluations = evaluations + 1
    // Greedy acceptance: the pathfinder slot keeps the better position
    if pf_val < best_val {
      animals[sorted_idx[0]] = pf_new
      values[sorted_idx[0]] = pf_val
      best_val = pf_val
      best_pos = pf_new
    }

    if pf_val < global_best_val {
      global_best_val = pf_val
      global_best_pos = pf_new
    }

    // Update followers (move towards pathfinder with noise)
    for i in 1..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r4 = RNG::next(current_rng)
        let r5 = RNG::next(current_rng)

        let pos = animals[i][d] +
          (best_pos[d] - animals[i][d]) * r4 +
          (r5 - 0.5) * (bounds.max - bounds.min) * 0.1 * (1.0 - iter)
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      animals[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }

      if new_val < global_best_val {
        global_best_val = new_val
        global_best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: global_best_pos,
    best_value: global_best_val,
    evaluations,
  }
}

// ===== Colony Predation Algorithm (CPA) =====

///|
pub struct CPAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cpa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CPAConfig {
  CPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn colony_predation_algorithm(
  config : CPAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let tigers : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    tigers.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = tigers[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let mut global_best_pos : Array[Double] = best_pos
  let mut global_best_val : Double = best_val

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    sorted_idx = sort_indices(values)
    best_pos = tigers[sorted_idx[0]]
    best_val = values[sorted_idx[0]]

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)

        // Predation: move towards best, with "escape" probability
        let pos = if r3 < 0.3 {
          // Random escape jump, shrinks to zero over time
          tigers[i][d] +
          (r1 - 0.5) * 2.0 * (bounds.max - bounds.min) * (1.0 - iter)
        } else {
          // Chase best prey, noise decays for fine convergence
          tigers[i][d] +
          r1 * (best_pos[d] - tigers[i][d]) +
          (r2 - 0.5) * 0.5 * (1.0 - iter)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      tigers[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
      if new_val < global_best_val {
        global_best_val = new_val
        global_best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: global_best_pos,
    best_value: global_best_val,
    evaluations,
  }
}

// ===== Bacteria Colony Optimization (BCO) =====

///|
pub struct BCOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bco_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BCOConfig {
  BCOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn bacteria_colony_optimization(
  config : BCOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let bacteria : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    bacteria.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = bacteria[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let mut global_best_pos : Array[Double] = best_pos
  let mut global_best_val : Double = best_val

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    sorted_idx = sort_indices(values)
    best_pos = bacteria[sorted_idx[0]]
    best_val = values[sorted_idx[0]]

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let chem = 0.2 * (1.0 - iter)

        // Chemotaxis: move towards best
        let pos = bacteria[i][d] +
          chem * r1 * (best_pos[d] - bacteria[i][d]) +
          (r2 - 0.5) * chem
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      bacteria[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
      if new_val < global_best_val {
        global_best_val = new_val
        global_best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: global_best_pos,
    best_value: global_best_val,
    evaluations,
  }
}

// ===== Flying Fox Optimization Algorithm (FFOA) =====

///|
pub struct FFOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ffoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FFOAConfig {
  FFOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn flying_fox_optimization(
  config : FFOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let foxes : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    foxes.push(sol)
    values.push(evaluate(f, sol))

    let vel : Array[Double] = []
    for _d in 0..<dim {
      let v = RNG::next(current_rng)
      vel.push((v - 0.5) * 2.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n
  let mut sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = foxes[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let mut global_best_pos : Array[Double] = best_pos
  let mut global_best_val : Double = best_val

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let inertia = 0.9 - iter * 0.5
    sorted_idx = sort_indices(values)
    best_pos = foxes[sorted_idx[0]]
    best_val = values[sorted_idx[0]]

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)

        velocities[i][d] = inertia * velocities[i][d] +
          r1 * (best_pos[d] - foxes[i][d]) +
          (r2 - 0.5) * (bounds.max - bounds.min) * 0.3 * (1.0 - iter)
        let pos = foxes[i][d] + velocities[i][d] * r3
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      foxes[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
      if new_val < global_best_val {
        global_best_val = new_val
        global_best_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: global_best_pos,
    best_value: global_best_val,
    evaluations,
  }
}

// ===== Artificial Fish Swarm Algorithm (AFSA) =====

///|
pub struct AFSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  visual : Double
  step : Double
  crowd_factor : Double
}

///|
pub fn make_afsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  visual : Double,
  step : Double,
  crowd : Double,
) -> AFSAConfig {
  AFSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    visual,
    step,
    crowd_factor: crowd,
  }
}

///|
pub fn artificial_fish_swarm(
  config : AFSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let visual = config.visual
  let step = config.step
  let crowd = config.crowd_factor
  let current_rng = rng

  let fish : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    fish.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = fish[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // Find neighbors within visual range
      let mut neighbor_count = 0
      let neighbor_sum : Array[Double] = []
      let mut neighbor_val_sum = 0.0
      let mut best_neighbor_idx = -1
      let mut best_neighbor_val = values[i]
      for _d in 0..<dim {
        neighbor_sum.push(0.0)
      }

      for j in 0..<n {
        if j == i {
          continue
        }
        let dist = euclidean_distance(fish[i], fish[j])
        if dist <= visual {
          neighbor_count = neighbor_count + 1
          for d in 0..<dim {
            neighbor_sum[d] = neighbor_sum[d] + fish[j][d]
          }
          neighbor_val_sum = neighbor_val_sum + values[j]
          if values[j] < best_neighbor_val {
            best_neighbor_val = values[j]
            best_neighbor_idx = j
          }
        }
      }

      let new_pos : Array[Double] = []
      if neighbor_count > 0 {
        let neighbor_avg_val = neighbor_val_sum / neighbor_count.to_double()
        // Swarm behavior: move towards center if not crowded and center is better
        if neighbor_avg_val < values[i] &&
          neighbor_count.to_double() < crowd * n.to_double() {
          for d in 0..<dim {
            let center = neighbor_sum[d] / neighbor_count.to_double()
            let r = RNG::next(current_rng)
            let pos = fish[i][d] + step * r * (center - fish[i][d]) / visual
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        } else if best_neighbor_idx >= 0 && best_neighbor_val < values[i] {
          // Follow behavior: move towards best neighbor
          for d in 0..<dim {
            let r = RNG::next(current_rng)
            let pos = fish[i][d] +
              step * r * (fish[best_neighbor_idx][d] - fish[i][d]) / visual
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        } else {
          // Prey behavior: random search in visual range
          for d in 0..<dim {
            let r = RNG::next(current_rng)
            let pos = fish[i][d] + (r - 0.5) * 2.0 * visual
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        }
      } else {
        // No neighbors: prey behavior (random search)
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          let pos = fish[i][d] + (r - 0.5) * 2.0 * visual
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < values[i] {
        fish[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Fireworks Algorithm (FWA) =====

///|
pub struct FWAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fwa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FWAConfig {
  FWAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn fireworks_algorithm(
  config : FWAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let fireworks : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    fireworks.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308
  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = fireworks[i]
    }
  }
  let range = bounds.max - bounds.min

  while evaluations < config.max_evaluations {
    // Refresh best/worst of the current population
    let cur_sorted = sort_indices(values)
    let best_v = values[cur_sorted[0]]
    let worst_v = values[cur_sorted[n - 1]]

    // Candidate pool: current fireworks plus all sparks generated this round
    let pool : Array[Array[Double]] = []
    let pool_values : Array[Double] = []
    for i in 0..<n {
      pool.push(fireworks[i])
      pool_values.push(values[i])
    }

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // ratio 0 = best firework, 1 = worst firework (minimization)
      let ratio = if worst_v == best_v {
        0.0
      } else {
        (values[i] - best_v) / (worst_v - best_v)
      }
      // Explosion amplitude shrinks across iterations so late sparks can
      // refine the best firework; ratio keeps better fireworks more local.
      let iter = evaluations.to_double() / config.max_evaluations.to_double()
      let iter_decay = max(0.001, 1.0 - iter)
      let amplitude = range * (0.05 + 0.95 * ratio) * iter_decay
      let num_sparks = max_int(2, min_int(8, 8 - (ratio * 6.0).to_int()))

      for _s in 0..<num_sparks {
        if evaluations >= config.max_evaluations {
          break
        }

        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          let gauss = box_muller(r, r * 0.5 + 0.1)
          let pos = fireworks[i][d] + amplitude * gauss
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        pool.push(new_pos)
        pool_values.push(new_val)

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      }
    }

    // Selection: keep the best n candidates as the next generation
    let pool_sorted = sort_indices(pool_values)
    let keep = min_int(n, pool.length())
    for i in 0..<keep {
      fireworks[i] = pool[pool_sorted[i]]
      values[i] = pool_values[pool_sorted[i]]
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Brain Storm Optimization (BSO) =====

///|
pub struct BSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BSOConfig {
  BSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn brain_storm_optimization(
  config : BSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let num_clusters = 3
  let current_rng = rng

  let ideas : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    ideas.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = ideas[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let sorted = sort_indices(values)

    // Assign to clusters (simple round-robin based on sorted order)
    let cluster_centers : Array[Array[Double]] = []
    for c in 0..<num_clusters {
      cluster_centers.push(ideas[sorted[c]])
    }

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r1 = RNG::next(current_rng)
      let _r2 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      // Select a cluster center and add perturbation
      let cluster_idx = (r1 * num_clusters.to_double()).to_int()
      let cluster_idx = if cluster_idx >= num_clusters {
        num_clusters - 1
      } else {
        cluster_idx
      }
      let center = cluster_centers[cluster_idx]

      for d in 0..<dim {
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        let gauss = box_muller(r3, r4)
        // Perturbation decreases over iterations
        let perturbation = gauss *
          (1.0 - iter) *
          (bounds.max - bounds.min) *
          0.1
        let pos = center[d] + perturbation
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < values[i] {
        ideas[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Backtracking Search Algorithm (BSA) =====

///|
pub struct BSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BSAConfig {
  BSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn backtracking_search_algorithm(
  config : BSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let population : Array[Array[Double]] = []
  let old_population : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    values.push(evaluate(f, sol))

    let old_sol = random_solution(dim, bounds, current_rng)
    old_population.push(old_sol)
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = population[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    // At interval, save current population as historical
    let r_flag = RNG::next(current_rng)
    if r_flag < 0.5 {
      for i in 0..<n {
        for d in 0..<dim {
          old_population[i][d] = population[i][d]
        }
      }
    }

    // Generate permutation for historical population indices
    let perm : Array[Int32] = []
    for _i in 0..<n {
      let r_p = RNG::next(current_rng)
      perm.push((r_p * n.to_double()).to_int() % n)
    }

    // Generate trial population via crossover with historical
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let old_idx = perm[i]
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)

        // Mix current with permuted historical, scaled by direction
        let direction = if r1 < 0.5 { 1.0 } else { -1.0 }
        let scale = 3.0 * r2
        let trial = population[i][d] +
          direction * scale * (old_population[old_idx][d] - population[i][d])
        new_pos.push(max(bounds.min, min(bounds.max, trial)))
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < values[i] {
        population[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Pigeon-Inspired Optimization (PIO) =====

///|
pub struct PIOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pio_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PIOConfig {
  PIOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn pigeon_inspired_optimization(
  config : PIOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let pigeons : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pigeons.push(sol)
    values.push(evaluate(f, sol))

    let vel : Array[Double] = []
    for _d in 0..<dim {
      let v = RNG::next(current_rng)
      vel.push((v - 0.5) * 2.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pigeons[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  // Phase 1: Map and compass operator (first 60% of evaluations)
  let phase1_end = config.max_evaluations * 6 / 10

  while evaluations < config.max_evaluations {
    if evaluations < phase1_end {
      // Map and compass phase
      let iter = evaluations.to_double() / phase1_end.to_double()
      let map_factor = exp_approx(-5.0 * iter)

      for i in 0..<n {
        if evaluations >= config.max_evaluations {
          break
        }

        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          velocities[i][d] = velocities[i][d] * map_factor +
            r * (best_pos[d] - pigeons[i][d])
          let pos = pigeons[i][d] + velocities[i][d]
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        pigeons[i] = new_pos
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        values[i] = new_val

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      }
    } else {
      // Landmark phase: move towards center of best pigeons
      let iter = (evaluations - phase1_end).to_double() /
        (config.max_evaluations - phase1_end).to_double()
      let num_landmark = max_int(
        1,
        n / 2 - (iter * n.to_double() / 2.0).to_int(),
      )

      let sorted = sort_indices(values)
      let center : Array[Double] = []
      for d in 0..<dim {
        let mut sum = 0.0
        for k in 0..<num_landmark {
          sum = sum + pigeons[sorted[k]][d]
        }
        center.push(sum / num_landmark.to_double())
      }

      for i in 0..<n {
        if evaluations >= config.max_evaluations {
          break
        }

        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          let pos = pigeons[i][d] + r * (center[d] - pigeons[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        pigeons[i] = new_pos
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        values[i] = new_val

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
    evaluations,
  }
}

// ===== Arithmetic Optimization Algorithm (AOA) =====

///|
pub struct AOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AOAConfig {
  AOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn arithmetic_optimization_algorithm(
  config : AOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let positions : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    positions.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = positions[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let alpha = 5.0
  let mu = 0.499

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Math Optimizer Acceleration function: decreases from 1 to 0
    let moa = 0.2 + iter * 0.8
    // Math Optimizer probability: increases from -2 to 2
    let mop = 1.0 - @math.pow(iter, 1.0 / alpha)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)

        let pos = if r1 > moa {
          // Exploration phase: division or multiplication
          if r2 < 0.5 {
            // D: best / (eps + mu * ((bounds.max - bounds.min) * r2 + bounds.min))
            let denom = 0.0001 +
              mu * ((bounds.max - bounds.min) * r2 + bounds.min)
            positions[i][d] * (1.0 - mop) + best_pos[d] / denom * mop
          } else {
            // M: best * (eps + mu * ((bounds.max - bounds.min) * r2 + bounds.min))
            positions[i][d] * (1.0 - mop) +
            best_pos[d] *
            (0.0001 + mu * ((bounds.max - bounds.min) * r2 + bounds.min)) *
            mop
          }
          // Exploitation phase: subtraction or addition
        } else if r2 < 0.5 {
          // S: best - mop * ((bounds.max - bounds.min) * r2 + bounds.min)
          positions[i][d] * mop +
          best_pos[d] * (1.0 - mop) -
          mop * ((bounds.max - bounds.min) * r2 + bounds.min) * 0.1
        } else {
          // A: best + mop * ((bounds.max - bounds.min) * r2 + bounds.min)
          positions[i][d] * mop +
          best_pos[d] * (1.0 - mop) +
          mop * ((bounds.max - bounds.min) * r2 + bounds.min) * 0.1
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      positions[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

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

// ===== Honey Badger Algorithm (HBA) =====

///|
pub struct HBAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hba_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HBAConfig {
  HBAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn honey_badger_algorithm(
  config : HBAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let positions : Array[Array[Double]] = []
  let values : Array[Double] = []
  let intensity : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    positions.push(sol)
    values.push(evaluate(f, sol))
    intensity.push(0.0)
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = positions[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let beta = 6.0
  let c = 2.0

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Density factor decreases over time
    let alpha = c - iter * c

    // Calculate smell intensity for each badger
    for i in 0..<n {
      let dist = euclidean_distance(positions[i], best_pos)
      let di = 0.01 + dist
      intensity[i] = exp_approx(-(di * di) / (2.0 * 0.5))
    }

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let pos = if r1 < 0.5 {
          // Digging phase: move towards best with intensity and random perturbation
          positions[i][d] +
          intensity[i] * beta * (best_pos[d] - positions[i][d]) +
          (r2 - 0.5) * alpha * (bounds.max - bounds.min) * 0.1
        } else {
          // Honey phase: guided by best with random factor
          positions[i][d] +
          (r3 - 0.5) * 2.0 * (best_pos[d] - positions[i][d]) * (1.0 - r4)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      positions[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

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

// ===== Pelican Optimization Algorithm (POA) =====

///|
pub struct POAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_poa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> POAConfig {
  POAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn pelican_optimization_algorithm(
  config : POAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let pelicans : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pelicans.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pelicans[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // Phase 1: Moving towards prey (exploration)
      let r1 = RNG::next(current_rng)
      let prey_idx = (r1 * n.to_double()).to_int() % n
      let prey = pelicans[prey_idx]

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r2 = RNG::next(current_rng)
        let pos = if values[prey_idx] < values[i] {
          pelicans[i][d] + r2 * (prey[d] - pelicans[i][d])
        } else {
          pelicans[i][d] + r2 * (pelicans[i][d] - prey[d])
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < values[i] {
        pelicans[i] = new_pos
        values[i] = new_val
        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      }

      // Phase 2: Surface expansion (exploitation)
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos2 : Array[Double] = []
      for d in 0..<dim {
        let r3 = RNG::next(current_rng)
        let pos2 = pelicans[i][d] +
          (r3 - 0.5) * 2.0 * (bounds.max - bounds.min) * (1.0 - iter) * 0.5
        new_pos2.push(max(bounds.min, min(bounds.max, pos2)))
      }

      let new_val2 = evaluate(f, new_pos2)
      evaluations = evaluations + 1

      if new_val2 < values[i] {
        pelicans[i] = new_pos2
        values[i] = new_val2
        if new_val2 < best_val {
          best_val = new_val2
          best_pos = new_pos2
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Tuna Swarm Optimization (TSA) =====

///|
pub struct TSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_tsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TSAConfig {
  TSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn tuna_swarm_optimization(
  config : TSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let tuna : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    tuna.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = tuna[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let z = 0.05

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let aa = 2.0 - 2.0 * iter

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)

        let t = if r1 < 0.5 {
          // Spiral movement towards best
          let dist = best_pos[d] - tuna[i][d]
          aa * dist + r2 * (best_pos[d] - tuna[i][d])
        } else {
          // Foraging behavior: random walk with best guidance
          tuna[i][d] + aa * z * (bounds.max - bounds.min) * (r3 - 0.5) * 2.0
        }
        new_pos.push(max(bounds.min, min(bounds.max, t)))
      }

      tuna[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

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

// ===== Golden Jackal Optimization (GJO) =====

///|
pub struct GJOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gjo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GJOConfig {
  GJOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn golden_jackal_optimization(
  config : GJOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let jackals : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    jackals.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  // Male and female jackals (top 2)
  let mut best_pos : Array[Double] = jackals[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let mut second_pos : Array[Double] = jackals[sorted_idx[1]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Energy decreases over iterations
    let e1 = 2.0 * (1.0 - iter)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)

        // Guided by male (best) and female (second) jackals
        let dist_male = best_pos[d] - jackals[i][d]
        let dist_female = second_pos[d] - jackals[i][d]
        // Average position with energy-based step size and random levy
        let pos = jackals[i][d] + e1 * r1 * dist_male + e1 * r2 * dist_female
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1

      if new_val < values[i] {
        jackals[i] = new_pos
        values[i] = new_val

        if new_val < best_val {
          second_pos = best_pos
          best_val = new_val
          best_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

///|
fn box_muller(u1 : Double, u2 : Double) -> Double {
  let safe_u1 = if u1 <= 0.0 { 1.0e-300 } else { u1 }
  let sqrt_val = (-2.0 * @math.ln(safe_u1)).sqrt()
  sqrt_val * @math.cos(2.0 * PI * u2)
}

///|
fn sort_indices(values : Array[Double]) -> Array[Int32] {
  let n = values.length()
  let indices : Array[Int32] = []
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

// ===== Grasshopper Optimization Algorithm (GOA) =====

///|
pub struct GOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_goa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GOAConfig {
  GOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
fn goa_s_func(r : Double) -> Double {
  let f = 0.5
  let l = 1.5
  f * exp_approx(-r / l) - exp_approx(-r)
}

///|
pub fn grasshopper_optimization_algorithm(
  config : GOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let grasshoppers : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    grasshoppers.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = grasshoppers[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Social influence decays quadratically so large budgets also converge
    let c = (1.0 - iter) * (1.0 - iter)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let mut sum_val = 0.0
        for j in 0..<n {
          if i != j {
            let dist = abs_val(grasshoppers[i][d] - grasshoppers[j][d])
            let s_val = goa_s_func(dist + 0.001)
            let dir = grasshoppers[j][d] - grasshoppers[i][d]
            sum_val = sum_val + s_val * dir / (dist + 0.001)
          }
        }
        let pos = c * sum_val + best_pos[d]
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        grasshoppers[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Seagull Optimization Algorithm (SOA) =====

///|
pub struct SOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_soa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SOAConfig {
  SOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn seagull_optimization_algorithm(
  config : SOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let seagulls : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    seagulls.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = seagulls[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let fc = 2.0 - iter * 2.0

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let b = r1 * fc
        let a = fc - r2 * b
        let spiral = a *
          cos_approx(2.0 * PI * r3) *
          (best_pos[d] - seagulls[i][d])
        let pos = b * (best_pos[d] - seagulls[i][d]) + spiral + seagulls[i][d]
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        seagulls[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Cat Swarm Optimization (CSO) =====

///|
pub struct CSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  mixture_ratio : Double
}

///|
pub fn make_cso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CSOConfig {
  CSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    mixture_ratio: 0.3,
  }
}

///|
pub fn cat_swarm_optimization(
  config : CSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let mr = config.mixture_ratio
  let current_rng = rng

  let cats : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    cats.push(sol)
    let vel : Array[Double] = []
    for _d in 0..<dim {
      let v = RNG::next(current_rng)
      vel.push((v - 0.5) * (bounds.max - bounds.min) * 0.1)
    }
    velocities.push(vel)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = cats[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let r = RNG::next(current_rng)
      let is_tracing = r < mr
      if is_tracing {
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let w = 0.8
          let c1 = 1.5
          let r1 = RNG::next(current_rng)
          let new_vel = w * velocities[i][d] +
            c1 * r1 * (best_pos[d] - cats[i][d])
          let pos = cats[i][d] + new_vel
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
          velocities[i][d] = max(
            bounds.min * 0.1,
            min(bounds.max * 0.1, new_vel),
          )
        }
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        cats[i] = new_pos
        values[i] = new_val
        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      } else {
        let srd = 0.2 * (bounds.max - bounds.min)
        let cdc = dim / 2
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r2 = RNG::next(current_rng)
          let mut p = cats[i][d]
          if d < cdc {
            p = p + (r2 - 0.5) * srd
          }
          new_pos.push(max(bounds.min, min(bounds.max, p)))
        }
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        if new_val < values[i] {
          cats[i] = new_pos
          values[i] = new_val
          if new_val < best_val {
            best_val = new_val
            best_pos = new_pos
          }
        }
      }
    }
  }
  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Gorilla Troops Optimizer (GTO) =====

///|
pub struct GTOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gto_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GTOConfig {
  GTOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn gorilla_troops_optimizer(
  config : GTOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let gorillas : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    gorillas.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = gorillas[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let c = 1.0 - iter

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r4 = RNG::next(current_rng)
        let rand_val = if r1 > 0.5 { r4 } else { -r4 }
        let pos = if r2 > 0.5 {
          gorillas[i][d] + c * rand_val * (best_pos[d] - gorillas[i][d])
        } else {
          let j = (r3 * n.to_double()).to_int() % n
          gorillas[i][d] + c * rand_val * (gorillas[j][d] - gorillas[i][d])
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        gorillas[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Slime Mould Algorithm (SMA) =====

///|
pub struct SMAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sma_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SMAConfig {
  SMAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

// ===== Aquila Optimizer (AO) =====

///|
pub struct AOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ao_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AOConfig {
  AOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

// ===== Chimp Optimization Algorithm (ChOA) =====

///|
pub struct ChOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_choa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ChOAConfig {
  ChOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

// ===== Dandelion Optimizer (DO) =====

///|
pub struct DOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_do_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DOConfig {
  DOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

// ===== Zebra Optimization Algorithm (ZOA) =====

///|
pub struct ZOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_zoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ZOAConfig {
  ZOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

// ===== Water Cycle Algorithm (WCA) =====

///|
pub struct WCAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  num_rivers : Int32
}

///|
pub fn make_wca_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  num_riv : Int32,
) -> WCAConfig {
  WCAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    num_rivers: num_riv,
  }
}

///|
pub fn water_cycle_algorithm(
  config : WCAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let nr = config.num_rivers
  let current_rng = rng

  let raindrops : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    raindrops.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = raindrops[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  let river_indices : Array[Int32] = []
  for i in 1..<=nr {
    river_indices.push(sorted_idx[i])
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let target = if i == sorted_idx[0] {
          best_pos[d]
        } else if river_indices.contains(i) {
          best_pos[d]
        } else {
          let riv_idx = i % nr
          raindrops[river_indices[riv_idx]][d]
        }
        let intensity = 1.0 - iter
        let pos = raindrops[i][d] +
          r1 * intensity * (target - raindrops[i][d]) +
          (r2 - 0.5) * 0.1 * (bounds.max - bounds.min)
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        raindrops[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

///|
pub fn slime_mould_algorithm(
  config : SMAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let slime : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    slime.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = slime[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  // Standard SMA: rare uniform restarts (z ~ 0.03) plus best-guided
  // exploitation with a decaying noise term. The previous version had two
  // identical "move toward best" branches and a z reset to [0.5, 1] every
  // generation, so the population homogenized and froze on local optima.
  let z = 0.03

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let pos = if r1 < z {
          // Rare random restart within bounds
          let rr = RNG::next(current_rng)
          bounds.min + rr * (bounds.max - bounds.min)
        } else {
          // Exploitation: converge toward best, noise fades out by mid-run
          slime[i][d] +
          r3 * (best_pos[d] - slime[i][d]) +
          (r2 - 0.5) *
          (bounds.max - bounds.min) *
          0.05 *
          max(0.0, 1.0 - 2.0 * iter)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        slime[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

///|
pub fn aquila_optimizer(
  config : AOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let eagles : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    eagles.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = eagles[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let r1 = RNG::next(current_rng)
    let r2 = RNG::next(current_rng)
    let r3 = RNG::next(current_rng)
    let _r4 = RNG::next(current_rng)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let pos = if iter < 1.0 / 3.0 {
          // exploration type 1
          let y = r1 * (bounds.max - bounds.min) + bounds.min
          eagles[i][d] - y
        } else if iter < 2.0 / 3.0 {
          // exploration type 2
          let y = r1 * (bounds.max - bounds.min) + bounds.min
          let rb1 = RNG::next(current_rng)
          let x = rb1 * (eagles[i][d] - y)
          x + r2 * (eagles[i][d] - best_pos[d])
        } else {
          // exploitation type 1
          let _rb2 = RNG::next(current_rng)
          let y = (best_pos[d] + eagles[i][d]) / 2.0
          eagles[i][d] - r3 * (y - eagles[i][d])
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        eagles[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

///|
pub fn chimp_optimization_algorithm(
  config : ChOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let chimps : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    chimps.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = chimps[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let _r1 = RNG::next(current_rng)
    let _r2 = RNG::next(current_rng)
    let _r3 = RNG::next(current_rng)
    let _r4 = RNG::next(current_rng)

    let r5 = RNG::next(current_rng)
    let _r6 = RNG::next(current_rng)
    let r7 = RNG::next(current_rng)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = if iter < 0.5 && r5 < 0.5 {
        // Driving: evaluate two full candidate solutions, keep the better one.
        // (The previous 1-d eval crashed on functions accessing x[1] and
        // skipped the evaluation budget.)
        let cand_a : Array[Double] = []
        let cand_b : Array[Double] = []
        for d in 0..<dim {
          let rd1 = RNG::next(current_rng)
          let rd2 = RNG::next(current_rng)
          cand_a.push(
            max(
              bounds.min,
              min(bounds.max, best_pos[d] + rd1 * (best_pos[d] - chimps[i][d])),
            ),
          )
          cand_b.push(
            max(
              bounds.min,
              min(bounds.max, best_pos[d] + rd2 * (best_pos[d] - chimps[i][d])),
            ),
          )
        }
        let val_a = evaluate(f, cand_a)
        let val_b = evaluate(f, cand_b)
        evaluations = evaluations + 2
        if val_a < val_b {
          cand_a
        } else {
          cand_b
        }
      } else {
        let arr : Array[Double] = []
        for d in 0..<dim {
          let pos = if iter < 0.5 {
            // chasing chimp
            let rc1 = RNG::next(current_rng)
            chimps[i][d] + rc1 * (best_pos[d] - chimps[i][d])
          } else {
            // exploitation
            chimps[i][d] + r7 * (best_pos[d] - chimps[i][d])
          }
          arr.push(max(bounds.min, min(bounds.max, pos)))
        }
        arr
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        chimps[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

///|
pub fn dandelion_optimizer(
  config : DOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let dandelions : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    dandelions.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = dandelions[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        // Two-phase: seeds drift randomly while ascending (first half), then
        // descend onto the incumbent best (second half, with a small fading
        // perturbation for diversity). The previous form used .abs() so the
        // drift was always positive and its step went negative after iter>0.5.
        // r2/r3 are drawn per component; sharing a single pair across the whole
        // population and every dimension made the seeds move rigidly together
        // and starved the search of diversity.
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let pos = if iter < 0.5 {
          dandelions[i][d] +
          (r2 - 0.5) *
          2.0 *
          (bounds.max - bounds.min) *
          0.15 *
          (1.0 - 2.0 * iter)
        } else {
          dandelions[i][d] +
          r3 * (best_pos[d] - dandelions[i][d]) +
          (r2 - 0.5) * 2.0 * (bounds.max - bounds.min) * 0.02 * (1.0 - iter)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        dandelions[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

///|
pub fn zebra_optimization_algorithm(
  config : ZOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let zebras : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    zebras.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = zebras[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let r1 = RNG::next(current_rng)
    let r2 = RNG::next(current_rng)
    let r3 = RNG::next(current_rng)
    let r4 = RNG::next(current_rng)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let pos = if r1 < 0.5 {
          // exploration
          let r5 = RNG::next(current_rng)
          let step = 2.0 * (bounds.max - bounds.min) * (0.5 - iter)
          zebras[i][d] + r2 * step * (r5 - 0.5)
        } else {
          // exploitation
          zebras[i][d] +
          r3 * (best_pos[d] - zebras[i][d]) +
          r4 * 0.1 * (bounds.max - bounds.min)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        zebras[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Manta Ray Foraging Optimization (MRFO) =====

///|
pub struct MRFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mrfo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MRFOConfig {
  MRFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Manta Ray Foraging Optimization (2020): each individual randomly performs
/// chain foraging or cyclone foraging, then every individual does a
/// somersault toward the best at the end of the generation.
pub fn manta_ray_foraging_optimization(
  config : MRFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng
  let range = bounds.max - bounds.min
  let somersault_factor = 2.0

  let rays : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    rays.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = rays[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Chain foraging or cyclone foraging, chosen per individual
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let choice = RNG::next(current_rng)
      let new_pos : Array[Double] = []
      if choice < 0.5 {
        // Cyclone foraging (spiral toward best, or toward a random point)
        let r1 = RNG::next(current_rng)
        let beta = 2.0 *
          exp_approx(r1 * (1.0 - iter) * sin_approx(2.0 * PI * r1))
        let phase = RNG::next(current_rng)
        let use_random = iter < phase
        for d in 0..<dim {
          let r2 = RNG::next(current_rng)
          let r3 = RNG::next(current_rng)
          let reference = if use_random {
            bounds.min + r2 * range
          } else {
            best_pos[d]
          }
          let pos = reference +
            r3 * (reference - rays[i][d]) +
            beta * (reference - rays[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Chain foraging: follow the previous ray toward the best
        let r_log = RNG::next(current_rng)
        let alpha = 2.0 * r_log * abs_val(@math.ln(max(r_log, 1.0e-10))).sqrt()
        for d in 0..<dim {
          let r2 = RNG::next(current_rng)
          let prev = if i == 0 { best_pos[d] } else { rays[i - 1][d] }
          let pos = rays[i][d] +
            r2 * (prev - rays[i][d]) +
            alpha * (best_pos[d] - rays[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }
      rays[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val
      if new_val < best_val {
        best_val = new_val
        best_pos = new_pos
      }
    }

    // Somersault foraging for every individual
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let pos = rays[i][d] +
          somersault_factor * (r2 * best_pos[d] - r3 * rays[i][d])
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      rays[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val
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

// ===== Hunger Games Search (HGS) =====

///|
pub struct HGSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hgs_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HGSConfig {
  HGSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Hunger Games Search (2021): movement is driven by a per-individual hunger
/// weight (adapted from the ranking), a contraction factor R that shrinks
/// with time, and a hunger effect E = sech(|f_i - gbest|).
pub fn hunger_games_search(
  config : HGSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng
  let range = bounds.max - bounds.min

  let individuals : Array[Array[Double]] = []
  let values : Array[Double] = []
  let hunger : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    individuals.push(sol)
    values.push(evaluate(f, sol))
    hunger.push(0.0)
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = individuals[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let shrink = 2.0 * (1.0 - iter)

    // Hunger accumulates except for the current best, which resets to zero
    let order = sort_indices(values)
    for i in 0..<n {
      if values[i] == values[order[0]] {
        hunger[i] = 0.0
      } else {
        hunger[i] = hunger[i] + 1.0
      }
    }
    let mut sum_hunger = 0.0
    for i in 0..<n {
      sum_hunger = sum_hunger + hunger[i]
    }
    let sum_denom = if sum_hunger > 0.0 { sum_hunger } else { 1.0 }

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let r_w = RNG::next(current_rng)
      let w_core = if hunger[i] <= 0.0 {
        1.0
      } else {
        (1.0 - hunger[i] / sum_denom) * r_w
      }
      let r_r = RNG::next(current_rng)
      let contraction = 2.0 * shrink * r_r - shrink
      let diff = abs_val(values[i] - best_val)
      let e_val = 2.0 / (exp_approx(diff) + exp_approx(-diff))

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let pos = if r1 < 0.5 {
          // Hunger-driven pull toward the best with a fading perturbation
          let g1 = RNG::next(current_rng)
          let g2 = RNG::next(current_rng)
          let gauss = box_muller(g1, g2)
          individuals[i][d] +
          (w_core + contraction) * (best_pos[d] - individuals[i][d]) +
          e_val * gauss * range * 0.05 * (1.0 - iter)
        } else {
          // Pull toward the best plus a hunger-weight pull from a peer
          let j = (RNG::next(current_rng) * n.to_double()).to_int() % n
          individuals[i][d] +
          w_core * (best_pos[d] - individuals[i][d]) +
          contraction * (individuals[j][d] - individuals[i][d])
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      individuals[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val
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

// ===== Reptile Search Algorithm (RSA) =====

///|
pub struct RSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RSAConfig {
  RSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Reptile Search Algorithm (2022): crocodile-inspired phases. The iteration
/// ratio selects high walking (encircling), belly walking, hunting or hunting
/// cooperation, with a shrinking factor driving the transition from
/// exploration to exploitation.
pub fn reptile_search_algorithm(
  config : RSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let current_rng = rng

  let reptiles : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    reptiles.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = reptiles[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let shrink = 2.0 * (1.0 - iter)

    // Per-dimension average position of the population
    let avg_pos : Array[Double] = []
    for _d in 0..<dim {
      avg_pos.push(0.0)
    }
    for i in 0..<n {
      for d in 0..<dim {
        avg_pos[d] = avg_pos[d] + reptiles[i][d]
      }
    }
    for d in 0..<dim {
      avg_pos[d] = avg_pos[d] / n.to_double()
    }

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      let ridx = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let eta = best_pos[d] - reptiles[i][d]
        let rr = (best_pos[d] - reptiles[ridx][d]) * r1 * shrink
        let beta = 1.0 + (best_pos[d] - avg_pos[d]) * r2
        // The later two phases multiply eta by a factor that vanishes as the
        // run progresses, so the reptile converges onto the best instead of
        // oscillating around it (which froze RSA on wide-range functions).
        let pos = if iter < 0.25 {
          // High walking (encircling, most exploratory)
          best_pos[d] - eta * beta - rr
        } else if iter < 0.5 {
          // Belly walking
          best_pos[d] + eta * beta - rr
        } else if iter < 0.75 {
          // Hunting
          best_pos[d] + eta * beta * (1.0 - iter) * 4.0
        } else {
          // Hunting cooperation (converges onto the best)
          best_pos[d] + eta * (1.0 - iter) * 4.0
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      reptiles[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val
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

// ===== Marine Predators Algorithm (MPA, 2020) =====

///|
pub struct MPAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mpa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MPAConfig {
  MPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Levy flight step (beta = 1.5), shared by the MPA foraging phases.
fn mpa_levy(rng : RNG) -> Double {
  let g1 = box_muller(RNG::next(rng), RNG::next(rng))
  let g2 = box_muller(RNG::next(rng), RNG::next(rng))
  0.05 * g1 / pow_approx(max(abs_val(g2), 1.0e-10), 1.0 / 1.5)
}

///|
/// Marine Predators Algorithm (2020): a top-predator "Elite" guides the prey
/// through three velocity regimes - Brownian exploration while the predator is
/// much faster than the prey, then a Levy/Brownian mix, then pure Levy
/// exploitation as the prey slows down.  The damping factor CF makes the last
/// phase decay to zero, so larger budgets still converge.
pub fn marine_predators_algorithm(
  config : MPAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let prey : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    prey.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = prey[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let half = n / 2

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    // Damping factor for the exploitation phase (fades to zero as iter -> 1)
    let cf = pow_approx(1.0 - iter, 2.0 * iter)

    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let new_pos : Array[Double] = []
      if iter < 1.0 / 3.0 {
        // Phase 1: Brownian motion, predator faster than prey
        for d in 0..<dim {
          let rb = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          let pos = prey[i][d] + 0.5 * rb * (best_pos[d] - rb * prey[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else if iter < 2.0 / 3.0 {
        if i < half {
          // Phase 2a: prey uses Levy flight
          for d in 0..<dim {
            let rl = mpa_levy(current_rng)
            let pos = prey[i][d] + 0.5 * rl * (best_pos[d] - rl * prey[i][d])
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        } else {
          // Phase 2b: predator uses Brownian motion
          for d in 0..<dim {
            let rb = box_muller(RNG::next(current_rng), RNG::next(current_rng))
            let pos = best_pos[d] + 0.5 * cf * (rb * best_pos[d] - prey[i][d])
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
        }
      } else {
        // Phase 3: predator uses Levy flight as the prey slows down
        for d in 0..<dim {
          let rl = mpa_levy(current_rng)
          let pos = best_pos[d] + 0.5 * cf * (rl * best_pos[d] - prey[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      }
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      // Greedy acceptance keeps the prey from being randomized by the moves.
      if new_val < values[i] {
        prey[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Golden Eagle Optimizer (GEO, 2021) =====

///|
pub struct GEOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_geo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GEOConfig {
  GEOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Golden Eagle Optimizer (2021): each eagle balances an attack move toward the
/// best (attack propensity A grows with time) against a cruise move along a
/// random direction orthogonal to the attack vector (cruise propensity
/// C = 1 - A shrinks).  Late in the run the cruise term vanishes, so the
/// population converges onto the incumbent best.
pub fn golden_eagle_optimizer(
  config : GEOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let eagles : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    eagles.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = eagles[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let tau = evaluations.to_double() / max_evals.to_double()
    let attack = tau
    let cruise = 1.0 - tau

    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }

      // Attack vector (toward the best) and its squared norm
      let mut sq = 0.0
      for d in 0..<dim {
        let diff = best_pos[d] - eagles[i][d]
        sq = sq + diff * diff
      }
      let norm_a = sq.sqrt()

      // Random cruise vector, orthogonalized against the attack vector
      let rc : Array[Double] = []
      for _d in 0..<dim {
        rc.push(RNG::next(current_rng) * 2.0 - 1.0)
      }
      let mut dot = 0.0
      for d in 0..<dim {
        dot = dot + rc[d] * (best_pos[d] - eagles[i][d])
      }
      let proj = dot / max(sq, 1.0e-12)
      let cruise_vec : Array[Double] = []
      let mut nc = 0.0
      for d in 0..<dim {
        let c = rc[d] - proj * (best_pos[d] - eagles[i][d])
        cruise_vec.push(c)
        nc = nc + c * c
      }
      nc = nc.sqrt()

      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let cd = if nc > 1.0e-12 { cruise_vec[d] / nc } else { 0.0 }
        let pos = eagles[i][d] +
          r1 * attack * (best_pos[d] - eagles[i][d]) +
          r2 * cruise * cd * norm_a
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      // Greedy acceptance keeps the swarm from being randomized
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        eagles[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Runge Kutta Optimizer (RUN, 2021) =====

///|
pub struct RUNConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_run_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RUNConfig {
  RUNConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Runge Kutta Optimizer (2021): the update is the classic RK4 weighted average
/// of four slope estimates, here formed from two guide individuals - a randomly
/// chosen peer (exploration) and the incumbent best (exploitation) - with the
/// guide switched at random each step.  Greedy acceptance and historical-best
/// tracking are added on top of the standard design.
pub fn runge_kutta_optimizer(
  config : RUNConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let individuals : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    individuals.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = individuals[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  let range = bounds.max - bounds.min

  while evaluations < max_evals {
    let tau = evaluations.to_double() / max_evals.to_double()
    // Exploration probability fades so larger budgets converge on the best
    let p_explore = 1.0 - 0.9 * tau
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }

      // Pick the two guide individuals
      let raw_a = (RNG::next(current_rng) * n.to_double()).floor().to_int()
      let idx_a = if raw_a >= n { n - 1 } else { raw_a }
      let raw_b = (RNG::next(current_rng) * n.to_double()).floor().to_int()
      let idx_b = if raw_b >= n { n - 1 } else { raw_b }
      let r_switch = RNG::next(current_rng)
      let guide_a : Array[Double] = if r_switch < p_explore {
        individuals[idx_a]
      } else {
        best_pos
      }
      let guide_b : Array[Double] = individuals[idx_b]

      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let k1 = r1 * (guide_a[d] - individuals[i][d])
        let k2 = r2 * (guide_b[d] - individuals[i][d])
        let k3 = r3 * (guide_a[d] - individuals[i][d])
        let k4 = r4 * (guide_b[d] - individuals[i][d])
        let dx = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        // Fading jitter keeps the swarm from freezing once it homogenizes
        let jitter = (RNG::next(current_rng) - 0.5) *
          2.0 *
          range *
          0.1 *
          (1.0 - tau)
        let pos = individuals[i][d] + dx + jitter
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      // Greedy acceptance keeps the swarm from being randomized
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        individuals[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== African Vultures Optimization Algorithm (AVOA, 2021) =====

///|
pub struct AVOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_avoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AVOAConfig {
  AVOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// African Vultures Optimization Algorithm (2021): a hunger factor F, which
/// fades with the budget, switches each vulture between exploration (|F| >= 1,
/// wide moves around a reference vulture) and exploitation (|F| < 1, moves that
/// close in on the reference).  The reference vulture is either the current best
/// or the runner-up, chosen at random.  Adds greedy acceptance and
/// historical-best tracking.
pub fn african_vultures_optimization(
  config : AVOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let vultures : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    vultures.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = vultures[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let tau = evaluations.to_double() / max_evals.to_double()
    // Current best and runner-up act as the two reference vultures
    let order = sort_indices(values)
    let cur_best = vultures[order[0]]
    let second_best = vultures[order[1]]

    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      // Hunger factor: fades with the budget so late phases exploit
      let z = 2.0 * r2 - 1.0
      let f_hunger = (2.0 * r1 + 1.0) * z * (1.0 - tau)
      let af = abs_val(f_hunger)
      let guide : Array[Double] = if r3 < 0.5 { cur_best } else { second_best }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let rnd = RNG::next(current_rng)
        let pos = if af >= 1.0 {
          // Exploration: wide moves around the reference vulture
          if r4 < 0.5 {
            guide[d] - f_hunger * (guide[d] - vultures[i][d])
          } else {
            guide[d] - f_hunger + (2.0 * rnd - 1.0) * range
          }
          // Exploitation: close in on the reference vulture
        } else if r4 < 0.5 {
          vultures[i][d] + af * (guide[d] - vultures[i][d])
        } else {
          guide[d] + (1.0 - af) * (2.0 * rnd - 1.0) * range * 0.1
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      // Greedy acceptance keeps the swarm from being randomized
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        vultures[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== White Shark Optimizer (WSO, 2022) =====

///|
pub struct WSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WSOConfig {
  WSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// White Shark Optimizer (2022): each shark carries a velocity that blends a
/// fading inertia term with an attraction toward the incumbent best (the prey).
/// A per-individual switch - whose "follow the prey" probability grows with the
/// budget - alternates between following the prey and a wider search move.
/// Adds greedy acceptance and historical-best tracking.
pub fn white_shark_optimizer(
  config : WSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let sharks : Array[Array[Double]] = []
  let velocity : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    sharks.push(sol)
    values.push(evaluate(f, sol))
    let v : Array[Double] = []
    for _d in 0..<dim {
      v.push(0.0)
    }
    velocity.push(v)
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = sharks[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let tau = evaluations.to_double() / max_evals.to_double()
    // Probability of following the prey grows so late phases converge
    let follow_prob = 0.3 + 0.7 * tau

    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let inertia = 0.5 * velocity[i][d] * (1.0 - tau)
        let attract = (0.5 * r1 + 0.5 * r2) * (best_pos[d] - sharks[i][d])
        let vel = inertia + attract
        velocity[i][d] = vel
        let pos = if r3 < follow_prob {
          // Follow the prey
          sharks[i][d] + vel
        } else {
          // Wider search move with a fading random component
          let rnd = RNG::next(current_rng)
          sharks[i][d] - vel + (2.0 * rnd - 1.0) * range * 0.1 * (1.0 - tau)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      // Greedy acceptance keeps the swarm from being randomized
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if new_val < values[i] {
        sharks[i] = new_pos
        values[i] = new_val
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
    evaluations,
  }
}

// ===== Coati Optimization Algorithm (COA, 2023) =====

///|
pub struct COAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_coa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> COAConfig {
  COAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Coati Optimization Algorithm (2023): the swarm hunts iguanas in two
/// sub-groups - half climb a tree after a randomly chosen iguana while the
/// others wait for one to fall to a random ground point - then every coati
/// performs a local escape move whose radius shrinks with the generation count.
/// Adds greedy acceptance and historical-best tracking.
pub fn coati_optimization_algorithm(
  config : COAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let half = n / 2

  let coatis : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    coatis.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = coatis[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }

      // Phase 1: hunt the iguana
      let raw = (RNG::next(current_rng) * n.to_double()).floor().to_int()
      let idx = if raw >= n { n - 1 } else { raw }
      let iguana = coatis[idx]
      let r = RNG::next(current_rng)
      let i_factor = if RNG::next(current_rng) < 0.5 { 1.0 } else { 2.0 }

      let new1 : Array[Double] = []
      for d in 0..<dim {
        let pos = if i < half {
          // Climb the tree and attack the iguana
          coatis[i][d] + r * (iguana[d] - i_factor * coatis[i][d])
        } else {
          // The iguana falls to a random position on the ground
          let g = bounds.min + RNG::next(current_rng) * range
          coatis[i][d] + r * (g - i_factor * coatis[i][d])
        }
        new1.push(max(bounds.min, min(bounds.max, pos)))
      }
      if evaluations < max_evals {
        let v1 = evaluate(f, new1)
        evaluations = evaluations + 1
        if v1 < values[i] {
          coatis[i] = new1
          values[i] = v1
          if v1 < best_val {
            best_val = v1
            best_pos = new1
          }
        }
      }

      // Phase 2: escape from a predator, with a shrinking local radius
      if evaluations < max_evals {
        let t = max(evaluations.to_double() / n.to_double(), 1.0)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let new2 : Array[Double] = []
        for d in 0..<dim {
          let lb_local = (bounds.min - coatis[i][d]) / t
          let ub_local = (bounds.max - coatis[i][d]) / t
          let pos = coatis[i][d] + r2 * (lb_local + r3 * (ub_local - lb_local))
          new2.push(max(bounds.min, min(bounds.max, pos)))
        }
        let v2 = evaluate(f, new2)
        evaluations = evaluations + 1
        if v2 < values[i] {
          coatis[i] = new2
          values[i] = v2
          if v2 < best_val {
            best_val = v2
            best_pos = new2
          }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== jDE: Self-adaptive Differential Evolution (2006) =====

///|
pub struct JDEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_jde_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> JDEConfig {
  JDEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// jDE (Brest et al., 2006): differential evolution in which every individual
/// carries its own scale factor F and crossover rate CR.  Before each trial the
/// parameters are occasionally re-sampled (probabilities tau1/tau2); the adapted
/// values are kept only if the trial succeeds and otherwise discarded - this
/// self-adaptation is what made jDE a landmark algorithm.
pub fn jde_algorithm(
  config : JDEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let pop_size = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let tau1 = 0.1
  let tau2 = 0.1
  let f_lower = 0.1
  let f_upper = 0.9

  let population : Array[Array[Double]] = []
  let fitness : Array[Double] = []
  let f_scale : Array[Double] = []
  let c_rate : Array[Double] = []
  for _i in 0..<pop_size {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    fitness.push(evaluate(f, sol))
    f_scale.push(0.5)
    c_rate.push(0.9)
  }
  let mut evaluations = pop_size
  let sorted_idx = sort_indices(fitness)
  let mut best_pos : Array[Double] = population[sorted_idx[0]]
  let mut best_val = fitness[sorted_idx[0]]

  while evaluations < max_evals {
    for target in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }

      // Self-adapt the control parameters
      let mut f_i = f_scale[target]
      let mut cr_i = c_rate[target]
      if RNG::next(current_rng) < tau1 {
        f_i = f_lower + RNG::next(current_rng) * (f_upper - f_lower)
      }
      if RNG::next(current_rng) < tau2 {
        cr_i = RNG::next(current_rng)
      }

      // Three distinct random indices different from `target`
      let mut a = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      while a == target {
        a = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      let mut b = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      while b == target || b == a {
        b = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      let mut c = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      while c == target || c == a || c == b {
        c = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }

      // DE/rand/1/bin trial
      let trial : Array[Double] = []
      for j in 0..<dim {
        let mutant = population[a][j] +
          f_i * (population[b][j] - population[c][j])
        let clamped = max(bounds.min, min(bounds.max, mutant))
        if RNG::next(current_rng) < cr_i {
          trial.push(clamped)
        } else {
          trial.push(population[target][j])
        }
      }

      let trial_fit = evaluate(f, trial)
      evaluations = evaluations + 1
      if trial_fit < fitness[target] {
        population[target] = trial
        fitness[target] = trial_fit
        // Keep the adapted parameters on success, discard them on failure
        f_scale[target] = f_i
        c_rate[target] = cr_i
        if trial_fit < best_val {
          best_val = trial_fit
          best_pos = trial
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SHADE: Success-History based Adaptive DE (2013) =====

///|
pub struct SHADEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_shade_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SHADEConfig {
  SHADEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Standard Cauchy sample (location 0, scale 1).
fn shade_cauchy(rng : RNG) -> Double {
  let u = RNG::next(rng)
  let t = 3.141592653589793 * (u - 0.5)
  let s = sin_approx(t)
  let c = cos_approx(t)
  let c2 = if abs_val(c) < 1.0e-8 {
    if c < 0.0 {
      -1.0e-8
    } else {
      1.0e-8
    }
  } else {
    c
  }
  s / c2
}

///|
/// SHADE (Tanabe & Fukunaga, 2013): current-to-pbest/1 mutation with
/// success-history memories M_F / M_CR, an external archive of replaced
/// parents and a p-best pool.  Successful parameters are folded back into one
/// memory slot per generation as weighted Lehmer/arithmetic means.  F is drawn
/// from a Cauchy and CR from a normal centred on the memory values.
pub fn shade_algorithm(
  config : SHADEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let pop_size = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let h_mem = 5

  let m_f : Array[Double] = []
  let m_cr : Array[Double] = []
  for _i in 0..<h_mem {
    m_f.push(0.5)
    m_cr.push(0.5)
  }
  let archive : Array[Array[Double]] = []
  let mut mem_idx = 0

  let population : Array[Array[Double]] = []
  let fitness : Array[Double] = []
  for _i in 0..<pop_size {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    fitness.push(evaluate(f, sol))
  }
  let mut evaluations = pop_size
  let sorted0 = sort_indices(fitness)
  let mut best_pos : Array[Double] = population[sorted0[0]]
  let mut best_val = fitness[sorted0[0]]

  while evaluations < max_evals {
    let order = sort_indices(fitness)

    // p-best pool size: p drawn from [2/N, 0.2]
    let p_min = 2.0 / pop_size.to_double()
    let p = p_min + RNG::next(current_rng) * (0.2 - p_min)
    let mut pnum = (p * pop_size.to_double()).to_int()
    if pnum < 2 {
      pnum = 2
    }
    if pnum > pop_size {
      pnum = pop_size
    }

    let s_f : Array[Double] = []
    let s_cr : Array[Double] = []
    let s_df : Array[Double] = []

    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }

      // Sample F (Cauchy) and CR (normal) from the success-history memory
      let mem_k = (RNG::next(current_rng) * h_mem.to_double()).to_int() % h_mem
      let mut cr_i = m_cr[mem_k] +
        0.1 * box_muller(RNG::next(current_rng), RNG::next(current_rng))
      if cr_i < 0.0 {
        cr_i = 0.0
      } else if cr_i > 1.0 {
        cr_i = 1.0
      }
      let mut f_i = 0.0
      let mut tries = 0
      while f_i <= 0.0 && tries < 10 {
        f_i = m_f[mem_k] + 0.1 * shade_cauchy(current_rng)
        tries = tries + 1
      }
      if f_i <= 0.0 {
        f_i = 0.5
      }
      if f_i > 1.0 {
        f_i = 1.0
      }

      // p-best individual
      let pbest_rank = (RNG::next(current_rng) * pnum.to_double()).to_int() %
        pnum
      let pbest = population[order[pbest_rank]]

      // r1 from the population
      let mut r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      while r1 == i {
        r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }

      // r2 from population ∪ archive
      let union_size = pop_size + archive.length()
      let mut r2 = (RNG::next(current_rng) * union_size.to_double()).to_int() %
        union_size
      while r2 == i || r2 == r1 {
        r2 = (RNG::next(current_rng) * union_size.to_double()).to_int() %
          union_size
      }
      let r2_pos : Array[Double] = if r2 < pop_size {
        population[r2]
      } else {
        archive[r2 - pop_size]
      }

      // current-to-pbest/1 with binomial crossover
      let jrand = (RNG::next(current_rng) * dim.to_double()).to_int() % dim
      let trial : Array[Double] = []
      for j in 0..<dim {
        let mutant = population[i][j] +
          f_i * (pbest[j] - population[i][j]) +
          f_i * (population[r1][j] - r2_pos[j])
        let clamped = max(bounds.min, min(bounds.max, mutant))
        if RNG::next(current_rng) < cr_i || j == jrand {
          trial.push(clamped)
        } else {
          trial.push(population[i][j])
        }
      }

      let trial_fit = evaluate(f, trial)
      evaluations = evaluations + 1
      if trial_fit < fitness[i] {
        s_f.push(f_i)
        s_cr.push(cr_i)
        s_df.push(abs_val(fitness[i] - trial_fit))
        // External archive of replaced parents (random replacement)
        if archive.length() < pop_size {
          archive.push(population[i])
        } else {
          let ridx = (RNG::next(current_rng) * archive.length().to_double()).to_int() %
            archive.length()
          archive[ridx] = population[i]
        }
        population[i] = trial
        fitness[i] = trial_fit
        if trial_fit < best_val {
          best_val = trial_fit
          best_pos = trial
        }
      }
    }

    // Fold this generation's successful parameters into one memory slot
    let m = s_df.length()
    if m > 0 {
      let mut sum_w = 0.0
      let mut sum_wf = 0.0
      let mut sum_wf2 = 0.0
      let mut sum_wcr = 0.0
      for t in 0..<m {
        let w = s_df[t]
        sum_w = sum_w + w
        sum_wf = sum_wf + w * s_f[t]
        sum_wf2 = sum_wf2 + w * s_f[t] * s_f[t]
        sum_wcr = sum_wcr + w * s_cr[t]
      }
      if sum_wf > 1.0e-12 {
        m_f[mem_idx] = sum_wf2 / sum_wf
        m_cr[mem_idx] = sum_wcr / sum_w
      } else {
        m_f[mem_idx] = 0.5
        m_cr[mem_idx] = 0.5
      }
      mem_idx = (mem_idx + 1) % h_mem
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EP: Evolutionary Programming (self-adaptive meta-EP) =====

///|
pub struct EPConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ep_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EPConfig {
  EPConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Evolutionary Programming with self-adaptive (meta-EP) mutation: every
/// individual carries a step size that is itself mutated by log-normal noise
/// before the position mutation, and (mu+mu) selection keeps the best parents
/// and children.  Adds historical-best tracking.
pub fn evolutionary_programming(
  config : EPConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let pop_size = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let dd = dim.to_double()
  let tau_prime = 1.0 / (2.0 * dd).sqrt()
  let tau = 1.0 / (2.0 * dd.sqrt()).sqrt()

  let mut population : Array[Array[Double]] = []
  let mut fitness : Array[Double] = []
  let mut sigma : Array[Double] = []
  for _i in 0..<pop_size {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    fitness.push(evaluate(f, sol))
    sigma.push(range * 0.2)
  }
  let mut evaluations = pop_size
  let sorted0 = sort_indices(fitness)
  let mut best_pos : Array[Double] = population[sorted0[0]]
  let mut best_val = fitness[sorted0[0]]

  while evaluations < max_evals {
    let children : Array[Array[Double]] = []
    let child_fit : Array[Double] = []
    let child_sigma : Array[Double] = []

    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      // Self-adapt the step size (log-normal mutation)
      let g1 = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      let g2 = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      let mut sig = sigma[i] * exp_approx(tau_prime * g1 + tau * g2)
      if sig > range {
        sig = range
      }
      if sig < 1.0e-12 {
        sig = 1.0e-12
      }
      // Mutate the position
      let trial : Array[Double] = []
      for j in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        let pos = population[i][j] + sig * g
        trial.push(max(bounds.min, min(bounds.max, pos)))
      }
      let val = evaluate(f, trial)
      evaluations = evaluations + 1
      children.push(trial)
      child_fit.push(val)
      child_sigma.push(sig)
      if val < best_val {
        best_val = val
        best_pos = trial
      }
    }

    // (mu + mu) selection over parents and children
    let all_x : Array[Array[Double]] = []
    let all_f : Array[Double] = []
    let all_s : Array[Double] = []
    for i in 0..<pop_size {
      all_x.push(population[i])
      all_f.push(fitness[i])
      all_s.push(sigma[i])
    }
    for k in 0..<children.length() {
      all_x.push(children[k])
      all_f.push(child_fit[k])
      all_s.push(child_sigma[k])
    }
    let order = sort_indices(all_f)
    let next_pop : Array[Array[Double]] = []
    let next_fit : Array[Double] = []
    let next_sig : Array[Double] = []
    for r in 0..<pop_size {
      let idx = order[r]
      next_pop.push(all_x[idx])
      next_fit.push(all_f[idx])
      next_sig.push(all_s[idx])
    }
    population = next_pop
    fitness = next_fit
    sigma = next_sig
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== OOA: Osprey Optimization Algorithm (2023) =====

///|
pub struct OOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ooa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> OOAConfig {
  OOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Osprey Optimization Algorithm (2023): in phase 1 each osprey dives at a
/// randomly chosen fish (another individual) and in phase 2 carries the fish to
/// a suitable position, modelled here as an attraction toward the incumbent best
/// with a shrinking random component.  Both updates are greedy.  Adds
/// historical-best tracking.
pub fn osprey_optimization_algorithm(
  config : OOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let ospreys : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    ospreys.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = ospreys[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }

      // Phase 1: identify and attack a random fish (another individual)
      let fish_idx = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let fish = ospreys[fish_idx]
      let r1 = RNG::next(current_rng)
      let i_factor = if RNG::next(current_rng) < 0.5 { 1.0 } else { 2.0 }
      let new1 : Array[Double] = []
      for d in 0..<dim {
        let pos = ospreys[i][d] + r1 * (fish[d] - i_factor * ospreys[i][d])
        new1.push(max(bounds.min, min(bounds.max, pos)))
      }
      if evaluations < max_evals {
        let v1 = evaluate(f, new1)
        evaluations = evaluations + 1
        if v1 < values[i] {
          ospreys[i] = new1
          values[i] = v1
          if v1 < best_val {
            best_val = v1
            best_pos = new1
          }
        }
      }

      // Phase 2: carry the fish to a suitable position (shrink toward the best)
      if evaluations < max_evals {
        let r2 = RNG::next(current_rng)
        let new2 : Array[Double] = []
        for d in 0..<dim {
          let jitter = (RNG::next(current_rng) - 0.5) *
            2.0 *
            range *
            0.1 *
            (1.0 - iter)
          let pos = ospreys[i][d] + r2 * (best_pos[d] - ospreys[i][d]) + jitter
          new2.push(max(bounds.min, min(bounds.max, pos)))
        }
        let v2 = evaluate(f, new2)
        evaluations = evaluations + 1
        if v2 < values[i] {
          ospreys[i] = new2
          values[i] = v2
          if v2 < best_val {
            best_val = v2
            best_pos = new2
          }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== TDO: Tasmanian Devil Optimization (2022) =====

///|
pub struct TDOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_tdo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TDOConfig {
  TDOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Tasmanian Devil Optimization (2022): each devil picks one of two feeding
/// strategies - scavenging (move toward a random point of the search space) or
/// hunting (move toward a randomly chosen prey/individual) - and then drags the
/// food to a safe position, modelled as an attraction toward the incumbent best
/// with a shrinking random component.  Greedy acceptance and historical-best
/// tracking are added.
pub fn tasmanian_devil_optimization(
  config : TDOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let devils : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    devils.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = devils[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }

      // Feeding strategy: scavenging vs hunting
      let prey_idx = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let prey = devils[prey_idx]
      let r_strategy = RNG::next(current_rng)
      let r1 = RNG::next(current_rng)
      let i_factor = if RNG::next(current_rng) < 0.5 { 1.0 } else { 2.0 }
      let new1 : Array[Double] = []
      for d in 0..<dim {
        let target = if r_strategy < 0.5 {
          // Scavenging: a random position in the search space
          bounds.min + RNG::next(current_rng) * range
        } else {
          // Hunting: the chosen prey
          prey[d]
        }
        let pos = devils[i][d] + r1 * (target - i_factor * devils[i][d])
        new1.push(max(bounds.min, min(bounds.max, pos)))
      }
      if evaluations < max_evals {
        let v1 = evaluate(f, new1)
        evaluations = evaluations + 1
        if v1 < values[i] {
          devils[i] = new1
          values[i] = v1
          if v1 < best_val {
            best_val = v1
            best_pos = new1
          }
        }
      }

      // Drag the food to a safe position (shrink toward the best)
      if evaluations < max_evals {
        let r2 = RNG::next(current_rng)
        let new2 : Array[Double] = []
        for d in 0..<dim {
          let jitter = (RNG::next(current_rng) - 0.5) *
            2.0 *
            range *
            0.1 *
            (1.0 - iter)
          let pos = devils[i][d] + r2 * (best_pos[d] - devils[i][d]) + jitter
          new2.push(max(bounds.min, min(bounds.max, pos)))
        }
        let v2 = evaluate(f, new2)
        evaluations = evaluations + 1
        if v2 < values[i] {
          devils[i] = new2
          values[i] = v2
          if v2 < best_val {
            best_val = v2
            best_pos = new2
          }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== PDO: Prairie Dog Optimization (2022) =====

///|
pub struct PDOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pdo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PDOConfig {
  PDOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Prairie Dog Optimization (2022): four behavioural phases selected by the
/// budget fraction - foraging toward a random peer, burrowing toward the
/// incumbent best, then two increasingly local searches around the best.  Each
/// phase keeps a shrinking random component so larger budgets converge.  Greedy
/// acceptance and historical-best tracking are added.
pub fn prairie_dog_optimization(
  config : PDOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let dogs : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    dogs.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = dogs[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let peer_idx = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let peer = dogs[peer_idx]
      let r1 = RNG::next(current_rng)
      let i_factor = if RNG::next(current_rng) < 0.5 { 1.0 } else { 2.0 }
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let rnd = RNG::next(current_rng)
        let pos = if iter < 0.25 {
          // Phase 1: foraging - move toward a random peer
          dogs[i][d] + r1 * (peer[d] - i_factor * dogs[i][d])
        } else if iter < 0.5 {
          // Phase 2: burrowing - move toward the incumbent best
          dogs[i][d] + r1 * (best_pos[d] - i_factor * dogs[i][d])
        } else if iter < 0.75 {
          // Phase 3: alarm - wider random search around the best
          best_pos[d] + (2.0 * rnd - 1.0) * range * 0.2 * (1.0 - iter)
        } else {
          // Phase 4: peaceful - fine search around the best
          best_pos[d] + (2.0 * rnd - 1.0) * range * 0.05 * (1.0 - iter)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if val < values[i] {
        dogs[i] = new_pos
        values[i] = val
        if val < best_val {
          best_val = val
          best_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== NGO: Northern Goshawk Optimization (2021) =====

///|
pub struct NGOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ngo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NGOConfig {
  NGOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Northern Goshawk Optimization (2021): each hawk alternates between an
/// exploration phase (dive at a randomly chosen prey) and an exploitation
/// phase (tail-chase with a shrinking radius).  Greedy acceptance and
/// historical-best tracking are added, and the chase radius decays with the
/// budget fraction so that larger budgets keep converging.
pub fn northern_goshawk_optimization(
  config : NGOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let hawks : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    hawks.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = hawks[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      // Phase 1: exploration - dive at a random prey
      let prey = hawks[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let r1 = RNG::next(current_rng)
      let i_factor = if RNG::next(current_rng) < 0.5 { 1.0 } else { 2.0 }
      let p1 : Array[Double] = []
      for d in 0..<dim {
        let pos = hawks[i][d] + r1 * (prey[d] - i_factor * hawks[i][d])
        p1.push(max(bounds.min, min(bounds.max, pos)))
      }
      let v1 = evaluate(f, p1)
      evaluations = evaluations + 1
      if v1 < values[i] {
        hawks[i] = p1
        values[i] = v1
        if v1 < best_val {
          best_val = v1
          best_pos = p1
        }
      }
      // Phase 2: exploitation - tail chase with a decaying radius
      if evaluations >= max_evals {
        break
      }
      let radius = 0.02 * (1.0 - iter)
      let p2 : Array[Double] = []
      for d in 0..<dim {
        let rnd = RNG::next(current_rng)
        let pos = hawks[i][d] + radius * (2.0 * rnd - 1.0) * hawks[i][d]
        p2.push(max(bounds.min, min(bounds.max, pos)))
      }
      let v2 = evaluate(f, p2)
      evaluations = evaluations + 1
      if v2 < values[i] {
        hawks[i] = p2
        values[i] = v2
        if v2 < best_val {
          best_val = v2
          best_pos = p2
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BWO: Beluga Whale Optimization (2022) =====

///|
pub struct BWOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bwo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BWOConfig {
  BWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Beluga Whale Optimization (2022): swimmers explore via a pairwise dive with
/// a sinusoidal weighting, feeders exploit around the incumbent best, and the
/// whale-fall probability decays with the budget fraction so that larger
/// budgets settle into local refinement.  Greedy acceptance and historical-best
/// tracking are added.
pub fn beluga_whale_optimization(
  config : BWOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let whales : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    whales.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = whales[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let fall_prob = 0.1 * (1.0 - iter)
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let r1 = RNG::next(current_rng)
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let curr = whales[i][d]
        let pos = if r1 < fall_prob {
          // whale fall: drift toward the best with a random recoil
          curr + r2 * (best_pos[d] - curr) - r3 * (whales[r][d] - curr)
        } else if r1 < 0.5 {
          // exploration: swim toward a peer with sinusoidal weighting
          curr + (whales[r][d] - curr) * (1.0 + r2) * sin_approx(TWO_PI * r3)
        } else {
          // exploitation: feed around the incumbent best
          best_pos[d] + r2 * (best_pos[d] - curr) + r3 * (whales[r][d] - curr)
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if val < values[i] {
        whales[i] = new_pos
        values[i] = val
        if val < best_val {
          best_val = val
          best_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GBO: Gradient-Based Optimizer (2020) =====

///|
pub struct GBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GBOConfig {
  GBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Gradient-Based Optimizer (2020): each individual is steered by a
/// gradient-search rule (GSR) built from the directional difference between
/// the incumbent best and worst members, with an occasional local escaping
/// operator (LEO) that samples around the best.  The GSR step size shrinks with
/// the budget fraction so larger budgets refine instead of re-randomising.
/// Greedy acceptance and historical-best tracking are added.
pub fn gradient_based_optimizer(
  config : GBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let eps = 1.0e-8
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let delta = 2.0 * (1.0 - iter)
    let order = sort_indices(values)
    let xbest = pop[order[0]]
    let xworst = pop[order[order.length() - 1]]
    let leo_prob = 0.1 * (1.0 - iter)
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let mu1 = 2.0 * delta * RNG::next(current_rng) - delta
      let mu2 = 2.0 * delta * RNG::next(current_rng) - delta
      let use_leo = RNG::next(current_rng) < leo_prob
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let denom = xworst[d] - xbest[d]
        let denom_safe = if abs_val(denom) < eps {
          if denom < 0.0 {
            -eps
          } else {
            eps
          }
        } else {
          denom
        }
        let pos = if use_leo {
          // local escaping: sample around the incumbent best
          best_pos[d] +
          (2.0 * RNG::next(current_rng) - 1.0) * delta * range * 0.1
        } else {
          // gradient-search rule
          let dx = 2.0 * RNG::next(current_rng) * abs_val(xbest[d] - xworst[d])
          pop[i][d] -
          mu1 * (2.0 * dx * pop[i][d]) / denom_safe +
          mu2 * (xbest[d] - pop[i][d])
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if val < values[i] {
        pop[i] = new_pos
        values[i] = val
        if val < best_val {
          best_val = val
          best_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GO: Growth Optimizer (2023) =====

///|
pub struct GOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_go_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GOConfig {
  GOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Growth Optimizer (2023): each member first runs a learning phase that blends
/// attraction to the incumbent best with a peer difference, then a reflection
/// phase that samples a point projected beyond the best.  The reflection
/// amplitude decays with the budget fraction so larger budgets settle into
/// local refinement.  Greedy acceptance and historical-best tracking are added.
pub fn growth_optimizer(
  config : GOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let reflect_scale = 1.0 - iter
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r = (RNG::next(current_rng) * n.to_double()).to_int() % n
      // learning phase: attraction to best, with a decaying peer term
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let learned : Array[Double] = []
      for d in 0..<dim {
        let pos = pop[i][d] +
          r1 * (best_pos[d] - pop[i][d]) +
          reflect_scale * r2 * (pop[r][d] - pop[i][d])
        learned.push(max(bounds.min, min(bounds.max, pos)))
      }
      let lv = evaluate(f, learned)
      evaluations = evaluations + 1
      if lv < values[i] {
        pop[i] = learned
        values[i] = lv
        if lv < best_val {
          best_val = lv
          best_pos = learned
        }
      }
      // reflection phase: shrinking Gaussian local search around the best
      if evaluations >= max_evals {
        break
      }
      let sigma = (bounds.max - bounds.min) * 0.1 * reflect_scale
      let reflected : Array[Double] = []
      for d in 0..<dim {
        let gauss = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        reflected.push(
          max(bounds.min, min(bounds.max, best_pos[d] + sigma * gauss)),
        )
      }
      let rv = evaluate(f, reflected)
      evaluations = evaluations + 1
      if rv < values[i] {
        pop[i] = reflected
        values[i] = rv
        if rv < best_val {
          best_val = rv
          best_pos = reflected
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WaOA: Walrus Optimizer (2023) =====

///|
pub struct WaOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_waoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WaOAConfig {
  WaOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Walrus Optimizer (2023): each walrus picks one of three behaviours per
/// generation - migrating toward a peer, feeding around the incumbent best, or
/// escaping a predator with a shrinking step.  The exploration amplitudes decay
/// with the budget fraction so larger budgets settle into local refinement.
/// Greedy acceptance and historical-best tracking are added.
pub fn walrus_optimizer(
  config : WaOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let safety = 1.0 - iter
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r3 = RNG::next(current_rng)
        let i_factor = 1.0 + r3
        let pos = if r1 < 0.33 {
          // migration: leap toward a peer
          pop[i][d] + r2 * safety * (pop[r][d] - i_factor * pop[i][d])
        } else if r1 < 0.66 {
          // feeding: approach the incumbent best
          pop[i][d] +
          r2 * (best_pos[d] - pop[i][d]) +
          r3 * safety * (pop[r][d] - pop[i][d])
        } else {
          // escaping a predator: withdraw with a shrinking step
          pop[i][d] - r2 * safety * (pop[r][d] - pop[i][d])
        }
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      let val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      if val < values[i] {
        pop[i] = new_pos
        values[i] = val
        if val < best_val {
          best_val = val
          best_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CLONALG: Clonal Selection Algorithm (2002) =====

///|
pub struct CLONALGConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_clonalg_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CLONALGConfig {
  CLONALGConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Clonal Selection Algorithm (2002), continuous variant: the best half of the
/// antibody population is selected, each selected antibody is cloned in
/// proportion to its affinity, and clones are mutated with a Gaussian step
/// whose amplitude grows with poorer rank and decays with the budget fraction.
/// A tiny fading jitter prevents late-phase freezing, and receptor editing
/// injects a random antibody each generation.  Greedy acceptance and
/// historical-best tracking are added.
pub fn clonal_selection_algorithm(
  config : CLONALGConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let n_select = if n / 2 < 1 { 1 } else { n / 2 }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let order = sort_indices(values)
    for k in 0..<n_select {
      if evaluations >= max_evals {
        break
      }
      let idx = order[k]
      let rank_frac = k.to_double() / n_select.to_double()
      let sigma = range * (0.001 + 0.3 * (1.0 - iter)) * (0.1 + rank_frac)
      let n_clones = 1 + (n_select - k) / 4
      for _c in 0..<n_clones {
        if evaluations >= max_evals {
          break
        }
        let clone : Array[Double] = []
        for d in 0..<dim {
          let gauss = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          let pos = pop[idx][d] + sigma * gauss
          clone.push(max(bounds.min, min(bounds.max, pos)))
        }
        let cv = evaluate(f, clone)
        evaluations = evaluations + 1
        if cv < values[idx] {
          pop[idx] = clone
          values[idx] = cv
          if cv < best_val {
            best_val = cv
            best_pos = clone
          }
        }
      }
    }
    // receptor editing: inject a random antibody over the worst member
    if evaluations < max_evals {
      let worst_idx = order[order.length() - 1]
      let fresh = random_solution(dim, bounds, current_rng)
      let fv = evaluate(f, fresh)
      evaluations = evaluations + 1
      if fv < values[worst_idx] {
        pop[worst_idx] = fresh
        values[worst_idx] = fv
        if fv < best_val {
          best_val = fv
          best_pos = fresh
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ABC: Artificial Bee Colony (2007) =====

///|
pub struct ABCConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_abc_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ABCConfig {
  ABCConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Artificial Bee Colony (2007), continuous variant: employed and onlooker bees
/// probe neighbours v_j = x_j + phi*(x_j - x_k,j) and greedily accept them,
/// while a scout re-seeds any solution whose improvement trial count exceeds
/// the limit.  Onlooker selection is rank-weighted.  Historical-best tracking
/// is added.
pub fn artificial_bee_colony(
  config : ABCConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let n_employed = if n / 2 < 1 { 1 } else { n / 2 }
  let limit = n_employed

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  let trials : Array[Int32] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
    trials.push(0)
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    // employed bees: each probes a neighbour of its own solution
    for i in 0..<n_employed {
      if evaluations >= max_evals {
        break
      }
      let mut k = (RNG::next(current_rng) * n.to_double()).to_int() % n
      if k == i {
        k = (k + 1) % n
      }
      let j = (RNG::next(current_rng) * dim.to_double()).to_int() % dim
      let phi = 2.0 * RNG::next(current_rng) - 1.0
      let pos = pop[i][j] + phi * (pop[i][j] - pop[k][j])
      let candidate : Array[Double] = []
      for d in 0..<dim {
        candidate.push(
          if d == j {
            max(bounds.min, min(bounds.max, pos))
          } else {
            pop[i][d]
          },
        )
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        trials[i] = 0
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      } else {
        trials[i] = trials[i] + 1
      }
    }
    // onlooker bees: rank-weighted roulette selection, then neighbour probe
    let order = sort_indices(values)
    let min_v = values[order[0]]
    let weights : Array[Double] = []
    let mut total = 0.0
    for i in 0..<n {
      let w = 1.0 / (values[i] - min_v + 1.0)
      weights.push(w)
      total = total + w
    }
    let n_onlookers = n - n_employed
    for _o in 0..<n_onlookers {
      if evaluations >= max_evals {
        break
      }
      let mut pick = RNG::next(current_rng) * total
      let mut sel = 0
      for i in 0..<n {
        pick = pick - weights[i]
        if pick <= 0.0 {
          sel = i
          break
        }
      }
      let mut k = (RNG::next(current_rng) * n.to_double()).to_int() % n
      if k == sel {
        k = (k + 1) % n
      }
      let j = (RNG::next(current_rng) * dim.to_double()).to_int() % dim
      let phi = 2.0 * RNG::next(current_rng) - 1.0
      let pos = pop[sel][j] + phi * (pop[sel][j] - pop[k][j])
      let candidate : Array[Double] = []
      for d in 0..<dim {
        candidate.push(
          if d == j {
            max(bounds.min, min(bounds.max, pos))
          } else {
            pop[sel][d]
          },
        )
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[sel] {
        pop[sel] = candidate
        values[sel] = cv
        trials[sel] = 0
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      } else {
        trials[sel] = trials[sel] + 1
      }
    }
    // scout bees: re-seed exhausted solutions
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      if trials[i] > limit {
        let fresh = random_solution(dim, bounds, current_rng)
        let fv = evaluate(f, fresh)
        evaluations = evaluations + 1
        pop[i] = fresh
        values[i] = fv
        trials[i] = 0
        if fv < best_val {
          best_val = fv
          best_pos = fresh
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EM: Electromagnetism-like Mechanism (2003) =====

///|
pub struct EMConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_em_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EMConfig {
  EMConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Electromagnetism-like Mechanism (2003): each particle carries a charge that
/// grows with its relative fitness, so it is attracted toward better particles
/// and repelled from worse ones.  The step along the resultant force is
/// randomised and shrinks with the budget fraction.  Greedy acceptance and
/// historical-best tracking are added.
pub fn electromagnetism_like_mechanism(
  config : EMConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let order = sort_indices(values)
    let best_v = values[order[0]]
    let worst_v = values[order[order.length() - 1]]
    let denom = worst_v - best_v
    let charges : Array[Double] = []
    for i in 0..<n {
      let q = if denom < 1.0e-12 { 1.0 } else { (worst_v - values[i]) / denom }
      charges.push(q + 1.0e-6)
    }
    let radius = range * 0.2 * (0.1 + 0.9 * (1.0 - iter))
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let force : Array[Double] = []
      for _d in 0..<dim {
        force.push(0.0)
      }
      for j in 0..<n {
        if j != i {
          let mut dist2 = 0.0
          for d in 0..<dim {
            let diff = pop[j][d] - pop[i][d]
            dist2 = dist2 + diff * diff
          }
          let dist = dist2.sqrt()
          if dist > 1.0e-12 {
            let sign = if values[j] < values[i] { 1.0 } else { -1.0 }
            let coeff = sign * charges[j] / dist
            for d in 0..<dim {
              force[d] = force[d] + coeff * (pop[j][d] - pop[i][d])
            }
          }
        }
      }
      let mut fnorm2 = 0.0
      for d in 0..<dim {
        fnorm2 = fnorm2 + force[d] * force[d]
      }
      let fnorm = fnorm2.sqrt()
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if fnorm > 1.0e-12 {
          pop[i][d] + RNG::next(current_rng) * radius * force[d] / fnorm
        } else {
          pop[i][d] + (2.0 * RNG::next(current_rng) - 1.0) * radius * 0.1
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== DBO: Dung Beetle Optimizer (2022) =====

///|
pub struct DBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_dbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DBOConfig {
  DBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Dung Beetle Optimizer (2022): the swarm is split into ball-rolling, brooding,
/// foraging and stealing beetles.  Rollers follow a sun-guided direction with
/// occasional obstacle dances, brooders sample inside a shrinking boundary
/// around the best, foragers refine toward the best, and thieves rush the best.
/// All amplitudes decay with the budget fraction.  Greedy acceptance and
/// historical-best tracking are added.
pub fn dung_beetle_optimizer(
  config : DBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let roll_end = n * 2 / 10
  let brood_end = n * 4 / 10
  let forage_end = n * 7 / 10

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let prev : Array[Array[Double]] = []
  for i in 0..<n {
    prev.push(pop[i])
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let old = pop[i]
      let candidate : Array[Double] = []
      if i < roll_end {
        let dir = if RNG::next(current_rng) < 0.5 { 1.0 } else { -1.0 }
        let dance = RNG::next(current_rng) < 0.1
        for d in 0..<dim {
          let pos = if dance {
            old[d] + (2.0 * RNG::next(current_rng) - 1.0) * range * 0.1 * rb
          } else {
            old[d] +
            dir * 0.1 * rb * (old[d] - prev[i][d]) +
            0.3 * (2.0 * RNG::next(current_rng) - 1.0) * range * 0.05 * rb
          }
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else if i < brood_end {
        let half = range * 0.5 * rb
        for d in 0..<dim {
          let pos = x_best[d] + (2.0 * RNG::next(current_rng) - 1.0) * half
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else if i < forage_end {
        for d in 0..<dim {
          let pos = old[d] +
            RNG::next(current_rng) * (x_best[d] - old[d]) +
            (2.0 * RNG::next(current_rng) - 1.0) * range * 0.2 * rb
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        for d in 0..<dim {
          let pos = x_best[d] +
            0.5 * (2.0 * RNG::next(current_rng) - 1.0) * range * 0.1 * rb
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
      prev[i] = old
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== TSO: Tunicate Swarm Algorithm (2020) =====

///|
pub struct TSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_tso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TSOConfig {
  TSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Tunicate Swarm Algorithm (2020): each tunicate propels toward the swarm best
/// with a jet coefficient built from three random control values, blended with
/// a water-current drift toward a random peer.  The step decays with the budget
/// fraction so larger budgets refine.  Greedy acceptance and historical-best
/// tracking are added.
pub fn tunicate_swarm_algorithm(
  config : TSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pmin = 1.0
  let pmax = 4.0

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let c1 = pmin + RNG::next(current_rng) * (pmax - pmin)
      let c2 = pmin + RNG::next(current_rng) * (pmax - pmin)
      let c3 = pmin + RNG::next(current_rng) * (pmax - pmin)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let jet = (c2 + c3 - 2.0 * c1) / c1
      let jet_c = max(-2.0, min(2.0, jet))
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = pop[i][d] +
          rb *
          (jet_c * (x_best[d] - pop[i][d]) + c2 * (peer[d] - pop[i][d]) * 0.5)
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== KOA: Kookaburra Optimization Algorithm (2024) =====

///|
pub struct KOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_koa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> KOAConfig {
  KOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Kookaburra Optimization Algorithm (2024): each bird either scans from its
/// perch (exploration), blending attraction to the best with a decaying peer
/// difference, or dives (exploitation) inside a shrinking neighbourhood of the
/// best prey.  Greedy acceptance and historical-best tracking are added.
pub fn kookaburra_optimization_algorithm(
  config : KOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let scanning = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if scanning {
          pop[i][d] +
          r1 * (x_best[d] - pop[i][d]) +
          r2 * rb * (peer[d] - pop[i][d])
        } else {
          x_best[d] + (2.0 * RNG::next(current_rng) - 1.0) * range * 0.3 * rb
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FHO: Fire Hawk Optimizer (2022) =====

///|
pub struct FHOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fho_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FHOConfig {
  FHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fire Hawk Optimizer (2022): each fire hawk either drops a burning stick to
/// flush prey (exploration) by chasing a peer while biased to the best, or
/// catches prey (exploitation) with a decaying random walk around the best
/// built from two prey differences.  Greedy acceptance and historical-best
/// tracking are added.
pub fn fire_hawk_optimizer(
  config : FHOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let peer1 = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let flushing = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if flushing {
          pop[i][d] +
          r1 * (x_best[d] - pop[i][d]) +
          r2 * rb * (peer1[d] - pop[i][d])
        } else {
          x_best[d] +
          rb *
          range *
          0.1 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RIME: Rime Optimization Algorithm (2023) =====

///|
pub struct RIMEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rime_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RIMEConfig {
  RIMEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Rime Optimization Algorithm (2023): a soft-rime search grows crystals around
/// the incumbent best with a decaying amplitude, while a hard-rime piercing
/// step fuses agents directly with the best, gated by a normalised rime factor.
/// Greedy acceptance and historical-best tracking are added.
pub fn rime_optimization_algorithm(
  config : RIMEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let beta = 5.0

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let p_soft = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let dist : Array[Double] = []
    let mut total = 0.0
    for i in 0..<n {
      let mut s = 0.0
      for d in 0..<dim {
        s = s + abs_val(pop[i][d] - x_best[d])
      }
      dist.push(s)
      total = total + s
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r2 = RNG::next(current_rng)
      let fr = if total < 1.0e-12 {
        1.0 / n.to_double()
      } else {
        dist[i] / total
      }
      let sigma = range * 0.1 * (1.0 - iter) * (0.5 + fr)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if r2 < p_soft {
          // soft-rime: per-dimension crystal growth around the best
          x_best[d] +
          (2.0 * RNG::next(current_rng) - 1.0) *
          cos_approx(PI * RNG::next(current_rng) / 10.0) *
          beta *
          range *
          0.05 *
          (1.0 - iter)
        } else {
          // hard-rime: shrinking local search around the best
          x_best[d] +
          sigma * box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== INFO: weIghted meaN oF vectOrs (2022) =====

///|
pub struct INFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_info_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> INFOConfig {
  INFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// weIghted meaN oF vectOrs (INFO, 2022): agents move toward a weighted mean of
/// the incumbent best, a randomly chosen better vector and the worst vector,
/// with a random scale that shrinks with the budget fraction.  Greedy
/// acceptance and historical-best tracking are added.
pub fn info_optimization_algorithm(
  config : INFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let n_better = if n / 2 < 1 { 1 } else { n / 2 }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let x_worst = pop[order[order.length() - 1]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let better_idx = order[(RNG::next(current_rng) * n_better.to_double()).to_int() %
        n_better]
      let x_better = pop[better_idx]
      let scale = 2.0 * RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let mean = 0.5 * x_best[d] + 0.3 * x_better[d] + 0.2 * x_worst[d]
        let pos = pop[i][d] +
          scale * rb * (mean - pop[i][d]) +
          (2.0 * RNG::next(current_rng) - 1.0) * range * 0.1 * rb
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SO: Snake Optimizer (2022) =====

///|
pub struct SOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_so_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SOConfig {
  SOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Snake Optimizer (2022): the population splits into males and females.  Low
/// food-quantity (worst-ranked) agents scatter, otherwise agents approach the
/// food (best) when the temperature is high or fight using a male-female
/// difference when it is low.  Steps decay with the budget fraction.  Greedy
/// acceptance and historical-best tracking are added.
pub fn snake_optimizer(
  config : SOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let n_male = if n / 2 < 1 { 1 } else { n / 2 }
  let n_female = n - n_male
  let c2 = 0.05
  let c3 = 2.0

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let temp = exp_approx(-iter)
    let order = sort_indices(values)
    let x_food = pop[order[0]]
    let rank : Array[Int32] = []
    for _i in 0..<n {
      rank.push(0)
    }
    for r in 0..<n {
      rank[order[r]] = r
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let q = 1.0 - rank[i].to_double() / n.to_double()
      let is_male = i < n_male
      let same_rand = if is_male {
        (RNG::next(current_rng) * n_male.to_double()).to_int() % n_male
      } else {
        n_male +
        (RNG::next(current_rng) * n_female.to_double()).to_int() % n_female
      }
      let other_rand = if is_male {
        n_male +
        (RNG::next(current_rng) * n_female.to_double()).to_int() % n_female
      } else {
        (RNG::next(current_rng) * n_male.to_double()).to_int() % n_male
      }
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if q < 0.25 {
          // scarce food: scatter around a random same-sex agent
          pop[same_rand][d] +
          c2 * (2.0 * RNG::next(current_rng) - 1.0) * range * rb
        } else if temp > 0.6 {
          // hot: approach the food
          x_food[d] +
          c3 *
          temp *
          RNG::next(current_rng) *
          rb *
          (x_food[d] - pop[same_rand][d])
        } else {
          // cold: fight using a male-female difference
          x_food[d] +
          c3 *
          RNG::next(current_rng) *
          rb *
          (pop[other_rand][d] - pop[same_rand][d])
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MGO: Mountain Gazelle Optimizer (2022) =====

///|
pub struct MGOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mgo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MGOConfig {
  MGOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Mountain Gazelle Optimizer (2022): gazelles either migrate toward the centre
/// of the adult-male group and the best territory, feed around the best with a
/// peer difference, or fight in a shrinking neighbourhood of the best.  All
/// amplitudes decay with the budget fraction.  Greedy acceptance and
/// historical-best tracking are added.
pub fn mountain_gazelle_optimizer(
  config : MGOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let n_top = if n / 2 < 1 { 1 } else { n / 2 }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let centre : Array[Double] = []
    for d in 0..<dim {
      let mut s = 0.0
      for t in 0..<n_top {
        s = s + pop[order[t]][d]
      }
      centre.push(s / n_top.to_double())
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let behavior = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.34 {
          // migrate toward the adult-male centre and the best territory
          pop[i][d] +
          r1 * rb * (centre[d] - pop[i][d]) +
          r2 * rb * (x_best[d] - pop[i][d])
        } else if behavior < 0.67 {
          // feed: attract toward the best with a decaying peer difference
          pop[i][d] +
          r1 * (x_best[d] - pop[i][d]) +
          r2 * rb * (peer[d] - pop[i][d])
        } else {
          // fight: shrinking Gaussian local search around the best
          x_best[d] +
          rb *
          range *
          0.15 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CHIO: Coronavirus Herd Immunity Optimizer (2021) =====

///|
pub struct CHIOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_chio_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CHIOConfig {
  CHIOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Coronavirus Herd Immunity Optimizer (2021): each agent either moves relative
/// to a "corona" source (a random agent, an immune agent, or the best) while
/// being attracted to the best, or takes an isolated random step.  Steps decay
/// with the budget fraction.  Greedy acceptance and historical-best tracking
/// are added.
pub fn coronavirus_herd_immunity_optimization(
  config : CHIOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let immunity = 0.5
  let n_immune = if (n.to_double() * immunity).to_int() < 1 {
    1
  } else {
    (n.to_double() * immunity).to_int()
  }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r = RNG::next(current_rng)
      let corona = if r < 0.3333 {
        pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      } else if r < 0.6667 {
        pop[order[n -
        1 -
        (RNG::next(current_rng) * n_immune.to_double()).to_int() % n_immune]]
      } else {
        x_best
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let isolate = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if isolate {
          // isolated movement: small random step
          pop[i][d] + rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // move relative to the corona, attracted to the best
          pop[i][d] +
          rb *
          (
            0.5 * r1 * (x_best[d] - pop[i][d]) +
            0.5 * r2 * (pop[i][d] - corona[d])
          )
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== PO: Political Optimizer (2020) =====

///|
pub struct POConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_po_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> POConfig {
  POConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Political Optimizer (2020): the population is split into parties, each with a
/// leader.  During the election campaign every member is attracted to its own
/// party leader, the global best, and another party's leader, with amplitudes
/// that decay with the budget fraction.  Greedy acceptance and historical-best
/// tracking are added.
pub fn political_optimizer(
  config : POConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let n_parties = if n / 4 < 2 { 1 } else { 4 }
  let mpp = if n / n_parties < 1 { 1 } else { n / n_parties }
  let np_minus = if n_parties > 1 { n_parties - 1 } else { 1 }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let leaders : Array[Int32] = []
    for p in 0..<n_parties {
      let mut best_member = -1
      let mut best_member_val = 1.0e308
      let mut i = p * mpp
      while i < n && i < (p + 1) * mpp {
        if values[i] < best_member_val {
          best_member_val = values[i]
          best_member = i
        }
        i = i + 1
      }
      if best_member < 0 {
        best_member = 0
      }
      leaders.push(best_member)
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let p = if i / mpp < n_parties { i / mpp } else { n_parties - 1 }
      let leader = pop[leaders[p]]
      let other_p = (
          p +
          1 +
          (RNG::next(current_rng) * np_minus.to_double()).to_int() % np_minus
        ) %
        n_parties
      let other = pop[leaders[other_p]]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = pop[i][d] +
          rb *
          (
            r1 * (leader[d] - pop[i][d]) +
            r2 * (x_best[d] - pop[i][d]) +
            r3 * (other[d] - pop[i][d])
          )
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SCSO: Sand Cat Swarm Optimization (2020) =====

///|
pub struct SCSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_scso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SCSOConfig {
  SCSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Sand Cat Swarm Optimization (2020): a linearly decaying sensitivity range R
/// makes each cat either search for prey (move toward a random cat) or attack
/// it (attract toward the best with a shrinking jitter).  Greedy acceptance and
/// historical-best tracking are added.
pub fn sand_cat_swarm_optimization(
  config : SCSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let r_g = 2.0 * rb
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r = r_g * RNG::next(current_rng)
      let rr = 2.0 * r_g * RNG::next(current_rng) - r_g
      let k = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let attacking = abs_val(rr) <= 1.0
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if attacking {
          // attack prey: attract to the best with a shrinking jitter
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          (2.0 * RNG::next(current_rng) - 1.0) * range * 0.1 * rb
        } else {
          // search for prey: move toward a random cat
          pop[i][d] + r * (pop[k][d] - pop[i][d])
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SABO: Subtraction-Average-Based Optimizer (2023) =====

///|
pub struct SABOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sabo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SABOConfig {
  SABOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Subtraction-Average-Based Optimizer (2023): each agent moves by a random
/// signed multiple of its average difference to the swarm (A_i = X_i - mean),
/// combined with an attraction to the best and a shrinking Gaussian local
/// search so larger budgets refine reliably.  Greedy acceptance and
/// historical-best tracking are added.
pub fn subtraction_average_based_optimizer(
  config : SABOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let mean : Array[Double] = []
    for d in 0..<dim {
      let mut s = 0.0
      for j in 0..<n {
        s = s + pop[j][d]
      }
      mean.push(s / n.to_double())
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let avg = pop[i][d] - mean[d]
        let r = 2.0 * RNG::next(current_rng) - 1.0
        let pos = pop[i][d] +
          r * avg +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb *
          range *
          0.05 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== HPO: Hunter-Prey Optimizer (2022) =====

///|
pub struct HPOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hpo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HPOConfig {
  HPOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Hunter-Prey Optimizer (2022): each agent acts either as a hunter, drawn by a
/// decaying balance factor toward the swarm mean and the best, or as prey doing
/// a shrinking Gaussian local search around the best.  Greedy acceptance and
/// historical-best tracking are added.
pub fn hunter_prey_optimizer(
  config : HPOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let z = 2.0 * rb
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let mean : Array[Double] = []
    for d in 0..<dim {
      let mut s = 0.0
      for j in 0..<n {
        s = s + pop[j][d]
      }
      mean.push(s / n.to_double())
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let hunter = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if hunter {
          // hunter: drawn toward the swarm mean and the best
          pop[i][d] +
          z *
          (
            RNG::next(current_rng) * (mean[d] - pop[i][d]) +
            RNG::next(current_rng) * (x_best[d] - pop[i][d])
          )
        } else {
          // prey: shrinking local search around the best
          x_best[d] +
          rb *
          range *
          0.15 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== AHA: Artificial Hummingbird Algorithm (2021) =====

///|
pub struct AHAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aha_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AHAConfig {
  AHAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Artificial Hummingbird Algorithm (2021): birds either perform guided
/// foraging toward the best with a random flight direction, territorial
/// foraging around their own position with a shrinking amplitude, or a rare
/// migration repositioning.  Greedy acceptance and historical-best tracking are
/// added.
pub fn artificial_hummingbird_algorithm(
  config : AHAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let candidate : Array[Double] = []
      if behavior < 0.6 {
        // guided foraging: fly toward the best with a random flight direction
        for d in 0..<dim {
          let dir = if RNG::next(current_rng) < 0.5 { 1.0 } else { -1.0 }
          let pos = pop[i][d] +
            dir * RNG::next(current_rng) * (x_best[d] - pop[i][d])
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else if behavior < 0.98 {
        // territorial foraging: shrinking local search around itself
        for d in 0..<dim {
          let pos = pop[i][d] +
            (2.0 * RNG::next(current_rng) - 1.0) * range * 0.1 * rb
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // migration: rare random repositioning
        for _d in 0..<dim {
          let pos = bounds.min + RNG::next(current_rng) * range
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BKA: Black-winged Kite Algorithm (2023) =====

///|
pub struct BKAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bka_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BKAConfig {
  BKAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Black-winged Kite Algorithm (2023): kites alternate between an attack
/// behaviour that pounces toward the best with a sinusoidal velocity and a
/// migration behaviour that leads toward the best habitat.  The attack
/// probability and amplitudes decay with the budget fraction.  Greedy
/// acceptance and historical-best tracking are added.
pub fn black_winged_kite_algorithm(
  config : BKAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let p_attack = 0.1 + 0.9 * rb
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let attack = RNG::next(current_rng) < p_attack
      let vel = 0.05 + 0.2 * rb
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        let pos = if attack {
          // attack: pounce toward the best with a sinusoidal velocity
          pop[i][d] +
          vel * (1.0 + sin_approx(TWO_PI * r)) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // migration: lead toward the best habitat
          pop[i][d] + r * (x_best[d] - pop[i][d])
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ARO: Artificial Rabbits Optimization (2022) =====

///|
pub struct AROConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aro_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AROConfig {
  AROConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Artificial Rabbits Optimization (2022): an energy factor that decays with the
/// budget fraction selects detour foraging (toward the best and a random rabbit
/// with a Gaussian disturbance), while low energy triggers random hiding (a
/// shrinking Gaussian local search around the best).  Greedy acceptance and
/// historical-best tracking are added.
pub fn artificial_rabbits_optimization(
  config : AROConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let energy = 4.0 * rb
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let j = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if energy > 1.0 {
          // detour foraging: explore toward the best and a random rabbit
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          RNG::next(current_rng) * rb * (pop[j][d] - pop[i][d]) +
          rb *
          range *
          0.05 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        } else {
          // random hiding: shrinking Gaussian local search around the best
          x_best[d] +
          rb *
          range *
          0.1 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FOX: Fox Optimizer (2023) =====

///|
pub struct FOXConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fox_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FOXConfig {
  FOXConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fox Optimizer (2023): a fox estimates the distance to its prey and either
/// explores with a decaying jump toward the best or pounces into a shrinking
/// neighbourhood of it.  Greedy acceptance and historical-best tracking are
/// added.
pub fn fox_optimizer(
  config : FOXConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let a = 2.0 * rb
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let mut dist2 = 0.0
      for d in 0..<dim {
        let diff = pop[i][d] - x_best[d]
        dist2 = dist2 + diff * diff
      }
      let dist = dist2.sqrt()
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        let pos = if r > 0.5 {
          // exploration: decaying jump toward the best
          pop[i][d] + a * RNG::next(current_rng) * (x_best[d] - pop[i][d])
        } else if dist > 1.0 {
          // pounce into the prey's vicinity
          x_best[d] + (2.0 * RNG::next(current_rng) - 1.0) * range * 0.1 * rb
        } else {
          // close prey: fine local refinement
          x_best[d] + (2.0 * RNG::next(current_rng) - 1.0) * range * 0.05 * rb
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== DMO: Dwarf Mongoose Optimization (2022) =====

///|
pub struct DMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_dmo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DMOConfig {
  DMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Dwarf Mongoose Optimization (2022): a pack is led by an alpha; mongooses
/// either forage by following the alpha with a peep-based attraction or scout
/// with a decaying local exploration.  Greedy acceptance and historical-best
/// tracking are added.
pub fn dwarf_mongoose_optimization(
  config : DMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_alpha = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.6 {
          // foraging: peep-guided attraction to the alpha
          pop[i][d] +
          RNG::next(current_rng) * (x_alpha[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // scouting: decaying local exploration
          pop[i][d] + (2.0 * RNG::next(current_rng) - 1.0) * range * 0.2 * rb
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RSO: Rat Swarm Optimizer (2021) =====

///|
pub struct RSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RSOConfig {
  RSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Rat Swarm Optimizer (2021): rats either chase a prey by pouncing toward the
/// best or fight over it with a shrinking Gaussian local search around the
/// best.  Greedy acceptance and historical-best tracking are added.
pub fn rat_swarm_optimizer(
  config : RSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let chasing = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if chasing {
          // chasing: convex pounce toward the best
          pop[i][d] + RNG::next(current_rng) * (x_best[d] - pop[i][d])
        } else {
          // fighting: shrinking Gaussian local search around the best
          x_best[d] +
          rb *
          range *
          0.1 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EVO: Egyptian Vulture Optimization (2021) =====

///|
pub struct EVOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_evo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EVOConfig {
  EVOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Egyptian Vulture Optimization (2021): vultures spend part of the budget
/// circling a shrinking neighbourhood of the best, and the rest attacking by
/// attraction to the best or a fine local search.  Greedy acceptance and
/// historical-best tracking are added.
pub fn egyptian_vulture_optimization(
  config : EVOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if r < 0.5 {
          // exploration: circling the best's neighbourhood
          x_best[d] + (2.0 * RNG::next(current_rng) - 1.0) * rb * range * 0.3
        } else if r < 0.75 {
          // exploitation: attack by attraction to the best
          pop[i][d] + RNG::next(current_rng) * (x_best[d] - pop[i][d])
        } else {
          // exploitation: fine local search around the best
          x_best[d] + rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EEFO: Electric Eel Foraging Optimization (2023) =====

///|
pub struct EEFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_eefo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EEFOConfig {
  EEFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Electric Eel Foraging Optimization (2023): eels either rest with a shrinking
/// local search, migrate toward the best, or hunt by combining attraction to the
/// best with a peer difference.  Greedy acceptance and historical-best tracking
/// are added.
pub fn electric_eel_foraging_optimization(
  config : EEFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let j = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.33 {
          // resting: shrinking local search around itself
          pop[i][d] + (2.0 * RNG::next(current_rng) - 1.0) * rb * range * 0.1
        } else if behavior < 0.66 {
          // migrating: move toward the best
          pop[i][d] + RNG::next(current_rng) * (x_best[d] - pop[i][d])
        } else {
          // hunting: toward the best and a peer
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          RNG::next(current_rng) * rb * (pop[j][d] - pop[i][d])
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Golden-SA: Golden Sine Algorithm (2017) =====

///|
pub struct GoldenSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_golden_sa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GoldenSAConfig {
  GoldenSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Golden Sine Algorithm (2017): positions are updated with a sine-modulated
/// step whose coefficients are drawn from golden-ratio-weighted intervals, plus
/// a shrinking Gaussian local search around the best for reliable convergence.
/// Greedy acceptance and historical-best tracking are added.
pub fn golden_sine_algorithm(
  config : GoldenSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let tau = 0.6180339887498949

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let r1 = TWO_PI * RNG::next(current_rng)
      let r2 = PI * RNG::next(current_rng)
      let r3 = -PI * tau + TWO_PI * tau * RNG::next(current_rng)
      let r4 = -PI * (1.0 - tau) + TWO_PI * (1.0 - tau) * RNG::next(current_rng)
      let coeff = r1 * sin_approx(r2) * rb
      let wander = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if wander {
          // golden-sine search step
          pop[i][d] + coeff * abs_val(r3 * x_best[d] - r4 * pop[i][d])
        } else {
          // shrinking Gaussian local search around the best
          x_best[d] +
          rb *
          range *
          0.1 *
          box_muller(RNG::next(current_rng), RNG::next(current_rng))
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== YSGA: Yellow Saddle Goatfish Algorithm (2019) =====

///|
pub struct YSGAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ysga_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> YSGAConfig {
  YSGAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Yellow Saddle Goatfish Algorithm (2019): the school splits into chasers that
/// pursue the prey (best) and blockers that surround it with a shrinking ring.
/// Greedy acceptance and historical-best tracking are added.
pub fn yellow_saddle_goatfish_algorithm(
  config : YSGAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let n_chasers = if n / 2 < 1 { 1 } else { n / 2 }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_prey = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let blocker = i >= n_chasers
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if blocker {
          // blocker: surround the prey with a shrinking ring
          x_prey[d] + (2.0 * RNG::next(current_rng) - 1.0) * rb * range * 0.3
        } else {
          // chaser: pursue the prey
          pop[i][d] +
          RNG::next(current_rng) * (x_prey[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FBI: Forensic-Based Investigation (2020) =====

///|
pub struct FBIConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fbi_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FBIConfig {
  FBIConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Forensic-Based Investigation (2020): an investigation phase moves agents
/// relative to a random clue while attracted to the best, and an interrogation
/// phase draws suspects toward the prosecutor (best) with a shrinking jitter.
/// Greedy acceptance and historical-best tracking are added.
pub fn forensic_based_investigation(
  config : FBIConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let investigate = RNG::next(current_rng) < 0.5
      let r = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if investigate {
          // investigation: relative to a random clue, attracted to the best
          pop[i][d] +
          RNG::next(current_rng) * rb * (pop[r][d] - pop[i][d]) +
          RNG::next(current_rng) * (x_best[d] - pop[i][d])
        } else {
          // interrogation: suspect drawn toward the prosecutor (best)
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CapSA: Capuchin Search Algorithm (2021) =====

///|
pub struct CapSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_capsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CapSAConfig {
  CapSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Capuchin Search Algorithm (2021): capuchins either jump (attraction to the
/// best with a shrinking jitter), swing (a shrinking search ring around the
/// best), or climb/forage (a local walk around themselves).  Greedy acceptance
/// and historical-best tracking are added.
pub fn capuchin_search_algorithm(
  config : CapSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.34 {
          // jumping: leap toward the best with a shrinking jitter
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        } else if behavior < 0.67 {
          // swinging: shrinking search ring around the best
          x_best[d] + rb * range * 0.2 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // climbing/foraging: local walk around itself
          pop[i][d] + rb * range * 0.2 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SWO: Spider Wasp Optimizer (2023) =====

///|
pub struct SWOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_swo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SWOConfig {
  SWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Spider Wasp Optimizer (2023): wasps search for spiders (toward a random
/// peer), follow/capture the best, or build a nest with a shrinking local
/// search around the best.  Greedy acceptance and historical-best tracking are
/// added.
pub fn spider_wasp_optimizer(
  config : SWOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let r = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.34 {
          // search: move toward a random peer
          pop[i][d] + RNG::next(current_rng) * (pop[r][d] - pop[i][d])
        } else if behavior < 0.67 {
          // follow/capture: toward the best spider
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // nest: shrinking local search around the best
          x_best[d] + rb * range * 0.15 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GAO: Green Anaconda Optimization (2023) =====

///|
pub struct GAOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gao_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GAOConfig {
  GAOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Green Anaconda Optimization (2023): anacondas either attack (approach the
/// prey/best) or coil around it with a shrinking ring.  Greedy acceptance and
/// historical-best tracking are added.
pub fn green_anaconda_optimization(
  config : GAOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_prey = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let attack = RNG::next(current_rng) < 0.5
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if attack {
          // attack: approach the prey
          pop[i][d] +
          RNG::next(current_rng) * (x_prey[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // coil: shrinking ring around the prey
          x_prey[d] + rb * range * 0.25 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MAO: Mexican Axolotl Optimization (2024) =====

///|
pub struct MAOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mao_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MAOConfig {
  MAOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Mexican Axolotl Optimization (2024): axolotls forage by blending attraction
/// to the best with a peer difference, or regenerate with a shrinking local
/// search around the best.  Greedy acceptance and historical-best tracking are
/// added.
pub fn mexican_axolotl_optimization(
  config : MAOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let r = (RNG::next(current_rng) * n.to_double()).to_int() % n
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.5 {
          // forage: attraction to the best with a peer difference
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          RNG::next(current_rng) * rb * (pop[r][d] - pop[i][d])
        } else {
          // regenerate: shrinking local search around the best
          x_best[d] + rb * range * 0.15 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== UMDA: Univariate Marginal Distribution Algorithm (1996) =====

///|
pub struct UMDAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_umda_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> UMDAConfig {
  UMDAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Univariate Marginal Distribution Algorithm (1996), continuous variant: a
/// Gaussian model (per-dimension mean and standard deviation) is sampled each
/// generation, then re-estimated from the best half of the sample.  A variance
/// floor prevents collapse and historical-best tracking is added.
pub fn univariate_marginal_distribution_algorithm(
  config : UMDAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let floor = range * 1.0e-4

  let mean : Array[Double] = []
  let sigma : Array[Double] = []
  for _d in 0..<dim {
    mean.push((bounds.min + bounds.max) * 0.5)
    sigma.push(range * 0.3)
  }
  let mut best_pos : Array[Double] = []
  for _d in 0..<dim {
    best_pos.push(bounds.min)
  }
  let mut best_val = 1.0e308
  let mut evaluations = 0

  while evaluations < max_evals {
    let pop : Array[Array[Double]] = []
    let values : Array[Double] = []
    for _i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let sol : Array[Double] = []
      for d in 0..<dim {
        let sample = mean[d] +
          sigma[d] * box_muller(RNG::next(current_rng), RNG::next(current_rng))
        sol.push(max(bounds.min, min(bounds.max, sample)))
      }
      let v = evaluate(f, sol)
      evaluations = evaluations + 1
      pop.push(sol)
      values.push(v)
      if v < best_val {
        best_val = v
        best_pos = sol
      }
    }
    let m = pop.length()
    if m == 0 {
      break
    }
    let order = sort_indices(values)
    let n_keep = if m / 2 < 1 { 1 } else { m / 2 }
    for d in 0..<dim {
      let mut s = 0.0
      for k in 0..<n_keep {
        s = s + pop[order[k]][d]
      }
      let mu = s / n_keep.to_double()
      let mut v2 = 0.0
      for k in 0..<n_keep {
        let diff = pop[order[k]][d] - mu
        v2 = v2 + diff * diff
      }
      let sd = (v2 / n_keep.to_double()).sqrt()
      mean[d] = mu
      sigma[d] = if sd < floor { floor } else { sd }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== cGA: Compact Genetic Algorithm (1999) =====

///|
pub struct CGAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cga_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CGAConfig {
  CGAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Compact Genetic Algorithm (1999), continuous analogue: a Gaussian model is
/// updated by a tournament of two samples, moving the mean toward the winner
/// with rate 1/N while the spread shrinks, so the model concentrates near the
/// optimum.  Historical-best tracking is added.
pub fn compact_genetic_algorithm(
  config : CGAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let floor = range * 1.0e-5

  let mean : Array[Double] = []
  let sigma : Array[Double] = []
  for _d in 0..<dim {
    mean.push((bounds.min + bounds.max) * 0.5)
    sigma.push(range * 0.3)
  }
  let mut best_pos : Array[Double] = []
  for _d in 0..<dim {
    best_pos.push(bounds.min)
  }
  let mut best_val = 1.0e308
  let mut evaluations = 0

  while evaluations < max_evals {
    let pop : Array[Array[Double]] = []
    let vals : Array[Double] = []
    for _i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let sol : Array[Double] = []
      for d in 0..<dim {
        let sample = mean[d] +
          sigma[d] * box_muller(RNG::next(current_rng), RNG::next(current_rng))
        sol.push(max(bounds.min, min(bounds.max, sample)))
      }
      let v = evaluate(f, sol)
      evaluations = evaluations + 1
      pop.push(sol)
      vals.push(v)
      if v < best_val {
        best_val = v
        best_pos = sol
      }
    }
    let m = pop.length()
    if m < 2 {
      break
    }
    let n_pairs = m / 2
    let winners : Array[Array[Double]] = []
    for k in 0..<n_pairs {
      let i = k * 2
      let j = i + 1
      let winner = if vals[i] <= vals[j] { pop[i] } else { pop[j] }
      winners.push(winner)
    }
    for d in 0..<dim {
      let mut wm = 0.0
      for k in 0..<n_pairs {
        wm = wm + winners[k][d]
      }
      wm = wm / n_pairs.to_double()
      let mut wv = 0.0
      for k in 0..<n_pairs {
        let diff = winners[k][d] - wm
        wv = wv + diff * diff
      }
      let sd = (wv / n_pairs.to_double()).sqrt()
      mean[d] = wm
      sigma[d] = if sd < floor { floor } else { sd }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CA: Cultural Algorithm (1994) =====

///|
pub struct CAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ca_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CAConfig {
  CAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Cultural Algorithm (1994), continuous variant: a population space evolves
/// under the guidance of a belief space whose situational knowledge is the best
/// individual and whose normative knowledge is a shrinking interval around it.
/// Individuals either sample that interval or are attracted to the best.
/// Greedy acceptance and historical-best tracking are added.
pub fn cultural_algorithm(
  config : CAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let half = range * 0.3 * rb
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.5 {
          // normative knowledge: sample the shrinking belief interval
          let lo = max(bounds.min, x_best[d] - half)
          let hi = min(bounds.max, x_best[d] + half)
          lo + RNG::next(current_rng) * (hi - lo)
        } else {
          // situational knowledge: move toward the best individual
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== COY: Coyote Optimization Algorithm (2018) =====

///|
pub struct COYConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_coy_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> COYConfig {
  COYConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Coyote Optimization Algorithm (2018): the population lives in packs; each
/// coyote grows toward its pack's alpha (best) and social tendency (pack mean),
/// and a pup born from two parents may replace the pack's worst member.  Greedy
/// acceptance and historical-best tracking are added.
pub fn coyote_optimization_algorithm(
  config : COYConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let n_packs = if n / 5 < 1 { 1 } else { n / 5 }
  let pack_size = if n / n_packs < 1 { 1 } else { n / n_packs }

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    for p in 0..<n_packs {
      let start = p * pack_size
      if start >= n {
        break
      }
      let end0 = start + pack_size
      let end = if end0 < n { end0 } else { n }
      let mut alpha = start
      let mut a_val = values[start]
      for c in (start + 1)..<end {
        if values[c] < a_val {
          a_val = values[c]
          alpha = c
        }
      }
      let tendency : Array[Double] = []
      for d in 0..<dim {
        let mut s = 0.0
        for c in start..<end {
          s = s + pop[c][d]
        }
        tendency.push(s / (end - start).to_double())
      }
      for c in start..<end {
        if evaluations >= max_evals {
          break
        }
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let candidate : Array[Double] = []
        for d in 0..<dim {
          let pos = pop[c][d] +
            r1 * (pop[alpha][d] - pop[c][d]) +
            r2 * rb * (tendency[d] - pop[c][d]) +
            rb * range * 0.03 * (2.0 * RNG::next(current_rng) - 1.0)
          candidate.push(max(bounds.min, min(bounds.max, pos)))
        }
        let cv = evaluate(f, candidate)
        evaluations = evaluations + 1
        if cv < values[c] {
          pop[c] = candidate
          values[c] = cv
          if cv < best_val {
            best_val = cv
            best_pos = candidate
          }
        }
      }
      if evaluations < max_evals {
        let mate = start +
          (RNG::next(current_rng) * (end - start).to_double()).to_int() %
          (end - start)
        let pup : Array[Double] = []
        for d in 0..<dim {
          let pos = pop[alpha][d] +
            RNG::next(current_rng) * (pop[mate][d] - pop[alpha][d]) +
            rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
          pup.push(max(bounds.min, min(bounds.max, pos)))
        }
        let pv = evaluate(f, pup)
        evaluations = evaluations + 1
        let mut worst = start
        for c in start..<end {
          if values[c] > values[worst] {
            worst = c
          }
        }
        if pv < values[worst] {
          pop[worst] = pup
          values[worst] = pv
          if pv < best_val {
            best_val = pv
            best_pos = pup
          }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BES: Bald Eagle Search (2019) =====

///|
pub struct BESConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bes_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BESConfig {
  BESConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Bald Eagle Search (2019): eagles select a space near the best using the
/// swarm mean, search with a shrinking local step, or swoop by attraction to
/// the best.  Greedy acceptance and historical-best tracking are added.
pub fn bald_eagle_search(
  config : BESConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    let mean : Array[Double] = []
    for d in 0..<dim {
      let mut s = 0.0
      for j in 0..<n {
        s = s + pop[j][d]
      }
      mean.push(s / n.to_double())
    }
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let phase = RNG::next(current_rng)
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if phase < 0.4 {
          // select: choose a space near the best using the swarm mean
          x_best[d] + 1.5 * RNG::next(current_rng) * (mean[d] - pop[i][d])
        } else if phase < 0.7 {
          // search: shrinking local step around itself
          pop[i][d] + rb * range * 0.25 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // swoop: attraction toward the best
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GOOSE: Goose Optimization Algorithm (2024) =====

///|
pub struct GOOSEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_goose_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GOOSEConfig {
  GOOSEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Goose Optimization Algorithm (2024): geese forage toward the best with a
/// peer difference, scatter with a shrinking step when alarmed, or migrate in
/// formation around the leader (best).  Greedy acceptance and historical-best
/// tracking are added.
pub fn goose_optimization_algorithm(
  config : GOOSEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.5 {
          // forage: toward the best with a peer difference
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          RNG::next(current_rng) * rb * (peer[d] - pop[i][d])
        } else if behavior < 0.8 {
          // alarm: scatter with a shrinking step
          pop[i][d] + rb * range * 0.3 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // migration: formation around the leader (best)
          x_best[d] + rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Powell: Powell's Conjugate Direction Method (1964) =====

///|
pub struct PowellConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_powell_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PowellConfig {
  PowellConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Powell's conjugate direction method (1964), the classic derivative-free
/// direct search: line searches along a set of directions, then the net move of
/// the sweep becomes a new conjugate direction.  This implementation is fully
/// deterministic (the step schedule does not depend on the evaluation budget),
/// starts from the best of a random sample, and keeps historical-best tracking.
pub fn powell_conjugate_direction_method(
  config : PowellConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let base_step = range * 0.05
  let min_step = range * 1.0e-9

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  let dirs : Array[Array[Double]] = []
  for i in 0..<dim {
    let u : Array[Double] = []
    for _d in 0..<dim {
      u.push(0.0)
    }
    u[i] = 1.0
    dirs.push(u)
  }
  let mut step = base_step

  while evaluations < max_evals {
    let start_pt : Array[Double] = []
    for d in 0..<dim {
      start_pt.push(best_pos[d])
    }
    let start_val = best_val
    let deltas : Array[Double] = []
    for i in 0..<dim {
      if evaluations >= max_evals {
        break
      }
      let mut guard_cnt = 0
      while evaluations < max_evals && guard_cnt < 32 {
        guard_cnt = guard_cnt + 1
        let cand : Array[Double] = []
        for d in 0..<dim {
          let pos = best_pos[d] + step * dirs[i][d]
          cand.push(max(bounds.min, min(bounds.max, pos)))
        }
        let cv = evaluate(f, cand)
        evaluations = evaluations + 1
        if cv < best_val {
          best_val = cv
          best_pos = cand
          step = step * 1.8
        } else {
          let cand2 : Array[Double] = []
          for d in 0..<dim {
            let pos = best_pos[d] - step * dirs[i][d]
            cand2.push(max(bounds.min, min(bounds.max, pos)))
          }
          let cv2 = evaluate(f, cand2)
          evaluations = evaluations + 1
          if cv2 < best_val {
            best_val = cv2
            best_pos = cand2
            step = step * 1.8
          } else {
            step = step * 0.5
          }
        }
        if step < min_step {
          break
        }
      }
      deltas.push(start_val - best_val)
    }
    if evaluations >= max_evals {
      break
    }
    // conjugate direction: the net move of the whole sweep
    let dirv : Array[Double] = []
    let mut norm = 0.0
    for d in 0..<dim {
      let diff = best_pos[d] - start_pt[d]
      dirv.push(diff)
      norm = norm + diff * diff
    }
    norm = norm.sqrt()
    if norm > 1.0e-12 {
      let mut guard_cnt = 0
      while evaluations < max_evals && guard_cnt < 32 {
        guard_cnt = guard_cnt + 1
        let cand : Array[Double] = []
        for d in 0..<dim {
          let pos = best_pos[d] + step * dirv[d] / norm
          cand.push(max(bounds.min, min(bounds.max, pos)))
        }
        let cv = evaluate(f, cand)
        evaluations = evaluations + 1
        if cv < best_val {
          best_val = cv
          best_pos = cand
          step = step * 1.8
        } else {
          let cand2 : Array[Double] = []
          for d in 0..<dim {
            let pos = best_pos[d] - step * dirv[d] / norm
            cand2.push(max(bounds.min, min(bounds.max, pos)))
          }
          let cv2 = evaluate(f, cand2)
          evaluations = evaluations + 1
          if cv2 < best_val {
            best_val = cv2
            best_pos = cand2
            step = step * 1.8
          } else {
            step = step * 0.5
          }
        }
        if step < min_step {
          break
        }
      }
      let delta_dir = start_val - best_val
      // replace the least productive direction with the conjugate one
      let mut min_i = 0
      for i in 1..<deltas.length() {
        if deltas[i] < deltas[min_i] {
          min_i = i
        }
      }
      if delta_dir > deltas[min_i] {
        let u : Array[Double] = []
        for d in 0..<dim {
          u.push(dirv[d] / norm)
        }
        dirs[min_i] = u
      }
    }
    // stagnation: the whole sweep made no progress - restart from a fresh
    // random point (historical best is preserved) and reset the step
    let round_gain = start_val - best_val
    if round_gain <= 0.0 {
      let fresh = random_solution(dim, bounds, current_rng)
      let fv = evaluate(f, fresh)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = fresh
      }
      step = base_step
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== TEO: Thermal Exchange Optimization (2017) =====

///|
pub struct TEOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_teo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TEOConfig {
  TEOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Thermal Exchange Optimization (2017): objects are randomly paired and heat
/// flows from the hotter (worse) to the cooler (better) one, while a fraction
/// of the population cools towards the environment (the incumbent best).  The
/// exchange amplitude decays with the budget fraction.  Greedy acceptance and
/// historical-best tracking are added.
pub fn thermal_exchange_optimization(
  config : TEOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    // heat exchange: random pairs, heat flows from the hotter to the cooler
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      if i % 2 == 0 && i + 1 < n {
        let j = i + 1
        let better = if values[i] <= values[j] { i } else { j }
        let worse = if better == i { j } else { i }
        let cand : Array[Double] = []
        for d in 0..<dim {
          let pos = pop[worse][d] +
            RNG::next(current_rng) * rb * (pop[better][d] - pop[worse][d]) +
            rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
          cand.push(max(bounds.min, min(bounds.max, pos)))
        }
        let cv = evaluate(f, cand)
        evaluations = evaluations + 1
        if cv < values[worse] {
          pop[worse] = cand
          values[worse] = cv
          if cv < best_val {
            best_val = cv
            best_pos = cand
          }
        }
        // the cooler object warms slightly towards thermal balance
        if evaluations < max_evals {
          let cand2 : Array[Double] = []
          for d in 0..<dim {
            let pos = pop[better][d] +
              RNG::next(current_rng) *
              rb *
              (pop[worse][d] - pop[better][d]) *
              0.3 +
              rb * range * 0.03 * (2.0 * RNG::next(current_rng) - 1.0)
            cand2.push(max(bounds.min, min(bounds.max, pos)))
          }
          let cv2 = evaluate(f, cand2)
          evaluations = evaluations + 1
          if cv2 < values[better] {
            pop[better] = cand2
            values[better] = cv2
            if cv2 < best_val {
              best_val = cv2
              best_pos = cand2
            }
          }
        }
      }
    }
    // environmental cooling: a fraction of the objects converge on the best
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      if RNG::next(current_rng) < 0.3 {
        let cand : Array[Double] = []
        for d in 0..<dim {
          let pos = x_best[d] +
            rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
          cand.push(max(bounds.min, min(bounds.max, pos)))
        }
        let cv = evaluate(f, cand)
        evaluations = evaluations + 1
        if cv < values[i] {
          pop[i] = cand
          values[i] = cv
          if cv < best_val {
            best_val = cv
            best_pos = cand
          }
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SHO: Spotted Hyena Optimizer (2017) =====

///|
pub struct SHOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sho_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SHOConfig {
  SHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Spotted Hyena Optimizer (2017): the pack encircles the prey (best) with a
/// linearly decaying attack coefficient and alternates between wide encircling
/// leaps and an average attack position.  Greedy acceptance and historical-best
/// tracking are added.
pub fn spotted_hyena_optimizer(
  config : SHOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_prey = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let b_coef = 2.0 * RNG::next(current_rng)
      // quadratic decay: keeps a long refinement phase on large budgets
      let h_coef = 5.0 * rb * rb
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if RNG::next(current_rng) < 0.5 {
          // encircling: a wide leap scaled by the decaying coefficients
          x_prey[d] - b_coef * h_coef * (x_prey[d] - pop[i][d])
        } else {
          // attacking: an average position between the prey and the hyena
          x_prey[d] + RNG::next(current_rng) * rb * (pop[i][d] - x_prey[d])
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ROA: Remora Optimization Algorithm (2021) =====

///|
pub struct ROAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_roa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ROAConfig {
  ROAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Remora Optimization Algorithm (2021): remoras either attach to the host
/// (best) and feed in its shrinking neighbourhood, swim with the host by
/// blending attraction to it with a peer difference, or forage alone with a
/// shrinking local walk.  Greedy acceptance and historical-best tracking are
/// added.
pub fn remora_optimization_algorithm(
  config : ROAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.4 {
          // attach: feed in the host's shrinking neighbourhood
          x_best[d] + rb * range * 0.15 * (2.0 * RNG::next(current_rng) - 1.0)
        } else if behavior < 0.7 {
          // swim with the host: attraction to the best plus a peer difference
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          RNG::next(current_rng) * rb * (peer[d] - pop[i][d])
        } else {
          // forage alone: shrinking local walk
          pop[i][d] + rb * range * 0.2 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== HOA: Horse Herd Optimization (2019) =====

///|
pub struct HOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HOAConfig {
  HOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Horse Herd Optimization (2019): horses either graze towards the dominant
/// stallion (best), breed by blending with a random companion, or fight for
/// dominance in a shrinking neighbourhood of the best.  Greedy acceptance and
/// historical-best tracking are added.
pub fn horse_herd_optimization(
  config : HOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.4 {
          // grazing: move towards the dominant stallion
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else if behavior < 0.7 {
          // breeding: blend with a random companion
          pop[i][d] +
          RNG::next(current_rng) * (peer[d] - pop[i][d]) +
          rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // fighting for dominance: shrinking neighbourhood of the best
          x_best[d] + rb * range * 0.15 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BWA: Black Widow Optimization Algorithm (2020) =====

///|
pub struct BWAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bwa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BWAConfig {
  BWAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Black Widow Optimization Algorithm (2020): males either approach the female
/// (best) across her web, move within the shrinking web around her, or mate
/// with a random companion.  Greedy acceptance and historical-best tracking
/// are added.
pub fn black_widow_optimization_algorithm(
  config : BWAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_female = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.4 {
          // approach the female: attraction to the best
          pop[i][d] +
          RNG::next(current_rng) * (x_female[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else if behavior < 0.7 {
          // web: shrink onto the female's web
          x_female[d] + rb * range * 0.15 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // mating: blend with a random companion
          pop[i][d] +
          RNG::next(current_rng) * (peer[d] - pop[i][d]) +
          rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CPSO: Constriction Particle Swarm Optimization (2002) =====

///|
pub struct CPSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  swarm_size : Int32
}

///|
pub fn make_cpso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  swarm_sz : Int32,
) -> CPSOConfig {
  CPSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    swarm_size: swarm_sz,
  }
}

///|
/// Constriction Particle Swarm Optimization (2002, Clerc & Kennedy): the
/// classic PSO update is multiplied by the constriction factor chi=0.7298 with
/// c1=c2=1.4962 (no inertia weight), and velocities are clamped to a fraction
/// of the search range.  Historical-best tracking is added.
pub fn constriction_pso(
  config : CPSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.swarm_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let chi = 0.729843788
  let c1 = 1.496180
  let c2 = 1.496180
  let v_max = 0.2 * range

  let positions : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let pbest : Array[Array[Double]] = []
  let pvals : Array[Double] = []
  for _i in 0..<n {
    let pos = random_solution(dim, bounds, current_rng)
    positions.push(pos)
    let vel : Array[Double] = []
    for _d in 0..<dim {
      vel.push((RNG::next(current_rng) * 2.0 - 1.0) * v_max)
    }
    velocities.push(vel)
    pbest.push(pos)
    pvals.push(evaluate(f, pos))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(pvals)
  let mut gbest : Array[Double] = pbest[sorted_idx[0]]
  let mut gval = pvals[sorted_idx[0]]

  while evaluations < max_evals {
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let nvel : Array[Double] = []
      let npos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let raw = velocities[i][d] +
          c1 * r1 * (pbest[i][d] - positions[i][d]) +
          c2 * r2 * (gbest[d] - positions[i][d])
        let v = chi * raw
        let v_clamped = max(-v_max, min(v_max, v))
        nvel.push(v_clamped)
        let pos = positions[i][d] + v_clamped
        npos.push(max(bounds.min, min(bounds.max, pos)))
      }
      velocities[i] = nvel
      positions[i] = npos
      let cv = evaluate(f, npos)
      evaluations = evaluations + 1
      if cv < pvals[i] {
        pvals[i] = cv
        pbest[i] = npos
        if cv < gval {
          gval = cv
          gbest = npos
        }
      }
    }
  }

  OptimizationResult::{ best_solution: gbest, best_value: gval, evaluations, }
}

// ===== LPSO: Local-best Particle Swarm Optimization (2002) =====

///|
pub struct LPSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  swarm_size : Int32
}

///|
pub fn make_lpso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  swarm_sz : Int32,
) -> LPSOConfig {
  LPSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    swarm_size: swarm_sz,
  }
}

///|
/// Local-best Particle Swarm Optimization (2002): the same constriction update
/// as CPSO, but the social term follows the best personal position inside the
/// ring neighbourhood {i-1, i, i+1}, which slows information spread and helps
/// on multimodal functions.  Historical-best tracking is added.
pub fn local_best_pso(
  config : LPSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.swarm_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let chi = 0.729843788
  let c1 = 1.496180
  let c2 = 1.496180
  let v_max = 0.2 * range

  let positions : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let pbest : Array[Array[Double]] = []
  let pvals : Array[Double] = []
  for _i in 0..<n {
    let pos = random_solution(dim, bounds, current_rng)
    positions.push(pos)
    let vel : Array[Double] = []
    for _d in 0..<dim {
      vel.push((RNG::next(current_rng) * 2.0 - 1.0) * v_max)
    }
    velocities.push(vel)
    pbest.push(pos)
    pvals.push(evaluate(f, pos))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(pvals)
  let mut gbest : Array[Double] = pbest[sorted_idx[0]]
  let mut gval = pvals[sorted_idx[0]]

  while evaluations < max_evals {
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      // ring neighbourhood {i-1, i, i+1}: local best among personal bests
      let left = if i == 0 { n - 1 } else { i - 1 }
      let right = if i + 1 == n { 0 } else { i + 1 }
      let mut lbest = i
      if pvals[left] < pvals[lbest] {
        lbest = left
      }
      if pvals[right] < pvals[lbest] {
        lbest = right
      }
      let nvel : Array[Double] = []
      let npos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let raw = velocities[i][d] +
          c1 * r1 * (pbest[i][d] - positions[i][d]) +
          c2 * r2 * (pbest[lbest][d] - positions[i][d])
        let v = chi * raw
        let v_clamped = max(-v_max, min(v_max, v))
        nvel.push(v_clamped)
        let pos = positions[i][d] + v_clamped
        npos.push(max(bounds.min, min(bounds.max, pos)))
      }
      velocities[i] = nvel
      positions[i] = npos
      let cv = evaluate(f, npos)
      evaluations = evaluations + 1
      if cv < pvals[i] {
        pvals[i] = cv
        pbest[i] = npos
        if cv < gval {
          gval = cv
          gbest = npos
        }
      }
    }
  }

  OptimizationResult::{ best_solution: gbest, best_value: gval, evaluations, }
}

// ===== WWO: Water Wave Optimization (2015) =====

///|
pub struct WWOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wwo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WWOConfig {
  WWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Water Wave Optimization (2015), adapted to the continuous black-box setting:
/// each wave either propagates with a wavelength that shrinks with the budget
/// fraction, refracts towards the midpoint with the incumbent best wave, or
/// breaks into fine ripples around itself.  Greedy acceptance and
/// historical-best tracking are added.
pub fn water_wave_optimization(
  config : WWOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let lambda = range * 0.3 * rb
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let cand : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.45 {
          // propagation: move with the current wavelength
          pop[i][d] +
          RNG::next(current_rng) * lambda * (2.0 * RNG::next(current_rng) - 1.0)
        } else if behavior < 0.75 {
          // refraction: collapse towards the midpoint with the best wave
          let mid = (pop[i][d] + x_best[d]) * 0.5
          mid +
          RNG::next(current_rng) *
          lambda *
          0.2 *
          (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // breaking: fine ripples around the wave itself
          pop[i][d] +
          RNG::next(current_rng) *
          lambda *
          0.1 *
          (2.0 * RNG::next(current_rng) - 1.0)
        }
        cand.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, cand)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = cand
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = cand
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WHO: Wild Horse Optimizer (2021) =====

///|
pub struct WHOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_who_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WHOConfig {
  WHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Wild Horse Optimizer (2021): foals graze towards the lead horse (best),
/// adults mate by blending with a random companion, and stallions drive
/// companions away from their territory.  Amplitudes decay with the budget
/// fraction.  Greedy acceptance and historical-best tracking are added.
pub fn wild_horse_optimizer(
  config : WHOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let values : Array[Double] = []
  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    pop.push(sol)
    values.push(evaluate(f, sol))
  }
  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = pop[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let order = sort_indices(values)
    let x_best = pop[order[0]]
    for i in 0..<n {
      if evaluations >= max_evals {
        break
      }
      let behavior = RNG::next(current_rng)
      let peer = pop[(RNG::next(current_rng) * n.to_double()).to_int() % n]
      let candidate : Array[Double] = []
      for d in 0..<dim {
        let pos = if behavior < 0.34 {
          // foal grazing: learn from the lead horse
          pop[i][d] +
          RNG::next(current_rng) * (x_best[d] - pop[i][d]) +
          rb * range * 0.05 * (2.0 * RNG::next(current_rng) - 1.0)
        } else if behavior < 0.67 {
          // mating: blend with a random companion
          pop[i][d] +
          RNG::next(current_rng) * (peer[d] - pop[i][d]) +
          rb * range * 0.1 * (2.0 * RNG::next(current_rng) - 1.0)
        } else {
          // stallion drive: push the companion away from the territory
          pop[i][d] + RNG::next(current_rng) * rb * (pop[i][d] - peer[d])
        }
        candidate.push(max(bounds.min, min(bounds.max, pos)))
      }
      let cv = evaluate(f, candidate)
      evaluations = evaluations + 1
      if cv < values[i] {
        pop[i] = candidate
        values[i] = cv
        if cv < best_val {
          best_val = cv
          best_pos = candidate
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}
