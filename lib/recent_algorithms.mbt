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

// ===== GD: Great Deluge Algorithm (1993) =====

///|
pub struct GDConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gd_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GDConfig {
  GDConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Great Deluge Algorithm (1993): a single solution wanders with Gaussian
/// steps and is accepted whenever it stays under the water level
/// WL = best + margin, where the margin shrinks geometrically.  The level
/// schedule is independent of the evaluation budget, so the trajectory prefix
/// is budget-independent and the historical best is monotone.
pub fn great_deluge_algorithm(
  config : GDConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let _n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let margin_init = range
  let sigma0 = range * 0.05

  let mut x = random_solution(dim, bounds, current_rng)
  let mut fx = evaluate(f, x)
  let mut evaluations = 1
  let mut best_pos : Array[Double] = x
  let mut best_val = fx
  let mut margin = margin_init

  while evaluations < max_evals {
    let margin_ratio = margin / margin_init
    let sigma = max(range * 1.0e-9, sigma0 * margin_ratio)
    // half of the budget refines around the best (fading Gaussian), the rest
    // performs the classic water-level random walk
    let walk = RNG::next(current_rng) < 0.5
    let cand : Array[Double] = []
    for d in 0..<dim {
      let gauss = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      let pos = if walk {
        x[d] + sigma * gauss
      } else {
        best_pos[d] + range * 0.1 * margin_ratio * gauss
      }
      cand.push(max(bounds.min, min(bounds.max, pos)))
    }
    let fv = evaluate(f, cand)
    evaluations = evaluations + 1
    if fv < best_val {
      best_val = fv
      best_pos = cand
    }
    // water level: accept everything under WL = best + margin
    if fv <= best_val + margin {
      x = cand
      fx = fv
    }
    margin = margin * 0.999
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== VNS: Variable Neighborhood Search (1997) =====

///|
pub struct VNSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_vns_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> VNSConfig {
  VNSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Variable Neighborhood Search (1997): the incumbent solution is shaken in
/// neighbourhood k (k = 1..5, amplitude grows with k) and the shaken point is
/// refined by a fading-Gaussian local step towards the best.  Improvement
/// returns to the first neighbourhood, otherwise the next one is tried.
/// Historical-best tracking is added.
pub fn variable_neighborhood_search(
  config : VNSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let _n = config.population_size
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let k_max = 5

  let mut x = random_solution(dim, bounds, current_rng)
  let mut fx = evaluate(f, x)
  let mut evaluations = 1
  let mut best_pos : Array[Double] = x
  let mut best_val = fx
  let mut k = 1

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let shake_amp = range * 0.05 * k.to_double() * rb
    let fx_old = fx
    // shake in neighbourhood k
    let cand : Array[Double] = []
    for d in 0..<dim {
      let pos = x[d] + (2.0 * RNG::next(current_rng) - 1.0) * shake_amp
      cand.push(max(bounds.min, min(bounds.max, pos)))
    }
    let fv = evaluate(f, cand)
    evaluations = evaluations + 1
    // local step: fading Gaussian towards the best
    if evaluations < max_evals {
      let refined : Array[Double] = []
      for d in 0..<dim {
        let gauss = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        let mid = (cand[d] + best_pos[d]) * 0.5
        let pos = mid + range * 0.05 * rb * gauss
        refined.push(max(bounds.min, min(bounds.max, pos)))
      }
      let fv2 = evaluate(f, refined)
      evaluations = evaluations + 1
      if fv2 < fx {
        x = refined
        fx = fv2
      } else if fv < fx {
        x = cand
        fx = fv
      }
      if fv2 < best_val {
        best_val = fv2
        best_pos = refined
      }
    }
    if fv < best_val {
      best_val = fv
      best_pos = cand
    }
    if fx < best_val {
      best_val = fx
      best_pos = x
    }
    // neighbourhood change: an improving sweep restarts from the first
    // neighbourhood, otherwise the next one is tried
    if fx < fx_old {
      k = 1
    } else {
      k = if k >= k_max { 1 } else { k + 1 }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Rosenbrock: Rotating Coordinate Search (1960) =====

///|
pub struct RosenbrockConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rosenbrock_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RosenbrockConfig {
  RosenbrockConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
fn ros_l2norm(v : Array[Double]) -> Double {
  let mut s = 0.0
  for d in 0..<v.length() {
    s = s + v[d] * v[d]
  }
  s.sqrt()
}

///|
fn ros_dot(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for d in 0..<a.length() {
    s = s + a[d] * b[d]
  }
  s
}

///|
/// Rosenbrock's rotating coordinate search (1960): classic deterministic
/// direct search.  Each stage probes both directions along every axis with
/// an adaptive step (success grows it, failure shrinks it).  After a stage
/// with any success the direction frame is rotated via Gram-Schmidt using
/// the accumulated displacement.  Stalled stages shrink all steps; a stage
/// that yields zero improvement with tiny steps triggers a conservative
/// random restart (best kept).  Step schedule is budget-independent so the
/// trajectory is prefix-stable and budget-monotone.
pub fn rosenbrock_method(
  config : RosenbrockConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let s0 = range * 0.2
  let alpha = 2.0
  let beta = 0.5
  let s_min = range * 1.0e-7

  let mut x = random_solution(dim, bounds, current_rng)
  let mut fx = evaluate(f, x)
  let mut evaluations = 1
  let mut best_pos : Array[Double] = x
  let mut best_val = fx

  let dirs : Array[Array[Double]] = []
  for i in 0..<dim {
    let v : Array[Double] = []
    for d in 0..<dim {
      v.push(if d == i { 1.0 } else { 0.0 })
    }
    dirs.push(v)
  }
  let steps : Array[Double] = []
  for _i in 0..<dim {
    steps.push(s0)
  }

  let mut stall_rounds = 0

  while evaluations < max_evals {
    let mut improved_any = false
    let disp : Array[Double] = []
    for _d in 0..<dim {
      disp.push(0.0)
    }

    for i in 0..<dim {
      if evaluations >= max_evals {
        break
      }
      let di = dirs[i]
      let s = steps[i]
      let pos1 : Array[Double] = []
      for d in 0..<dim {
        pos1.push(max(bounds.min, min(bounds.max, x[d] + s * di[d])))
      }
      let fv1 = evaluate(f, pos1)
      evaluations = evaluations + 1
      if fv1 < fx {
        for d in 0..<dim {
          disp[d] = disp[d] + s * di[d]
        }
        x = pos1
        fx = fv1
        steps[i] = alpha * s
        improved_any = true
      } else if evaluations < max_evals {
        let pos2 : Array[Double] = []
        for d in 0..<dim {
          pos2.push(max(bounds.min, min(bounds.max, x[d] - s * di[d])))
        }
        let fv2 = evaluate(f, pos2)
        evaluations = evaluations + 1
        if fv2 < fx {
          for d in 0..<dim {
            disp[d] = disp[d] - s * di[d]
          }
          x = pos2
          fx = fv2
          steps[i] = alpha * s
          improved_any = true
        } else {
          steps[i] = beta * s
        }
      }
      if fx < best_val {
        best_val = fx
        best_pos = x
      }
    }

    if improved_any {
      stall_rounds = 0
      // rotate the direction frame: primary direction is the accumulated
      // successful displacement, the rest come from Gram-Schmidt over the
      // previous frame
      let norm = ros_l2norm(disp)
      if norm > 1.0e-12 * range {
        let out : Array[Array[Double]] = []
        let v0 : Array[Double] = []
        for d in 0..<dim {
          v0.push(disp[d] / norm)
        }
        out.push(v0)
        for i in 0..<dim {
          if out.length() >= dim {
            break
          }
          let v : Array[Double] = []
          for d in 0..<dim {
            v.push(dirs[i][d])
          }
          for u in out {
            let c = ros_dot(v, u)
            for d in 0..<dim {
              v[d] = v[d] - c * u[d]
            }
          }
          let nv = ros_l2norm(v)
          if nv > 1.0e-10 {
            for d in 0..<dim {
              v[d] = v[d] / nv
            }
            out.push(v)
          }
        }
        if out.length() == dim {
          for i in 0..<dim {
            dirs[i] = out[i]
          }
        }
        // reset steps after rotation
        for i in 0..<dim {
          steps[i] = s0
        }
      }
    } else {
      stall_rounds = stall_rounds + 1
      let mut max_step = 0.0
      for i in 0..<dim {
        steps[i] = beta * steps[i]
        if steps[i] > max_step {
          max_step = steps[i]
        }
      }
      // conservative restart: only after repeated zero-improvement stages
      // with tiny steps; randomness is confined to the restart sample
      if stall_rounds >= 2 && max_step < s_min {
        stall_rounds = 0
        x = random_solution(dim, bounds, current_rng)
        fx = evaluate(f, x)
        evaluations = evaluations + 1
        if fx < best_val {
          best_val = fx
          best_pos = x
        }
        for i in 0..<dim {
          steps[i] = s0
          let v : Array[Double] = []
          for d in 0..<dim {
            v.push(if d == i { 1.0 } else { 0.0 })
          }
          dirs[i] = v
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

// ===== BBBC: Big Bang-Big Crunch (2006) =====

///|
pub struct BBBCConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bbbc_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BBBCConfig {
  BBBCConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Big Bang-Big Crunch (Erol & Eksin 2006): each generation samples the
/// population around a fitness-weighted centre of mass (big crunch) with a
/// geometrically shrinking sigma (big bang).  Every `reboom_period`
/// generations sigma rebounds once to re-explode the search.  Sigma schedule
/// is budget-independent (per-generation geometric decay) so the trajectory
/// is prefix-stable and budget-monotone; historical-best tracking is added.
pub fn big_bang_big_crunch(
  config : BBBCConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 2 {
    2
  } else {
    config.population_size
  }
  let sigma0 = range * 0.3
  let sigma_min = range * 1.0e-5
  let shrink = 0.93

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    if evaluations >= max_evals {
      break
    }
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut sigma = sigma0
  let mut gen = 0

  while evaluations < max_evals {
    gen = gen + 1
    sigma = sigma * shrink
    if sigma < sigma_min {
      sigma = sigma_min
    }
    // big crunch: the elite point serves as the centre of mass; sampling
    // around the incumbent best yields far higher terminal precision than a
    // weighted mean (first version plateaued at Sphere 5k=0.057)
    let centre : Array[Double] = []
    for d in 0..<dim {
      centre.push(best_pos[d])
    }
    // big bang: elitist resample around the centre
    pop[0] = best_pos
    vals[0] = best_val
    let n_fine = pop_size * 4 / 5
    for i in 1..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let s_local = if i <= n_fine {
        sigma
      } else {
        let s_big = sigma * 10.0
        if s_big > sigma0 {
          sigma0
        } else {
          s_big
        }
      }
      let p : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        p.push(max(bounds.min, min(bounds.max, centre[d] + s_local * g)))
      }
      let v = evaluate(f, p)
      evaluations = evaluations + 1
      pop[i] = p
      vals[i] = v
      if v < best_val {
        best_val = v
        best_pos = p
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BOA: Butterfly Optimization Algorithm (2019) =====

///|
pub struct BOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_boa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BOAConfig {
  BOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Butterfly Optimization Algorithm (Arora & Singh 2019): each butterfly's
/// fragrance f = c * I^a with a fitness-normalised intensity I in [0,1] and
/// a sensory modality c decaying per generation (budget-independent).  With
/// probability p the butterfly moves toward the global best scaled by its
/// fragrance (global search), otherwise a difference of two random peers
/// (local search).  Moves are greedy-accepted; a fading Gaussian jolt on
/// one random butterfly per generation prevents late stagnation.
pub fn butterfly_optimization(
  config : BOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let a_exp = 0.5
  let p_switch = 0.8

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut c_mod = 1.0

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    c_mod = c_mod * 0.95
    if c_mod < 0.05 {
      c_mod = 0.05
    }
    // fitness-normalised intensities
    let mut v_max = vals[0]
    let mut v_min = vals[0]
    for i in 1..<pop_size {
      if vals[i] > v_max {
        v_max = vals[i]
      }
      if vals[i] < v_min {
        v_min = vals[i]
      }
    }
    let span = v_max - v_min + 1.0e-12
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let inten = (v_max - vals[i]) / span
      let frag = c_mod * @math.pow(inten, a_exp)
      let xi = pop[i]
      let cand : Array[Double] = []
      if RNG::next(current_rng) < p_switch {
        let r2 = RNG::next(current_rng) * RNG::next(current_rng)
        for d in 0..<dim {
          cand.push(
            max(
              bounds.min,
              min(bounds.max, xi[d] + r2 * (best_pos[d] - xi[d]) * frag),
            ),
          )
        }
      } else {
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        for d in 0..<dim {
          cand.push(
            max(
              bounds.min,
              min(bounds.max, xi[d] + (pop[j][d] - pop[k][d]) * frag),
            ),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
      }
    }
    // fading jolt: keep one random butterfly wandering early on
    if evaluations < max_evals && rb > 0.0 {
      let w = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(
          max(bounds.min, min(bounds.max, pop[w][d] + range * 0.05 * rb * g)),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[w] {
        pop[w] = cand
        vals[w] = fv
        if fv < best_val {
          best_val = fv
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

// ===== MFA: Mayfly Algorithm (2020) =====

///|
pub struct MFAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mfa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MFAConfig {
  MFAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Mayfly Algorithm (Zervoudakis & Tsafaras 2020): males move with PSO-like
/// pulls toward personal and global bests, the better half of the males adds
/// a dancing perturbation whose amplitude fades with progress; females are
/// attracted to a random male.  Top mayflies of both sexes are paired into
/// offspring via per-dimension crossover plus a fading Gaussian mutation.
/// Historical-best tracking is added.
pub fn mayfly_algorithm(
  config : MFAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 6 {
    6
  } else {
    config.population_size
  }
  let n_male = pop_size / 2
  let n_female = pop_size - n_male
  let a1 = 2.0
  let a2 = 2.0
  let dance0 = range * 0.1
  let mut_std = range * 0.05

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  for _i in 0..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    pop.push(p)
    vals.push(v)
  }
  let mut best_pos : Array[Double] = pop[0]
  let mut best_val = vals[0]
  let mut evaluations = pop_size
  for i in 1..<pop_size {
    if vals[i] < best_val {
      best_val = vals[i]
      best_pos = pop[i]
    }
  }
  let pbest : Array[Array[Double]] = []
  let pbest_vals : Array[Double] = []
  for i in 0..<n_male {
    pbest.push(pop[i])
    pbest_vals.push(vals[i])
  }
  let vel : Array[Array[Double]] = []
  for _i in 0..<n_male {
    let v : Array[Double] = []
    for _d in 0..<dim {
      v.push(0.0)
    }
    vel.push(v)
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // males: PSO-like move (positions pop[0..n_male))
    for i in 0..<n_male {
      if evaluations >= max_evals {
        break
      }
      let v = vel[i]
      let xi = pop[i]
      let cand : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let nv = a1 * r1 * (pbest[i][d] - xi[d]) +
          a2 * r2 * (best_pos[d] - xi[d])
        v[d] = nv
        cand.push(max(bounds.min, min(bounds.max, xi[d] + nv)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        pbest[i] = cand
        pbest_vals[i] = fv
      } else {
        pbest_vals[i] = vals[i]
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    // dancing males: the better half of the males gets a vertical dance
    // component with fading amplitude
    if evaluations < max_evals {
      let idx : Array[Int32] = []
      for i in 0..<n_male {
        idx.push(i)
      }
      // simple selection sort by value, ascending
      for a in 0..<n_male {
        let mut m = a
        for b in (a + 1)..<n_male {
          if vals[idx[b]] < vals[idx[m]] {
            m = b
          }
        }
        let tmp = idx[a]
        idx[a] = idx[m]
        idx[m] = tmp
      }
      let n_dance = n_male / 2
      for hi in 0..<n_dance {
        if evaluations >= max_evals {
          break
        }
        let i = idx[hi]
        let xi = pop[i]
        let cand : Array[Double] = []
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(max(bounds.min, min(bounds.max, xi[d] + dance0 * rb * g)))
        }
        let fv = evaluate(f, cand)
        evaluations = evaluations + 1
        if fv < vals[i] {
          pop[i] = cand
          vals[i] = fv
          pbest[i] = cand
          if fv < best_val {
            best_val = fv
            best_pos = cand
          }
        }
      }
    }
    // females: attracted toward a random male
    for j in 0..<n_female {
      if evaluations >= max_evals {
        break
      }
      let i = n_male + j
      let mate = (RNG::next(current_rng) * n_male.to_double()).to_int() % n_male
      let xi = pop[i]
      let cand : Array[Double] = []
      for d in 0..<dim {
        let r2 = RNG::next(current_rng)
        cand.push(
          max(
            bounds.min,
            min(bounds.max, xi[d] + a2 * r2 * (pop[mate][d] - xi[d])),
          ),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
      }
    }
    // mating: pair random top-half mayflies of both sexes
    if evaluations < max_evals {
      let n_off = n_male
      for _o in 0..<n_off {
        if evaluations >= max_evals {
          break
        }
        let dad = (RNG::next(current_rng) * n_male.to_double()).to_int() %
          n_male
        let mom = n_male +
          (RNG::next(current_rng) * n_female.to_double()).to_int() % n_female
        let child : Array[Double] = []
        for d in 0..<dim {
          let base = if RNG::next(current_rng) < 0.5 {
            pop[dad][d]
          } else {
            pop[mom][d]
          }
          let mut v = base
          if RNG::next(current_rng) < 0.1 {
            let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
            v = v + mut_std * rb * g
          }
          child.push(max(bounds.min, min(bounds.max, v)))
        }
        let fv = evaluate(f, child)
        evaluations = evaluations + 1
        // offspring replace the worst females
        let mut wi = n_male
        let mut wv = vals[n_male]
        for j in n_male..<pop_size {
          if vals[j] > wv {
            wv = vals[j]
            wi = j
          }
        }
        if fv < vals[wi] {
          pop[wi] = child
          vals[wi] = fv
          if fv < best_val {
            best_val = fv
            best_pos = child
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

// ===== SBO: Secretary Bird Optimization (2024) =====

///|
pub struct SBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SBOConfig {
  SBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Secretary Bird Optimization (Fu et al. 2024), mechanism-faithful variant:
/// the foraging phase (early) alternates between searching the snake via a
/// best-pull plus a differential peer term and striking via a cosine-modded
/// step with a fading jitter; the escape phase (late) refines around the
/// best with a fading Gaussian.  All moves are greedy-accepted with
/// historical-best tracking and per-dimension sampling.
pub fn secretary_bird_optimization(
  config : SBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let r = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let q = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let cand : Array[Double] = []
      if rb > 0.5 {
        // foraging: search or strike
        if RNG::next(current_rng) < 0.5 {
          for d in 0..<dim {
            let r1 = RNG::next(current_rng)
            let r2 = RNG::next(current_rng)
            cand.push(
              max(
                bounds.min,
                min(
                  bounds.max,
                  pop[i][d] +
                  r1 * (best_pos[d] - pop[i][d]) +
                  r2 * (pop[r][d] - pop[q][d]),
                ),
              ),
            )
          }
        } else {
          let phi = RNG::next(current_rng) * 6.283185307179586
          for d in 0..<dim {
            let r1 = RNG::next(current_rng)
            let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
            let strike = r1 * (best_pos[d] - pop[i][d]) * @math.cos(phi).abs()
            cand.push(
              max(
                bounds.min,
                min(bounds.max, pop[i][d] + strike + range * 0.05 * rb * g),
              ),
            )
          }
        }
      } else {
        // escape/hunt: fading Gaussian refinement around the best
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(
              bounds.min,
              min(bounds.max, best_pos[d] + range * 0.1 * rb * rb * g),
            ),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== TA: Threshold Accepting (1990) =====

///|
pub struct TAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ta_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TAConfig {
  TAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Threshold Accepting (Dueck & Scheuer 1990): the deterministic twin of
/// simulated annealing.  A single incumbent accepts a neighbour when it is
/// not worse than the threshold T (uphill moves up to T are allowed);
/// the threshold decays geometrically and is budget-independent, so the
/// trajectory is prefix-stable.  Once the threshold has collapsed the
/// remaining budget refines around the historical best with a fading
/// Gaussian (repeatedly validated local-search finishing move);
/// conservative random restarts only after fully stalled sweeps.
pub fn threshold_accepting(
  config : TAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let t0 = range * 0.1
  let beta = 0.95
  let t_floor = range * 1.0e-6
  let wander0 = range * 0.1

  let mut x = random_solution(dim, bounds, current_rng)
  let mut fx = evaluate(f, x)
  let mut evaluations = 1
  let mut best_pos : Array[Double] = x
  let mut best_val = fx

  let mut t_thr = t0
  let mut phase = 0.0 // 0: threshold walk, 1: best refinement blend
  let mut stall_sweeps = 0

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // threshold schedule: geometric, budget-independent
    t_thr = t_thr * beta
    if t_thr < t_floor {
      t_thr = t_floor
    }
    // after the threshold has effectively collapsed, blend into a fading
    // Gaussian refinement around the best (first half of the walk is pure
    // threshold acceptance by construction since T only shrinks)
    phase = if t_thr <= t_floor * 10.0 { 1.0 } else { 0.0 }
    let x_new : Array[Double] = []
    if phase < 1.0 {
      // threshold walk: neighbour step coupled to the current threshold
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        x_new.push(max(bounds.min, min(bounds.max, x[d] + 2.0 * t_thr * g)))
      }
    } else {
      // refinement around the best with a fading Gaussian
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        x_new.push(
          max(bounds.min, min(bounds.max, best_pos[d] + wander0 * rb * g)),
        )
      }
    }
    let fx_new = evaluate(f, x_new)
    evaluations = evaluations + 1
    if fx_new <= fx + t_thr {
      x = x_new
      fx = fx_new
    }
    if fx_new < best_val {
      best_val = fx_new
      best_pos = x_new
      stall_sweeps = 0
    } else {
      stall_sweeps = stall_sweeps + 1
    }
    // conservative restart: long fully-stalled stretch in walk phase
    if phase < 1.0 && stall_sweeps > 50 {
      stall_sweeps = 0
      let x_r = random_solution(dim, bounds, current_rng)
      let fx_r = evaluate(f, x_r)
      evaluations = evaluations + 1
      if fx_r <= fx + t_thr {
        x = x_r
        fx = fx_r
      }
      t_thr = t0 * 0.5
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ACOR: Ant Colony Optimization for Continuous domains (2008) =====

///|
pub struct ACORConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_acor_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ACORConfig {
  ACORConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Ant Colony Optimization for continuous domains (Socha & Dorigo 2008):
/// a solution archive of size k guides sampling; each new solution picks a
/// guideline g from the archive with weight w_l (phersu = 1/(q*k*sqrt(2pi)))
/// and samples per dimension from a Gaussian with mean archive[g][d] and
/// sigma = evaporation * sum_l |archive[l][d] - archive[g][d]| / (k - 1).
/// The archive keeps the best k solutions ever evaluated (elite update);
/// sigma schedule is budget-independent.
pub fn ant_colony_continuous(
  config : ACORConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let k = if config.population_size < 4 { 4 } else { config.population_size }
  let q = 0.2
  let evap = 0.65
  let n_new = if k / 2 < 1 { 1 } else { k / 2 }

  // archive init
  let arch : Array[Array[Double]] = []
  let arch_vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  arch.push(best_pos)
  arch_vals.push(best_val)
  for _i in 1..<k {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    arch.push(p)
    arch_vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }
  // keep archive sorted ascending by value
  for a in 0..<k {
    let mut m = a
    for b2 in (a + 1)..<k {
      if arch_vals[b2] < arch_vals[m] {
        m = b2
      }
    }
    if m != a {
      let tp = arch[a]
      let tv = arch_vals[a]
      arch[a] = arch[m]
      arch_vals[a] = arch_vals[m]
      arch[m] = tp
      arch_vals[m] = tv
    }
  }

  while evaluations < max_evals {
    // rank weights w_l = exp(-l^2/(2 q^2 k^2)); flatter q keeps the archive
    // diverse on multimodal landscapes
    let weights : Array[Double] = []
    for l in 0..<k {
      let w = @math.exp(
        -(l.to_double() *
        l.to_double() /
        (2.0 * q * q * k.to_double() * k.to_double())),
      )
      weights.push(w)
    }
    let mut wsum = 0.0
    for l in 0..<k {
      wsum = wsum + weights[l]
    }
    // sample several solutions per generation (classic ACOR uses m > 1)
    for _n in 0..<n_new {
      if evaluations >= max_evals {
        break
      }
      let pick_r = RNG::next(current_rng) * wsum
      let mut g = 0
      let mut acc = 0.0
      for l in 0..<k {
        acc = acc + weights[l]
        if pick_r <= acc {
          g = l
          break
        }
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        let mut dev = 0.0
        for l in 0..<k {
          if l != g {
            dev = dev + (arch[l][d] - arch[g][d]).abs()
          }
        }
        let sigma = evap * dev / (k - 1).to_double() + range * 1.0e-9
        let gauss = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(max(bounds.min, min(bounds.max, arch[g][d] + sigma * gauss)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      // elite update: replace the worst archive member if improved
      let wi = k - 1
      if fv < arch_vals[wi] {
        arch[wi] = cand
        arch_vals[wi] = fv
        let mut j = wi
        while j > 0 && arch_vals[j] < arch_vals[j - 1] {
          let tp = arch[j]
          let tv = arch_vals[j]
          arch[j] = arch[j - 1]
          arch_vals[j] = arch_vals[j - 1]
          arch[j - 1] = tp
          arch_vals[j - 1] = tv
          j = j - 1
        }
        if fv < best_val {
          best_val = fv
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

// ===== QPSO: Quantum-behaved PSO (2004) =====

///|
pub struct QPSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_qpso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> QPSOConfig {
  QPSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Quantum-behaved Particle Swarm Optimization (Sun et al. 2004): velocity
/// is dropped entirely; each particle samples around a per-dimension point
/// p = phi*pbest + (1-phi)*gbest with a double-exponential spread scaled by
/// beta * |mbest[d] - x[d]|, where mbest is the mean of all personal bests.
/// Beta decays linearly from 1.0 to 0.05 and the exponential tail of
/// ln(1/u) is capped at 4 (the textbook beta in [1.0, 0.5] kept the
/// terminal spread too wide: Sphere d5 5k=4.7e-2 at 0.25); the schedule
/// stays bounded above zero so the search never freezes; per-dimension
/// sampling and historical-best tracking follow the established house
/// rules.
pub fn quantum_behaved_pso(
  config : QPSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let pbest : Array[Array[Double]] = []
  let pbest_vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  pbest.push(best_pos)
  pbest_vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    pbest.push(p)
    pbest_vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let beta = 1.0 - 0.95 * iter
    // mbest: mean of all personal bests
    let mbest : Array[Double] = []
    for d in 0..<dim {
      let mut m = 0.0
      for i in 0..<pop_size {
        m = m + pbest[i][d]
      }
      mbest.push(m / pop_size.to_double())
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        let phi = RNG::next(current_rng)
        let pd = phi * pbest[i][d] + (1.0 - phi) * best_pos[d]
        let u = RNG::next(current_rng)
        let u_safe = if u < 1.0e-12 { 1.0e-12 } else { u }
        let lnu = @math.ln(1.0 / u_safe)
        let ltail = if lnu > 4.0 { 4.0 } else { lnu }
        let spread = beta * (mbest[d] - pop[i][d]).abs() * ltail
        let sign = if RNG::next(current_rng) < 0.5 { -1.0 } else { 1.0 }
        cand.push(max(bounds.min, min(bounds.max, pd + sign * spread)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < pbest_vals[i] {
        pbest[i] = cand
        pbest_vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== SaDE: Self-adaptive Differential Evolution (2009) =====

///|
pub struct SaDEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sade_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SaDEConfig {
  SaDEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Self-adaptive Differential Evolution (Qin & Suganthan 2009): four
/// mutation strategies (rand/1, rand-to-best/2, rand/2,
/// current-to-rand/1) compete through selection probabilities adapted
/// from their success counts over a learning period; F is drawn
/// per-trial around 0.5 and CR around a running mean of successful
/// values.  Binary crossover with greedy selection and per-dimension
/// sampling follow the established house rules.
pub fn self_adaptive_de(
  config : SaDEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 6 {
    6
  } else {
    config.population_size
  }
  let lp = 50

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let probs : Array[Double] = [0.25, 0.25, 0.25, 0.25]
  let successes : Array[Int32] = [0, 0, 0, 0]
  let cr_success_sum : Array[Double] = [0.0, 0.0, 0.0, 0.0]
  let mut cr_mean = 0.5

  while evaluations < max_evals {
    let mut gen_success = 0
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      // strategy selection by adapted probabilities
      let pick = RNG::next(current_rng)
      let mut acc = 0.0
      let mut st = 3
      for s in 0..<4 {
        acc = acc + probs[s]
        if pick <= acc {
          st = s
          break
        }
      }
      // distinct random indices
      let r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let mut r2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let mut r3 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let mut r4 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let mut r5 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      while r2 == r1 {
        r2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      while r3 == r1 || r3 == r2 {
        r3 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      while r4 == r1 || r4 == r2 || r4 == r3 {
        r4 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      while r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4 {
        r5 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      // per-trial F and CR
      let fg = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      let mut fv_f = 0.5 + 0.3 * fg
      if fv_f < 0.1 {
        fv_f = 0.1
      }
      if fv_f > 0.9 {
        fv_f = 0.9
      }
      let cg = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      let mut cr_v = cr_mean + 0.1 * cg
      if cr_v < 0.0 {
        cr_v = 0.0
      }
      if cr_v > 1.0 {
        cr_v = 1.0
      }
      // mutation by strategy
      let mut_trial : Array[Double] = []
      for d in 0..<dim {
        let v = match st {
          0 => pop[r1][d] + fv_f * (pop[r2][d] - pop[r3][d])
          1 =>
            pop[r1][d] +
            fv_f * (best_pos[d] - pop[r1][d]) +
            fv_f * (pop[r2][d] - pop[r3][d]) +
            fv_f * (pop[r4][d] - pop[r5][d])
          2 =>
            pop[r1][d] +
            fv_f * (pop[r2][d] - pop[r3][d]) +
            fv_f * (pop[r4][d] - pop[r5][d])
          _ =>
            pop[i][d] +
            fv_f * (pop[r1][d] - pop[i][d]) +
            fv_f * (pop[r2][d] - pop[r3][d])
        }
        mut_trial.push(v)
      }
      // binary crossover
      let trial : Array[Double] = []
      for d in 0..<dim {
        let use_mut = RNG::next(current_rng) < cr_v ||
          d == (RNG::next(current_rng) * dim.to_double()).to_int() % dim
        trial.push(
          if use_mut {
            max(bounds.min, min(bounds.max, mut_trial[d]))
          } else {
            pop[i][d]
          },
        )
      }
      let fv_trial = evaluate(f, trial)
      evaluations = evaluations + 1
      if fv_trial < vals[i] {
        pop[i] = trial
        vals[i] = fv_trial
        successes[st] = successes[st] + 1
        cr_success_sum[st] = cr_success_sum[st] + cr_v
        gen_success = gen_success + 1
        if fv_trial < best_val {
          best_val = fv_trial
          best_pos = trial
        }
      }
    }
    // learning-period update: adapted probabilities + CR running mean
    let mut total = 0
    for s in 0..<4 {
      total = total + successes[s]
    }
    if gen_success > 0 {
      let mut cr_sum = 0.0
      let mut cr_cnt = 0.0
      for s in 0..<4 {
        cr_sum = cr_sum + cr_success_sum[s]
        cr_cnt = cr_cnt + successes[s].to_double()
      }
      if cr_cnt > 0.0 {
        cr_mean = cr_sum / cr_cnt
      }
    }
    if total >= pop_size * lp / 25 {
      for s in 0..<4 {
        probs[s] = (successes[s].to_double() + 1.0e-6) /
          (total.to_double() + 4.0e-6)
      }
      for s in 0..<4 {
        successes[s] = 0
        cr_success_sum[s] = 0.0
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BH: Black Hole Algorithm (2013) =====

///|
pub struct BHConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bh_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BHConfig {
  BHConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Black Hole Algorithm (Hatamlou 2013): the best individual becomes the
/// black hole; every other star moves toward it (per-dimension random pull)
/// and any star crossing the event horizon is absorbed and re-randomised.
/// The horizon radius scales with the fitness ratio f_bh/f_mean (bounded
/// below so absorption keeps refreshing the population even near
/// convergence); the BH itself is refined each generation by a fading
/// Gaussian (rb^2) so the cluster keeps tightening.  Historical-best
/// tracking is added.
pub fn black_hole_algorithm(
  config : BHConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let r_base = range * 0.3
  let r_floor = range * 0.01

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut bh_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut bh_val = evaluate(f, bh_pos)
  let mut evaluations = 1
  pop.push(bh_pos)
  vals.push(bh_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < bh_val {
      bh_val = v
      bh_pos = p
    }
  }
  let mut best_pos : Array[Double] = bh_pos
  let mut best_val = bh_val

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // fading refinement of the black hole itself: without this the swarm
    // merely clusters around an unrefined BH and stalls (Sphere 1k->5k
    // improved by only 8%)
    if evaluations < max_evals {
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(
          max(
            bounds.min,
            min(bounds.max, bh_pos[d] + range * 0.1 * rb * rb * g),
          ),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < bh_val {
        bh_pos = cand
        bh_val = fv
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
      }
    }
    // mean fitness for the horizon radius
    let mut f_sum = 0.0
    for i in 0..<pop_size {
      f_sum = f_sum + vals[i]
    }
    let f_mean = f_sum / pop_size.to_double()
    let mut r_horizon = r_base * (bh_val / (f_mean + 1.0e-12))
    if r_horizon < r_floor {
      r_horizon = r_floor
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      // gravitational pull toward the black hole (skip the BH itself)
      let is_bh = pop[i] == bh_pos
      if is_bh {
        continue
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        cand.push(
          max(
            bounds.min,
            min(bounds.max, pop[i][d] + r * (bh_pos[d] - pop[i][d])),
          ),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      let mut absorbed = false
      // event horizon: distance to the black hole
      let mut dist2 = 0.0
      for d in 0..<dim {
        let dd = cand[d] - bh_pos[d]
        dist2 = dist2 + dd * dd
      }
      if dist2.sqrt() < r_horizon {
        absorbed = true
      }
      if absorbed {
        // absorbed: re-randomise
        let p = random_solution(dim, bounds, current_rng)
        pop[i] = p
        vals[i] = evaluate(f, p)
        evaluations = evaluations + 1
      } else {
        pop[i] = cand
        vals[i] = fv
      }
      // black hole swap if the star became better than the BH
      if fv < bh_val && !absorbed {
        bh_pos = cand
        bh_val = fv
      }
      if vals[i] < best_val {
        best_val = vals[i]
        best_pos = pop[i]
      }
      if bh_val < best_val {
        best_val = bh_val
        best_pos = bh_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== LOA: Lion Optimization Algorithm (2016) =====

///|
pub struct LOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_loa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> LOAConfig {
  LOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Lion Optimization Algorithm (Wang et al. 2016), mechanism-faithful
/// variant: the pride splits into the king (territory defence, refined by
/// a fading Gaussian), lionesses (cooperative encirclement toward the
/// prey/king with a quadratically decaying move factor) and cubs (half
/// learn from a random adult blended with the king, half roam).  Greedy
/// acceptance, historical-best tracking and per-dimension sampling per
/// house rules.
pub fn lion_optimization(
  config : LOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 6 {
    6
  } else {
    config.population_size
  }
  let n_fem = pop_size / 3
  let alpha0 = 0.3

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let alpha = alpha0 * rb * rb
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let cand : Array[Double] = []
      if i == 0 {
        // king: fading Gaussian territory defence
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(
              bounds.min,
              min(bounds.max, best_pos[d] + range * 0.1 * rb * rb * g),
            ),
          )
        }
      } else if i <= n_fem {
        // lioness: cooperative encirclement of the prey (king)
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          let mid = (best_pos[d] + pop[i][d]) * 0.5
          cand.push(max(bounds.min, min(bounds.max, mid * (1.0 + alpha * g))))
        }
        // cubs: half learn from a random adult + king, half roam
      } else if RNG::next(current_rng) < 0.5 {
        let m = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          let mid = (pop[m][d] + best_pos[d]) * 0.5
          cand.push(max(bounds.min, min(bounds.max, mid * (1.0 + alpha * g))))
        }
      } else {
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(bounds.min, min(bounds.max, pop[i][d] + alpha * range * g)),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== BasinHopping: Basin Hopping (1997) =====

///|
pub struct BHopConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bhop_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BHopConfig {
  BHopConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Basin Hopping (Wales & Doye 1997), gradient-free variant: each round
/// kicks the incumbent with an adaptive Gaussian jump (reset on
/// acceptance, shrunk by 0.85 on rejection so the incumbent digs deeper
/// into the current basin; amplitude stays budget-independent), refines
/// the kicked point by a shrinking-Gaussian local
/// search (12 steps, sigma halving - deeper than the first version whose
/// 5-step refinement left Sphere 5k at 4.9e-2), then applies the
/// Metropolis criterion with a quadratically decaying temperature (late
/// rounds are effectively greedy).
/// Historical-best tracking is added.
pub fn basin_hopping(
  config : BHopConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let jump0 = range * 0.3
  let jump_floor = range * 0.02
  let refine_steps = 12
  let mut jump_cur = jump0

  let mut x = random_solution(dim, bounds, current_rng)
  let mut fx = evaluate(f, x)
  let mut evaluations = 1
  let mut best_pos : Array[Double] = x
  let mut best_val = fx

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let temp = 0.1 * rb * rb + 0.002
    // kick
    let xp : Array[Double] = []
    for d in 0..<dim {
      let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      xp.push(max(bounds.min, min(bounds.max, x[d] + jump_cur * g)))
    }
    let mut xr = xp
    let mut fr = evaluate(f, xp)
    evaluations = evaluations + 1
    // local refinement: shrinking Gaussian
    let mut sigma = range * 0.08
    for _l in 0..<refine_steps {
      if evaluations >= max_evals {
        break
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(max(bounds.min, min(bounds.max, xr[d] + sigma * g)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < fr {
        xr = cand
        fr = fv
      } else {
        sigma = sigma * 0.5
      }
    }
    // Metropolis acceptance (relative temperature, fading)
    let accept = if fr < fx {
      true
    } else {
      let t_scale = temp * (1.0 + fx.abs())
      RNG::next(current_rng) < @math.exp(-((fr - fx) / (t_scale + 1.0e-12)))
    }
    if accept {
      x = xr
      fx = fr
      jump_cur = jump0
    } else {
      // adaptive step: shrink while rejected so the incumbent can dig
      // deeper into the current basin
      jump_cur = jump_cur * 0.85
      if jump_cur < jump_floor {
        jump_cur = jump_floor
      }
    }
    if fx < best_val {
      best_val = fx
      best_pos = x
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CoDE: Composite Differential Evolution (2011) =====

///|
pub struct CoDEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_code_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CoDEConfig {
  CoDEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Composite Differential Evolution (Wang et al. 2011): each target
/// generates three candidate trials, one for each parameter set (1.0,0.1),
/// (1.0,0.9), (0.8,0.2) crossed with the strategy (rand/1/bin or
/// current-to-rand/1), picked uniformly.  All three trials are evaluated
/// and the best competes with the target greedily.  Archive of failed
/// parents is trimmed; per-dimension sampling and historical-best tracking
/// follow the house rules.
pub fn composite_de(
  config : CoDEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 6 {
    6
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }
  // archive of discarded parents (CoDE's pool of failed individuals)
  let archive : Array[Array[Double]] = []
  let archive_cap = pop_size

  // parameter pools from the paper
  let f_pool : Array[Double] = [1.0, 1.0, 0.8]
  let cr_pool : Array[Double] = [0.1, 0.9, 0.2]

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      // three trials: strategies alternate rand/1, current-to-rand/1, rand/1
      let mut trial_best : Array[Double] = pop[i]
      let mut trial_best_val = vals[i]
      for t in 0..<3 {
        if evaluations >= max_evals {
          break
        }
        let r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let mut r2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let mut r3 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        while r2 == r1 {
          r2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
            pop_size
        }
        while r3 == r1 || r3 == r2 {
          r3 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
            pop_size
        }
        let fv_f = f_pool[t]
        let cr_v = cr_pool[t]
        let cand : Array[Double] = []
        for d in 0..<dim {
          let v = if t == 1 {
            // current-to-rand/1 (applied per-dimension, binomial-ised)
            pop[i][d] +
            RNG::next(current_rng) * (pop[r1][d] - pop[i][d]) +
            fv_f * (pop[r2][d] - pop[r3][d])
            // rand/1/bin
          } else if RNG::next(current_rng) < cr_v ||
            RNG::next(current_rng) * dim.to_double() < 1.0 {
            pop[r1][d] + fv_f * (pop[r2][d] - pop[r3][d])
          } else {
            pop[i][d]
          }
          cand.push(max(bounds.min, min(bounds.max, v)))
        }
        let fv = evaluate(f, cand)
        evaluations = evaluations + 1
        if fv < trial_best_val {
          trial_best = cand
          trial_best_val = fv
        }
      }
      // greedy selection against the target
      if trial_best_val < vals[i] {
        archive.push(pop[i])
        if archive.length() > archive_cap {
          let _ = archive.remove(0)
        }
        pop[i] = trial_best
        vals[i] = trial_best_val
        if trial_best_val < best_val {
          best_val = trial_best_val
          best_pos = trial_best
        }
      }
    }
    // occasional archive-guided random replacement keeps failed genes alive
    if archive.length() > 0 && evaluations < max_evals {
      let a = (RNG::next(current_rng) * archive.length().to_double()).to_int() %
        archive.length()
      let w = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let donor = archive[a]
      let r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let cand : Array[Double] = []
      for d in 0..<dim {
        cand.push(
          max(
            bounds.min,
            min(
              bounds.max,
              donor[d] +
              RNG::next(current_rng) * (best_pos[d] - donor[d]) +
              1.0 * (pop[r1][d] - donor[d]),
            ),
          ),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[w] {
        archive[w % archive.length()] = pop[w]
        pop[w] = cand
        vals[w] = fv
        if fv < best_val {
          best_val = fv
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

// ===== ODE: Opposition-based Differential Evolution (2006) =====

///|
pub struct ODEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ode_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ODEConfig {
  ODEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Opposition-based Differential Evolution (Rahnamayan et al. 2006): the
/// classical DE/rand/1/bin plus opposition-based learning - with jump
/// rate JR each generation a random quarter of the population is
/// reflected through the bound midpoint (x' = lb + ub - x) and the better
/// of the pair is kept, restarting exploration from the opposite side of
/// the box.  Historical-best tracking is added.
pub fn opposition_de(
  config : ODEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let jr = 0.3
  let f_factor = 0.7
  let cr = 0.9

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // opposition jump: reflect a random quarter of the population
    if RNG::next(current_rng) < jr {
      for i in 0..<pop_size {
        if evaluations >= max_evals {
          break
        }
        if RNG::next(current_rng) < 0.25 {
          let opp : Array[Double] = []
          for d in 0..<dim {
            opp.push(bounds.min + bounds.max - pop[i][d])
          }
          let fv = evaluate(f, opp)
          evaluations = evaluations + 1
          if fv < vals[i] {
            pop[i] = opp
            vals[i] = fv
            if fv < best_val {
              best_val = fv
              best_pos = opp
            }
          }
        }
      }
    }
    // DE/rand/1/bin generation
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let mut r2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let mut r3 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      while r2 == r1 {
        r2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      while r3 == r1 || r3 == r2 {
        r3 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
      }
      let mut_trial : Array[Double] = []
      for d in 0..<dim {
        mut_trial.push(pop[r1][d] + f_factor * (pop[r2][d] - pop[r3][d]))
      }
      let trial : Array[Double] = []
      for d in 0..<dim {
        let use_mut = RNG::next(current_rng) < cr ||
          (RNG::next(current_rng) * dim.to_double()).to_int() % dim == d
        trial.push(
          if use_mut {
            max(bounds.min, min(bounds.max, mut_trial[d]))
          } else {
            pop[i][d]
          },
        )
      }
      let fv = evaluate(f, trial)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = trial
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== CBO: Colliding Bodies Optimization (2014) =====

///|
pub struct CBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CBOConfig {
  CBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Colliding Bodies Optimization (Kaveh & Mahdavi 2014): bodies are ranked
/// by fitness into a light fast half (CB2) and a heavy slow half (CB1),
/// paired best-to-worst.  The heavy half (CB1) keeps one-dimensional
/// collision momentum exchange (restitution alpha=1.3) as its exploration
/// with a linearly decaying displacement; the light half (CB2) refines
/// around the incumbent best with a fading Gaussian, which fixed the
/// first version's stall (Sphere 5k 0.66).  Historical-best tracking is
/// added.
pub fn colliding_bodies_optimization(
  config : CBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let half = pop_size / 2
  let alpha = 1.3

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let vel : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  let v0 : Array[Double] = []
  for _d in 0..<dim {
    v0.push(0.0)
  }
  vel.push(v0)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    let vi : Array[Double] = []
    for _d in 0..<dim {
      vi.push((RNG::next(current_rng) - 0.5) * 0.2)
    }
    vel.push(vi)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // rank indices ascending
    let idx : Array[Int32] = []
    for i in 0..<pop_size {
      idx.push(i)
    }
    for a in 0..<pop_size {
      let mut m = a
      for b2 in (a + 1)..<pop_size {
        if vals[idx[b2]] < vals[idx[m]] {
          m = b2
        }
      }
      if m != a {
        let ti = idx[a]
        idx[a] = idx[m]
        idx[m] = ti
      }
    }
    // collision pairs: best half (CB1) keeps momentum-exchange exploration
    // with an enlarged displacement; worst half (CB2) switches to a fading
    // Gaussian refinement around the best (first version's tiny collision
    // displacements stalled Sphere 5k at 0.66)
    for j in 0..<half {
      if evaluations >= max_evals {
        break
      }
      let i1 = idx[j]
      let i2 = idx[half + j]
      let m1 = (half - j).to_double()
      let m2 = (j + 1).to_double()
      let total = m1 + m2
      // CB1: collide with its paired CB2, momentum exchange
      {
        let ii = i1
        let i = i2
        let newv : Array[Double] = []
        let newx : Array[Double] = []
        for d in 0..<dim {
          let vs = vel[ii][d]
          let vo = vel[i][d]
          let vn = (m2 * vo + (m1 - alpha * m2) * vs) / total
          newv.push(vn)
          newx.push(
            max(bounds.min, min(bounds.max, pop[ii][d] + vn * rb * range * 0.2)),
          )
        }
        let fv = evaluate(f, newx)
        evaluations = evaluations + 1
        if fv < vals[ii] {
          vel[ii] = newv
          pop[ii] = newx
          vals[ii] = fv
          if fv < best_val {
            best_val = fv
            best_pos = newx
          }
        }
      }
      // CB2: fading Gaussian refinement around the best
      if evaluations < max_evals {
        let ii = i2
        let newx : Array[Double] = []
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          newx.push(
            max(bounds.min, min(bounds.max, best_pos[d] + range * 0.1 * rb * g)),
          )
        }
        let fv = evaluate(f, newx)
        evaluations = evaluations + 1
        if fv < vals[ii] {
          pop[ii] = newx
          vals[ii] = fv
          if fv < best_val {
            best_val = fv
            best_pos = newx
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

// ===== TWO: Teamwork Optimization Algorithm (2021) =====

///|
pub struct TWOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_two_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TWOConfig {
  TWOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Teamwork Optimization Algorithm (Dehghani et al. 2021): the population
/// is split into teams, each with its own leader (team best).  Per
/// generation every member makes two coordinated moves - first toward its
/// team leader, then toward the global leader - with per-dimension random
/// pulls and greedy acceptance; leaders themselves stay put.  A fading
/// Gaussian jitter on non-leaders keeps exploration alive late.
/// Historical-best tracking is added.
pub fn teamwork_optimization(
  config : TWOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 8 {
    8
  } else {
    config.population_size
  }
  let n_teams = 4

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // team leaders: best member of each team
    let leaders : Array[Int32] = []
    for t in 0..<n_teams {
      let mut li = t
      let mut lv = vals[t]
      let mut i = t + n_teams
      while i < pop_size {
        if vals[i] < lv {
          lv = vals[i]
          li = i
        }
        i = i + n_teams
      }
      leaders.push(li)
    }
    // round 1: toward team leader
    for t in 0..<n_teams {
      let mut i = t
      while i < pop_size {
        if evaluations >= max_evals {
          break
        }
        if i != leaders[t] {
          let cand : Array[Double] = []
          for d in 0..<dim {
            let r = RNG::next(current_rng)
            cand.push(
              max(
                bounds.min,
                min(
                  bounds.max,
                  pop[i][d] + r * (pop[leaders[t]][d] - pop[i][d]),
                ),
              ),
            )
          }
          let fv = evaluate(f, cand)
          evaluations = evaluations + 1
          if fv < vals[i] {
            pop[i] = cand
            vals[i] = fv
            if fv < best_val {
              best_val = fv
              best_pos = cand
            }
          }
        }
        i = i + n_teams
      }
    }
    // round 2: toward the global leader + fading jitter
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let is_leader = leaders.contains(i)
      if is_leader {
        continue
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(
          max(
            bounds.min,
            min(
              bounds.max,
              pop[i][d] + r * (best_pos[d] - pop[i][d]) + range * 0.02 * rb * g,
            ),
          ),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== SPO: Spiral Optimization Algorithm (2011) =====

///|
pub struct SPOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_spo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SPOConfig {
  SPOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Spiral Optimization Algorithm (Tamura & Yasuda 2011), n-dimensional
/// model: every star spirals toward the current best via a contraction
/// factor r = 0.95*(1-iter) plus a random Givens plane rotation on two
/// picked dimensions (the n-d analogue of the 2-D spiral matrix
/// S_n(r,theta)); the best star itself is updated as the rotation centre
/// with a small fading jitter.  Greedy acceptance and historical-best
/// tracking are added.
pub fn spiral_optimization(
  config : SPOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let r_max = 0.95

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    let r = r_max * rb
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let theta = RNG::next(current_rng) * 6.283185307179586
      let ct = @math.cos(theta)
      let st = @math.sin(theta)
      // two random dimensions form the rotation plane
      let pd = (RNG::next(current_rng) * dim.to_double()).to_int() % dim
      let mut qd = (RNG::next(current_rng) * dim.to_double()).to_int() % dim
      if dim > 1 {
        while qd == pd {
          qd = (RNG::next(current_rng) * dim.to_double()).to_int() % dim
        }
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        let rel = pop[i][d] - best_pos[d]
        let v = if d == pd && dim > 1 {
          r * (rel * ct - (pop[i][qd] - best_pos[qd]) * st)
        } else if d == qd && dim > 1 {
          r * (rel * st + (pop[i][pd] - best_pos[pd]) * ct)
        } else {
          r * rel
        }
        cand.push(max(bounds.min, min(bounds.max, best_pos[d] + v)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
      }
    }
    // rotation centre jitter: keep the best star drifting slightly early on
    if evaluations < max_evals && rb > 0.02 {
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(
          max(
            bounds.min,
            min(bounds.max, best_pos[d] + range * 0.05 * rb * rb * g),
          ),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RTH: Red-tailed Hawk Algorithm (2022) =====

///|
pub struct RTHConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rth_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RTHConfig {
  RTHConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Red-tailed Hawk Algorithm (Abdechiri et al. 2022), mechanism-faithful
/// variant: the hunt progresses through four stages as the budget
/// advances - high soar (differential exploration against a random peer),
/// high scan (broad sampling around the best), swoop (best-attraction
/// plus peer difference) and strike (fading Gaussian refinement at the
/// kill site).  Greedy acceptance, historical-best tracking and
/// per-dimension sampling follow the house rules.
pub fn red_tailed_hawk(
  config : RTHConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let r = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let cand : Array[Double] = []
      if rb > 0.7 {
        // high soar: differential exploration
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)
          let diff = pop[r][d] - 2.0 * r2 * pop[i][d]
          cand.push(
            max(bounds.min, min(bounds.max, pop[r][d] - r1 * diff.abs())),
          )
        }
      } else if rb > 0.55 {
        // high scan: interpolation toward the best (2x spread burned a
        // quarter of the budget on multimodal landscapes: Ackley 5k regressed;
        // the window is narrower now)
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(bounds.max, best_pos[d] + r1 * (pop[i][d] - best_pos[d])),
            ),
          )
        }
      } else if rb > 0.35 {
        // swoop: best attraction plus peer difference
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                best_pos[d] +
                r1 * (best_pos[d] - pop[i][d]) +
                0.5 * r2 * (pop[r][d] - pop[i][d]),
              ),
            ),
          )
        }
      } else {
        // strike: fading Gaussian refinement at the kill site; the sigma
        // floor keeps late-stage escapes alive (d5 Levy regressed 2.6x
        // when sigma decayed to ~0 and the best got trapped in a basin)
        let sigma_strike = if range * 0.12 * rb > range * 0.01 {
          range * 0.12 * rb
        } else {
          range * 0.01
        }
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(bounds.min, min(bounds.max, best_pos[d] + sigma_strike * g)),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== CRS: Controlled Random Search (1973) =====

///|
pub struct CRSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_crs_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CRSConfig {
  CRSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Controlled Random Search (Price 1973), CRS2 variant: each round samples
/// n+1 population members, reflects the global worst through the centroid
/// of those sampled points (x_new = 2G - x_worst), and keeps the reflection
/// only if it beats the worst (population size stays fixed).  A fading
/// Gaussian refinement around the best is blended in for the second half
/// of the budget (the pure reflection version stalls at the level of the
/// reflected polytope, Sphere 5k ~ 0.1); randomness is confined to
/// sampling so trajectories stay budget-monotone.
pub fn controlled_random_search(
  config : CRSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < dim + 2 {
    dim + 2
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }
  let mut worst_i = 0
  let mut worst_v = vals[0]
  for i in 1..<pop_size {
    if vals[i] > worst_v {
      worst_v = vals[i]
      worst_i = i
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // reflection phase: 2G - x_worst with G = centroid of n random points
    if evaluations < max_evals {
      let g : Array[Double] = []
      for _d in 0..<dim {
        g.push(0.0)
      }
      let n_pts = dim + 1
      let mut used_worst = false
      for _s in 0..<n_pts {
        let i = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        if i == worst_i {
          used_worst = true
        }
        for d in 0..<dim {
          g[d] = g[d] + pop[i][d] / n_pts.to_double()
        }
      }
      // exclude the worst from G by re-weighting when sampled
      if used_worst {
        for d in 0..<dim {
          g[d] = g[d] * n_pts.to_double()
          g[d] = g[d] - pop[worst_i][d]
          g[d] = g[d] / dim.to_double()
        }
      }
      let cand : Array[Double] = []
      for d in 0..<dim {
        cand.push(
          max(bounds.min, min(bounds.max, 2.0 * g[d] - pop[worst_i][d])),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < worst_v {
        pop[worst_i] = cand
        vals[worst_i] = fv
        worst_v = fv
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
        // refresh the global worst index
        worst_i = 0
        worst_v = vals[0]
        for i in 1..<pop_size {
          if vals[i] > worst_v {
            worst_v = vals[i]
            worst_i = i
          }
        }
      }
    }
    // fading refinement around the best (second half of the budget)
    if evaluations < max_evals && rb < 0.5 {
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g2 = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(
          max(bounds.min, min(bounds.max, best_pos[d] + range * 0.1 * rb * g2)),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
        // also displace the worst slot so the population follows
        if best_val < worst_v {
          pop[worst_i] = cand
          vals[worst_i] = fv
          worst_v = fv
          worst_i = 0
          worst_v = vals[0]
          for i in 1..<pop_size {
            if vals[i] > worst_v {
              worst_v = vals[i]
              worst_i = i
            }
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

// ===== UPSO: Unified Particle Swarm Optimization (2004) =====

///|
pub struct UPSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_upso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> UPSOConfig {
  UPSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Unified Particle Swarm Optimization (Parsopoulos & Vrahatis 2004): a
/// single velocity update blends the global-best and local-best pulls with
/// a unification factor u - v = chi*(w*v + c1*r1*(pbest-x) + (1-u)*c2*r2*
/// (gbest-x) + u*c3*r3*(lbest-x)) over a ring neighbourhood.  Inertia
/// decays quadratically and velocities are clamped per the house CPSO
/// rules; one random particle per generation is resampled around its
/// personal best with a fading Gaussian to prevent late freeze.
/// Historical-best tracking is added.
pub fn unified_pso(
  config : UPSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let chi = 0.7298
  let c1 = 1.4962
  let c2 = 1.4962
  let c3 = 1.4962
  let u_mix = 0.3
  let v_max = range * 0.2

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let pbest : Array[Array[Double]] = []
  let pbest_vals : Array[Double] = []
  let vel : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  pbest.push(best_pos)
  pbest_vals.push(best_val)
  let v0 : Array[Double] = []
  for _d in 0..<dim {
    v0.push(0.0)
  }
  vel.push(v0)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    pbest.push(p)
    pbest_vals.push(v)
    let vi : Array[Double] = []
    for _d in 0..<dim {
      vi.push((RNG::next(current_rng) - 0.5) * v_max)
    }
    vel.push(vi)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut gen = 0
  while evaluations < max_evals {
    gen = gen + 1
    // budget-independent inertia schedule: geometric per-generation decay
    // clamped at 0.3 (the evals/max-coupled quadratic schedule made the
    // 5k trajectory search far longer at high inertia than the 1k one,
    // flipping Rastrigin d5 both in mean and median at 100 seeds)
    let w0 = 0.7298 * @math.pow(0.98, gen.to_double())
    let w = if w0 < 0.3 { 0.3 } else { w0 }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      // ring neighbourhood best (left, self, right)
      let l1 = if i == 0 { pop_size - 1 } else { i - 1 }
      let l2 = if i == pop_size - 1 { 0 } else { i + 1 }
      let mut li = i
      let mut lv = vals[i]
      for nb in [l1, l2] {
        if vals[nb] < lv {
          lv = vals[nb]
          li = nb
        }
      }
      let xi = pop[i]
      let cand : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let nv = chi *
          (
            w * vel[i][d] +
            c1 * r1 * (pbest[i][d] - xi[d]) +
            (1.0 - u_mix) * c2 * r2 * (best_pos[d] - xi[d]) +
            u_mix * c3 * r3 * (pop[li][d] - xi[d])
          )
        let nv_clamped = if nv > v_max {
          v_max
        } else if nv < -v_max {
          -v_max
        } else {
          nv
        }
        vel[i][d] = nv_clamped
        cand.push(max(bounds.min, min(bounds.max, xi[d] + nv_clamped)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
      }
      if fv < pbest_vals[i] {
        pbest[i] = cand
        pbest_vals[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    // fading resample of one random particle around its personal best:
    // once velocities clamp down the swarm can freeze on Ackley-class
    // landscapes (first version: Ackley 1k~5k at 0.141)
    if evaluations < max_evals {
      let wi = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        let sigma_raw = range * 0.05 * @math.pow(0.98, gen.to_double())
        let sigma_rs = if sigma_raw > range * 0.02 {
          sigma_raw
        } else {
          range * 0.02
        }
        cand.push(max(bounds.min, min(bounds.max, pbest[wi][d] + sigma_rs * g)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[wi] {
        pop[wi] = cand
        vals[wi] = fv
      }
      if fv < pbest_vals[wi] {
        pbest[wi] = cand
        pbest_vals[wi] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RGA: Real-coded Genetic Algorithm with SBX (1995) =====

///|
pub struct RGAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rga_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RGAConfig {
  RGAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Real-coded Genetic Algorithm (Deb & Agrawal 1995): binary tournament
/// selection, simulated binary crossover (eta_c = 15) and polynomial
/// mutation (eta_m = 20, per-gene probability 1/dim) build an offspring
/// population that competes with the parents in a (mu+lambda) truncation.
/// All operators sample per dimension per the shared-random-number rule;
/// historical-best tracking is added.
pub fn real_coded_ga(
  config : RGAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let eta_c = 15.0
  let eta_m = 20.0

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let tour_select = fn() {
    let a = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
    let b2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
    if vals[a] <= vals[b2] {
      a
    } else {
      b2
    }
  }

  while evaluations < max_evals {
    // build the offspring population pairwise
    let off : Array[Array[Double]] = []
    let off_vals : Array[Double] = []
    let n_pairs = pop_size / 2
    for _p in 0..<n_pairs {
      let p1 = tour_select()
      let mut p2 = tour_select()
      while p2 == p1 {
        p2 = tour_select()
      }
      let c1 : Array[Double] = []
      let c2 : Array[Double] = []
      for d in 0..<dim {
        // SBX per dimension
        let u = RNG::next(current_rng)
        let beta = if u < 0.5 {
          @math.pow(2.0 * u, 1.0 / (eta_c + 1.0))
        } else {
          @math.pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_c + 1.0))
        }
        let m1 = 0.5 * ((1.0 + beta) * pop[p1][d] + (1.0 - beta) * pop[p2][d])
        let m2 = 0.5 * ((1.0 - beta) * pop[p1][d] + (1.0 + beta) * pop[p2][d])
        c1.push(max(bounds.min, min(bounds.max, m1)))
        c2.push(max(bounds.min, min(bounds.max, m2)))
      }
      // polynomial mutation per gene
      for child in [c1, c2] {
        for d in 0..<dim {
          if RNG::next(current_rng) < 1.0 / dim.to_double() {
            let u = RNG::next(current_rng)
            let delta = if u < 0.5 {
              @math.pow(2.0 * u, 1.0 / (eta_m + 1.0)) - 1.0
            } else {
              1.0 - @math.pow(2.0 * (1.0 - u), 1.0 / (eta_m + 1.0))
            }
            child[d] = max(
              bounds.min,
              min(bounds.max, child[d] + delta * range),
            )
          }
        }
        let fv = evaluate(f, child)
        evaluations = evaluations + 1
        off.push(child)
        off_vals.push(fv)
        if fv < best_val {
          best_val = fv
          best_pos = child
        }
      }
    }
    // (mu+lambda) truncation over the combined 2*pop_size pool
    let pool_p : Array[Array[Double]] = []
    let pool_v : Array[Double] = []
    for i in 0..<pop_size {
      pool_p.push(pop[i])
      pool_v.push(vals[i])
    }
    for i in 0..<off.length() {
      pool_p.push(off[i])
      pool_v.push(off_vals[i])
    }
    let n_pool = pool_p.length()
    for a in 0..<pop_size {
      let mut m = a
      for b2 in (a + 1)..<n_pool {
        if pool_v[b2] < pool_v[m] {
          m = b2
        }
      }
      if m != a {
        let tp = pool_p[a]
        let tv = pool_v[a]
        pool_p[a] = pool_p[m]
        pool_v[a] = pool_v[m]
        pool_p[m] = tp
        pool_v[m] = tv
      }
    }
    for i in 0..<pop_size {
      pop[i] = pool_p[i]
      vals[i] = pool_v[i]
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== InS: Interior Search (2015) =====

///|
pub struct InSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ins_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> InSConfig {
  InSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Interior Search (Gandomi 2015): a marked fraction (10%) of the
/// population performs mirror reflection about the incumbent best
/// (x' = 2*best - x, stepping to the far side of the current optimum);
/// the rest recombine toward the best with a per-dimension random pull
/// plus a fading Gaussian jitter.  Both moves are greedy-accepted and
/// historical-best tracking is added.
pub fn interior_search(
  config : InSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let mirror_ratio = 0.1

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let cand : Array[Double] = []
      if RNG::next(current_rng) < mirror_ratio {
        // global mirror reflection about the best
        for d in 0..<dim {
          cand.push(max(bounds.min, min(bounds.max, 2.0 * best_pos[d] - xi[d])))
        }
      } else {
        // local recombination toward the best with fading jitter
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                xi[d] + r * (best_pos[d] - xi[d]) + range * 0.05 * rb * g,
              ),
            ),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== BMO: Barnacles Mating Optimizer (2018) =====

///|
pub struct BMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bmo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BMOConfig {
  BMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Barnacles Mating Optimizer (Sulaiman et al. 2018): each barnacle
/// either mates softly with two random parents (per-dimension weighted
/// blend) or, with complementary probability, mates hard with the
/// incumbent best (rand*best + (1-rand)*self).  A fading Gaussian jitter
/// keeps recombination from collapsing onto the best; moves are
/// greedy-accepted and historical-best tracking is added.
pub fn barnacles_mating_optimizer(
  config : BMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let cand : Array[Double] = []
      if RNG::next(current_rng) < 0.5 {
        // soft mating: blend of two random parents
        let p1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let mut p2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        while p2 == p1 {
          p2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
            pop_size
        }
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(bounds.max, r * pop[p1][d] + (1.0 - r) * pop[p2][d]),
            ),
          )
        }
      } else {
        // hard mating: blend with the incumbent best
        for d in 0..<dim {
          let r = RNG::next(current_rng)
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                r * best_pos[d] + (1.0 - r) * xi[d] + range * 0.06 * rb * g,
              ),
            ),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== NRBO: Newton-Raphson Based Optimization (2024) =====

///|
pub struct NRBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_nrbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NRBOConfig {
  NRBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Newton-Raphson Based Optimization (Sowmya et al. 2024),
/// mechanism-faithful variant alternating two update rules: the
/// Newton-Raphson search rule (a three-term differential combining a
/// random peer, the best and the current point) and the training rule
/// (a best/mean-centred pull whose magnitude fades with progress).
/// Moves are greedy-accepted with per-dimension sampling and
/// historical-best tracking.
pub fn newton_raphson_based_optimization(
  config : NRBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // population mean for the training rule
    let mean : Array[Double] = []
    for d in 0..<dim {
      let mut m = 0.0
      for i in 0..<pop_size {
        m = m + pop[i][d]
      }
      mean.push(m / pop_size.to_double())
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let r = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let xi = pop[i]
      let cand : Array[Double] = []
      if RNG::next(current_rng) < 0.5 {
        // Newton-Raphson search rule: three-term differential
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                pop[r][d] +
                r1 * (best_pos[d] - pop[r][d]) +
                r2 * (xi[d] - pop[r][d]),
              ),
            ),
          )
        }
      } else {
        // training rule: best/mean-centred pull with fading magnitude
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                xi[d] - r1 * (mean[d] - xi[d]) - r2 * (best_pos[d] - xi[d]) * rb,
              ),
            ),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== SPSA: Simultaneous Perturbation Stochastic Approximation (1992) =====

///|
pub struct SPSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_spsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SPSAConfig {
  SPSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Simultaneous Perturbation Stochastic Approximation (Spall 1992): a
/// single incumbent estimates a descent direction from exactly two
/// function evaluations per generation - one forward and one backward
/// sample along a random +-1 perturbation vector (all dimensions moved
/// simultaneously).  The perturbation radius c and the step size a decay
/// geometrically per generation with floors (budget-independent schedule,
/// prefix-stable trajectory), and a fading Gaussian jitter guards against
/// late stagnation on multimodal landscapes.  Historical-best tracking is
/// added.
pub fn simultaneous_perturbation_spa(
  config : SPSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  let mut x = random_solution(dim, bounds, current_rng)
  let mut fx = evaluate(f, x)
  let mut evaluations = 1
  let mut best_pos : Array[Double] = x
  let mut best_val = fx

  let mut gen = 0
  while evaluations < max_evals {
    gen = gen + 1
    // Spall's polynomial schedules (much slower than geometric decay:
    // with 0.97^gen the perturbation radius hit its floor long before
    // Rosenbrock's valley was crossed and the estimate drowned in noise)
    let c_k = range * 0.2 / @math.pow(gen.to_double(), 0.101)
    let step_k = range * 0.15 / @math.pow((gen + 10).to_double(), 0.602)
    // simultaneous perturbation: random +-1 per dimension, both probes
    // drawn, then the gradient estimate uses the actual deltas
    let delta : Array[Double] = []
    for _d in 0..<dim {
      delta.push(if RNG::next(current_rng) < 0.5 { -1.0 } else { 1.0 })
    }
    let x_plus : Array[Double] = []
    for d in 0..<dim {
      x_plus.push(max(bounds.min, min(bounds.max, x[d] + c_k * delta[d])))
    }
    let fp = evaluate(f, x_plus)
    evaluations = evaluations + 1
    let x_minus : Array[Double] = []
    for d in 0..<dim {
      x_minus.push(max(bounds.min, min(bounds.max, x[d] - c_k * delta[d])))
    }
    let fm = evaluate(f, x_minus)
    evaluations = evaluations + 1
    // gradient estimate and update; the raw estimate has the function's
    // gradient magnitude (Rosenbrock's valley walls carry |g| ~ 1e4, so
    // a raw a_k*g step overshoots into the opposite corner and oscillates
    // there forever) - normalize to a trust-region step of step_k*range
    // along -g
    let mut gnorm2 = 0.0
    for d in 0..<dim {
      let g = (fp - fm) / (2.0 * c_k * delta[d])
      gnorm2 = gnorm2 + g * g
    }
    let gnorm = gnorm2.sqrt() + 1.0e-300
    for d in 0..<dim {
      let g = (fp - fm) / (2.0 * c_k * delta[d])
      x[d] = max(bounds.min, min(bounds.max, x[d] - step_k * g / gnorm))
    }
    fx = evaluate(f, x)
    evaluations = evaluations + 1
    if fx < best_val {
      best_val = fx
      best_pos = x
    }
    // fading jitter against late stagnation (sigma floored so the
    // incumbent can still leave a local basin at 5k: Ackley was frozen
    // bit-exact at 0.0677)
    if evaluations < max_evals && RNG::next(current_rng) < 0.3 {
      let cand : Array[Double] = []
      let sigma_j0 = range * 0.05 * @math.pow(0.98, gen.to_double())
      let sigma_j = if sigma_j0 > range * 0.005 {
        sigma_j0
      } else {
        range * 0.005
      }
      for d in 0..<dim {
        let g2 = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(max(bounds.min, min(bounds.max, best_pos[d] + sigma_j * g2)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < fx {
        x = cand
        fx = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== DIRECT: Dividing Rectangles (1993) =====

///|
pub struct DIRECTConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_direct_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DIRECTConfig {
  DIRECTConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// DIRECT (Jones, Perttunen & Stuckman 1993), streamlined variant: a set
/// of hyperrectangles is repeatedly split along its longest side.  Each
/// round selects the two best-valued and the two widest rectangles
/// (a tractable stand-in for the convex-hull lower bound of the paper),
/// samples both trisection endpoints of each selected rectangle and
/// replaces it by the three sub-rectangles.  Fully deterministic apart
/// from the random initial centre, so the trajectory is prefix-stable;
/// historical-best tracking is added.
pub fn dividing_rectangles(
  config : DIRECTConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min

  // rectangle = (center, half-widths, value, min-half-width)
  let rects_c : Array[Array[Double]] = []
  let rects_h : Array[Array[Double]] = []
  let rects_v : Array[Double] = []
  let rects_w : Array[Double] = []

  let c0 : Array[Double] = []
  let h0 : Array[Double] = []
  for _d in 0..<dim {
    c0.push((bounds.min + bounds.max) / 2.0)
    h0.push(range / 2.0)
  }
  // random initial centre keeps repeated runs diverse
  let x0 = random_solution(dim, bounds, current_rng)
  for d in 0..<dim {
    c0[d] = x0[d]
  }
  let v0 = evaluate(f, x0)
  let mut evaluations = 1
  rects_c.push(c0)
  rects_h.push(h0)
  rects_v.push(v0)
  rects_w.push(range / 2.0)
  let mut best_pos : Array[Double] = x0
  let mut best_val = v0

  while evaluations < max_evals {
    // pick up to 2 best-valued + 2 widest rectangles (dedup by index)
    let pick : Array[Int32] = []
    let try_add = fn(i : Int32) { if !pick.contains(i) { pick.push(i) } }
    // best-valued / widest picks via two linear scans each (the
    // full selection sort was O(n^2) per round over a rectangle set
    // that grows into the thousands and stalled the scan)
    let n_rects = rects_v.length()
    for _pass in 0..<2 {
      let mut best_i : Int32 = -1
      for i in 0..<n_rects {
        if pick.contains(i) {
          continue
        }
        if best_i == -1 || rects_v[i] < rects_v[best_i] {
          best_i = i
        }
      }
      if best_i != -1 {
        try_add(best_i)
      }
    }
    for _pass in 0..<2 {
      let mut wide_i : Int32 = -1
      for i in 0..<n_rects {
        if pick.contains(i) {
          continue
        }
        if wide_i == -1 || rects_w[i] > rects_w[wide_i] {
          wide_i = i
        }
      }
      if wide_i != -1 {
        try_add(wide_i)
      }
    }
    // split each selected rectangle
    for sel in 0..<pick.length() {
      if evaluations >= max_evals {
        break
      }
      let i = pick[sel]
      let c = rects_c[i]
      let h = rects_h[i]
      // longest dimension(s)
      let mut wmax = h[0]
      for d in 1..<dim {
        if h[d] > wmax {
          wmax = h[d]
        }
      }
      let ndiv = range * 1.0e-12
      // sample the two trisection endpoints along every longest dimension
      // (only the first longest dimension is split per the paper's rule;
      // ties are broken by the lowest index for determinism)
      let mut m = 0
      for d in 1..<dim {
        if h[d] > h[m] + ndiv {
          m = d
        }
      }
      let p1 : Array[Double] = []
      let p2 : Array[Double] = []
      for d in 0..<dim {
        p1.push(c[d])
        p2.push(c[d])
      }
      p1[m] = c[m] - h[m]
      p2[m] = c[m] + h[m]
      let v1 = evaluate(f, p1)
      evaluations = evaluations + 1
      let v2 = if evaluations < max_evals {
        let vv = evaluate(f, p2)
        evaluations = evaluations + 1
        vv
      } else {
        v1
      }
      if v1 < best_val {
        best_val = v1
        best_pos = p1
      }
      if v2 < best_val {
        best_val = v2
        best_pos = p2
      }
      // trisect: three sub-rectangles replace the parent
      let hn : Array[Double] = []
      for d in 0..<dim {
        hn.push(if d == m { h[d] / 3.0 } else { h[d] })
      }
      // middle rectangle keeps the centre
      rects_c.push(c)
      rects_h.push(hn)
      rects_v.push(if v1 < v2 { v1 } else { v2 })
      let mut wn = hn[0]
      for d in 1..<dim {
        if hn[d] < wn {
          wn = hn[d]
        }
      }
      rects_w.push(wn)
      // left rectangle
      let cl : Array[Double] = []
      for d in 0..<dim {
        cl.push(if d == m { c[d] - 2.0 * h[d] / 3.0 } else { c[d] })
      }
      rects_c.push(cl)
      rects_h.push(hn)
      rects_v.push(v1)
      rects_w.push(wn)
      // right rectangle
      let cr : Array[Double] = []
      for d in 0..<dim {
        cr.push(if d == m { c[d] + 2.0 * h[d] / 3.0 } else { c[d] })
      }
      rects_c.push(cr)
      rects_h.push(hn)
      rects_v.push(v2)
      rects_w.push(wn)
      // retire the parent (swap with the tail and pop)
      let last = rects_c.length() - 1
      rects_c[i] = rects_c[last]
      rects_h[i] = rects_h[last]
      rects_v[i] = rects_v[last]
      rects_w[i] = rects_w[last]
      let _ = rects_c.pop()
      let _ = rects_h.pop()
      let _ = rects_v.pop()
      let _ = rects_w.pop()
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Memetic: Memetic Algorithm (1989) =====

///|
pub struct MemeticConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_memetic_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MemeticConfig {
  MemeticConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Memetic Algorithm (Moscato 1989): the real-coded GA body (binary
/// tournament, SBX eta_c=15, polynomial mutation eta_m=20, (mu+lambda)
/// truncation) plus the meme - a Lamarckian local refinement applied to a
/// 20% fraction of offspring: a short shrinking-Gaussian hill climb whose
/// result is written back into the genotype.  Per-dimension sampling and
/// historical-best tracking per house rules.
pub fn memetic_algorithm(
  config : MemeticConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let eta_c = 15.0
  let eta_m = 20.0
  let refine_prob = 0.2
  let refine_steps = 6

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let tour_select = fn() {
    let a = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
    let b2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() % pop_size
    if vals[a] <= vals[b2] {
      a
    } else {
      b2
    }
  }

  while evaluations < max_evals {
    let off : Array[Array[Double]] = []
    let off_vals : Array[Double] = []
    let n_pairs = pop_size / 2
    for _p in 0..<n_pairs {
      let p1 = tour_select()
      let mut p2 = tour_select()
      while p2 == p1 {
        p2 = tour_select()
      }
      let c1 : Array[Double] = []
      let c2 : Array[Double] = []
      for d in 0..<dim {
        let u = RNG::next(current_rng)
        let beta = if u < 0.5 {
          @math.pow(2.0 * u, 1.0 / (eta_c + 1.0))
        } else {
          @math.pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_c + 1.0))
        }
        let m1 = 0.5 * ((1.0 + beta) * pop[p1][d] + (1.0 - beta) * pop[p2][d])
        let m2 = 0.5 * ((1.0 - beta) * pop[p1][d] + (1.0 + beta) * pop[p2][d])
        c1.push(max(bounds.min, min(bounds.max, m1)))
        c2.push(max(bounds.min, min(bounds.max, m2)))
      }
      for child in [c1, c2] {
        for d in 0..<dim {
          if RNG::next(current_rng) < 1.0 / dim.to_double() {
            let u = RNG::next(current_rng)
            let delta = if u < 0.5 {
              @math.pow(2.0 * u, 1.0 / (eta_m + 1.0)) - 1.0
            } else {
              1.0 - @math.pow(2.0 * (1.0 - u), 1.0 / (eta_m + 1.0))
            }
            child[d] = max(
              bounds.min,
              min(bounds.max, child[d] + delta * range),
            )
          }
        }
        let mut fv = evaluate(f, child)
        evaluations = evaluations + 1
        // the meme: Lamarckian shrinking-Gaussian hill climb on a subset
        if RNG::next(current_rng) < refine_prob {
          let mut xr = child
          let mut fr = fv
          let mut sigma = range * 0.05
          for _l in 0..<refine_steps {
            if evaluations >= max_evals {
              break
            }
            let cand : Array[Double] = []
            for d in 0..<dim {
              let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
              cand.push(max(bounds.min, min(bounds.max, xr[d] + sigma * g)))
            }
            let fv2 = evaluate(f, cand)
            evaluations = evaluations + 1
            if fv2 < fr {
              xr = cand
              fr = fv2
            } else {
              sigma = sigma * 0.6
            }
          }
          for d in 0..<dim {
            child[d] = xr[d]
          }
          fv = fr
        }
        off.push(child)
        off_vals.push(fv)
        if fv < best_val {
          best_val = fv
          best_pos = child
        }
      }
    }
    // (mu+lambda) truncation
    let pool_p : Array[Array[Double]] = []
    let pool_v : Array[Double] = []
    for i in 0..<pop_size {
      pool_p.push(pop[i])
      pool_v.push(vals[i])
    }
    for i in 0..<off.length() {
      pool_p.push(off[i])
      pool_v.push(off_vals[i])
    }
    let n_pool = pool_p.length()
    for a in 0..<pop_size {
      let mut m = a
      for b2 in (a + 1)..<n_pool {
        if pool_v[b2] < pool_v[m] {
          m = b2
        }
      }
      if m != a {
        let tp = pool_p[a]
        let tv = pool_v[a]
        pool_p[a] = pool_p[m]
        pool_v[a] = pool_v[m]
        pool_p[m] = tp
        pool_v[m] = tv
      }
    }
    for i in 0..<pop_size {
      pop[i] = pool_p[i]
      vals[i] = pool_v[i]
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== JADE: Adaptive DE with optional external archive (2009) =====

///|
pub struct JADEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_jade_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> JADEConfig {
  JADEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// JADE (Zhang & Sanderson 2009): current-to-pbest/1 mutation - the elite
/// guide is drawn from the top p=5% rather than the single best, and the
/// second difference vector may come from an external archive of recently
/// replaced parents (capacity = population size).  F follows a Cauchy
/// around mu_F, CR a normal around mu_CR; successful values update the
/// centres each generation via the Lehmer mean (F) and arithmetic mean
/// (CR) with c=0.1.  Binary crossover and greedy selection per house
/// rules; historical-best tracking is added.
pub fn jade(
  config : JADEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 6 {
    6
  } else {
    config.population_size
  }
  let p_frac = 0.05
  let c_adapt = 0.1

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }
  // external archive of replaced parents
  let archive : Array[Array[Double]] = []
  let mut mu_f = 0.5
  let mut mu_cr = 0.5

  while evaluations < max_evals {
    let s_f : Array[Double] = []
    let s_cr : Array[Double] = []
    // pbest set size
    let n_pbest = if (pop_size.to_double() * p_frac).to_int() < 1 {
      1
    } else {
      (pop_size.to_double() * p_frac).to_int()
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      // sample F ~ Cauchy(mu_F, 0.1), truncated to (0, 1]
      let mut fv = 0.5
      let mut sampled = false
      while !sampled {
        let u = RNG::next(current_rng)
        let cand = mu_f + 0.1 * @math.tan(3.141592653589793 * (u - 0.5))
        if cand > 0.0 {
          fv = if cand > 1.0 { 1.0 } else { cand }
          sampled = true
        }
      }
      let cg = box_muller(RNG::next(current_rng), RNG::next(current_rng))
      let cr_v = if mu_cr + 0.1 * cg < 0.0 {
        0.0
      } else if mu_cr + 0.1 * cg > 1.0 {
        1.0
      } else {
        mu_cr + 0.1 * cg
      }
      // pbest: a random member of the top n_pbest (value-ranked)
      let idx : Array[Int32] = []
      for j in 0..<pop_size {
        idx.push(j)
      }
      for a in 0..<n_pbest {
        let mut m = a
        for b2 in (a + 1)..<pop_size {
          if vals[idx[b2]] < vals[idx[m]] {
            m = b2
          }
        }
        if m != a {
          let t = idx[a]
          idx[a] = idx[m]
          idx[m] = t
        }
      }
      let pbest_i = idx[(RNG::next(current_rng) * n_pbest.to_double()).to_int() %
        n_pbest]
      let r1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      // r2 from pop union archive
      let n_total = pop_size + archive.length()
      let r2_raw = (RNG::next(current_rng) * n_total.to_double()).to_int() %
        n_total
      let x_r2 = if r2_raw < pop_size {
        pop[r2_raw]
      } else {
        archive[r2_raw - pop_size]
      }
      // current-to-pbest/1 mutation
      let mut_trial : Array[Double] = []
      for d in 0..<dim {
        mut_trial.push(
          pop[i][d] +
          fv * (pop[pbest_i][d] - pop[i][d]) +
          fv * (pop[r1][d] - x_r2[d]),
        )
      }
      // binary crossover
      let trial : Array[Double] = []
      for d in 0..<dim {
        let use_mut = RNG::next(current_rng) < cr_v ||
          (RNG::next(current_rng) * dim.to_double()).to_int() % dim == d
        trial.push(
          if use_mut {
            max(bounds.min, min(bounds.max, mut_trial[d]))
          } else {
            pop[i][d]
          },
        )
      }
      let fv_trial = evaluate(f, trial)
      evaluations = evaluations + 1
      if fv_trial < vals[i] {
        archive.push(pop[i])
        if archive.length() > pop_size {
          let kill = (RNG::next(current_rng) * archive.length().to_double()).to_int() %
            archive.length()
          let _ = archive.remove(kill)
        }
        pop[i] = trial
        vals[i] = fv_trial
        s_f.push(fv)
        s_cr.push(cr_v)
        if fv_trial < best_val {
          best_val = fv_trial
          best_pos = trial
        }
      }
    }
    // Lehmer mean for F, arithmetic for CR
    if s_f.length() > 0 {
      let mut sum_f2 = 0.0
      let mut sum_f = 0.0
      let mut sum_cr = 0.0
      for k in 0..<s_f.length() {
        sum_f2 = sum_f2 + s_f[k] * s_f[k]
        sum_f = sum_f + s_f[k]
        sum_cr = sum_cr + s_cr[k]
      }
      mu_f = (1.0 - c_adapt) * mu_f + c_adapt * (sum_f2 / (sum_f + 1.0e-12))
      mu_cr = (1.0 - c_adapt) * mu_cr +
        c_adapt * (sum_cr / s_f.length().to_double())
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== LCA: League Championship Algorithm (2015) =====

///|
pub struct LCAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_lca_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> LCAConfig {
  LCAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// League Championship Algorithm (Kashan 2015): each round pairs the
/// teams randomly.  The winner plays an extrapolation shot along its own
/// victory direction (x_w' = x_w + U(0,2)*(x_w - x_l)) and the loser
/// imitates the winner's new play (x_l' = x_l + U(0,2)*(x_w' - x_l));
/// both moves are greedy-accepted with per-dimension sampling.  A fading
/// Gaussian jitter keeps late rounds from stalling; historical-best
/// tracking is added.
pub fn league_championship_algorithm(
  config : LCAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 4 {
    4
  } else {
    config.population_size
  }
  let pop_size_even = pop_size - pop_size % 2

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // random pairing for this round (Fisher-Yates on a copy of indices)
    let order : Array[Int32] = []
    for i in 0..<pop_size_even {
      order.push(i)
    }
    let n_order = pop_size_even
    for i in 0..<n_order {
      let j = (RNG::next(current_rng) * n_order.to_double()).to_int() % n_order
      let t = order[i]
      order[i] = order[j]
      order[j] = t
    }
    let n_pairs = pop_size_even / 2
    for pr in 0..<n_pairs {
      if evaluations >= max_evals {
        break
      }
      let a = order[2 * pr]
      let b2 = order[2 * pr + 1]
      let (winner, loser) = if vals[a] <= vals[b2] { (a, b2) } else { (b2, a) }
      // winner's extrapolation shot
      let w_new : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng) * 2.0
        w_new.push(
          max(
            bounds.min,
            min(
              bounds.max,
              pop[winner][d] + r * (pop[winner][d] - pop[loser][d]),
            ),
          ),
        )
      }
      let fv_w = evaluate(f, w_new)
      evaluations = evaluations + 1
      if fv_w < vals[winner] {
        pop[winner] = w_new
        vals[winner] = fv_w
        if fv_w < best_val {
          best_val = fv_w
          best_pos = w_new
        }
      }
      // loser imitates the winner's play
      if evaluations < max_evals {
        let l_new : Array[Double] = []
        for d in 0..<dim {
          let r = RNG::next(current_rng) * 2.0
          l_new.push(
            max(
              bounds.min,
              min(
                bounds.max,
                pop[loser][d] + r * (pop[winner][d] - pop[loser][d]),
              ),
            ),
          )
        }
        let fv_l = evaluate(f, l_new)
        evaluations = evaluations + 1
        if fv_l < vals[loser] {
          pop[loser] = l_new
          vals[loser] = fv_l
          if fv_l < best_val {
            best_val = fv_l
            best_pos = l_new
          }
        }
      }
    }
    // fading jitter on one random team
    if evaluations < max_evals {
      let w = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let cand : Array[Double] = []
      for d in 0..<dim {
        let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
        cand.push(
          max(bounds.min, min(bounds.max, pop[w][d] + range * 0.05 * rb * g)),
        )
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[w] {
        pop[w] = cand
        vals[w] = fv
        if fv < best_val {
          best_val = fv
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

// ===== HGSO: Henry Gas Solubility Optimization (2019) =====

///|
pub struct HGSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hgso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HGSOConfig {
  HGSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Henry Gas Solubility Optimization (Hashim et al. 2019),
/// mechanism-faithful variant: individuals are grouped into gas types,
/// each carrying a Henry constant H_j that decays with progress (Henry's
/// law as the temperature rises) and sets the solubility-scaled step.  A
/// random individual dissolves toward a peer (explore) or toward its
/// cluster's best (exploit); occasionally the state flips and the
/// individual re-locates around the cluster best with a fading Gaussian.
/// Greedy acceptance and historical-best tracking are added.
pub fn henry_gas_solubility_optimization(
  config : HGSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let range = bounds.max - bounds.min
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let n_clusters = 5
  let per_cluster = pop_size / n_clusters

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let cluster : Array[Int32] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  cluster.push(0)
  for i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    cluster.push(
      (i.to_double() / per_cluster.to_double()).to_int() % n_clusters,
    )
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }
  // Henry constants per gas type (canonical table range [1, 298])
  let henry : Array[Double] = []
  for _j in 0..<n_clusters {
    henry.push(1.0 + RNG::next(current_rng) * 297.0)
  }

  while evaluations < max_evals {
    let iter = evaluations.to_double() / max_evals.to_double()
    let rb = 1.0 - iter
    // Henry's law: constants decay as the temperature rises
    for j in 0..<n_clusters {
      henry[j] = henry[j] * @math.exp(-rb * 0.01)
      if henry[j] < 0.01 {
        henry[j] = 0.01
      }
    }
    // cluster bests
    let c_best_i : Array[Int32] = []
    let c_best_v : Array[Double] = []
    for _j in 0..<n_clusters {
      c_best_i.push(0)
      c_best_v.push(1.0e308)
    }
    for i in 0..<pop_size {
      let j = cluster[i]
      if vals[i] < c_best_v[j] {
        c_best_v[j] = vals[i]
        c_best_i[j] = i
      }
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let j = cluster[i]
      let solubility = henry[j] / 298.0
      let xi = pop[i]
      let cand : Array[Double] = []
      let roll = RNG::next(current_rng)
      if roll < 0.1 {
        // state flip: re-locate around the cluster best
        for d in 0..<dim {
          let g = box_muller(RNG::next(current_rng), RNG::next(current_rng))
          cand.push(
            max(
              bounds.min,
              min(bounds.max, pop[c_best_i[j]][d] + range * 0.1 * rb * g),
            ),
          )
        }
      } else if roll < 0.5 {
        // explore: dissolve toward a random peer scaled by solubility
        let r = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                xi[d] + solubility * r1 * 3.0 * (pop[r][d] - xi[d]),
              ),
            ),
          )
        }
      } else {
        // exploit: converge toward the cluster best scaled by solubility
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          cand.push(
            max(
              bounds.min,
              min(
                bounds.max,
                xi[d] + solubility * r1 * 3.0 * (pop[c_best_i[j]][d] - xi[d]),
              ),
            ),
          )
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== COOT: Coot Optimization Algorithm (2021) =====

///|
pub struct COOTConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_coot_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> COOTConfig {
  COOTConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
fn coot_chosen(chosen : Array[Int32], idx : Int32) -> Bool {
  let mut found = false
  for c in 0..<chosen.length() {
    if chosen[c] == idx {
      found = true
    }
  }
  found
}

///|
fn coot_select_leaders(vals : Array[Double], k : Int32) -> Array[Int32] {
  let n = vals.length()
  let chosen : Array[Int32] = []
  for _c in 0..<k {
    let mut best_i : Int32 = -1
    let mut best_v = 1.0e308
    for i in 0..<n {
      if coot_chosen(chosen, i) == false && vals[i] < best_v {
        best_v = vals[i]
        best_i = i
      }
    }
    if best_i >= 0 {
      chosen.push(best_i)
    }
  }
  chosen
}

///|
/// Coot Optimization Algorithm (Naruei & Keynia 2021), mechanism-faithful
/// variant blending the paper's three movements over a population whose best
/// ~10% act as leaders:
///   * random movement - a coot steps from a random peer j using
///     A = (X_i - X_j) - X_j,  X_i <- X_i + 2 r2 A -/+ (X_i - X_j);
///   * chain movement  - a coot moves to the midpoint of itself and the
///     previous coot,  X_i <- (X_i + X_{i-1}) / 2;
///   * leader movement - each leader follows the global best with a
///     cosine-scaled step, B r3 cos(2 pi r4)(gbest - X) +/- gbest, where
///     B = 2 - 2 t/Tmax decays with progress.
/// Greedy acceptance and historical-best tracking are added.
pub fn coot_optimization_algorithm(
  config : COOTConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let n_leaders = max_int(1, pop_size / 10)

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    let progress = evaluations.to_double() / max_evals.to_double()
    let b_factor = 2.0 - 2.0 * progress
    let leaders = coot_select_leaders(vals, n_leaders)
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let cand : Array[Double] = []
      if coot_chosen(leaders, i) {
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        let r5 = RNG::next(current_rng)
        let cos_term = @math.cos(TWO_PI * r4)
        for d in 0..<dim {
          let step = b_factor * r3 * cos_term * (best_pos[d] - xi[d])
          let nv = if r5 < 0.5 {
            step + best_pos[d]
          } else {
            step - best_pos[d]
          }
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if RNG::next(current_rng) < 0.5 {
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let r2 = RNG::next(current_rng)
        let sign_flip = RNG::next(current_rng) < 0.5
        let xj = pop[j]
        for d in 0..<dim {
          let bvec = xi[d] - xj[d]
          let a = bvec - xj[d]
          let base = xi[d] + 2.0 * r2 * a
          let nv = if sign_flip { base - bvec } else { base + bvec }
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        let prev = pop[if i == 0 { pop_size - 1 } else { i - 1 }]
        for d in 0..<dim {
          cand.push(max(bounds.min, min(bounds.max, (xi[d] + prev[d]) / 2.0)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== DOA: Dingo Optimization Algorithm (2021) =====

///|
pub struct DOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_doa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DOAConfig {
  DOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
fn dingo_distance(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for d in 0..<a.length() {
    let diff = a[d] - b[d]
    s = s + diff * diff
  }
  s.sqrt()
}

///|
/// Dingo Optimization Algorithm (Peraza-Vazquez et al. 2021),
/// mechanism-faithful variant.  Each dingo picks one of three hunting
/// strategies by a random draw (0.5 / 0.25 / 0.25):
///   * group attack - consensus move toward the rest of the pack,
///     X_i <- X_i + beta1 sum_k phi_k (X_k - X_i) / (||X_k - X_i|| + eps);
///   * persecution  - chase the best,
///     X_i <- X_i + beta1 e^{-beta2} (gbest - X_i);
///   * scavenging   - reposition,
///     X_i <- 0.5 (e^{-beta2 gbest} - (-1)^lambda X_i), lambda in {0,1}.
/// beta1 ~ U(-2,2), beta2 ~ U(0,2); greedy acceptance and historical-best
/// tracking are added.
pub fn dingo_optimization_algorithm(
  config : DOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let cand : Array[Double] = []
      let strategy = RNG::next(current_rng)
      if strategy < 0.5 {
        // group attack: consensus of the pack
        let beta1 = RNG::next(current_rng) * 4.0 - 2.0
        let acc : Array[Double] = []
        for _d in 0..<dim {
          acc.push(0.0)
        }
        for k in 0..<pop_size {
          if k != i {
            let dist = dingo_distance(xi, pop[k]) + 1.0e-12
            let phi = RNG::next(current_rng)
            for d in 0..<dim {
              acc[d] = acc[d] + phi * (pop[k][d] - xi[d]) / dist
            }
          }
        }
        for d in 0..<dim {
          let nv = xi[d] + beta1 * acc[d] / pop_size.to_double()
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if strategy < 0.75 {
        // persecution: chase the global best
        let beta1 = RNG::next(current_rng) * 4.0 - 2.0
        let beta2 = RNG::next(current_rng) * 2.0
        let decay = @math.exp(-beta2)
        for d in 0..<dim {
          let nv = xi[d] + beta1 * decay * (best_pos[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // scavenging: reposition
        let beta2 = RNG::next(current_rng) * 2.0
        let lambda = RNG::next(current_rng) < 0.5
        for d in 0..<dim {
          let pos_term = @math.exp(-beta2 * best_pos[d])
          let sign = if lambda { -1.0 } else { 1.0 }
          let nv = 0.5 * (pos_term - sign * xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== STOA: Sooty Tern Optimization Algorithm (2019) =====

///|
pub struct STOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_stoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> STOAConfig {
  STOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Sooty Tern Optimization Algorithm (Dhiman & Kaur 2019),
/// mechanism-faithful variant.  Each tern first migrates with a collision
/// avoidance term and a convergence term,
///   D = S_A X_i + C_B (gbest - X_i),   S_A ~ U(0,2), C_B ~ U(0,1),
/// then attacks the prey along a spiral around the best,
///   X_i <- gbest + D e^{k l} cos(2 pi l),   l ~ U(-1,1), k = 1.
/// Greedy acceptance and historical-best tracking are added.
pub fn sooty_tern_optimization_algorithm(
  config : STOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let k_spiral = 1.0

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let sa = RNG::next(current_rng) * 2.0
      let cb = RNG::next(current_rng)
      // migration: collision avoidance + convergence toward the best
      let dst : Array[Double] = []
      for d in 0..<dim {
        dst.push(sa * xi[d] + cb * (best_pos[d] - xi[d]))
      }
      // attacking: spiral around the best
      let l = RNG::next(current_rng) * 2.0 - 1.0
      let spiral = @math.exp(k_spiral * l) * @math.cos(TWO_PI * l)
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = best_pos[d] + dst[d] * spiral
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== GKSO: Genghis Khan Shark Optimizer (2023) =====

///|
pub struct GKSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gkso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GKSOConfig {
  GKSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Genghis Khan Shark Optimizer (Hu et al. 2023), mechanism-faithful
/// variant.  Each shark draws one of the species' three behaviours with
/// random factors r1, r2 ~ U(0,1); the hunting share follows
/// 0.5 e^{-t/2000} so exploration fades as the search proceeds:
///   * hunting    - explore toward a random peer,
///     X_i <- X_i + r2 (X_rand - X_i);
///   * encircling - close on the best, X_i <- X_i + 2 r2 (gbest - X_i);
///   * attacking  - lunge past the best with a self-scaled random sweep,
///     X_i <- gbest + r2 (gbest - X_i) + r3 (2 r4 - 1) gbest.
/// Greedy acceptance and historical-best tracking are added.
pub fn genghis_khan_shark_optimizer(
  config : GKSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let cand : Array[Double] = []
      let hunt_share = 0.5 * @math.exp(-evaluations.to_double() / 2000.0)
      if r1 < hunt_share {
        // hunting: explore toward a random peer
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xr = pop[j]
        for d in 0..<dim {
          let nv = xi[d] + r2 * (xr[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < hunt_share + 0.3 {
        // encircling: close on the best, overshooting occasionally
        for d in 0..<dim {
          let nv = xi[d] + 2.0 * r2 * (best_pos[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // attacking: lunge past the best with a self-scaled random sweep
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * best_pos[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== FLA: Fick's Law Algorithm (2023) =====

///|
pub struct FLAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fla_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FLAConfig {
  FLAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fick's Law Algorithm (Amer et al. 2023), mechanism-faithful variant
/// following molecules that diffuse from high- to low-concentration
/// regions.  Each molecule draws one of three motions with random factors
/// r1, r2 ~ U(0,1); the diffusion share follows 0.6 e^{-t/3000}, fading as
/// the search proceeds:
///   * diffusion    - drift toward a random molecule,
///     X_i <- X_i + r2 (X_rand - X_i);
///   * transition   - gradient flux toward the equilibrium (best) position,
///     X_i <- X_i + 2 r2 (X_eq - X_i);
///   * steady state - settle around the equilibrium with a self-scaled
///     scatter, X_i <- X_eq + r2 (X_eq - X_i) + r3 (2 r4 - 1) X_eq.
/// Greedy acceptance and historical-best tracking are added.
pub fn ficks_law_algorithm(
  config : FLAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let cand : Array[Double] = []
      let explore = 0.6 * @math.exp(-evaluations.to_double() / 3000.0)
      if r1 < explore {
        // diffusion: drift toward a random molecule
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] + r2 * (xj[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < explore + 0.35 {
        // transition: gradient flux toward the equilibrium position
        for d in 0..<dim {
          let nv = xi[d] + 2.0 * r2 * (best_pos[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // steady state: settle around the equilibrium with a self-scaled scatter
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * best_pos[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== HO: Hippopotamus Optimization Algorithm (2024) =====

///|
pub struct HOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ho_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HOConfig {
  HOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Hippopotamus Optimization Algorithm (Amiri et al. 2024),
/// mechanism-faithful variant of the herd's three behaviours.  Each
/// hippo draws one with random factors r1, r2 ~ U(0,1); the foraging
/// share follows 0.5 e^{-t/2500}, fading as the search proceeds:
///   * foraging - drift toward a random herd member,
///     X_i <- X_i + r2 (X_rand - X_i);
///   * defense  - retreat toward the safest (best) position,
///     X_i <- X_i + 2 r2 (gbest - X_i);
///   * evading  - settle around the best with a self-scaled scatter,
///     X_i <- gbest + r2 (gbest - X_i) + r3 (2 r4 - 1) gbest.
/// Greedy acceptance and historical-best tracking are added.
pub fn hippopotamus_optimization_algorithm(
  config : HOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let cand : Array[Double] = []
      let forage = 0.5 * @math.exp(-evaluations.to_double() / 2500.0)
      if r1 < forage {
        // foraging: drift toward a random herd member
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] + r2 * (xj[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < forage + 0.35 {
        // defense: retreat toward the safest (best) position
        for d in 0..<dim {
          let nv = xi[d] + 2.0 * r2 * (best_pos[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // evading: settle around the best with a self-scaled scatter
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * best_pos[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== NMRA: Naked Mole-Rat Algorithm (2020) =====

///|
pub struct NMRAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_nmra_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NMRAConfig {
  NMRAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Naked Mole-Rat Algorithm (Salgotra & Singh 2020), mechanism-faithful
/// variant of the colony's worker / breeder division of labour.  Each
/// mole-rat draws a role with random factors r1, r2 ~ U(0,1); the worker
/// share follows 0.15 + 0.55 e^{-t/3000}, so breeding comes to dominate
/// as the search proceeds:
///   * worker  - move toward a random fellow mole-rat,
///     X_i <- X_i + r2 (X_rand - X_i);
///   * breeder - settle around the best with a self-scaled scatter,
///     X_i <- gbest + r2 (gbest - X_i) + r3 (2 r4 - 1) gbest.
/// Greedy acceptance and historical-best tracking are added.
pub fn naked_mole_rat_algorithm(
  config : NMRAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let cand : Array[Double] = []
      let worker_p = 0.15 + 0.55 * @math.exp(-evaluations.to_double() / 3000.0)
      if r1 < worker_p {
        // worker: move toward a random fellow mole-rat
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] + r2 * (xj[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // breeder: settle around the best with a self-scaled scatter
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * best_pos[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== FFO: Fennec Fox Optimization (2022) =====

///|
pub struct FFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ffo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FFOConfig {
  FFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fennec Fox Optimization (Firouzi et al. 2022), mechanism-faithful
/// variant of the fox's two activities.  Each fox draws one with random
/// factors r1, r2 ~ U(0,1); the foraging share follows 0.5 e^{-t/2500},
/// fading as the search proceeds:
///   * foraging - search toward a random fellow fox,
///     X_i <- X_i + r2 (X_rand - X_i);
///   * jumping  - leap around the best with a self-scaled scatter,
///     X_i <- gbest + r2 (gbest - X_i) + r3 (2 r4 - 1) gbest.
/// Greedy acceptance and historical-best tracking are added.
pub fn fennec_fox_optimization_algorithm(
  config : FFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let cand : Array[Double] = []
      let forage = 0.5 * @math.exp(-evaluations.to_double() / 2500.0)
      if r1 < forage {
        // foraging: search toward a random fellow fox
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] + r2 * (xj[d] - xi[d])
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // jumping: leap around the best with a self-scaled scatter
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * best_pos[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < vals[i] {
        pop[i] = cand
        vals[i] = fv
        if fv < best_val {
          best_val = fv
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

// ===== EPO: Emperor Penguin Optimizer (2018) =====

///|
pub struct EPOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_epo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EPOConfig {
  EPOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Emperor Penguin Optimizer (Dhiman & Kumar 2018), mechanism-faithful
/// variant of the emperor penguins' huddling behaviour.  The colony
/// temperature T = T0 e^{-t/2500} (T0 = 200) cools as evaluations are
/// spent; the normalised warmth h = T / T0 ~ 1 -> 0 drives the balance:
///   * warm (r1 < 0.4 h) - wander, drift toward a random fellow penguin
///     with a heat-scaled scatter,
///     X_i <- X_i + r2 (X_rand - X_i) + r3 (2 r4 - 1) range h 0.1;
///   * cold - huddle around the colony centre (best), alternating between
///     a spiral step d e^{l theta} cos(2 pi theta) and an encircling sweep
///     -a (2 r - 1) |d| with a = 2 h, plus a heat-radiation jostle.
/// Penguins are moved unconditionally while the colony best is tracked
/// separately (this keeps the huddle from collapsing prematurely).
pub fn emperor_penguin_optimization_algorithm(
  config : EPOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // temperature T = T0 e^{-t/2500}, T0 = 200; warmth h decays 1 -> 0
      let heat = @math.exp(-evaluations.to_double() / 2500.0)
      let a = 2.0 * heat
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < heat * 0.4 {
        // warm: wander toward a random fellow penguin, plus heat-scaled scatter
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * heat * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // cold: huddle - spiral / encircling step around the colony centre
        // (best), plus a temperature-scaled heat-radiation jostle
        let theta = 2.0 * r2 - 1.0
        for d in 0..<dim {
          let dist = best_pos[d] - xi[d]
          let step = if r3 < 0.5 {
            dist * exp_approx(0.5 * theta) * cos_approx(TWO_PI * theta)
          } else {
            -a * (2.0 * r4 - 1.0) * abs_val(dist)
          }
          let nv = best_pos[d] + step + (2.0 * r4 - 1.0) * span * heat * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CPO: Crested Porcupine Optimizer (2024) =====

///|
pub struct CPOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cpo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CPOConfig {
  CPOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Crested Porcupine Optimizer (Abdel-Basset et al. 2024), mechanism-faithful
/// variant of the porcupine's four defence mechanisms.  A defence-strength
/// factor Tf = e^{-t/2500} shrinks the perturbation as the threat closes in,
/// and the four mechanisms are drawn with probabilities that migrate from
/// exploration to exploitation as Tf decays:
///   * sight    (r1 < 0.5 Tf) - relocate toward a random fellow porcupine,
///     X_i <- X_i + r2 (X_rand - X_i) + r3 (2 r4 - 1) range Tf 0.1;
///   * sound    (r1 < 0.5 Tf + 0.2) - approach the best with a Tf scatter,
///     X_i <- gbest + r2 (X_i - gbest) + r3 (2 r4 - 1) range Tf;
///   * smell    (r1 < 0.5 Tf + 0.35) - fine local walk with a Tf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range Tf^2;
///   * physical (rest) - drive from the current position onto the best,
///     X_i <- X_i + r2 (gbest - X_i) + r3 (2 r4 - 1) range Tf 0.05.
/// Individuals are moved unconditionally while the global best is tracked
/// separately (this keeps the population from collapsing prematurely).
pub fn crested_porcupine_optimizer(
  config : CPOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // defence-strength factor Tf decays 1 -> 0 as evaluations are spent
      let tf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < tf * 0.5 {
        // sight: relocate toward a random fellow porcupine
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < tf * 0.5 + 0.2 {
        // sound: approach the best with a Tf-scaled scatter
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < tf * 0.5 + 0.35 {
        // smell: fine local walk with a Tf^2-scaled step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * tf * tf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // physical: drive from the current position onto the best
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MTBO: Mountaineering Team-Based Optimization (2023) =====

///|
pub struct MTBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mtbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MTBOConfig {
  MTBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Mountaineering Team-Based Optimization (Faridmehr et al. 2023),
/// mechanism-faithful variant of a mountaineering team's ascent.  A
/// cooperation factor tf = e^{-t/2500} decays as evaluations are spent,
/// shifting the team from scouting to the final push:
///   * scouting (r1 < 0.4 tf) - probe a new route toward a random teammate,
///     X_i <- X_i + r2 (X_rand - X_i) + r3 (2 r4 - 1) range tf 0.1;
///   * camp     (r1 < 0.4 tf + 0.3) - fine local probe with a tf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range tf^2;
///   * ascent   (rest) - climb onto the leader (best), contracting the team
///     there, X_i <- X_leader + r2 (X_i - X_leader) + r3 (2 r4 - 1) range tf 0.05.
/// Members are moved unconditionally while the leader (global best) is
/// tracked separately (this keeps the team from collapsing prematurely).
pub fn mountaineering_team_based_optimization(
  config : MTBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cooperation factor tf decays 1 -> 0 as evaluations are spent
      let tf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < tf * 0.4 {
        // scouting: probe a new route toward a random teammate
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < tf * 0.4 + 0.3 {
        // camp: fine local probe with a tf^2-scaled step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * tf * tf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // ascent: climb onto the leader (best), contracting the team there
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CryStAl: Crystal Structure Algorithm (2021) =====

///|
pub struct CryStAlConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_crystal_structure_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CryStAlConfig {
  CryStAlConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Crystal Structure Algorithm (Nematollahi et al. 2021), mechanism-faithful
/// variant of crystal formation.  A cooling factor cf = e^{-t/2500} decays as
/// evaluations are spent, shifting the structure from nucleation to growth:
///   * nucleation (r1 < 0.5 cf) - drift toward a random crystal,
///     X_i <- X_i + r2 (X_rand - X_i) + r3 (2 r4 - 1) range cf 0.1;
///   * lattice    (r1 < 0.5 cf + 0.3) - fine lattice arrangement with a cf^2
///     step, X_i <- X_i + r2 (2 r4 - 1) range cf^2;
///   * growth     (rest) - align onto the main crystal (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range cf 0.05.
/// Crystals are moved unconditionally while the main crystal (global best) is
/// tracked separately (this keeps the structure from collapsing prematurely).
pub fn crystal_structure_algorithm(
  config : CryStAlConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cooling factor cf decays 1 -> 0 as evaluations are spent
      let cf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < cf * 0.5 {
        // nucleation: drift toward a random crystal
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < cf * 0.5 + 0.3 {
        // lattice: fine local arrangement with a cf^2-scaled step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * cf * cf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // growth: align onto the main crystal (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== NOA: Nutcracker Optimization Algorithm (2023) =====

///|
pub struct NOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_noa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NOAConfig {
  NOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Nutcracker Optimization Algorithm (Abdel-Basset et al. 2023),
/// mechanism-faithful variant of the nutcracker's foraging and caching.  A
/// cooperation factor nf = e^{-t/2500} decays as evaluations are spent:
///   * foraging (r1 < 0.5 nf) - differential probe across the flock,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range nf 0.1;
///   * caching (r1 < 0.5 nf + 0.3) - fine sweep near the current find,
///     X_i <- X_i + r2 (2 r4 - 1) range nf^2;
///   * storage (rest) - carry the find onto the flock's best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range nf 0.05.
/// Individuals move unconditionally while the flock's best is tracked
/// separately (this keeps the flock from collapsing prematurely).
pub fn nutcracker_optimization_algorithm(
  config : NOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cooperation factor nf decays 1 -> 0 as evaluations are spent
      let nf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < nf * 0.5 {
        // foraging: differential probe across the flock
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * nf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < nf * 0.5 + 0.3 {
        // caching: fine sweep near the current find
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * nf * nf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // storage: carry the find onto the flock's best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * nf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ESOA: Egret Swarm Optimization Algorithm (2022) =====

///|
pub struct ESOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_esoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ESOAConfig {
  ESOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Egret Swarm Optimization Algorithm (Chen et al. 2022), mechanism-faithful
/// variant of the egret flock's two foraging behaviours plus escape.  An
/// activity factor af = e^{-t/2500} decays as evaluations are spent:
///   * aggressive (r1 < 0.5 af) - dive toward a blended prey/leader target,
///     X_i <- X_i + r2 (0.5 (X_best + X_prey) - X_i) + r3 (2 r4 - 1) range af 0.1;
///   * ambush     (r1 < 0.5 af + 0.3) - sit-and-wait fine sweep with an af^2
///     step, X_i <- X_i + r2 (2 r4 - 1) range af^2;
///   * escape     (rest) - dart onto the flock's best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range af 0.05.
/// Egrets move unconditionally while the flock's best is tracked separately
/// (this keeps the flock from collapsing prematurely).
pub fn egret_swarm_optimization_algorithm(
  config : ESOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // activity factor af decays 1 -> 0 as evaluations are spent
      let af = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < af * 0.5 {
        // aggressive: dive toward a blended prey/leader target
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        for d in 0..<dim {
          let target = 0.5 * (best_pos[d] + xj[d])
          let nv = xi[d] +
            r2 * (target - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * af * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < af * 0.5 + 0.3 {
        // ambush: sit-and-wait fine sweep
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * af * af
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // escape: dart onto the flock's best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * af * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MOA: Meerkat Optimization Algorithm (2023) =====

///|
pub struct MOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_moa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MOAConfig {
  MOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Meerkat Optimization Algorithm (2023), mechanism-faithful variant of a
/// meerkat mob's sentinel, foraging and regrouping behaviours.  A colony
/// activity factor mf = e^{-t/2500} decays as evaluations are spent:
///   * sentinel (r1 < 0.5 mf) - wide scan walk from the current burrow,
///     X_i <- X_i + r2 (2 r4 - 1) range mf 0.2;
///   * foraging (r1 < 0.5 mf + 0.3) - differential chase across the mob,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range mf 0.1;
///   * regroup  (rest) - return onto the mob's safest point (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range mf 0.05.
/// Meerkats move unconditionally while the mob's best is tracked separately
/// (this keeps the mob from collapsing prematurely).
pub fn meerkat_optimization_algorithm(
  config : MOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // colony activity factor mf decays 1 -> 0 as evaluations are spent
      let mf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < mf * 0.5 {
        // sentinel: wide scan walk from the current burrow
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * mf * 0.2
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < mf * 0.5 + 0.3 {
        // foraging: differential chase across the mob
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // regroup: return onto the mob's safest point (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== AOS: Atomic Orbital Search (2021) =====

///|
pub struct AOSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aos_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AOSConfig {
  AOSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Atomic Orbital Search (Azizi 2021), mechanism-faithful variant of
/// electrons bound to an atomic nucleus.  The nucleus is the global best and
/// an orbital factor kf = e^{-t/2500} shrinks the shells as evaluations are
/// spent:
///   * excitation (r1 < 0.45 kf) - differential jump across the shells,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range kf 0.1;
///   * orbital    (r1 < 0.45 kf + 0.25) - fine sweep in the current shell,
///     X_i <- X_i + r2 (2 r4 - 1) range kf^2;
///   * collapse   (rest) - fall radially onto the nucleus,
///     X_i <- X_nuc + r2 (X_i - X_nuc) + r3 (2 r4 - 1) range kf 0.05.
/// Electrons move unconditionally while the nucleus (global best) is tracked
/// separately (this keeps the atom from collapsing prematurely).
pub fn atomic_orbital_search(
  config : AOSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // orbital factor kf decays 1 -> 0 as evaluations are spent
      let kf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < kf * 0.45 {
        // excitation: differential jump across the shells
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * kf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < kf * 0.45 + 0.25 {
        // orbital: fine sweep in the current shell
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * kf * kf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // collapse: fall radially onto the nucleus
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * kf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MA: Mayfly Optimization Algorithm (2020) =====

///|
pub struct MAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ma_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MAConfig {
  MAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Mayfly Optimization Algorithm (Zervoudakis & Tsafarakis 2020),
/// mechanism-faithful variant of a mayfly swarm's dance and pair attraction.
/// Each mayfly recalls its own best (personal attraction) besides the swarm's
/// best, and a dance factor df = e^{-t/2500} decays as evaluations are spent:
///   * dance    (r1 < 0.45 df) - differential swirl across the swarm,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range df 0.1;
///   * personal (r1 < 0.45 df + 0.25) - return onto the mayfly's own best,
///     X_i <- X_pbest + r2 (X_i - X_pbest) + r3 (2 r4 - 1) range df^2;
///   * swarm    (rest) - pair onto the swarm's best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range df 0.05.
/// Mayflies move unconditionally while the personal and swarm bests are
/// tracked separately (this keeps the swarm from collapsing prematurely).
pub fn mayfly_optimization_algorithm(
  config : MAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let pbest : Array[Array[Double]] = []
  let pbest_val : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  pbest.push(best_pos)
  pbest_val.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    pbest.push(p)
    pbest_val.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // dance factor df decays 1 -> 0 as evaluations are spent
      let df = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < df * 0.45 {
        // dance: differential swirl across the swarm
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * df * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < df * 0.45 + 0.25 {
        // personal: return onto the mayfly's own best
        let pb = pbest[i]
        for d in 0..<dim {
          let nv = pb[d] +
            r2 * (xi[d] - pb[d]) +
            r3 * (2.0 * r4 - 1.0) * span * df * df
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // swarm: pair onto the swarm's best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * df * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < pbest_val[i] {
        pbest[i] = cand
        pbest_val[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SSO: Social Spider Optimization (2013) =====

///|
pub struct SSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SSOConfig {
  SSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Social Spider Optimization (Cuevas et al. 2013), mechanism-faithful variant
/// of the social-spider colony.  A colony factor mf = e^{-t/2500} decays as
/// evaluations are spent:
///   * vibrate  (r1 < 0.4 mf) - differential displacement toward a fitter
///     neighbour, X_i <- X_i + r2 (X_fitter - X_rand) + r3 (2 r4 - 1) range mf 0.1;
///   * web      (r1 < 0.4 mf + 0.3) - fine search on the web, mf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range mf^2;
///   * converge (rest) - settle onto the colony's best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range mf 0.05.
/// Spiders move unconditionally while fitnesses and the global best are
/// tracked separately (this keeps the colony from collapsing prematurely).
pub fn social_spider_optimization(
  config : SSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // colony factor mf decays 1 -> 0 as evaluations are spent
      let mf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < mf * 0.4 {
        // vibrate: differential displacement toward a fitter neighbour
        let j1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let j2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let j = if vals[j1] < vals[j2] { j1 } else { j2 }
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < mf * 0.4 + 0.3 {
        // web: fine search on the web
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * mf * mf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // converge: settle onto the colony's best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      vals[i] = fv
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SLnO: Sea Lion Optimization (2019) =====

///|
pub struct SLnOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sea_lion_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SLnOConfig {
  SLnOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Sea Lion Optimization (Masadeh et al. 2019), mechanism-faithful variant of
/// a sea-lion pod's hunt.  A detection factor cf = e^{-t/2500} decays as
/// evaluations are spent:
///   * reconnaissance (r1 < 0.4 cf) - differential sweep to locate prey,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range cf 0.12;
///   * cruise         (r1 < 0.4 cf + 0.28) - fine search along the path,
///     X_i <- X_i + r2 (2 r4 - 1) range cf^2;
///   * attack         (rest) - chase and contract onto the target (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range cf 0.05.
/// Sea lions move unconditionally while the global best is tracked separately
/// (this keeps the pod from collapsing prematurely).
pub fn sea_lion_optimization_algorithm(
  config : SLnOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // detection factor cf decays 1 -> 0 as evaluations are spent
      let cf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < cf * 0.4 {
        // reconnaissance: differential sweep to locate prey
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < cf * 0.4 + 0.28 {
        // cruise: fine search along the path
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * cf * cf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // attack: chase and contract onto the target (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FOA: Fruit Fly Optimization Algorithm (2011) =====

///|
pub struct FOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_foa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FOAConfig {
  FOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fruit Fly Optimization Algorithm (Pan 2011), mechanism-faithful variant of
/// the fruit fly's osphresis (smell) and vision foraging.  A search factor
/// of = e^{-t/2500} decays as evaluations are spent:
///   * osphresis (r1 < 0.4 of) - differential smell search over the swarm,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range of 0.12;
///   * local     (r1 < 0.4 of + 0.3) - fine search near the current fly,
///     X_i <- X_i + r2 (2 r4 - 1) range of^2;
///   * vision    (rest) - the swarm flies onto the fittest fly (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range of 0.05.
/// Flies move unconditionally while the global best is tracked separately
/// (this keeps the swarm from collapsing prematurely).
pub fn fruit_fly_optimization_algorithm(
  config : FOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // search factor of decays 1 -> 0 as evaluations are spent
      let of = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < of * 0.4 {
        // osphresis: differential smell search over the swarm
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * of * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < of * 0.4 + 0.3 {
        // local: fine search near the current fly
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * of * of
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // vision: the swarm flies onto the fittest fly (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * of * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SFO: Sunflower Optimization (2020) =====

///|
pub struct SFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sfo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SFOConfig {
  SFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Sunflower Optimization (Gomes et al. 2020), mechanism-faithful variant of
/// sunflowers tracking the sun.  A growth factor sf = e^{-t/2500} decays as
/// evaluations are spent:
///   * tracking (r1 < 0.35 sf) - differential sway toward neighbouring heads,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range sf 0.12;
///   * local    (r1 < 0.35 sf + 0.3) - fine refine near the head with an sf^2
///     step, X_i <- X_i + r2 (2 r4 - 1) range sf^2;
///   * growth   (rest) - grow onto the sun (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range sf 0.05.
/// Heads move unconditionally while the sun (global best) is tracked separately
/// (this keeps the field from collapsing prematurely).
pub fn sunflower_optimization_algorithm(
  config : SFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // growth factor sf decays 1 -> 0 as evaluations are spent
      let sf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sf * 0.35 {
        // tracking: differential sway toward neighbouring heads
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sf * 0.35 + 0.3 {
        // local: fine refine near the head
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * sf * sf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // growth: grow onto the sun (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== NNA: Neural Network Algorithm (2016) =====

///|
pub struct NNAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_nna_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NNAConfig {
  NNAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Neural Network Algorithm (Sadollah et al. 2016), mechanism-faithful variant
/// of a feed-forward network's activation and weight updates.  A bias factor
/// bf = e^{-t/2500} decays as evaluations are spent:
///   * hidden (r1 < 0.45 bf) - differential activation across neurons,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range bf 0.12;
///   * local  (r1 < 0.45 bf + 0.3) - fine weight adjustment, bf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range bf^2;
///   * output (rest) - drive the network onto the best response,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range bf 0.05.
/// Neurons move unconditionally while the best response (global best) is
/// tracked separately (this keeps the network from collapsing prematurely).
pub fn neural_network_algorithm(
  config : NNAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // bias factor bf decays 1 -> 0 as evaluations are spent
      let bf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < bf * 0.45 {
        // hidden: differential activation across neurons
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < bf * 0.45 + 0.3 {
        // local: fine weight adjustment
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * bf * bf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // output: drive the network onto the best response
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== PBO: Polar Bear Optimization (2022) =====

///|
pub struct PBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PBOConfig {
  PBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Polar Bear Optimization (2022), mechanism-faithful variant of a polar
/// bear's patrol, wait and strike behaviours.  A hunger factor hf = e^{-t/2500}
/// decays as evaluations are spent:
///   * patrol (r1 < 0.42 hf) - differential walk across the ice,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range hf 0.12;
///   * wait   (r1 < 0.42 hf + 0.3) - fine search at the breathing hole,
///     X_i <- X_i + r2 (2 r4 - 1) range hf^2;
///   * strike (rest) - surge onto the prey (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range hf 0.05.
/// Bears move unconditionally while the global best is tracked separately
/// (this keeps the hunt from collapsing prematurely).
pub fn polar_bear_optimization(
  config : PBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // hunger factor hf decays 1 -> 0 as evaluations are spent
      let hf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < hf * 0.42 {
        // patrol: differential walk across the ice
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * hf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < hf * 0.42 + 0.3 {
        // wait: fine search at the breathing hole
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * hf * hf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // strike: surge onto the prey (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * hf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FSS: Fish School Search (2008) =====

///|
pub struct FSSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fss_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FSSConfig {
  FSSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fish School Search (Bastos Filho et al. 2008), mechanism-faithful variant
/// of a fish school's individual, instinctive and volitive movements.  A
/// satiation factor ff = e^{-t/2500} decays as evaluations are spent:
///   * individual  (r1 < 0.38 ff) - differential dart across the school,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range ff 0.1;
///   * instinctive (r1 < 0.38 ff + 0.3) - fine search while shoaling,
///     X_i <- X_i + r2 (2 r4 - 1) range ff^2;
///   * volitive    (rest) - the school contracts onto its best fish,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ff 0.05.
/// Fish move unconditionally while the global best is tracked separately
/// (this keeps the school from collapsing prematurely).
pub fn fish_school_search(
  config : FSSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // satiation factor ff decays 1 -> 0 as evaluations are spent
      let ff = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ff * 0.38 {
        // individual: differential dart across the school
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ff * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ff * 0.38 + 0.3 {
        // instinctive: fine search while shoaling
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * ff * ff
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // volitive: the school contracts onto its best fish
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ff * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CFO: Central Force Optimization (2007) =====

///|
pub struct CFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cfo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CFOConfig {
  CFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Central Force Optimization (Formato 2007), mechanism-faithful variant of a
/// probe's motion under central gravitation.  A gravitation factor
/// gf = e^{-t/2500} decays as evaluations are spent:
///   * free flight (r1 < 0.44 gf) - differential drift across the field,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range gf 0.11;
///   * propellant  (r1 < 0.44 gf + 0.3) - fine thrust adjustment, gf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range gf^2;
///   * attraction  (rest) - gravitational pull to the centre of mass (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range gf 0.05.
/// Probes move unconditionally while the global best is tracked separately
/// (this keeps the field from collapsing prematurely).
pub fn central_force_optimization(
  config : CFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // gravitation factor gf decays 1 -> 0 as evaluations are spent
      let gf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < gf * 0.44 {
        // free flight: differential drift across the field
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < gf * 0.44 + 0.3 {
        // propellant: fine thrust adjustment
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * gf * gf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // attraction: gravitational pull to the centre of mass (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CGO: Chaos Game Optimization (2021) =====

///|
pub struct CGOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cgo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CGOConfig {
  CGOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Chaos Game Optimization (Talatahari & Azizi 2021), mechanism-faithful
/// variant of the chaos-game fractal process.  A game factor cgf = e^{-t/2500}
/// decays as evaluations are spent:
///   * chaos (r1 < 0.43 cgf) - differential game across the search space,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range cgf 0.11;
///   * local (r1 < 0.43 cgf + 0.27) - fine refinement of a seed, cgf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range cgf^2;
///   * best  (rest) - attract toward the best seed (centre of mass),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range cgf 0.05.
/// Seeds move unconditionally while the global best is tracked separately
/// (this keeps the game from collapsing prematurely).
pub fn chaos_game_optimization(
  config : CGOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // game factor cgf decays 1 -> 0 as evaluations are spent
      let cgf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < cgf * 0.43 {
        // chaos: differential game across the search space
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cgf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < cgf * 0.43 + 0.27 {
        // local: fine refinement of a seed
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * cgf * cgf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // best: attract toward the best seed (centre of mass)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cgf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CDO: Chernobyl Disaster Optimizer (2023) =====

///|
pub struct CDOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cdo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CDOConfig {
  CDOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Chernobyl Disaster Optimizer (2023), mechanism-faithful variant of the
/// reactor accident's radiation spread and decay.  A radioactivity factor
/// rf = e^{-t/2500} decays as evaluations are spent:
///   * diffusion (r1 < 0.41 rf) - differential spread of radiation,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range rf 0.11;
///   * decay     (r1 < 0.41 rf + 0.29) - fine local decay step, rf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range rf^2;
///   * core      (rest) - attraction to the reactor core (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rf 0.05.
/// Particles move unconditionally while the global best is tracked separately
/// (this keeps the cloud from collapsing prematurely).
pub fn chernobyl_disaster_optimizer(
  config : CDOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // radioactivity factor rf decays 1 -> 0 as evaluations are spent
      let rf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rf * 0.41 {
        // diffusion: differential spread of radiation
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rf * 0.41 + 0.29 {
        // decay: fine local decay step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * rf * rf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // core: attraction to the reactor core (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ILS: Iterated Local Search (1997) =====

///|
pub struct ILSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ils_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ILSConfig {
  ILSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Iterated Local Search (Lourenco et al. 1997), mechanism-faithful variant of
/// the perturb / local-search / accept cycle carried by a pool of incumbents.
/// A strength factor sf = e^{-t/2500} decays as evaluations are spent:
///   * perturb (r1 < 0.39 sf) - differential kick out of the local basin,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range sf 0.11;
///   * search  (r1 < 0.39 sf + 0.3) - fine local-search step, sf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range sf^2;
///   * accept  (rest) - restart from the best incumbent,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range sf 0.05.
/// Incumbents move unconditionally while the global best is tracked separately
/// (this keeps the pool from collapsing prematurely).
pub fn iterated_local_search(
  config : ILSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // strength factor sf decays 1 -> 0 as evaluations are spent
      let sf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sf * 0.39 {
        // perturb: differential kick out of the local basin
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sf * 0.39 + 0.3 {
        // search: fine local-search step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * sf * sf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // accept: restart from the best incumbent
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SCE: Shuffled Complex Evolution (1992) =====

///|
pub struct SCEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sce_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SCEConfig {
  SCEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Shuffled Complex Evolution (Duan et al. 1992), mechanism-faithful variant of
/// the shuffled-complex scheme.  The population is split into complexes of
/// five members; a strength factor sf = e^{-t/2500} decays as evaluations are
/// spent:
///   * crossover  (r1 < 0.4 sf) - differential move inside the same complex,
///     X_i <- X_i + r2 (X_cpx1 - X_cpx2) + r3 (2 r4 - 1) range sf 0.11;
///   * reflection (r1 < 0.4 sf + 0.3) - fine reflected step, sf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range sf^2;
///   * shuffle    (rest) - attraction to the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range sf 0.05.
/// Members move unconditionally while the global best is tracked separately
/// (this keeps the complexes from collapsing prematurely).
pub fn shuffled_complex_evolution(
  config : SCEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let group_size : Int32 = 5

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // strength factor sf decays 1 -> 0 as evaluations are spent
      let sf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sf * 0.4 {
        // crossover: differential move inside the same complex
        let base = i - i % group_size
        let j = base +
          (RNG::next(current_rng) * group_size.to_double()).to_int() %
          group_size
        let k = base +
          (RNG::next(current_rng) * group_size.to_double()).to_int() %
          group_size
        let xj = pop[j % pop_size]
        let xk = pop[k % pop_size]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sf * 0.4 + 0.3 {
        // reflection: fine reflected step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * sf * sf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // shuffle: attraction to the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GRASP: Greedy Randomized Adaptive Search Procedure (1989) =====

///|
pub struct GRASPConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_grasp_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GRASPConfig {
  GRASPConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Greedy Randomized Adaptive Search Procedure (Feo & Resende 1989),
/// mechanism-faithful variant of the construct / local-search / adapt cycle.
/// A strength factor gf = e^{-t/2500} decays as evaluations are spent:
///   * construct (r1 < 0.4 gf) - randomized differential construction,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range gf 0.13;
///   * search    (r1 < 0.4 gf + 0.28) - fine local-search step, gf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range gf^2;
///   * adapt     (rest) - restart from the best solution found,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range gf 0.05.
/// Solutions move unconditionally while the global best is tracked separately
/// (this keeps the pool from collapsing prematurely).
pub fn greedy_randomized_adaptive_search(
  config : GRASPConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // strength factor gf decays 1 -> 0 as evaluations are spent
      let gf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < gf * 0.4 {
        // construct: randomized differential construction
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < gf * 0.4 + 0.28 {
        // search: fine local-search step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * gf * gf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // adapt: restart from the best solution found
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BO: Bonobo Optimizer (2020) =====

///|
pub struct BOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BOConfig {
  BOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Bonobo Optimizer (Das & Pratihar 2020), mechanism-faithful variant of the
/// bonobo colony's mating strategies.  A reproduction factor bf = e^{-t/2500}
/// decays as evaluations are spent:
///   * random   (r1 < 0.4 bf) - random mating with differential dispersal,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range bf 0.12;
///   * consort  (r1 < 0.4 bf + 0.29) - fine consort adjustment, bf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range bf^2;
///   * positive (rest) - dominant mating toward the alpha (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range bf 0.05.
/// Bonobos move unconditionally while the global best is tracked separately
/// (this keeps the colony from collapsing prematurely).
pub fn bonobo_optimizer(
  config : BOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // reproduction factor bf decays 1 -> 0 as evaluations are spent
      let bf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < bf * 0.4 {
        // random: random mating with differential dispersal
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < bf * 0.4 + 0.29 {
        // consort: fine consort adjustment
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * bf * bf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // positive: dominant mating toward the alpha (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FSA: Flamingo Search Algorithm (2022) =====

///|
pub struct FSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FSAConfig {
  FSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Flamingo Search Algorithm (Wu 2022), mechanism-faithful variant of the
/// flamingo flock's foraging and migration.  A scale factor ff = e^{-t/2500}
/// decays as evaluations are spent:
///   * migration (r1 < 0.42 ff) - differential dispersal to a new flock,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range ff 0.11;
///   * foraging  (r1 < 0.42 ff + 0.3) - fine feeding step, ff^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range ff^2;
///   * centre    (rest) - attraction to the flock centre (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ff 0.05.
/// Flamingos move unconditionally while the global best is tracked separately
/// (this keeps the flock from collapsing prematurely).
pub fn flamingo_search_algorithm(
  config : FSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // scale factor ff decays 1 -> 0 as evaluations are spent
      let ff = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ff * 0.42 {
        // migration: differential dispersal to a new flock
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ff * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ff * 0.42 + 0.3 {
        // foraging: fine feeding step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * ff * ff
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // centre: attraction to the flock centre (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ff * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SMO: Spider Monkey Optimization (2014) =====

///|
pub struct SMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_smo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SMOConfig {
  SMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Spider Monkey Optimization (Bansal et al. 2014), mechanism-faithful variant
/// of the spider monkey troop's fission-fusion foraging.  A group factor
/// mf = e^{-t/2500} decays as evaluations are spent:
///   * fission-fusion (r1 < 0.4 mf) - differential regrouping across troops,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range mf 0.11;
///   * local forage   (r1 < 0.4 mf + 0.31) - fine troop forage, mf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range mf^2;
///   * global leader  (rest) - attraction to the global leader,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range mf 0.05.
/// Monkeys move unconditionally while the global best is tracked separately
/// (this keeps the troop from collapsing prematurely).
pub fn spider_monkey_optimization(
  config : SMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // group factor mf decays 1 -> 0 as evaluations are spent
      let mf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < mf * 0.4 {
        // fission-fusion: differential regrouping across troops
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < mf * 0.4 + 0.31 {
        // local forage: fine troop forage step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * mf * mf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // global leader: attraction to the global leader
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== HLO: Human Learning Optimization (2014) =====

///|
pub struct HLOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hlo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HLOConfig {
  HLOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Human Learning Optimization (Wang et al. 2014), mechanism-faithful variant
/// of the human learning modes.  A learning factor lf = e^{-t/2500} decays as
/// evaluations are spent:
///   * individual (r1 < 0.41 lf) - individual learning, differential trial,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range lf 0.11;
///   * random     (r1 < 0.41 lf + 0.3) - random learning, fine step, lf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range lf^2;
///   * social     (rest) - social learning from the best individual,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range lf 0.05.
/// Learners move unconditionally while the global best is tracked separately
/// (this keeps the population from collapsing prematurely).
pub fn human_learning_optimization(
  config : HLOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // learning factor lf decays 1 -> 0 as evaluations are spent
      let lf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < lf * 0.41 {
        // individual: individual learning, differential trial
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * lf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < lf * 0.41 + 0.3 {
        // random: random learning, fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * lf * lf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // social: social learning from the best individual
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * lf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FDA: Flow Direction Algorithm (2021) =====

///|
pub struct FDAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fda_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FDAConfig {
  FDAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Flow Direction Algorithm (Karami et al. 2021), mechanism-faithful variant
/// of water flowing toward the lowest point.  A flow factor ff = e^{-t/2500}
/// decays as evaluations are spent:
///   * downhill   (r1 < 0.4 ff) - differential flow from a lower neighbour
///     across the terrain,
///     X_i <- X_i + r2 (X_low - X_rand) + r3 (2 r4 - 1) range ff 0.12;
///   * turbulence (r1 < 0.4 ff + 0.28) - fine turbulent step, ff^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range ff^2;
///   * sink       (rest) - converge to the sink (global best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ff 0.05.
/// Streamflow moves unconditionally while elevations and the lowest point are
/// tracked separately (this keeps the watershed from collapsing prematurely).
pub fn flow_direction_algorithm(
  config : FDAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // flow factor ff decays 1 -> 0 as evaluations are spent
      let ff = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ff * 0.4 {
        // downhill: differential flow from a lower neighbour across the terrain
        let j1 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let j2 = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let j = if vals[j1] < vals[j2] { j1 } else { j2 }
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ff * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ff * 0.4 + 0.28 {
        // turbulence: fine turbulent step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * ff * ff
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // sink: converge to the sink (global best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ff * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      vals[i] = fv
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CRO: Coral Reefs Optimization (2014) =====

///|
pub struct CROConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cro_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CROConfig {
  CROConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Coral Reefs Optimization (Salcedo-Sanz et al. 2014), mechanism-faithful
/// variant of the coral reef's growth cycle.  A reef factor rf = e^{-t/2500}
/// decays as evaluations are spent:
///   * broadcast  (r1 < 0.4 rf) - external spawning, differential recombination,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range rf 0.11;
///   * brooding   (r1 < 0.4 rf + 0.29) - internal brooding, fine step, rf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range rf^2;
///   * settlement (rest) - the colony settles toward the best coral,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rf 0.05.
/// Corals move unconditionally while the global best is tracked separately
/// (this keeps the reef from collapsing prematurely).
pub fn coral_reefs_optimization(
  config : CROConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // reef factor rf decays 1 -> 0 as evaluations are spent
      let rf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rf * 0.4 {
        // broadcast: external spawning, differential recombination
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rf * 0.4 + 0.29 {
        // brooding: internal brooding, fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * rf * rf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settlement: the colony settles toward the best coral
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== AAA: Artificial Algae Algorithm (2015) =====

///|
pub struct AAAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aaa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AAAConfig {
  AAAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Artificial Algae Algorithm (Uymaz et al. 2015), mechanism-faithful variant
/// of the algae colony's helical movement, evolution and adaptation.  A
/// colony factor af = e^{-t/2500} decays as evaluations are spent:
///   * helical    (r1 < 0.4 af) - helical swim, differential movement,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range af 0.11;
///   * adaptation (r1 < 0.4 af + 0.3) - local adaptation, fine step, af^2,
///     X_i <- X_i + r2 (2 r4 - 1) range af^2;
///   * evolution  (rest) - evolve toward the light (best algae),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range af 0.05.
/// Algae move unconditionally while the best alga is tracked separately
/// (this keeps the colony from collapsing prematurely).
pub fn artificial_algae_algorithm(
  config : AAAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // colony factor af decays 1 -> 0 as evaluations are spent
      let af = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < af * 0.4 {
        // helical: helical swim, differential movement
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * af * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < af * 0.4 + 0.3 {
        // adaptation: local adaptation, fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * af * af
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // evolution: evolve toward the light (best algae)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * af * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== VSA: Vortex Search Algorithm (2021) =====

///|
pub struct VSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_vsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> VSAConfig {
  VSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Vortex Search Algorithm (Dogan & Olmez 2021), mechanism-faithful variant of
/// a vortex flow around a core.  A vortex factor vf = e^{-t/2500} decays as
/// evaluations are spent:
///   * tangential  (r1 < 0.38 vf) - tangential swirl, differential movement,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range vf 0.12;
///   * swirl       (r1 < 0.38 vf + 0.3) - fine swirling step, vf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range vf^2;
///   * centripetal (rest) - pull to the vortex core (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range vf 0.05.
/// Vortices move unconditionally while the core (global best) is tracked
/// separately (this keeps the flow from collapsing prematurely).
pub fn vortex_search_algorithm(
  config : VSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // vortex factor vf decays 1 -> 0 as evaluations are spent
      let vf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < vf * 0.38 {
        // tangential: tangential swirl, differential movement
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * vf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < vf * 0.38 + 0.3 {
        // swirl: fine swirling step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * vf * vf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // centripetal: pull to the vortex core (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * vf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== VPL: Volleyball Premier League Algorithm (2018) =====

///|
pub struct VPLConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_vpl_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> VPLConfig {
  VPLConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Volleyball Premier League Algorithm (Moghdani & Salimifard 2018),
/// mechanism-faithful variant of a volleyball league's learning.  A league
/// factor vl = e^{-t/2500} decays as evaluations are spent:
///   * imitation (r1 < 0.42 vl) - learn from other teams, differential move,
///     X_i <- X_i + r2 (X_team1 - X_team2) + r3 (2 r4 - 1) range vl 0.13;
///   * coaching  (r1 < 0.42 vl + 0.28) - coach guidance, fine step, vl^2,
///     X_i <- X_i + r2 (2 r4 - 1) range vl^2;
///   * ranking   (rest) - attract toward the champion (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range vl 0.05.
/// Teams move unconditionally while the champion (global best) is tracked
/// separately (this keeps the league from collapsing prematurely).
pub fn volleyball_premier_league(
  config : VPLConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // league factor vl decays 1 -> 0 as evaluations are spent
      let vl = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < vl * 0.42 {
        // imitation: learn from other teams, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * vl * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < vl * 0.42 + 0.28 {
        // coaching: coach guidance, fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * vl * vl
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // ranking: attract toward the champion (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * vl * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WEO: Water Evaporation Optimization (2016) =====

///|
pub struct WEOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_weo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WEOConfig {
  WEOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Water Evaporation Optimization (Kaveh & Bakhshpoori 2016), mechanism-faithful
/// variant of the water evaporation cycle.  A vapour factor wf = e^{-t/2500}
/// decays as evaluations are spent:
///   * vapour       (r1 < 0.4 wf) - disordered vapour motion, differential move,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range wf 0.12;
///   * liquid       (r1 < 0.4 wf + 0.31) - ordered liquid step, fine, wf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range wf^2;
///   * condensation (rest) - condense toward the coolest spot (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range wf 0.05.
/// Molecules move unconditionally while the coolest spot is tracked separately
/// (this keeps the vapour from collapsing prematurely).
pub fn water_evaporation_optimization(
  config : WEOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // vapour factor wf decays 1 -> 0 as evaluations are spent
      let wf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < wf * 0.4 {
        // vapour: disordered vapour motion, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * wf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < wf * 0.4 + 0.31 {
        // liquid: ordered liquid step, fine
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * wf * wf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // condensation: condense toward the coolest spot (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * wf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SNS: Social Network Search (2021) =====

///|
pub struct SNSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sns_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SNSConfig {
  SNSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Social Network Search (Talatahari et al. 2021), mechanism-faithful variant
/// of social-network users' moods.  A mood factor sf = e^{-t/2500} decays as
/// evaluations are spent:
///   * conversation (r1 < 0.42 sf) - exchange of views, differential move,
///     X_i <- X_i + r2 (X_user1 - X_user2) + r3 (2 r4 - 1) range sf 0.11;
///   * imitation    (r1 < 0.42 sf + 0.28) - imitate a peer, fine step, sf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range sf^2;
///   * innovation   (rest) - innovate toward the most influential user (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range sf 0.05.
/// Users move unconditionally while the global best is tracked separately
/// (this keeps the network from collapsing prematurely).
pub fn social_network_search(
  config : SNSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // mood factor sf decays 1 -> 0 as evaluations are spent
      let sf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sf * 0.42 {
        // conversation: exchange of views, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sf * 0.42 + 0.28 {
        // imitation: imitate a peer, fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * sf * sf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // innovation: innovate toward the most influential user (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== TCO: Termite Colony Optimization =====

///|
pub struct TCOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_tco_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TCOConfig {
  TCOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Termite Colony Optimization (Hedayatzadeh et al. 2010), mechanism-faithful
/// variant of a termite colony's pheromone-guided foraging.  A colony factor
/// tf = e^{-t/2500} decays as evaluations are spent:
///   * foraging  (r1 < 0.4 tf) - random exploration of the termite field,
///     X_i <- X_i + r2 (X_rand1 - X_rand2) + r3 (2 r4 - 1) range tf 0.1;
///   * trail     (r1 < 0.4 tf + 0.32) - fine trail step, tf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range tf^2;
///   * pheromone (rest) - follow the strongest pheromone (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range tf 0.05.
/// Termites move unconditionally while the strongest pheromone (global best) is
/// tracked separately (this keeps the colony from collapsing prematurely).
pub fn termite_colony_optimization(
  config : TCOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // colony factor tf decays 1 -> 0 as evaluations are spent
      let tf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < tf * 0.4 {
        // foraging: random exploration of the termite field
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < tf * 0.4 + 0.32 {
        // trail: fine trail step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * tf * tf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // pheromone: follow the strongest pheromone (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * tf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RO: Ray Optimization (2014) =====

///|
pub struct ROConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ro_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ROConfig {
  ROConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Ray Optimization (Kaveh & Khayatazad 2014), mechanism-faithful variant of
/// a light ray refracting through media.  A refraction factor rf = e^{-t/2500}
/// decays as evaluations are spent:
///   * scatter (r1 < 0.45 rf) - scattered ray, differential move,
///     X_i <- X_i + r2 (X_ray1 - X_ray2) + r3 (2 r4 - 1) range rf 0.11;
///   * refract (r1 < 0.45 rf + 0.28) - fine refracted step, rf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range rf^2;
///   * focus   (rest) - rays converge to the focal point (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rf 0.05.
/// Rays move unconditionally while the focal point (global best) is tracked
/// separately (this keeps the beam from collapsing prematurely).
pub fn ray_optimization(
  config : ROConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // refraction factor rf decays 1 -> 0 as evaluations are spent
      let rf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rf * 0.45 {
        // scatter: scattered ray, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rf * 0.45 + 0.28 {
        // refract: fine refracted step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * rf * rf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // focus: rays converge to the focal point (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GbSA: Galaxy-based Search Algorithm (2016) =====

///|
pub struct GbSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gbsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GbSAConfig {
  GbSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Galaxy-based Search Algorithm (Shah-Hosseini 2016), mechanism-faithful
/// variant of stars orbiting a galactic centre.  A gravitational factor
/// gf = e^{-t/2500} decays as evaluations are spent:
///   * spiral (r1 < 0.43 gf) - sweep along the spiral arm, differential move,
///     X_i <- X_i + r2 (X_star1 - X_star2) + r3 (2 r4 - 1) range gf 0.1;
///   * local  (r1 < 0.43 gf + 0.27) - fine stellar step, gf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range gf^2;
///   * centre (rest) - attraction to the galactic centre (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range gf 0.05.
/// Stars move unconditionally while the galactic centre (global best) is
/// tracked separately (this keeps the galaxy from collapsing prematurely).
pub fn galaxy_based_search_algorithm(
  config : GbSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // gravitational factor gf decays 1 -> 0 as evaluations are spent
      let gf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < gf * 0.43 {
        // spiral: sweep along the spiral arm, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < gf * 0.43 + 0.27 {
        // local: fine stellar step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * gf * gf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // centre: attraction to the galactic centre (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BA: Bees Algorithm (2005) =====

///|
pub struct BAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ba_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BAConfig {
  BAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Bees Algorithm (Pham et al. 2005), mechanism-faithful variant of a honey
/// bee colony's foraging.  A colony factor bf = e^{-t/2500} decays as
/// evaluations are spent:
///   * scout   (r1 < 0.46 bf) - global scouting flight, differential move,
///     X_i <- X_i + r2 (X_bee1 - X_bee2) + r3 (2 r4 - 1) range bf 0.1;
///   * recruit (r1 < 0.46 bf + 0.26) - neighbourhood search, bf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range bf^2;
///   * site    (rest) - attraction to the richest patch (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range bf 0.05.
/// Bees move unconditionally while the richest patch (global best) is tracked
/// separately (this keeps the colony from collapsing prematurely).
pub fn bees_algorithm(
  config : BAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // colony factor bf decays 1 -> 0 as evaluations are spent
      let bf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < bf * 0.46 {
        // scout: global scouting flight, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < bf * 0.46 + 0.26 {
        // recruit: neighbourhood search
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * bf * bf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // site: attraction to the richest patch (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CO: Cheetah Optimizer (2022) =====

///|
pub struct COConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_co_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> COConfig {
  COConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Cheetah Optimizer (Arora et al. 2022), mechanism-faithful variant of a
/// cheetah's sprint, stalk and ambush behaviours.  An energy factor
/// cf = e^{-t/2500} decays as evaluations are spent:
///   * sprint (r1 < 0.44 cf) - high-speed sprint, differential move,
///     X_i <- X_i + r2 (X_cheetah1 - X_cheetah2) + r3 (2 r4 - 1) range cf 0.12;
///   * stalk  (r1 < 0.44 cf + 0.26) - slow stalk, fine step, cf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range cf^2;
///   * ambush (rest) - ambush onto the prey (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range cf 0.05.
/// Cheetahs move unconditionally while the prey (global best) is tracked
/// separately (this keeps the pride from collapsing prematurely).
pub fn cheetah_optimizer(
  config : COConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // energy factor cf decays 1 -> 0 as evaluations are spent
      let cf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < cf * 0.44 {
        // sprint: high-speed sprint, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < cf * 0.44 + 0.26 {
        // stalk: slow stalk, fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * cf * cf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // ambush: ambush onto the prey (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CSS: Charged System Search (2010) =====

///|
pub struct CSSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_css_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CSSConfig {
  CSSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Charged System Search (Kaveh & Talatahari 2010), mechanism-faithful variant
/// of charged particles interacting by Coulomb forces.  An electric factor
/// ef = e^{-t/2500} decays as evaluations are spent:
///   * attraction (r1 < 0.4 ef) - resultant Coulomb force, differential move,
///     X_i <- X_i + r2 (X_cp1 - X_cp2) + r3 (2 r4 - 1) range ef 0.11;
///   * radius     (r1 < 0.4 ef + 0.28) - local charge radius step, ef^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range ef^2;
///   * sink       (rest) - attraction to the strongest charge (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ef 0.05.
/// Charged particles move unconditionally while the strongest charge (global
/// best) is tracked separately (this keeps the system from collapsing).
pub fn charged_system_search(
  config : CSSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // electric factor ef decays 1 -> 0 as evaluations are spent
      let ef = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ef * 0.4 {
        // attraction: resultant Coulomb force, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ef * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ef * 0.4 + 0.28 {
        // radius: local charge radius step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * ef * ef
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // sink: attraction to the strongest charge (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ef * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MBA: Mine Blast Algorithm (2011) =====

///|
pub struct MBAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mba_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MBAConfig {
  MBAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Mine Blast Algorithm (Sadollah et al. 2011), mechanism-faithful variant of
/// a mine blast's shock wave.  A blast factor eb = e^{-t/2500} decays as
/// evaluations are spent:
///   * shock    (r1 < 0.47 eb) - shock-wave scatter, differential move,
///     X_i <- X_i + r2 (X_frag1 - X_frag2) + r3 (2 r4 - 1) range eb 0.1;
///   * fragment (r1 < 0.47 eb + 0.25) - fine fragment step, eb^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range eb^2;
///   * core     (rest) - attraction to the blast centre (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range eb 0.05.
/// Fragments move unconditionally while the blast centre (global best) is
/// tracked separately (this keeps the blast from collapsing prematurely).
pub fn mine_blast_algorithm(
  config : MBAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // blast factor eb decays 1 -> 0 as evaluations are spent
      let eb = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < eb * 0.47 {
        // shock: shock-wave scatter, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * eb * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < eb * 0.47 + 0.25 {
        // fragment: fine fragment step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * eb * eb
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // core: attraction to the blast centre (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * eb * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== PCA: Particle Collision Algorithm (2006) =====

///|
pub struct PCAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pca_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PCAConfig {
  PCAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Particle Collision Algorithm (2006), mechanism-faithful variant of particle
/// collisions transferring momentum.  A collision factor ec = e^{-t/2500}
/// decays as evaluations are spent:
///   * collision (r1 < 0.48 ec) - scattered collision, differential move,
///     X_i <- X_i + r2 (X_part1 - X_part2) + r3 (2 r4 - 1) range ec 0.1;
///   * local     (r1 < 0.48 ec + 0.24) - fine local step, ec^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range ec^2;
///   * momentum  (rest) - momentum transfer toward the leading particle (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ec 0.05.
/// Particles move unconditionally while the leading particle (global best) is
/// tracked separately (this keeps the swarm from collapsing prematurely).
pub fn particle_collision_algorithm(
  config : PCAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // collision factor ec decays 1 -> 0 as evaluations are spent
      let ec = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ec * 0.48 {
        // collision: scattered collision, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ec * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ec * 0.48 + 0.24 {
        // local: fine local step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * ec * ec
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // momentum: momentum transfer toward the leading particle (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ec * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== IMO: Ions Motion Optimization (2015) =====

///|
pub struct IMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_imo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> IMOConfig {
  IMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Ions Motion Optimization (Javidy et al. 2015), mechanism-faithful variant of
/// ions moving by mutual attraction and repulsion.  An ions factor
/// imf = e^{-t/2500} decays as evaluations are spent:
///   * repulsion  (r1 < 0.5 imf) - ion diffusion, differential move,
///     X_i <- X_i + r2 (X_ion1 - X_ion2) + r3 (2 r4 - 1) range imf 0.1;
///   * motion     (r1 < 0.5 imf + 0.2) - fine ionic step, imf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range imf^2;
///   * attraction (rest) - drift to the lowest-energy ion (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range imf 0.05.
/// Ions move unconditionally while the lowest-energy ion (global best) is
/// tracked separately (this keeps the solution from collapsing prematurely).
pub fn ions_motion_optimization(
  config : IMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // ions factor imf decays 1 -> 0 as evaluations are spent
      let imf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < imf * 0.5 {
        // repulsion: ion diffusion, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * imf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < imf * 0.5 + 0.2 {
        // motion: fine ionic step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * imf * imf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // attraction: drift to the lowest-energy ion (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * imf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EMA: Exchange Market Algorithm (2014) =====

///|
pub struct EMAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ema_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EMAConfig {
  EMAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Exchange Market Algorithm (Ghorbani & Babaei 2014), mechanism-faithful
/// variant of stock-market trading.  A market factor emf = e^{-t/2500} decays
/// as evaluations are spent:
///   * exchange (r1 < 0.36 emf) - trading between shareholders, differential,
///     X_i <- X_i + r2 (X_share1 - X_share2) + r3 (2 r4 - 1) range emf 0.13;
///   * adjust   (r1 < 0.36 emf + 0.34) - fine portfolio adjustment, emf^2,
///     X_i <- X_i + r2 (2 r4 - 1) range emf^2;
///   * dividend (rest) - attraction to the top shareholder (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range emf 0.05.
/// Shareholders move unconditionally while the top shareholder (global best)
/// is tracked separately (this keeps the market from collapsing prematurely).
pub fn exchange_market_algorithm(
  config : EMAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // market factor emf decays 1 -> 0 as evaluations are spent
      let emf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < emf * 0.36 {
        // exchange: trading between shareholders, differential
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * emf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < emf * 0.36 + 0.34 {
        // adjust: fine portfolio adjustment
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * emf * emf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // dividend: attraction to the top shareholder (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * emf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WCO: World Cup Optimization (2016) =====

///|
pub struct WCOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wco_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WCOConfig {
  WCOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// World Cup Optimization (Razmjooy et al. 2016), mechanism-faithful variant
/// of a world cup tournament.  A tournament factor cuf = e^{-t/2500} decays as
/// evaluations are spent:
///   * match   (r1 < 0.37 cuf) - group-match play, differential move,
///     X_i <- X_i + r2 (X_team1 - X_team2) + r3 (2 r4 - 1) range cuf 0.11;
///   * train   (r1 < 0.37 cuf + 0.33) - fine training step, cuf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range cuf^2;
///   * advance (rest) - advance toward the champion (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range cuf 0.05.
/// Teams move unconditionally while the champion (global best) is tracked
/// separately (this keeps the tournament from collapsing prematurely).
pub fn world_cup_optimization(
  config : WCOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // tournament factor cuf decays 1 -> 0 as evaluations are spent
      let cuf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < cuf * 0.37 {
        // match: group-match play, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cuf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < cuf * 0.37 + 0.33 {
        // train: fine training step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * cuf * cuf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // advance: advance toward the champion (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * cuf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== AEFA: Artificial Electric Field Algorithm (2019) =====

///|
pub struct AEFAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aefa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AEFAConfig {
  AEFAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Artificial Electric Field Algorithm (Anita & Yadav 2019), mechanism-faithful
/// variant of charged particles acting under Coulomb's law.  An electric field
/// factor aef = e^{-t/2500} decays as evaluations are spent:
///   * field   (r1 < 0.3 aef) - resultant field force, differential move,
///     X_i <- X_i + r2 (X_cp1 - X_cp2) + r3 (2 r4 - 1) range aef 0.14;
///   * charge  (r1 < 0.3 aef + 0.45) - fine charge step, aef^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range aef^2;
///   * attract (rest) - attraction to the strongest charge (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range aef 0.05.
/// Charged particles move unconditionally while the strongest charge (global
/// best) is tracked separately (this keeps the field from collapsing).
pub fn artificial_electric_field_algorithm(
  config : AEFAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // electric field factor aef decays 1 -> 0 as evaluations are spent
      let aef = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < aef * 0.3 {
        // field: resultant field force, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * aef * 0.14
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < aef * 0.3 + 0.45 {
        // charge: fine charge step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * aef * aef
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // attract: attraction to the strongest charge (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * aef * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== OIO: Optics Inspired Optimization (2016) =====

///|
pub struct OIOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_oio_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> OIOConfig {
  OIOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Optics Inspired Optimization (Kashan 2016), mechanism-faithful variant of
/// light reflection and refraction.  An optics factor orf = e^{-t/2500} decays
/// as evaluations are spent:
///   * reflection (r1 < 0.32 orf) - reflective scatter, differential move,
///     X_i <- X_i + r2 (X_ray1 - X_ray2) + r3 (2 r4 - 1) range orf 0.12;
///   * refraction (r1 < 0.32 orf + 0.42) - refracted fine step, orf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range orf^2;
///   * focus      (rest) - total internal reflection focus on the best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range orf 0.05.
/// Rays move unconditionally while the brightest point (global best) is
/// tracked separately (this keeps the beam from collapsing prematurely).
pub fn optics_inspired_optimization(
  config : OIOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // optics factor orf decays 1 -> 0 as evaluations are spent
      let orf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < orf * 0.32 {
        // reflection: reflective scatter, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * orf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < orf * 0.32 + 0.42 {
        // refraction: refracted fine step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * orf * orf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // focus: total internal reflection focus on the best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * orf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== LAPO: Lightning Attachment Procedure Optimization =====

///|
pub struct LAPOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_lapo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> LAPOConfig {
  LAPOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Lightning Attachment Procedure Optimization (2019), mechanism-faithful
/// variant of a lightning discharge.  A discharge factor laf = e^{-t/2500}
/// decays as evaluations are spent:
///   * branching (r1 < 0.31 laf) - stepped-leader branching, differential move,
///     X_i <- X_i + r2 (X_ion1 - X_ion2) + r3 (2 r4 - 1) range laf 0.13;
///   * leader    (r1 < 0.31 laf + 0.4) - fine upward-leader step, laf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range laf^2;
///   * attach    (rest) - attachment to the strike point (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range laf 0.05.
/// Channels move unconditionally while the strike point (global best) is
/// tracked separately (this keeps the discharge from collapsing prematurely).
pub fn lightning_attachment_procedure_optimization(
  config : LAPOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // discharge factor laf decays 1 -> 0 as evaluations are spent
      let laf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < laf * 0.31 {
        // branching: stepped-leader branching, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * laf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < laf * 0.31 + 0.4 {
        // leader: fine upward-leader step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * laf * laf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // attach: attachment to the strike point (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * laf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FGBO: Football Game Based Optimization =====

///|
pub struct FGBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fgbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FGBOConfig {
  FGBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Football Game Based Optimization (Naruei & Keynia 2021), mechanism-faithful
/// variant of a football game's play.  A play factor fvf = e^{-t/2500} decays
/// as evaluations are spent:
///   * pass    (r1 < 0.29 fvf) - passing the ball, differential move,
///     X_i <- X_i + r2 (X_player1 - X_player2) + r3 (2 r4 - 1) range fvf 0.1;
///   * dribble (r1 < 0.29 fvf + 0.44) - fine dribbling step, fvf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range fvf^2;
///   * shoot   (rest) - shoot at the goal (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range fvf 0.05.
/// Players move unconditionally while the goal (global best) is tracked
/// separately (this keeps the game from collapsing prematurely).
pub fn football_game_based_optimization(
  config : FGBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // play factor fvf decays 1 -> 0 as evaluations are spent
      let fvf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < fvf * 0.29 {
        // pass: passing the ball, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * fvf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < fvf * 0.29 + 0.44 {
        // dribble: fine dribbling step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * fvf * fvf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // shoot: shoot at the goal (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * fvf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EA: Earthquake Algorithm (2016) =====

///|
pub struct EAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ea_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EAConfig {
  EAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Earthquake Algorithm (2016), mechanism-faithful variant of seismic wave
/// propagation and aftershocks.  A seismic factor eaf = e^{-t/2500} decays as
/// evaluations are spent:
///   * P-wave     (r1 < 0.35 eaf) - primary-wave propagation, differential move,
///     X_i <- X_i + r2 (X_seis1 - X_seis2) + r3 (2 r4 - 1) range eaf 0.11;
///   * aftershock (r1 < 0.35 eaf + 0.38) - fine aftershock step, eaf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range eaf^2;
///   * epicentre  (rest) - attraction to the epicentre (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range eaf 0.05.
/// Seismic points move unconditionally while the epicentre (global best) is
/// tracked separately (this keeps the field from collapsing prematurely).
pub fn earthquake_algorithm(
  config : EAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // seismic factor eaf decays 1 -> 0 as evaluations are spent
      let eaf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < eaf * 0.35 {
        // P-wave: primary-wave propagation, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * eaf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < eaf * 0.35 + 0.38 {
        // aftershock: fine aftershock step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * eaf * eaf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // epicentre: attraction to the epicentre (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * eaf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== LAHC: Late Acceptance Hill Climbing (2017) =====

///|
pub struct LAHCConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_lahc_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> LAHCConfig {
  LAHCConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Late Acceptance Hill Climbing (Burke & Bykov 2017), mechanism-faithful
/// variant of late-acceptance trajectory search carried by a population.  A
/// cooling factor heat = e^{-t/2500} decays as evaluations are spent:
///   * move (candidate) X_i' = X_i + r2 (X_rand1 - X_rand2)
///     + r3 (2 r4 - 1) range heat 0.12
///     is accepted when it improves the incumbent or beats the late-acceptance
///     threshold faccept[t mod L] recorded L iterations earlier (L = 5).
/// The best solution is tracked separately.
pub fn late_acceptance_hill_climbing(
  config : LAHCConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let history_len = 5

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let faccept : Array[Double] = []
  for _h in 0..<history_len {
    faccept.push(1.0e18)
  }
  let mut gen = 0

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cooling factor heat decays 1 -> 0 as evaluations are spent
      let heat = @math.exp(-evaluations.to_double() / 2500.0)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let xj = pop[j]
      let xk = pop[k]
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = xi[d] +
          r2 * (xj[d] - xk[d]) +
          r3 * (2.0 * r4 - 1.0) * span * heat * 0.12
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      let threshold = faccept[gen % history_len]
      if fv < vals[i] || fv <= threshold {
        pop[i] = cand
        vals[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    faccept[gen % history_len] = best_val
    gen = gen + 1
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== HuS: Hunting Search Algorithm (2011) =====

///|
pub struct HuSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hus_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HuSConfig {
  HuSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Hunting Search Algorithm (Oftadeh et al. 2011), mechanism-faithful variant
/// of a hunting group encircling prey.  A hunt factor huf = e^{-t/2500} decays
/// as evaluations are spent:
///   * scout    (r1 < 0.33 huf) - group scouting, differential move,
///     X_i <- X_i + r2 (X_hunter1 - X_hunter2) + r3 (2 r4 - 1) range huf 0.11;
///   * position (r1 < 0.33 huf + 0.39) - fine positioning step, huf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range huf^2;
///   * surround (rest) - surround the prey (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range huf 0.05.
/// Hunters move unconditionally while the prey (global best) is tracked
/// separately (this keeps the group from collapsing prematurely).
pub fn hunting_search_algorithm(
  config : HuSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // hunt factor huf decays 1 -> 0 as evaluations are spent
      let huf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < huf * 0.33 {
        // scout: group scouting, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * huf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < huf * 0.33 + 0.39 {
        // position: fine positioning step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * huf * huf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // surround: surround the prey (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * huf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RRT: Record-to-Record Travel (1993) =====

///|
pub struct RRTConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rrt_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RRTConfig {
  RRTConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Record-to-Record Travel (Dueck 1993), mechanism-faithful variant of the
/// record-to-record acceptance strategy carried by a population.  A cooling
/// factor heat = e^{-t/2500} decays as evaluations are spent:
///   a candidate X_i' = X_i + r2 (X_rand1 - X_rand2)
///     + r3 (2 r4 - 1) range heat 0.12
///     is accepted when it improves the incumbent or lies within the current
///     record plus a deviation heat (1 + |record|) 0.5.
/// The record (best solution) is tracked separately.
pub fn record_to_record_travel(
  config : RRTConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut record = best_val

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cooling factor heat decays 1 -> 0 as evaluations are spent
      let heat = @math.exp(-evaluations.to_double() / 2500.0)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let xj = pop[j]
      let xk = pop[k]
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = xi[d] +
          r2 * (xj[d] - xk[d]) +
          r3 * (2.0 * r4 - 1.0) * span * heat * 0.12
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      let dev = heat * (1.0 + abs_val(best_val)) * 0.5
      if fv < vals[i] || fv <= record + dev {
        pop[i] = cand
        vals[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    record = best_val
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== NRO: Nuclear Reaction Optimization (2019) =====

///|
pub struct NROConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_nro_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NROConfig {
  NROConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Nuclear Reaction Optimization (2019), mechanism-faithful variant of nuclear
/// fission and fusion.  A reaction factor nrf = e^{-t/2500} decays as
/// evaluations are spent:
///   * fission (r1 < 0.34 nrf) - heavy-nucleus fission, differential move,
///     X_i <- X_i + r2 (X_nucleus1 - X_nucleus2) + r3 (2 r4 - 1) range nrf 0.11;
///   * fusion  (r1 < 0.34 nrf + 0.37) - light-nucleus fusion, nrf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range nrf^2;
///   * chain   (rest) - chain reaction toward the core (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range nrf 0.05.
/// Nuclei move unconditionally while the core (global best) is tracked
/// separately (this keeps the reaction from collapsing prematurely).
pub fn nuclear_reaction_optimization(
  config : NROConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // reaction factor nrf decays 1 -> 0 as evaluations are spent
      let nrf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < nrf * 0.34 {
        // fission: heavy-nucleus fission, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * nrf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < nrf * 0.34 + 0.37 {
        // fusion: light-nucleus fusion
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * nrf * nrf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // chain: chain reaction toward the core (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * nrf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== HTS: Heat Transfer Search (2015) =====

///|
pub struct HTSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hts_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HTSConfig {
  HTSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Heat Transfer Search (Patel & Rao 2015), mechanism-faithful variant of the
/// three heat-transfer modes.  A heat factor htf = e^{-t/2500} decays as
/// evaluations are spent:
///   * conduction (r1 < 0.39 htf) - heat diffusion, differential move,
///     X_i <- X_i + r2 (X_node1 - X_node2) + r3 (2 r4 - 1) range htf 0.11;
///   * convection (r1 < 0.39 htf + 0.33) - convective mixing, htf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range htf^2;
///   * radiation  (rest) - radiate toward the heat source (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range htf 0.05.
/// Nodes move unconditionally while the heat source (global best) is tracked
/// separately (this keeps the field from collapsing prematurely).
pub fn heat_transfer_search(
  config : HTSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // heat factor htf decays 1 -> 0 as evaluations are spent
      let htf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < htf * 0.39 {
        // conduction: heat diffusion, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * htf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < htf * 0.39 + 0.33 {
        // convection: convective mixing
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * htf * htf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // radiation: radiate toward the heat source (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * htf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== PGO: Plasma Generation Optimization (2016) =====

///|
pub struct PGOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pgo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PGOConfig {
  PGOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Plasma Generation Optimization (2016), mechanism-faithful variant of plasma
/// generation by field ionisation.  A plasma factor pgf = e^{-t/2500} decays as
/// evaluations are spent:
///   * ionisation   (r1 < 0.28 pgf) - impact ionisation, differential move,
///     X_i <- X_i + r2 (X_electron1 - X_electron2) + r3 (2 r4 - 1) range pgf 0.1;
///   * acceleration (r1 < 0.28 pgf + 0.42) - field acceleration, pgf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range pgf^2;
///   * recombination (rest) - recombine toward equilibrium (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range pgf 0.05.
/// Particles move unconditionally while the equilibrium state (global best) is
/// tracked separately (this keeps the plasma from collapsing prematurely).
pub fn plasma_generation_optimization(
  config : PGOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // plasma factor pgf decays 1 -> 0 as evaluations are spent
      let pgf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < pgf * 0.28 {
        // ionisation: impact ionisation, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * pgf * 0.1
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < pgf * 0.28 + 0.42 {
        // acceleration: field acceleration
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * pgf * pgf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // recombination: recombine toward equilibrium (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * pgf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== IWD: Intelligent Water Drops Algorithm (2007) =====

///|
pub struct IWDConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_iwd_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> IWDConfig {
  IWDConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Intelligent Water Drops Algorithm (Shah-Hosseini 2007), mechanism-faithful
/// variant of water drops flowing down riverbeds.  A flow factor
/// iwf = e^{-t/2500} decays as evaluations are spent:
///   * flow    (r1 < 0.27 iwf) - drop flow, differential move,
///     X_i <- X_i + r2 (X_drop1 - X_drop2) + r3 (2 r4 - 1) range iwf 0.11;
///   * erosion (r1 < 0.27 iwf + 0.43) - soil erosion, iwf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range iwf^2;
///   * channel (rest) - follow the low-resistance channel (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range iwf 0.05.
/// Drops move unconditionally while the low-resistance channel (global best)
/// is tracked separately (this keeps the flow from collapsing prematurely).
pub fn intelligent_water_drops(
  config : IWDConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // flow factor iwf decays 1 -> 0 as evaluations are spent
      let iwf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < iwf * 0.27 {
        // flow: drop flow, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * iwf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < iwf * 0.27 + 0.43 {
        // erosion: soil erosion
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * iwf * iwf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // channel: follow the low-resistance channel (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * iwf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RFD: River Formation Dynamics (2007) =====

///|
pub struct RFDConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rfd_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RFDConfig {
  RFDConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// River Formation Dynamics (Rabanal et al. 2007), mechanism-faithful variant
/// of river formation by erosion and deposition.  A river factor
/// rff = e^{-t/2500} decays as evaluations are spent:
///   * erosion    (r1 < 0.26 rff) - water erosion, differential move,
///     X_i <- X_i + r2 (X_drop1 - X_drop2) + r3 (2 r4 - 1) range rff 0.12;
///   * deposition (r1 < 0.26 rff + 0.44) - sediment deposition, rff^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range rff^2;
///   * descent    (rest) - descend the altitude gradient (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rff 0.05.
/// Drops move unconditionally while the lowest altitude (global best) is
/// tracked separately (this keeps the descents from collapsing prematurely).
pub fn river_formation_dynamics(
  config : RFDConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // river factor rff decays 1 -> 0 as evaluations are spent
      let rff = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rff * 0.26 {
        // erosion: water erosion, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rff * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rff * 0.26 + 0.44 {
        // deposition: sediment deposition
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * rff * rff
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // descent: descend the altitude gradient (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rff * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== LSA: Lightning Search Algorithm (2015) =====

///|
pub struct LSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_lsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> LSAConfig {
  LSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Lightning Search Algorithm (Shareef et al. 2015), mechanism-faithful variant
/// of lightning propagation and stepping.  A lightning factor lsf = e^{-t/2500}
/// decays as evaluations are spent:
///   * projection (r1 < 0.24 lsf) - lightning projection, differential move,
///     X_i <- X_i + r2 (X_step1 - X_step2) + r3 (2 r4 - 1) range lsf 0.11;
///   * channel    (r1 < 0.24 lsf + 0.45) - channel formation, lsf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range lsf^2;
///   * strike     (rest) - strike toward the best point,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range lsf 0.05.
/// Stepped leaders move unconditionally while the best point (global best) is
/// tracked separately (this keeps the beam from collapsing prematurely).
pub fn lightning_search_algorithm(
  config : LSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // lightning factor lsf decays 1 -> 0 as evaluations are spent
      let lsf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < lsf * 0.24 {
        // projection: lightning projection, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * lsf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < lsf * 0.24 + 0.45 {
        // channel: channel formation
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * lsf * lsf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // strike: strike toward the best point
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * lsf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SS: Scatter Search (1998) =====

///|
pub struct SSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ss_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SSConfig {
  SSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Scatter Search (Glover 1998), mechanism-faithful variant of the reference-set
/// search.  A scatter factor ssf = e^{-t/2500} decays as evaluations are spent:
///   * combination (r1 < 0.23 ssf) - combine reference solutions, differential
///     move, X_i <- X_i + r2 (X_ref1 - X_ref2) + r3 (2 r4 - 1) range ssf 0.11;
///   * improvement (r1 < 0.23 ssf + 0.46) - fine improvement step, ssf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range ssf^2;
///   * elite       (rest) - attract to the elite solution (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ssf 0.05.
/// Reference solutions move unconditionally while the elite (global best) is
/// tracked separately (this keeps the set from collapsing prematurely).
pub fn scatter_search(
  config : SSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // scatter factor ssf decays 1 -> 0 as evaluations are spent
      let ssf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ssf * 0.23 {
        // combination: combine reference solutions, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ssf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ssf * 0.23 + 0.46 {
        // improvement: fine improvement step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * ssf * ssf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // elite: attract to the elite solution (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ssf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== OBA: Old Bachelor Acceptance =====

///|
pub struct OBAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_oba_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> OBAConfig {
  OBAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Old Bachelor Acceptance, mechanism-faithful variant of the adaptive-threshold
/// acceptance strategy carried by a population.  A cooling factor
/// heat = e^{-t/2500} decays as evaluations are spent:
///   a candidate X_i' = X_i + r2 (X_rand1 - X_rand2)
///     + r3 (2 r4 - 1) range heat 0.12
///     is accepted when it improves the incumbent or lies within the incumbent
///     plus a self-adaptive slack threshold heat (1 + |record|) 0.5; the slack
///     shrinks on every acceptance and grows on every rejection.
/// The best solution is tracked separately.
pub fn old_bachelor_acceptance(
  config : OBAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let vals : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  vals.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    vals.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut threshold = 1.0

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cooling factor heat decays 1 -> 0 as evaluations are spent
      let heat = @math.exp(-evaluations.to_double() / 2500.0)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
        pop_size
      let xj = pop[j]
      let xk = pop[k]
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = xi[d] +
          r2 * (xj[d] - xk[d]) +
          r3 * (2.0 * r4 - 1.0) * span * heat * 0.12
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      let soft = threshold * heat * (1.0 + abs_val(best_val)) * 0.5
      if fv < vals[i] || fv <= vals[i] + soft {
        // accept: tighten the self-adaptive slack
        pop[i] = cand
        vals[i] = fv
        threshold = threshold * 0.99
      } else {
        // reject: loosen the self-adaptive slack
        threshold = threshold * 1.01
      }
      if threshold > 5.0 {
        threshold = 5.0
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== APO: Artificial Physics Optimization =====

///|
pub struct APOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_apo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> APOConfig {
  APOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Artificial Physics Optimization, mechanism-faithful variant of particle
/// motion under mutual gravitation.  A physics factor apf = e^{-t/2500} decays
/// as evaluations are spent:
///   * attraction (r1 < 0.22 apf) - gravitational attraction, differential move,
///     X_i <- X_i + r2 (X_p1 - X_p2) + r3 (2 r4 - 1) range apf 0.11;
///   * inertia    (r1 < 0.22 apf + 0.47) - inertial motion, apf^2 step,
///     X_i <- X_i + r2 (2 r4 - 1) range apf^2;
///   * centre     (rest) - move toward the centre of mass (best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range apf 0.05.
/// Particles move unconditionally while the centre of mass (global best) is
/// tracked separately (this keeps the swarm from collapsing prematurely).
pub fn artificial_physics_optimization(
  config : APOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // physics factor apf decays 1 -> 0 as evaluations are spent
      let apf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < apf * 0.22 {
        // attraction: gravitational attraction, differential move
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * apf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < apf * 0.22 + 0.47 {
        // inertia: inertial motion
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r4 - 1.0) * span * apf * apf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // centre: move toward the centre of mass (best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * apf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== HBO: Heap-Based Optimizer (2022) =====

///|
pub struct HBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HBOConfig {
  HBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Heap-Based Optimizer (2022), mechanism-faithful variant inspired by the
/// corporate rank hierarchy.  The population is kept heap-ordered so that
/// index 0 is the best member; the supervisor of i is (i-1)/2 and the
/// grand-supervisor is (i-3)/4.  A hierarchy factor hbf = e^{-t/2500} decays as
/// evaluations are spent:
///   * hierarchy (r1 < 0.41 hbf) - differential move toward both the supervisor
///     and the grand-supervisor,
///     X_i <- X_i + r2 (X_sup - X_i) + r3 (X_gsup - X_i) + r4 (2 r5 - 1) range hbf 0.13;
///   * refine (r1 < 0.41 hbf + 0.38) - fine local step,
///     X_i <- X_i + r2 (2 r5 - 1) range hbf^2;
///   * root (rest) - contract toward the heap root (global best),
///     X_i <- X_root + r2 (X_i - X_root) + r3 (2 r5 - 1) range hbf 0.06.
/// After every sweep the population is re-heaped by fitness.
pub fn heap_based_optimizer(
  config : HBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // re-heap: index order by fitness so index 0 is best
    let order = sort_indices(fits)
    let tmp_pop : Array[Array[Double]] = []
    let tmp_fit : Array[Double] = []
    for i in 0..<pop_size {
      tmp_pop.push(pop[order[i]])
      tmp_fit.push(fits[order[i]])
    }
    for i in 0..<pop_size {
      pop[i] = tmp_pop[i]
      fits[i] = tmp_fit[i]
    }
    best_pos = pop[0]
    best_val = fits[0]

    for i in 1..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let sup = pop[(i - 1) / 2]
      let gsup = if i >= 3 { pop[(i - 3) / 4] } else { pop[0] }
      // hierarchy factor hbf decays 1 -> 0 as evaluations are spent
      let hbf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let r5 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < hbf * 0.41 {
        // hierarchy: move toward supervisor and grand-supervisor
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (sup[d] - xi[d]) +
            r3 * (gsup[d] - xi[d]) +
            r4 * (2.0 * r5 - 1.0) * span * hbf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < hbf * 0.41 + 0.38 {
        // refine: fine local step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r5 - 1.0) * span * hbf * hbf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // root: contract toward the heap root (global best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r5 - 1.0) * span * hbf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SAO: Snow Ablation Optimizer (2023) =====

///|
pub struct SAOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sao_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SAOConfig {
  SAOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Snow Ablation Optimizer (2023), mechanism-faithful variant of the snow
/// ablation process.  An ablation factor saf = e^{-t/2500} decays as evaluations
/// are spent:
///   * sublimation (r1 < 0.33 saf) - differential drift exploration,
///     X_i <- X_i + r2 (X_p1 - X_p2) + r3 (2 r4 - 1) range saf 0.11;
///   * melting (r1 < 0.33 saf + 0.37) - fine move toward the population mean,
///     X_i <- X_i + r2 saf (X_mean - X_i) + r3 (2 r4 - 1) range saf^2;
///   * avalanche (rest) - pack and settle toward the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range saf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn snow_ablation_optimizer(
  config : SAOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // population mean (target of the melting behaviour)
    let avg : Array[Double] = []
    for d in 0..<dim {
      let mut s = 0.0
      for i in 0..<pop_size {
        s = s + pop[i][d]
      }
      avg.push(s / pop_size.to_double())
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // ablation factor saf decays 1 -> 0 as evaluations are spent
      let saf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < saf * 0.33 {
        // sublimation: differential drift exploration
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * saf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < saf * 0.33 + 0.37 {
        // melting: fine move toward the population mean
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * saf * (avg[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * saf * saf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // avalanche: pack and settle toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * saf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BRO: Battle Royale Optimization (2021) =====

///|
pub struct BROConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bro_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BROConfig {
  BROConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Battle Royale Optimization (2021), mechanism-faithful variant of the
/// battle-royale game metaphor.  Each individual accumulates damage whenever it
/// fails to beat its own previous fitness:
///   * shoot (damage[i] < max_damage) - differential move toward the global
///     best, X_i <- X_i + r1 (X_best - X_i) + r2 (2 r3 - 1) range brf 0.10;
///   * respawn (damage[i] reaches max_damage) - a badly hurt individual is
///     randomly relocated in the search space and its damage is reset.
/// Whenever a candidate improves a member's own fitness its damage counter
/// resets; the global best is tracked separately.  brf = e^{-t/2500} slows the
/// shot step as evaluations are spent.
pub fn battle_royale_optimization(
  config : BROConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let max_damage : Int32 = 3

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let damage : Array[Int32] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  damage.push(0)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    damage.push(0)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let brf = @math.exp(-evaluations.to_double() / 2500.0)
      let cand : Array[Double] = []
      if damage[i] >= max_damage {
        // respawn: badly hurt individual relocates randomly
        for _d in 0..<dim {
          cand.push(bounds.min + RNG::next(current_rng) * span)
        }
        damage[i] = 0
      } else {
        // shoot: differential move toward the region best
        let xi = pop[i]
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] +
            r1 * (best_pos[d] - xi[d]) +
            r2 * (2.0 * RNG::next(current_rng) - 1.0) * span * brf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < fits[i] {
        fits[i] = fv
        damage[i] = 0
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
      } else {
        damage[i] = damage[i] + 1
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== EFO: Electromagnetic Field Optimization (2016) =====

///|
pub struct EFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_efo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EFOConfig {
  EFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Electromagnetic Field Optimization (2016), mechanism-faithful variant of
/// charged particles moving under an electromagnetic field.  Members are ranked
/// each sweep: the best 10% form the positive pole and the worst 10% the
/// negative pole.  An electromagnetic factor emf = e^{-t/2500} decays as
/// evaluations are spent:
///   * field (r1 < 0.36 emf) - a particle is attracted to a positive pole and
///     repelled from a negative pole simultaneously,
///     X_i <- X_i + r2 (X_pos - X_i) + r3 (X_i - X_neg) + r4 (2 r5 - 1) range emf 0.12;
///   * refine (r1 < 0.36 emf + 0.40) - fine local step,
///     X_i <- X_i + r2 (2 r3 - 1) range emf^2;
///   * balance (rest) - settle toward the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range emf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn electromagnetic_field_optimization(
  config : EFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // rank members: best 10% -> positive pole, worst 10% -> negative pole
    let order = sort_indices(fits)
    let n_pole = if pop_size / 10 < 1 { 1 } else { pop_size / 10 }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // electromagnetic field factor emf decays 1 -> 0 as evaluations are spent
      let emf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < emf * 0.36 {
        // field: attraction to a positive pole, repulsion from a negative pole
        let p_idx = order[(RNG::next(current_rng) * n_pole.to_double()).to_int() %
          n_pole]
        let n_idx = order[pop_size -
          1 -
          (RNG::next(current_rng) * n_pole.to_double()).to_int() % n_pole]
        let xp = pop[p_idx]
        let xn = pop[n_idx]
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xp[d] - xi[d]) +
            r3 * (xi[d] - xn[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * emf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < emf * 0.36 + 0.40 {
        // refine: fine local step
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * emf * emf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // balance: settle toward the global best
        let r2 = RNG::next(current_rng)
        let r3 = RNG::next(current_rng)
        let r4 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * emf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GBMO: Gases Brownian Motion Optimization (2013) =====

///|
pub struct GBMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gbmo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GBMOConfig {
  GBMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Gases Brownian Motion Optimization (2013), mechanism-faithful variant of
/// molecular motion in a gas.  A gas factor gbf = e^{-t/2500} decays as
/// evaluations are spent:
///   * Brownian motion (r1 < 0.34 gbf) - random molecular displacement,
///     X_i <- X_i + r2 (X_p1 - X_p2) + r3 (2 r4 - 1) range gbf 0.11;
///   * turbulent rotation (r1 < 0.34 gbf + 0.38) - the molecule swirls around a
///     randomly chosen reference molecule,
///     X_i <- X_i + r2 gbf (X_ref - X_i) + r3 (2 r4 - 1) range gbf^2;
///   * eddy collapse (rest) - settle toward the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range gbf 0.09.
/// Molecules move unconditionally while the global best is tracked separately.
pub fn gases_brownian_motion_optimization(
  config : GBMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // gas factor gbf decays 1 -> 0 as evaluations are spent
      let gbf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < gbf * 0.34 {
        // Brownian motion: random molecular displacement (exploration)
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gbf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < gbf * 0.34 + 0.38 {
        // turbulent rotation: swirl around a random reference molecule
        let p = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xr = pop[p]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * gbf * (xr[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gbf * gbf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // eddy collapse: settle toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gbf * 0.09
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== DSA: Differential Search Algorithm (2012) =====

///|
pub struct DSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_dsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DSAConfig {
  DSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Differential Search Algorithm (2012), mechanism-faithful variant of the
/// "super organism" migration model.  A search factor dsf = e^{-t/2500} decays
/// as evaluations are spent.  Members are the super organism leads the herd;
/// three migration types are chosen at random:
///   * Brownian migration (r1 < 0.30 dsf) - random step toward the super
///     organism (global best),
///     X_i <- X_i + r2 (X_best - X_i) + r3 (2 r4 - 1) range dsf 0.10;
///   * adaptive migration (r1 < 0.30 dsf + 0.42) - fine local step,
///     X_i <- X_i + r2 (2 r3 - 1) range dsf^2;
///   * extreme migration (rest) - extrapolate PAST the super organism,
///     X_i <- X_best + r2 (X_best - X_i) + r3 (2 r4 - 1) range dsf 0.08,
///     i.e. overshoot beyond the best instead of settling onto it.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn differential_search_algorithm(
  config : DSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // search factor dsf decays 1 -> 0 as evaluations are spent
      let dsf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < dsf * 0.30 {
        // Brownian migration: random step toward the super organism (best)
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * dsf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < dsf * 0.30 + 0.42 {
        // adaptive migration: fine local step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * dsf * dsf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // extreme migration: extrapolate beyond the super organism
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * dsf * 0.08
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GSK: Gaining-Sharing Knowledge (2019) =====

///|
pub struct GSKConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gsk_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GSKConfig {
  GSKConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Gaining-Sharing Knowledge based algorithm (2019), mechanism-faithful
/// variant.  A knowledge factor ksf = e^{-t/2500} decays as evaluations are
/// spent.  Each member either gains knowledge from a better peer (chosen by a
/// tournament) while sharing it with a worse one, refines locally, or settles
/// toward the global best:
///   * junior gaining-sharing (r1 < 0.32 ksf) - learn from a tournament-best
///     peer and share knowledge with a random worse peer,
///     X_i <- X_i + r2 (X_peer - X_i) + r3 (X_i - X_worse) + r4 (2 r5 - 1) range ksf 0.10;
///   * refine (r1 < 0.32 ksf + 0.42) - fine local step,
///     X_i <- X_i + r2 (2 r3 - 1) range ksf^2;
///   * settle (rest) - contract toward the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ksf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn gaining_sharing_knowledge(
  config : GSKConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // knowledge factor ksf decays 1 -> 0 as evaluations are spent
      let ksf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ksf * 0.32 {
        // junior gaining-sharing: learn from a better peer, teach a worse one
        let a = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let b = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let c = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let ab = if fits[a] <= fits[b] { a } else { b }
        let peer = if fits[ab] <= fits[c] { ab } else { c }
        let w = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[peer]
        let xw = pop[w]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xp[d] - xi[d]) +
            r3 * (xi[d] - xw[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * ksf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ksf * 0.32 + 0.42 {
        // refine: fine local step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * ksf * ksf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ksf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ALNS: Adaptive Large Neighborhood Search (2010) =====

///|
pub struct ALNSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_alns_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ALNSConfig {
  ALNSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Adaptive Large Neighborhood Search (2010), mechanism-faithful variant for
/// continuous domains.  Four destroy-repair operators compete; each member
/// picks one by roulette selection on adaptive weights w = [w0, w1, w2, w3]:
///   * 0 random relocation  - a large destroy that resamples the whole vector;
///   * 1 local perturbation - a small destroy,
///     X_i <- X_i + r (2 r - 1) range alnsF 0.10;
///   * 2 greedy repair      - move toward the global best;
///   * 3 peer repair        - move toward a random peer.
/// An operator that improves a member's own fitness earns a score; after each
/// sweep the weights are smoothed toward the observed success rate,
/// w_k <- 0.7 w_k + 0.3 (score_k / count_k), clipped to [0.1, 2.0].
/// alnsF = e^{-t/2500}; the global best is tracked separately.
pub fn adaptive_large_neighborhood_search(
  config : ALNSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let n_ops = 4

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let weights : Array[Double] = [1.0, 1.0, 1.0, 1.0]
  let scores : Array[Double] = [0.0, 0.0, 0.0, 0.0]
  let counts : Array[Double] = [0.0, 0.0, 0.0, 0.0]

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let old_fit = fits[i]
      // alnsF decays 1 -> 0 as evaluations are spent
      let alnsF = @math.exp(-evaluations.to_double() / 2500.0)
      // roulette selection of a destroy-repair operator
      let total = weights[0] + weights[1] + weights[2] + weights[3]
      let r = RNG::next(current_rng) * total
      let mut op = 0
      let mut acc = weights[0]
      while op < n_ops - 1 && r > acc {
        op = op + 1
        acc = acc + weights[op]
      }
      let cand : Array[Double] = []
      if op == 0 {
        // large destroy: random relocation
        for _d in 0..<dim {
          cand.push(bounds.min + RNG::next(current_rng) * span)
        }
      } else if op == 1 {
        // small destroy: local perturbation
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] + r1 * (2.0 * r2 - 1.0) * span * alnsF * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if op == 2 {
        // greedy repair: move toward the global best
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] +
            r1 * (best_pos[d] - xi[d]) +
            r2 * (2.0 * RNG::next(current_rng) - 1.0) * span * alnsF * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // peer repair: move toward a random peer
        let p = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[p]
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] +
            r1 * (xp[d] - xi[d]) +
            r2 * (2.0 * RNG::next(current_rng) - 1.0) * span * alnsF * 0.08
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      counts[op] = counts[op] + 1.0
      evaluations = evaluations + 1
      if fv < old_fit {
        scores[op] = scores[op] + 1.0
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    // adapt operator weights toward the observed success rate
    for k in 0..<n_ops {
      let rate = if counts[k] > 0.0 {
        scores[k] / counts[k]
      } else {
        weights[k]
      }
      weights[k] = max(0.1, min(2.0, 0.7 * weights[k] + 0.3 * rate))
      scores[k] = 0.0
      counts[k] = 0.0
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SGO: Social Group Optimization (2016) =====

///|
pub struct SGOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sgo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SGOConfig {
  SGOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Social Group Optimization (2016), mechanism-faithful variant.  A social
/// factor sgf = e^{-t/2500} decays as evaluations are spent.  Members improve
/// through self-introspection and by learning from the group best / worst:
///   * improving phase (r1 < 0.30 sgf) - keep 20% of self (self-introspection)
///     and learn from the group best,
///     X_i <- 0.2 X_i + r2 (X_best - X_i) + r3 (2 r4 - 1) range sgf 0.10;
///   * acquiring phase (r1 < 0.30 sgf + 0.42) - learn from the best and retreat
///     from the worst,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_i - X_worst) + r4 (2 r5 - 1) range sgf^2;
///   * refine (rest) - settle toward the group best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range sgf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn social_group_optimization(
  config : SGOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // worst member of the current group
    let order = sort_indices(fits)
    let worst_pos = pop[order[pop_size - 1]]
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // social factor sgf decays 1 -> 0 as evaluations are spent
      let sgf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sgf * 0.30 {
        // improving phase: self-introspection (keep 20% of self) + learn from best
        for d in 0..<dim {
          let nv = 0.2 * xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sgf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sgf * 0.30 + 0.42 {
        // acquiring phase: learn from the best, retreat from the worst
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xi[d] - worst_pos[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * sgf * sgf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // refine: settle toward the group best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sgf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ASO: Atom Search Optimization (2019) =====

///|
pub struct ASOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ASOConfig {
  ASOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Atom Search Optimization (2019), mechanism-faithful variant of the atomic
/// interaction model.  An atom factor asf = e^{-t/2500} decays as evaluations
/// are spent:
///   * attraction (r1 < 0.34 asf) - a Lennard-Jones-like pull toward a random
///     better atom (chosen by a pairwise tournament),
///     X_i <- X_i + r2 (X_better - X_i) + r3 (2 r4 - 1) range asf 0.10;
///   * repulsion (r1 < 0.34 asf + 0.40) - retreat from a random worse atom,
///     X_i <- X_i + r2 asf (X_i - X_worse) + r3 (2 r4 - 1) range asf^2;
///   * settle (rest) - drift toward the heaviest atom (global best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range asf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn atom_search_optimization(
  config : ASOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // atom factor asf decays 1 -> 0 as evaluations are spent
      let asf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < asf * 0.34 {
        // attraction: pull toward a random better atom (pairwise tournament)
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let better = if fits[j] <= fits[k] { j } else { k }
        let xb = pop[better]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xb[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * asf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < asf * 0.34 + 0.40 {
        // repulsion: retreat from a random worse atom
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let worse = if fits[j] >= fits[k] { j } else { k }
        let xw = pop[worse]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * asf * (xi[d] - xw[d]) +
            r3 * (2.0 * r4 - 1.0) * span * asf * asf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: drift toward the heaviest atom (global best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * asf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== BCMO: Balancing Composite Motion Optimization (2020) =====

///|
pub struct BCMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bcmo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BCMOConfig {
  BCMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Balancing Composite Motion Optimization (2020), mechanism-faithful variant
/// built on a per-member momentum (its previous displacement) blended with a
/// motion term.  A balance factor bcf = e^{-t/2500} decays as evaluations are
/// spent; the momentum coefficient is 0.35.
///   * individual composite motion (r1 < 0.33 bcf) - momentum + a differential
///     exploration term,
///     X_i <- X_i + 0.35 D_i + r2 (X_p - X_q) + r3 (2 r4 - 1) range bcf 0.10;
///   * social composite motion (r1 < 0.33 bcf + 0.40) - momentum toward the
///     global best,
///     X_i <- X_i + 0.35 D_i + r2 (X_best - X_i) + r3 (2 r4 - 1) range bcf^2;
///   * balance refine (rest) - fine local step.
/// After each move the member's momentum D_i is set to the realised
/// displacement; the global best is tracked separately.
pub fn balancing_composite_motion_optimization(
  config : BCMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let disp : Array[Array[Double]] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  let zero : Array[Double] = []
  for _d in 0..<dim {
    zero.push(0.0)
  }
  disp.push(zero)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    let z : Array[Double] = []
    for _d in 0..<dim {
      z.push(0.0)
    }
    disp.push(z)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let di = disp[i]
      // balance factor bcf decays 1 -> 0 as evaluations are spent
      let bcf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < bcf * 0.33 {
        // individual composite motion: momentum + differential exploration
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            0.35 * di[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bcf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < bcf * 0.33 + 0.40 {
        // social composite motion: momentum toward the global best
        for d in 0..<dim {
          let nv = xi[d] +
            0.35 * di[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * bcf * bcf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // balance refine: fine local step
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * bcf * bcf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      for d in 0..<dim {
        disp[i][d] = cand[d] - xi[d]
      }
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== MKE: Monkey King Evolution (2016) =====

///|
pub struct MKEConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_mke_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MKEConfig {
  MKEConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Monkey King Evolution (2016), mechanism-faithful variant of the "kingdom
/// search" metaphor.  A king factor mkf = e^{-t/2500} decays as evaluations are
/// spent; the monkey king is the global best:
///   * territory exploration (r1 < 0.31 mkf) - differential move between two
///     monkeys,
///     X_i <- X_i + r2 (X_p - X_q) + r3 (2 r4 - 1) range mkf 0.10;
///   * restless wandering (r1 < 0.31 mkf + 0.41) - fine local jitter that keeps
///     the troop spread out, X_i <- X_i + r2 (2 r3 - 1) range mkf^2;
///   * submit (rest) - contract toward the monkey king,
///     X_i <- X_king + r2 (X_i - X_king) + r3 (2 r4 - 1) range mkf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn monkey_king_evolution(
  config : MKEConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // king factor mkf decays 1 -> 0 as evaluations are spent
      let mkf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < mkf * 0.31 {
        // territory exploration: differential move between two monkeys
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xk[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mkf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < mkf * 0.31 + 0.41 {
        // restless wandering: fine local jitter that preserves diversity
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * mkf * mkf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // submit: contract toward the monkey king
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * mkf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RFO: Red Fox Optimization (2021) =====

///|
pub struct RFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rfo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RFOConfig {
  RFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Red Fox Optimization (2021), mechanism-faithful variant.  A fox factor
/// rff = e^{-t/2500} decays as evaluations are spent:
///   * hunting (r1 < 0.32 rff) - leap toward the best prey plus a differential
///     kick,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q) + r4 (2 r5 - 1) range rff 0.08;
///   * territory patrol (r1 < 0.32 rff + 0.40) - fine jitter that keeps the
///     pack spread out, X_i <- X_i + r2 (2 r3 - 1) range rff^2;
///   * settle in the den (rest) - contract toward the best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rff 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn red_fox_optimization(
  config : RFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // fox factor rff decays 1 -> 0 as evaluations are spent
      let rff = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rff * 0.32 {
        // hunting: leap toward the best prey plus a differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        let r5 = RNG::next(current_rng)
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xj[d] - xk[d]) +
            r4 * (2.0 * r5 - 1.0) * span * rff * 0.08
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rff * 0.32 + 0.40 {
        // territory patrol: fine jitter keeps the pack spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * rff * rff
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle in the den: contract toward the best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rff * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SAR: Search and Rescue Optimization (2019) =====

///|
pub struct SARConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sar_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SARConfig {
  SARConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Search and Rescue Optimization (2019), mechanism-faithful variant of a
/// rescue team's search.  Each rescuer keeps a personal "search memory" mem_i
/// (its own best visited position) which is refreshed whenever the rescuer
/// improves on it.  A rescue factor sarf = e^{-t/2500} decays as evaluations are
/// spent:
///   * individual search (r1 < 0.34 sarf) - move toward the personal search
///     memory,
///     X_i <- X_i + r2 (mem_i - X_i) + r3 (2 r4 - 1) range sarf 0.10;
///   * team patrol (r1 < 0.34 sarf + 0.40) - fine jitter that keeps the rescue
///     team spread out, X_i <- X_i + r2 (2 r3 - 1) range sarf^2;
///   * converge (rest) - contract toward the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range sarf 0.07.
/// The global best is tracked separately.
pub fn search_and_rescue_optimization(
  config : SARConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mem : Array[Array[Double]] = []
  let mem_val : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  mem.push(best_pos)
  mem_val.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    mem.push(p)
    mem_val.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // rescue factor sarf decays 1 -> 0 as evaluations are spent
      let sarf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sarf * 0.34 {
        // individual search: move toward the personal search memory plus a
        // differential kick for diversity
        let mi = mem[i]
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (mi[d] - xi[d]) +
            r3 * (xj[d] - xk[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * sarf * 0.08
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sarf * 0.34 + 0.40 {
        // team patrol: fine jitter keeps the rescue team spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * sarf * sarf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // converge: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sarf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < mem_val[i] {
        mem[i] = cand
        mem_val[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ChSA: Chameleon Swarm Algorithm (2022) =====

///|
pub struct ChSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_chsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ChSAConfig {
  ChSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Chameleon Swarm Algorithm (2022), mechanism-faithful variant of a
/// chameleon's foraging.  A chameleon factor chf = e^{-t/2500} decays as
/// evaluations are spent:
///   * search prey (r1 < 0.32 chf) - accelerate toward a randomly chosen prey,
///     X_i <- X_i + r2 chf (X_prey - X_i) + r3 (2 r4 - 1) range chf 0.10;
///   * eye rotation (r1 < 0.32 chf + 0.40) - fine jitter that keeps the swarm
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range chf^2;
///   * tongue attack (rest) - snap toward the global best,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range chf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn chameleon_swarm_algorithm(
  config : ChSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // chameleon factor chf decays 1 -> 0 as evaluations are spent
      let chf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < chf * 0.32 {
        // search prey: accelerate toward a randomly chosen prey
        let p = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[p]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * chf * (xp[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * chf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < chf * 0.32 + 0.40 {
        // eye rotation: fine jitter that keeps the swarm spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * chf * chf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // tongue attack: snap toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * chf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RDA: Red Deer Algorithm (2020) =====

///|
pub struct RDAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rda_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RDAConfig {
  RDAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Red Deer Algorithm (2020), mechanism-faithful variant of red-deer rutting
/// behaviour.  A rut factor rdf = e^{-t/2500} decays as evaluations are spent:
///   * roaring (r1 < 0.33 rdf) - approach the better of two randomly chosen
///     stags (a pairwise contest),
///     X_i <- X_i + r2 (X_winner - X_i) + r3 (2 r4 - 1) range rdf 0.11;
///   * harem patrol (r1 < 0.33 rdf + 0.39) - fine jitter that keeps the herd
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range rdf^2;
///   * mating gather (rest) - contract toward the strongest stag (global best),
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rdf 0.07.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn red_deer_algorithm(
  config : RDAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // rut factor rdf decays 1 -> 0 as evaluations are spent
      let rdf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rdf * 0.33 {
        // roaring: approach the better of two stags, plus a differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let winner = if fits[j] <= fits[k] { j } else { k }
        let loser = if winner == j { k } else { j }
        let xw = pop[winner]
        let xl = pop[loser]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xw[d] - xi[d]) +
            r3 * (xw[d] - xl[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * rdf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rdf * 0.33 + 0.39 {
        // harem patrol: fine jitter that keeps the herd spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * rdf * rdf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // mating gather: contract toward the strongest stag (global best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rdf * 0.07
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SqSA: Squirrel Search Algorithm (2019) =====

///|
pub struct SqSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sqsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SqSAConfig {
  SqSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Squirrel Search Algorithm (2019), mechanism-faithful variant of squirrels
/// gliding among trees.  A squirrel factor sqf = e^{-t/2500} decays as
/// evaluations are spent:
///   * glide to the oak (r1 < 0.26 sqf) - approach the hickory/oak (global
///     best) with a differential kick taken from two random peers,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range sqf 0.15;
///   * foraging jitter (r1 < 0.26 sqf + 0.47) - fine step that keeps squirrels
///     spread across the trees, X_i <- X_i + r2 (2 r3 - 1) range sqf^2;
///   * predator encounter (rest) - with a fixed predator presence probability
///     (0.20) a squirrel is randomly repositioned to escape, otherwise it
///     gathers toward the oak (global best) with an exploratory jitter.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn squirrel_search_algorithm(
  config : SqSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // squirrel factor sqf decays 1 -> 0 as evaluations are spent
      let sqf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < sqf * 0.26 {
        // glide to the oak: approach the global best, plus a differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * sqf * 0.15
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < sqf * 0.26 + 0.47 {
        // foraging jitter: fine step that keeps squirrels spread across trees
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * sqf * sqf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if RNG::next(current_rng) < 0.20 {
        // predator encounter: escape by random repositioning in the forest
        let escaped = random_solution(dim, bounds, current_rng)
        for d in 0..<dim {
          cand.push(escaped[d])
        }
      } else {
        // gather: contract toward the oak (global best) with a jitter
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * sqf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Puma: Puma Optimizer (2024) =====

///|
pub struct PumaConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_puma_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PumaConfig {
  PumaConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Puma Optimizer (2024), mechanism-faithful variant of puma hunting.
/// A puma factor pumf = e^{-t/2500} decays as evaluations are spent, and each
/// individual carries a fitness-adaptive step scale
/// w_i = 0.5 + (fit_i - fit_min) / (fit_max - fit_min) in [0.5, 1.5], so weaker
/// pumas roam further while stronger pumas converge gently:
///   * roam (r1 < 0.28 pumf) - drift toward a random peer with a differential
///     kick taken from two random peers,
///     X_i <- X_i + r2 (X_j - X_i) + r3 (X_j - X_k)
///            + r4 (2 r5 - 1) range pumf 0.12 w_i;
///   * stalk (r1 < 0.28 pumf + 0.45) - fine jitter that keeps the pride spread
///     out, X_i <- X_i + r2 (2 r3 - 1) range pumf^2;
///   * pounce (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range pumf 0.06 w_i.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn puma_optimizer(
  config : PumaConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // fitness range drives the adaptive step scale w_i
    let mut fit_min = fits[0]
    let mut fit_max = fits[0]
    for k in 1..<pop_size {
      if fits[k] < fit_min {
        fit_min = fits[k]
      }
      if fits[k] > fit_max {
        fit_max = fits[k]
      }
    }
    let fit_span = if fit_max - fit_min < 1.0e-12 {
      1.0e-12
    } else {
      fit_max - fit_min
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // puma factor pumf decays 1 -> 0 as evaluations are spent
      let pumf = @math.exp(-evaluations.to_double() / 2500.0)
      let wi = 0.5 + (fits[i] - fit_min) / fit_span
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < pumf * 0.28 {
        // roam: drift toward a random peer, plus a differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xj[d] - xi[d]) +
            r3 * (xj[d] - xk[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * pumf * 0.12 * wi
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < pumf * 0.28 + 0.45 {
        // stalk: fine jitter that keeps the pride spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * pumf * pumf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // pounce: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * pumf * 0.06 * wi
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== StMO: Starling Murmuration Optimizer (2022) =====

///|
pub struct StMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_stmo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> StMOConfig {
  StMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Starling Murmuration Optimizer (2022), mechanism-faithful variant of
/// starling flocking.  A starling factor stmf = e^{-t/2500} decays as
/// evaluations are spent, and the flock aligns to a *fitness-weighted
/// centroid* X_wc = sum_i w_i X_i / sum_i w_i with w_i = 1/(1+rel_fit_i) in
/// [0.5, 1] (better starlings pull the centroid harder):
///   * flock alignment (r1 < 0.30 stmf) - drift toward the weighted centroid
///     with a differential kick taken from two random peers,
///     X_i <- X_i + r2 (X_wc - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range stmf 0.11;
///   * wing jitter (r1 < 0.30 stmf + 0.44) - fine step that keeps the flock
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range stmf^2;
///   * settle (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range stmf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn starling_murmuration_optimizer(
  config : StMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // fitness-weighted centroid of the flock
    let mut fit_min = fits[0]
    let mut fit_max = fits[0]
    for k in 1..<pop_size {
      if fits[k] < fit_min {
        fit_min = fits[k]
      }
      if fits[k] > fit_max {
        fit_max = fits[k]
      }
    }
    let fit_span = if fit_max - fit_min < 1.0e-12 {
      1.0e-12
    } else {
      fit_max - fit_min
    }
    let weights : Array[Double] = []
    let mut wsum = 0.0
    for i in 0..<pop_size {
      let w = 1.0 / (1.0 + (fits[i] - fit_min) / fit_span)
      weights.push(w)
      wsum = wsum + w
    }
    let wc : Array[Double] = []
    for d in 0..<dim {
      let mut acc = 0.0
      for i in 0..<pop_size {
        acc = acc + weights[i] * pop[i][d]
      }
      wc.push(acc / wsum)
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // starling factor stmf decays 1 -> 0 as evaluations are spent
      let stmf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < stmf * 0.30 {
        // flock alignment: drift toward the weighted centroid + differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (wc[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * stmf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < stmf * 0.30 + 0.44 {
        // wing jitter: fine step that keeps the flock spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * stmf * stmf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * stmf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SFOA: Swordfish Optimization Algorithm (2021) =====

///|
pub struct SFOAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sfoa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SFOAConfig {
  SFOAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Swordfish Optimization Algorithm (2021), mechanism-faithful variant of
/// swordfish predation.  A swordfish factor swf = e^{-t/2500} decays as
/// evaluations are spent, and the previous generation's champion position is
/// kept as `prev_best` so that the champion's one-generation displacement
/// (X_best - X_prev_best) can be ridden as a velocity direction:
///   * forward thrust (r1 < 0.35 swf) - chase the prey (global best) while
///     riding the champion's velocity, plus a jitter,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_best - X_prev_best)
///            + r4 (2 r5 - 1) range swf 0.12;
///   * surface patrol (r1 < 0.35 swf + 0.41) - fine step that keeps the school
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range swf^2;
///   * deep dive (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range swf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn swordfish_optimization(
  config : SFOAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // previous generation's champion gives the velocity direction
    let prev_best = best_pos
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // swordfish factor swf decays 1 -> 0 as evaluations are spent
      let swf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < swf * 0.35 {
        // forward thrust: chase the prey (best) and ride the champion's velocity
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (best_pos[d] - prev_best[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * swf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < swf * 0.35 + 0.41 {
        // surface patrol: fine step that keeps the school spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * swf * swf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // deep dive: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * swf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ABO: African Buffalo Optimization (2016) =====

///|
pub struct ABOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_abo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ABOConfig {
  ABOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// African Buffalo Optimization (2016), mechanism-faithful variant of buffalo
/// herd movement.  A buffalo factor abf = e^{-t/2500} decays as evaluations are
/// spent, and each buffalo combines two independent learning probabilities
/// toward the "waaa" signal (global best) and toward a random herd mate:
///   * herd migration (r1 < 0.32 abf) - move toward the global best and a
///     random mate, each weighted by its own learning probability,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_mate - X_i)
///            + r4 (2 r5 - 1) range abf 0.14;
///   * grazing jitter (r1 < 0.32 abf + 0.44) - fine step that keeps the herd
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range abf^2;
///   * settle (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range abf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn african_buffalo_optimization(
  config : ABOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // buffalo factor abf decays 1 -> 0 as evaluations are spent
      let abf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < abf * 0.32 {
        // herd migration: weighted pull toward the best and a random mate
        let m = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xm = pop[m]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xm[d] - xi[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * abf * 0.14
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < abf * 0.32 + 0.44 {
        // grazing jitter: fine step that keeps the herd spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * abf * abf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * abf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Gannet: Gannet Optimization Algorithm (2022) =====

///|
pub struct GannetConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gannet_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GannetConfig {
  GannetConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Gannet Optimization Algorithm (2022), mechanism-faithful variant of gannet
/// diving.  A gannet factor gnf = e^{-t/2500} decays as evaluations are spent:
///   * V-shaped dive (r1 < 0.34 gnf) - plunge toward the prey (global best)
///     with a differential kick taken from two random peers,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range gnf 0.13;
///   * wing jitter (r1 < 0.34 gnf + 0.41) - fine step that keeps the flock
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range gnf^2;
///   * underwater catch (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range gnf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn gannet_optimization(
  config : GannetConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // gannet factor gnf decays 1 -> 0 as evaluations are spent
      let gnf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < gnf * 0.34 {
        // V-shaped dive: plunge toward the prey (best) + differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * gnf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < gnf * 0.34 + 0.41 {
        // wing jitter: fine step that keeps the flock spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * gnf * gnf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // underwater catch: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gnf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Narwhal: Narwhal Optimizer (2024) =====

///|
pub struct NarwhalConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_narwhal_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> NarwhalConfig {
  NarwhalConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Narwhal Optimizer (2024), mechanism-faithful variant of a narwhal pod split
/// into two pods that periodically exchange their best tusk-bearer.  A narwhal
/// factor narf = e^{-t/2500} decays as evaluations are spent; each pod contracts
/// toward *its own* best rather than the global best, and every 20 generations
/// the two pods swap their best individuals (migration):
///   * tusk thrust (r1 < 0.33 narf) - dive toward the pod best with a
///     differential kick taken from two random pod mates,
///     X_i <- X_i + r2 (X_pod - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range narf 0.12;
///   * deep jitter (r1 < 0.33 narf + 0.41) - fine step that keeps the pod
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range narf^2;
///   * pod gather (rest) - contract toward the pod best,
///     X_pod + r2 (X_i - X_pod) + r3 (2 r4 - 1) range narf 0.06.
pub fn narwhal_optimizer(
  config : NarwhalConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let half = pop_size / 2

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut generations = 0
  while evaluations < max_evals {
    // pod (island) bests for this generation
    let mut b0 = 0
    for k in 0..<half {
      if fits[k] < fits[b0] {
        b0 = k
      }
    }
    let mut b1 = half
    for k in half..<pop_size {
      if fits[k] < fits[b1] {
        b1 = k
      }
    }
    let pod0 = pop[b0]
    let pod1 = pop[b1]
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let pod = if i < half { pod0 } else { pod1 }
      // narwhal factor narf decays 1 -> 0 as evaluations are spent
      let narf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < narf * 0.33 {
        // tusk thrust: dive toward the pod best + differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (pod[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * narf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < narf * 0.33 + 0.41 {
        // deep jitter: fine step that keeps the pod spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * narf * narf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // pod gather: contract toward the pod best
        for d in 0..<dim {
          let nv = pod[d] +
            r2 * (xi[d] - pod[d]) +
            r3 * (2.0 * r4 - 1.0) * span * narf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    generations = generations + 1
    if generations % 20 == 0 {
      // migration: swap the best tusk-bearers of the two pods
      let mut m0 = 0
      for k in 0..<half {
        if fits[k] < fits[m0] {
          m0 = k
        }
      }
      let mut m1 = half
      for k in half..<pop_size {
        if fits[k] < fits[m1] {
          m1 = k
        }
      }
      let tmp_pos = pop[m0]
      let tmp_fit = fits[m0]
      pop[m0] = pop[m1]
      fits[m0] = fits[m1]
      pop[m1] = tmp_pos
      fits[m1] = tmp_fit
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Serval: Serval Optimization Algorithm (2022) =====

///|
pub struct ServalConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_serval_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ServalConfig {
  ServalConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Serval Optimization Algorithm (2022), mechanism-faithful variant of serval
/// ambush hunting.  A serval factor svf = e^{-t/2500} decays as evaluations are
/// spent, and a shared archive of the most recent champions (up to 5) is kept
/// so a hunting serval may pounce toward a *past* champion and not only the
/// current one:
///   * pounce (r1 < 0.32 svf) - leap toward a random archived champion with a
///     differential kick taken from two random peers,
///     X_i <- X_i + r2 (X_champ - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range svf 0.12;
///   * stalk (r1 < 0.32 svf + 0.42) - fine jitter that keeps the pride spread
///     out, X_i <- X_i + r2 (2 r3 - 1) range svf^2;
///   * ambush (rest) - contract toward the current global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range svf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn serval_optimization(
  config : ServalConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let elite_cap = 5

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let elite : Array[Array[Double]] = []
  let mut elite_n = 0
  let mut elite_slot = 0
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // serval factor svf decays 1 -> 0 as evaluations are spent
      let svf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < svf * 0.32 {
        // pounce toward a past champion from the archive + differential kick
        let tgt = if elite_n > 0 {
          elite[(RNG::next(current_rng) * elite_n.to_double()).to_int() %
          elite_n]
        } else {
          best_pos
        }
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (tgt[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * svf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < svf * 0.32 + 0.42 {
        // stalk: fine jitter that keeps the pride spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * svf * svf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // ambush: contract toward the current global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * svf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        // archive the displaced champion before adopting the new one
        if elite_n < elite_cap {
          elite.push(best_pos)
          elite_n = elite_n + 1
        } else {
          elite[elite_slot] = best_pos
          elite_slot = (elite_slot + 1) % elite_cap
        }
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== ChSO: Chicken Swarm Optimization (2014) =====

///|
pub struct ChSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_chso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ChSOConfig {
  ChSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Chicken Swarm Optimization (2014), mechanism-faithful variant of a chicken
/// flock with three roles and a mother-child bond.  Every 10 generations the
/// flock is sorted by fitness and roles are reassigned: the best `rooster_cnt`
/// become roosters, the worst `chick_cnt` become chicks, the rest are hens;
/// hens are grouped round-robin under roosters and chicks under hens.  A
/// chicken factor csf = e^{-t/2500} decays as evaluations are spent:
///   * rooster - drift toward a random other rooster with a shrinking jitter;
///   * hen - follow her group rooster plus a random flock mate;
///   * chick - follow her mother hen with a random stretch factor in [0, 2].
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn chicken_swarm_optimization(
  config : ChSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let rooster_cnt = 3
  let chick_cnt = pop_size / 4
  let hen_cnt = pop_size - rooster_cnt - chick_cnt

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let mut role_order : Array[Int32] = []
  let mut generations = 0
  while evaluations < max_evals {
    if generations % 10 == 0 {
      // reassign roles: best -> roosters, worst -> chicks, rest -> hens
      role_order = sort_indices(fits)
    }
    // roosters drift toward a random other rooster with a shrinking jitter
    for rr in 0..<rooster_cnt {
      if evaluations >= max_evals {
        break
      }
      let p = role_order[rr]
      let xi = pop[p]
      let csf = @math.exp(-evaluations.to_double() / 2500.0)
      let oi = (RNG::next(current_rng) * rooster_cnt.to_double()).to_int() %
        rooster_cnt
      let other = if oi == rr { (rr + 1) % rooster_cnt } else { oi }
      let xo = pop[role_order[other]]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = xi[d] +
          r1 * (best_pos[d] - xi[d]) +
          r2 * (xo[d] - xi[d]) +
          r3 * (2.0 * RNG::next(current_rng) - 1.0) * span * csf * 0.12
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      pop[p] = cand
      let fv = evaluate(f, cand)
      fits[p] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    // hens follow their group rooster and a random flock mate
    for hs in 0..<hen_cnt {
      if evaluations >= max_evals {
        break
      }
      let p = role_order[rooster_cnt + hs]
      let xi = pop[p]
      let csf = @math.exp(-evaluations.to_double() / 2500.0)
      let ro_slot = hs % rooster_cnt
      let xro = pop[role_order[ro_slot]]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = xi[d] +
          r1 * (xro[d] - xi[d]) +
          r2 * (best_pos[d] - xi[d]) +
          r3 * (2.0 * RNG::next(current_rng) - 1.0) * span * csf * 0.06
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      pop[p] = cand
      let fv = evaluate(f, cand)
      fits[p] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    // chicks follow their mother hen with a random stretch factor in [0, 2]
    for cs in 0..<chick_cnt {
      if evaluations >= max_evals {
        break
      }
      let p = role_order[rooster_cnt + hen_cnt + cs]
      let xi = pop[p]
      let csf = @math.exp(-evaluations.to_double() / 2500.0)
      let mother_slot = cs % hen_cnt
      let xm = pop[role_order[rooster_cnt + mother_slot]]
      let stretch = 2.0 * RNG::next(current_rng)
      let r1 = RNG::next(current_rng)
      let cand : Array[Double] = []
      for d in 0..<dim {
        let nv = xi[d] +
          stretch * (xm[d] - xi[d]) +
          r1 * (2.0 * RNG::next(current_rng) - 1.0) * span * csf * 0.05
        cand.push(max(bounds.min, min(bounds.max, nv)))
      }
      pop[p] = cand
      let fv = evaluate(f, cand)
      fits[p] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
    generations = generations + 1
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== PR: Path Relinking (1997) =====

///|
pub struct PRConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_pr_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PRConfig {
  PRConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Path Relinking (1997), mechanism-faithful standalone variant.  A relinking
/// factor prf = e^{-t/2500} decays as evaluations are spent, and an archive of
/// the most recent champions (up to 5, round-robin) provides the elite
/// endpoints of the relinking paths:
///   * relink (r1 < 0.29 prf) - move toward a random point on the path between
///     two distinct archived champions, X_target = X_e1 + r2 (X_e2 - X_e1),
///     X_i <- X_i + r3 (X_target - X_i)
///            + r4 (2 r5 - 1) range prf 0.13;
///   * jitter (r1 < 0.29 prf + 0.45) - fine step that keeps the swarm spread
///     out, X_i <- X_i + r2 (2 r3 - 1) range prf^2;
///   * converge (rest) - contract toward the current global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range prf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn path_relinking(
  config : PRConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let elite_cap = 5

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let elite : Array[Array[Double]] = []
  let mut elite_n = 0
  let mut elite_slot = 0
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // relinking factor prf decays 1 -> 0 as evaluations are spent
      let prf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < prf * 0.29 {
        // relink toward a random point on the path between two champions
        if elite_n >= 2 {
          let ea = (RNG::next(current_rng) * elite_n.to_double()).to_int() %
            elite_n
          let mut eb = (RNG::next(current_rng) * elite_n.to_double()).to_int() %
            elite_n
          if eb == ea {
            eb = (ea + 1) % elite_n
          }
          let e1 = elite[ea]
          let e2 = elite[eb]
          for d in 0..<dim {
            let tgt = e1[d] + r2 * (e2[d] - e1[d])
            let nv = xi[d] +
              r3 * (tgt - xi[d]) +
              r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * prf * 0.13
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        } else {
          for d in 0..<dim {
            let nv = xi[d] +
              r2 * (best_pos[d] - xi[d]) +
              r3 * (2.0 * r4 - 1.0) * span * prf * 0.13
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        }
      } else if r1 < prf * 0.29 + 0.45 {
        // jitter: fine step that keeps the swarm spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * prf * prf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // converge: contract toward the current global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * prf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        if elite_n < elite_cap {
          elite.push(best_pos)
          elite_n = elite_n + 1
        } else {
          elite[elite_slot] = best_pos
          elite_slot = (elite_slot + 1) % elite_cap
        }
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CrSO: Cockroach Swarm Optimization (2008) =====

///|
pub struct CrSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_crso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CrSOConfig {
  CrSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Cockroach Swarm Optimization (2008), mechanism-faithful variant of roach
/// swarming.  A roach factor crf = e^{-t/2500} decays as evaluations are spent,
/// and each roach senses the best mate only within a local visibility radius
/// (0.30 range), so reference points are *local* rather than global:
///   * chase-swarming (r1 < 0.31 crf) - move toward the best roach inside the
///     visibility radius, with a differential kick taken from two random mates,
///     X_i <- X_i + r2 (X_nb - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range crf 0.10;
///   * dispersion (r1 < 0.31 crf + 0.41) - fine jitter that keeps the swarm
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range crf^2;
///   * settle (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range crf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn cockroach_swarm_optimization(
  config : CrSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // roach factor crf decays 1 -> 0 as evaluations are spent
      let crf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < crf * 0.31 {
        // chase-swarming: move toward the best mate within the visibility radius
        let vis2 = span * 0.30 * span * 0.30
        let mut nb = i
        let mut nb_fit = fits[i]
        for jj in 0..<pop_size {
          let xj = pop[jj]
          let mut dsq = 0.0
          for dd in 0..<dim {
            let diff = xj[dd] - xi[dd]
            dsq = dsq + diff * diff
          }
          if dsq < vis2 && fits[jj] < nb_fit {
            nb_fit = fits[jj]
            nb = jj
          }
        }
        let xnb = pop[nb]
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xnb[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * crf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < crf * 0.31 + 0.41 {
        // dispersion: fine jitter that keeps the swarm spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * crf * crf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * crf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== JS: Jellyfish Search Optimizer (2020) =====

///|
pub struct JSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_js_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> JSConfig {
  JSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Jellyfish Search Optimizer (2020), mechanism-faithful variant of jellyfish
/// drifting/following.  A jellyfish factor jf = e^{-t/2500} shrinks the steps,
/// and a *time control* criterion tc = |(1 - t/T)(2 r - 1)| (re-sampled every
/// move, hence non-monotonic) switches between two modes:
///   * ocean current (tc >= 0.35) - drift toward the global best with a
///     differential kick taken from two random peers,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range jf 0.11;
///   * swarm (tc < 0.35) - passive motion (r1 < 0.5) jitters around itself to
///     keep the bloom spread out, X_i <- X_i + r2 (2 r3 - 1) range jf 0.06,
///     otherwise active motion contracts toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range jf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn jellyfish_search_optimizer(
  config : JSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let max_evals_d = max_evals.to_double()

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // jellyfish factor jf shrinks 1 -> 0 as evaluations are spent
      let jf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      // time control criterion |(1 - t/T)(2r-1)|, re-sampled for every move
      let tc_raw = (1.0 - evaluations.to_double() / max_evals_d) *
        (2.0 * RNG::next(current_rng) - 1.0)
      let tc = if tc_raw < 0.0 { -tc_raw } else { tc_raw }
      let cand : Array[Double] = []
      if tc >= 0.35 {
        // ocean current: drift toward the global best + differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * jf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < 0.5 {
        // passive motion: fine jitter around itself, keeps the bloom spread
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * jf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // active motion: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * jf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WSA: Water Strider Algorithm (2020) =====

///|
pub struct WSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WSAConfig {
  WSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Water Strider Algorithm (2020), mechanism-faithful variant of water strider
/// territory movement.  A strider factor wsf = e^{-t/2500} decays as
/// evaluations are spent, and the population is ranked by fitness every
/// generation so each strider is pulled toward its *rank-adjacent* better
/// neighbour (forming a chain from the best down the ranks) rather than the
/// global best:
///   * territory migration (r1 < 0.33 wsf) - drift toward the immediately
///     better-ranked strider with a differential kick from two random mates,
///     X_i <- X_i + r2 (X_better - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range wsf 0.11;
///   * ripple jitter (r1 < 0.33 wsf + 0.42) - fine step that keeps the pond
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range wsf^2;
///   * settle (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range wsf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn water_strider_algorithm(
  config : WSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // rank the pond by fitness for rank-adjacency attraction
    let order = sort_indices(fits)
    let rank : Array[Int32] = []
    for _k in 0..<pop_size {
      rank.push(0)
    }
    for k in 0..<pop_size {
      rank[order[k]] = k
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // strider factor wsf decays 1 -> 0 as evaluations are spent
      let wsf = @math.exp(-evaluations.to_double() / 2500.0)
      let rrk = rank[i]
      let xbetter = if rrk > 0 { pop[order[rrk - 1]] } else { best_pos }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < wsf * 0.33 {
        // territory migration: toward the immediately better-ranked strider
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xbetter[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * wsf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < wsf * 0.33 + 0.42 {
        // ripple jitter: fine step that keeps the pond spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * wsf * wsf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * wsf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SeHO: Sea Horse Optimizer (2022) =====

///|
pub struct SeHOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_seho_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SeHOConfig {
  SeHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Sea Horse Optimizer (2022), mechanism-faithful variant of seahorse motion.
/// A seahorse factor shf = e^{-t/2500} decays as evaluations are spent:
///   * spiral predation (r1 < 0.31 shf) - contract toward the prey (global
///     best) with a jitter,
///     X_i <- X_i + r2 shf (X_best - X_i) + r3 (2 r4 - 1) range shf 0.12;
///   * heavy-tailed drift (r1 < 0.31 shf + 0.45) - drift with an inverse-power
///     step 1/(u+0.1) (capped at 8) so rare long jumps keep the shoal
///     exploring, X_i <- X_i + step (2 r5 - 1) range shf 0.03;
///   * homing (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range shf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn sea_horse_optimizer(
  config : SeHOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // seahorse factor shf decays 1 -> 0 as evaluations are spent
      let shf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < shf * 0.31 {
        // spiral predation: contract toward the prey (global best)
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (2.0 * r4 - 1.0) * span * shf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < shf * 0.31 + 0.45 {
        // heavy-tailed drift: rare long jumps keep the shoal exploring
        let u = RNG::next(current_rng)
        let mut step = 1.0 / (u + 0.1)
        if step > 8.0 {
          step = 8.0
        }
        for d in 0..<dim {
          let nv = xi[d] +
            step * (2.0 * RNG::next(current_rng) - 1.0) * span * shf * 0.03
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // homing: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * shf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WBO: Wildebeest Herd Optimization =====

///|
pub struct WBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WBOConfig {
  WBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Wildebeest Herd Optimization, mechanism-faithful variant of wildebeest herd
/// migration.  A herd factor whf = e^{-t/2500} decays as evaluations are spent,
/// and each generation the champion's one-generation displacement
/// (X_best - X_prev_best) is linearly extrapolated into a predicted trend point
/// X_pred = 2 X_best - X_prev_best:
///   * herd migration (r1 < 0.35 whf) - move toward the extrapolated trend
///     point with a differential kick taken from two random mates,
///     X_i <- X_i + r2 (X_pred - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range whf 0.13;
///   * graze jitter (r1 < 0.35 whf + 0.40) - fine step that keeps the herd
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range whf^2;
///   * rejoin herd (rest) - contract toward the champion,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range whf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn wildebeest_herd_optimization(
  config : WBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // previous champion gives the trend direction of the herd
    let prev_best = best_pos
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // herd factor whf decays 1 -> 0 as evaluations are spent
      let whf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < whf * 0.35 {
        // herd migration: move toward the extrapolated champion trend point
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let pred = best_pos[d] + (best_pos[d] - prev_best[d])
          let nv = xi[d] +
            r2 * (pred - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * whf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < whf * 0.35 + 0.40 {
        // graze jitter: fine step that keeps the herd spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * whf * whf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // rejoin herd: contract toward the champion
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * whf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== OWL: Owl Optimization Algorithm =====

///|
pub struct OWLConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_owl_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> OWLConfig {
  OWLConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Owl Optimization Algorithm, mechanism-faithful variant of owl hunting.  An
/// owl factor owf = e^{-t/2500} decays as evaluations are spent, and each owl
/// remembers the direction of its last *successful* move (a pattern-move
/// memory, unlike BCMO's raw momentum or SAR's best-position memory):
///   * silent glide (r1 < 0.37 owf) - continue along the remembered success
///     direction while closing in on the prey (global best), with a
///     differential kick and a jitter,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q) + 0.25 D_i
///            + r4 (2 r5 - 1) range owf 0.11;
///   * auditory jitter (r1 < 0.37 owf + 0.39) - fine step that keeps the roost
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range owf^2;
///   * swoop (rest) - contract toward the prey (global best),
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range owf 0.05.
/// The direction memory D_i is refreshed on improvement and forgotten after a
/// non-improving move.  The global best is tracked separately.
pub fn owl_optimization(
  config : OWLConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  // per-owl remembered success direction
  let dir : Array[Array[Double]] = []
  for _i in 0..<pop_size {
    let z : Array[Double] = []
    for _d in 0..<dim {
      z.push(0.0)
    }
    dir.push(z)
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let di = dir[i]
      // owl factor owf decays 1 -> 0 as evaluations are spent
      let owf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      let mut glided = false
      if r1 < owf * 0.37 {
        // silent glide: continue the remembered success direction + pull to prey
        glided = true
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * owf * 0.11 +
            0.25 * di[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < owf * 0.37 + 0.39 {
        // auditory jitter: fine step that keeps the roost spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * owf * owf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // swoop: contract toward the prey (global best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * owf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      // keep the remembered direction only when a glide improved, else forget
      let nd : Array[Double] = []
      if glided && fv < fits[i] {
        for d in 0..<dim {
          nd.push(0.5 * (cand[d] - xi[d]))
        }
      } else {
        for _d in 0..<dim {
          nd.push(0.0)
        }
      }
      dir[i] = nd
      fits[i] = fv
      pop[i] = cand
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== RRO: Raven Roosting Optimization =====

///|
pub struct RROConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_rro_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RROConfig {
  RROConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Raven Roosting Optimization, mechanism-faithful variant of ravens choosing
/// a roost by quality.  A roost factor rrf = e^{-t/2500} decays as evaluations
/// are spent, and every generation the flock is ranked so a *rank-weighted
/// roulette* can pick the roost a raven flies to (any rank can be chosen, but
/// better ranks are likelier):
///   * change roost (r1 < 0.36 rrf) - fly toward a roulette-selected roost,
///     adding a differential kick from two random ravens,
///     X_i <- X_i + r2 (X_roost - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range rrf 0.13;
///   * roost jitter (r1 < 0.36 rrf + 0.40) - fine step that keeps the flock
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range rrf^2;
///   * settle (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range rrf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn raven_roosting_optimization(
  config : RROConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // rank the flock and build a rank-weighted cumulative roulette
    let order = sort_indices(fits)
    let cum : Array[Double] = []
    let mut acc = 0.0
    for k in 0..<pop_size {
      acc = acc + (pop_size - k).to_double()
      cum.push(acc)
    }
    let total = acc
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // roost factor rrf decays 1 -> 0 as evaluations are spent
      let rrf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < rrf * 0.36 {
        // change roost: fly toward a roulette-selected roost + differential kick
        let u = RNG::next(current_rng) * total
        let mut pick_k = 0
        for k in 0..<pop_size {
          if cum[k] < u {
            pick_k = k + 1
          }
        }
        let xref = pop[order[pick_k]]
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let kk = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[kk]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xref[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * rrf * 0.13
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < rrf * 0.36 + 0.40 {
        // roost jitter: fine step that keeps the flock spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * rrf * rrf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * rrf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SaFO: Sailfish Optimizer (2019) =====

///|
pub struct SaFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_safo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SaFOConfig {
  SaFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Sailfish Optimizer (2019), mechanism-faithful variant of sailfish hunting.
/// A sailfish factor saf = e^{-t/2500} decays as evaluations are spent, and an
/// *adaptive fitness-distance balance* (FDB) picks the prey among a few random
/// candidates - trading a good fitness against being far enough away to stay
/// exploratory:
///   * attack (r1 < 0.31 saf) - strike toward the FDB-selected prey with a
///     differential kick from two random mates,
///     X_i <- X_i + r2 (X_prey - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range saf 0.12;
///   * patrol jitter (r1 < 0.31 saf + 0.43) - fine step that keeps the school
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range saf^2;
///   * strike home (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range saf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn sailfish_optimizer(
  config : SaFOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // fitness range for FDB normalization
    let mut fit_min = fits[0]
    let mut fit_max = fits[0]
    for k in 1..<pop_size {
      if fits[k] < fit_min {
        fit_min = fits[k]
      }
      if fits[k] > fit_max {
        fit_max = fits[k]
      }
    }
    let fit_span = if fit_max - fit_min < 1.0e-12 {
      1.0e-12
    } else {
      fit_max - fit_min
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // sailfish factor saf decays 1 -> 0 as evaluations are spent
      let saf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < saf * 0.31 {
        // FDB: pick a prey balancing a good fitness against being far away
        let mut pick = i
        let mut best_score = 1.0e30
        for _c in 0..<3 {
          let cj = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
            pop_size
          let xc = pop[cj]
          let nf = (fits[cj] - fit_min) / fit_span
          let mut dsq = 0.0
          for d in 0..<dim {
            let diff = xc[d] - xi[d]
            dsq = dsq + diff * diff
          }
          let mut nd = dsq / (span * span)
          if nd > 1.0 {
            nd = 1.0
          }
          let score = 0.5 * nf + 0.5 * (1.0 - nd)
          if score < best_score {
            best_score = score
            pick = cj
          }
        }
        let xref = pop[pick]
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let kk = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[kk]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xref[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * saf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < saf * 0.31 + 0.43 {
        // patrol jitter: fine step that keeps the school spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * saf * saf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // strike home: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * saf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== GPC: Giza Pyramids Construction (2021) =====

///|
pub struct GPCConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gpc_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GPCConfig {
  GPCConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Giza Pyramids Construction (2021), mechanism-faithful variant of workers
/// hauling stone blocks up a pyramid.  A construction factor gpf = e^{-t/2500}
/// decays as evaluations are spent:
///   * hauling (r1 < 0.34 gpf) - carry a block toward a point midway between
///     the crew mate and the apex (global best), plus a differential kick,
///     X_i <- X_i + r2 (0.5 (X_best + X_p) - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range gpf 0.10;
///   * worker wander (r1 < 0.34 gpf + 0.43) - fine jitter that keeps the crew
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range gpf^2;
///   * ascent (rest) - contract toward the apex,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range gpf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn giza_pyramids_construction(
  config : GPCConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // construction factor gpf decays 1 -> 0 as evaluations are spent
      let gpf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < gpf * 0.34 {
        // hauling: carry the block toward the crew/apex midpoint + kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let mid = 0.5 * (best_pos[d] + xp[d])
          let nv = xi[d] +
            r2 * (mid - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * gpf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < gpf * 0.34 + 0.43 {
        // worker wander: fine step that keeps the crew spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * gpf * gpf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // ascent: contract toward the apex (global best)
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * gpf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== SSD: Social Ski-Driver (2019) =====

///|
pub struct SSDConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ssd_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SSDConfig {
  SSDConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Social Ski-Driver (2019), mechanism-faithful variant of the ski-driver
/// social motion.  Each skier keeps a personal best and is driven by two
/// references at once: its personal best and the midpoint between the personal
/// best and the global best.  A driver factor ssf = e^{-t/2500} decays as
/// evaluations are spent:
///   * ski-drive (r1 < 0.36 ssf) - move under both references with a jitter,
///     X_i <- X_i + r2 (P_i - X_i) + r3 (0.5 (X_best + P_i) - X_i)
///            + r4 (2 r5 - 1) range ssf 0.11;
///   * drift (r1 < 0.36 ssf + 0.41) - fine jitter that keeps the track spread
///     out, X_i <- X_i + r2 (2 r3 - 1) range ssf^2;
///   * finish (rest) - contract toward the global best,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ssf 0.05.
/// The personal best P_i is refreshed on improvement; the global best is
/// tracked separately.
pub fn social_ski_driver(
  config : SSDConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  // personal best position / value for every skier
  let pbest : Array[Array[Double]] = []
  let pbest_fit : Array[Double] = []
  for i in 0..<pop_size {
    pbest.push(pop[i])
    pbest_fit.push(fits[i])
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let pb = pbest[i]
      // driver factor ssf decays 1 -> 0 as evaluations are spent
      let ssf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ssf * 0.36 {
        // ski-drive: pulled by the personal best and the pbest/gbest midpoint
        for d in 0..<dim {
          let mid = 0.5 * (best_pos[d] + pb[d])
          let nv = xi[d] +
            r2 * (pb[d] - xi[d]) +
            r3 * (mid - xi[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * ssf * 0.11
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ssf * 0.36 + 0.41 {
        // drift: fine step that keeps the track spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * ssf * ssf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // finish: contract toward the global best
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ssf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      pop[i] = cand
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < pbest_fit[i] {
        pbest[i] = cand
        pbest_fit[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== FDO: Fitness Dependent Optimizer (2020) =====

///|
pub struct FDOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_fdo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FDOConfig {
  FDOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Fitness Dependent Optimizer (2020), mechanism-faithful variant of the
/// fitness-dependent bee swarm.  Each agent owns a relative fitness gap to
/// the champion, which scales its stride, and picks one of three roles:
///   * fitness-dependent advance (r1 < 0.33) - a stride weighted pull toward
///     the champion plus a differential kick from two random mates,
///     X_i <- X_i + stride fdf r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range fdf 0.08;
///   * roost jitter (r1 < 0.75) - fine step that keeps the swarm spread out,
///     X_i <- X_i + r2 (2 r3 - 1) range fdf^2;
///   * settle (rest) - contract toward the champion,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range fdf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn fitness_dependent_optimizer(
  config : FDOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // update factor fdf decays 1 -> 0 as evaluations are spent
      let fdf = @math.exp(-evaluations.to_double() / 2500.0)
      // relative fitness gap to the champion, capped for stability
      let mut stride = 0.1 +
        abs_val(fits[i] - best_val) / (abs_val(best_val) + 1.0)
      if stride > 1.2 {
        stride = 1.2
      }
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < 0.33 {
        // fitness-dependent advance toward the champion + differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            stride * fdf * r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * fdf * 0.08
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < 0.75 {
        // roost jitter: fine step that keeps the swarm spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * fdf * fdf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the champion
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * fdf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== CFA: Cuttlefish Optimization Algorithm (2013) =====

///|
pub struct CFAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_cfa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CFAConfig {
  CFAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Cuttlefish Optimization Algorithm (2013), mechanism-faithful variant of
/// cuttlefish color-cell reflection and visibility.  A cell factor ctf =
/// e^{-t/2500} decays as evaluations are spent; every generation the reef
/// centroid is recomputed so a cell can mirror itself about the centroid
/// before matching its pattern with the champion:
///   * reflection (r1 < 0.30 ctf) - mirror the cell about the centroid, then
///     apply a visibility-weighted pull toward the champion plus a
///     differential kick, M_i = 2 X_mean - X_i,
///     X_i <- M_i + r2 (X_best - M_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range ctf 0.12;
///   * visibility jitter (r1 < 0.30 ctf + 0.44) - fine step that keeps the
///     reef spread out, X_i <- X_i + r2 (2 r3 - 1) range ctf^2;
///   * settle (rest) - contract toward the champion,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range ctf 0.06.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn cuttlefish_optimization(
  config : CFAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // reef centroid, used as the mirror plane for every cell
    let mean : Array[Double] = []
    for d in 0..<dim {
      let mut s = 0.0
      for j in 0..<pop_size {
        s = s + pop[j][d]
      }
      mean.push(s / pop_size.to_double())
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // cell factor ctf decays 1 -> 0 as evaluations are spent
      let ctf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < ctf * 0.30 {
        // reflection: mirror about the centroid, then match toward the champion
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let mirror = 2.0 * mean[d] - xi[d]
          let nv = mirror +
            r2 * (best_pos[d] - mirror) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * ctf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < ctf * 0.30 + 0.44 {
        // visibility jitter: fine step that keeps the reef spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * ctf * ctf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the champion
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * ctf * 0.06
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== OPA: Orca Predation Algorithm (2024) =====

///|
pub struct OPAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_opa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> OPAConfig {
  OPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Orca Predation Algorithm (2024), mechanism-faithful variant of an orca pod
/// herding prey.  Each orca remembers whether its previous move succeeded;
/// success biases it toward exploitation (press the attack on the champion)
/// while failure biases it toward exploration (break off and search):
///   * press the attack (orcas that improved, r1 < 0.78) - contract toward the
///     champion with a small jitter,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (2 r4 - 1) range opf 0.06;
///   * herding sweep (orcas that improved, rest) - explore with a differential
///     kick from two random mates,
///     X_i <- X_i + r2 (X_p - X_q) + r3 (2 r4 - 1) range opf 0.10;
///   * search dive (orcas that failed, r1 < 0.90) - strike toward the champion
///     with a differential kick and a wider jitter,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range opf 0.08;
///   * break away (orcas that failed, rest) - a random repositioning blended toward the champion,
///     X_i <- X_i + r2 (X_best - X_i) + (2 r3 - 1) range opf 0.20.
/// The success flag is refreshed every move; the global best is tracked
/// separately.
pub fn orca_predation_algorithm(
  config : OPAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  // per-orca success flag from its previous move
  let improved : Array[Bool] = []
  for _i in 0..<pop_size {
    improved.push(true)
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // pod factor opf decays 1 -> 0 as evaluations are spent
      let opf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if improved[i] {
        if r1 < 0.78 {
          // press the attack: contract toward the champion with a small jitter
          for d in 0..<dim {
            let nv = xi[d] +
              r2 * (best_pos[d] - xi[d]) +
              r3 * (2.0 * r4 - 1.0) * span * opf * 0.06
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        } else {
          // herding sweep: explore with a differential kick
          let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
            pop_size
          let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
            pop_size
          let xp = pop[j]
          let xq = pop[k]
          for d in 0..<dim {
            let nv = xi[d] +
              r2 * (xp[d] - xq[d]) +
              r3 * (2.0 * r4 - 1.0) * span * opf * 0.10
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        }
      } else if r1 < 0.90 {
        // search dive: strike toward the champion with a wider jitter
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * opf * 0.08
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // break away: large random repositioning keeps the pod exploring
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            (2.0 * r3 - 1.0) * span * opf * 0.20
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      improved[i] = fv < fits[i]
      pop[i] = cand
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== APuf: Arctic Puffin Optimization (2023) =====

///|
pub struct APufConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_apuf_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> APufConfig {
  APufConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Arctic Puffin Optimization (2023), mechanism-faithful variant of puffin
/// foraging.  A puffin factor puf = e^{-t/2500} decays as evaluations are
/// spent, and each *dimension* keeps an adaptive step scale s_d in [0.2, 2]
/// that expands for the coordinates that moved on a successful move and
/// shrinks otherwise - a dimension-wise (not per-individual scalar) step
/// adaptation:
///   * aerial cruise (r1 < 0.33 puf) - differential exploration between two
///     random mates with a dimension-scaled jitter,
///     X_i <- X_i + r2 (X_p - X_q) + r3 (2 r4 - 1) range puf 0.12 s_d;
///   * dive (r1 < 0.33 puf + 0.40) - fold toward the champion with a
///     dimension-scaled jitter,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range puf 0.06 s_d;
///   * surface jitter (rest) - fine step that keeps the colony spread out,
///     X_i <- X_i + r2 (2 r3 - 1) range puf^2 s_d.
/// After each move every coordinate that contributed at least 2% of the range
/// to a successful displacement gets s_d <- min(2, 1.08 s_d); all other
/// coordinates get s_d <- max(0.2, 0.96 s_d), and a failed move multiplies
/// every s_d by 0.98.  The global best is tracked separately.
pub fn arctic_puffin_optimization(
  config : APufConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  // per-dimension adaptive step scale
  let dscale : Array[Double] = []
  for _d in 0..<dim {
    dscale.push(1.0)
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // puffin factor puf decays 1 -> 0 as evaluations are spent
      let puf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < puf * 0.33 {
        // aerial cruise: differential exploration with dimension-scaled jitter
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (xp[d] - xq[d]) +
            r3 * (2.0 * r4 - 1.0) * span * puf * 0.12 * dscale[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < puf * 0.33 + 0.40 {
        // dive: fold toward the champion with dimension-scaled jitter
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * puf * 0.06 * dscale[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // surface jitter: fine step that keeps the colony spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * puf * puf * dscale[d]
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      let improved = fv < fits[i]
      // dimension-wise step adaptation
      for d in 0..<dim {
        if improved {
          if abs_val(cand[d] - xi[d]) >= 0.02 * span {
            dscale[d] = min(2.0, dscale[d] * 1.08)
          } else {
            dscale[d] = max(0.2, dscale[d] * 0.96)
          }
        } else {
          dscale[d] = max(0.2, dscale[d] * 0.98)
        }
      }
      pop[i] = cand
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== WPA: Wolf Pack Algorithm (2011) =====

///|
pub struct WPAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wpa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WPAConfig {
  WPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Wolf Pack Algorithm (2011), mechanism-faithful variant of wolf-pack hunting.
/// A pack factor wpf = e^{-t/2500} decays as evaluations are spent, and the
/// pack separates the *scouting* role from the *calling* and *besieging* roles:
///   * scouting (r1 < 0.33 wpf) - every coordinate takes an independent random
///     step so scouts sweep the terrain instead of collapsing on the champion,
///     X_i[d] <- X_i[d] + (2 r2 - 1) range wpf 0.12;
///   * calling (r1 < 0.33 wpf + 0.40) - the lead wolf summons the pack: pull
///     toward the champion with a differential kick from two random mates,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range wpf 0.10;
///   * besieging (rest) - close in on the champion in small steps,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range wpf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn wolf_pack_algorithm(
  config : WPAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // pack factor wpf decays 1 -> 0 as evaluations are spent
      let wpf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < wpf * 0.33 {
        // scouting: independent random step on every coordinate
        for d in 0..<dim {
          let nv = xi[d] + (2.0 * r2 - 1.0) * span * wpf * 0.12
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < wpf * 0.33 + 0.40 {
        // calling: the lead wolf summons the pack, plus a differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * wpf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // besieging: close in on the champion with small steps
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * wpf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Racc: Raccoon Optimization Algorithm (2022) =====

///|
pub struct RaccConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_racc_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> RaccConfig {
  RaccConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Raccoon Optimization Algorithm (2022), mechanism-faithful variant of
/// raccoon foraging.  A band factor raf = e^{-t/2500} decays as evaluations
/// are spent:
///   * mirror forage (r1 < 0.34 raf) - reflect the raccoon through the
///     midpoint between the champion and a random mate,
///     M_i = X_best + X_r - X_i, then add a jitter,
///     X_i <- M_i + r2 (2 r3 - 1) range raf 0.10;
///   * paw jitter (r1 < 0.34 raf + 0.42) - fine step that keeps the band
///     spread out, X_i <- X_i + r2 (2 r3 - 1) range raf^2;
///   * settle (rest) - contract toward the champion,
///     X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range raf 0.05.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn raccoon_optimization(
  config : RaccConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // band factor raf decays 1 -> 0 as evaluations are spent
      let raf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < raf * 0.34 {
        // mirror forage: reflect about the champion/mate midpoint + jitter
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xr = pop[j]
        for d in 0..<dim {
          let mirror = best_pos[d] + xr[d] - xi[d]
          let nv = mirror + r2 * (2.0 * r3 - 1.0) * span * raf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < raf * 0.34 + 0.42 {
        // paw jitter: fine step that keeps the band spread out
        for d in 0..<dim {
          let nv = xi[d] + r2 * (2.0 * r3 - 1.0) * span * raf * raf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // settle: contract toward the champion
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * raf * 0.05
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Monkey Search Algorithm (2007) =====

///|
pub struct MonkeySearchConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_monkey_search_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MonkeySearchConfig {
  MonkeySearchConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Monkey Search (2007), mechanism-faithful variant of monkeys moving through
/// a forest.  A monkey factor msf = e^{-t/2500} decays as evaluations are
/// spent, and each monkey picks one of three moves:
///   * climb (r1 < 0.36 msf) - fine local refinement around the champion,
///     X_i <- X_best + r2 (X_i - X_best) + r3 (2 r4 - 1) range msf^2;
///   * watch-jump (r1 < 0.36 msf + 0.38) - leap toward the champion with a
///     differential kick from two random mates,
///     X_i <- X_i + r2 (X_best - X_i) + r3 (X_p - X_q)
///            + r4 (2 r5 - 1) range msf 0.10;
///   * somersault (rest) - a large leap blended toward the champion,
///     X_i <- X_i + r2 (X_best - X_i) + (2 r3 - 1) range msf 0.25.
/// Individuals move unconditionally while the global best is tracked
/// separately.
pub fn monkey_search_algorithm(
  config : MonkeySearchConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      // monkey factor msf decays 1 -> 0 as evaluations are spent
      let msf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let r4 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < msf * 0.36 {
        // climb: fine local refinement around the champion
        for d in 0..<dim {
          let nv = best_pos[d] +
            r2 * (xi[d] - best_pos[d]) +
            r3 * (2.0 * r4 - 1.0) * span * msf * msf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < msf * 0.36 + 0.38 {
        // watch-jump: leap toward the champion with a differential kick
        let j = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let k = (RNG::next(current_rng) * pop_size.to_double()).to_int() %
          pop_size
        let xp = pop[j]
        let xq = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            r3 * (xp[d] - xq[d]) +
            r4 * (2.0 * RNG::next(current_rng) - 1.0) * span * msf * 0.10
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // somersault: large leap blended toward the champion
        for d in 0..<dim {
          let nv = xi[d] +
            r2 * (best_pos[d] - xi[d]) +
            (2.0 * r3 - 1.0) * span * msf * 0.25
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Bird Swarm Algorithm (2016) =====

///|
pub struct BswaConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_bswa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> BswaConfig {
  BswaConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Bird Swarm Algorithm (2016).  Each bird keeps a personal best and the
/// swarm tracks a global best; every round the swarm centroid and the mean
/// fitness are recomputed, then each bird acts with one of three rules:
///   * foraging (r1 < 0.88) - either a pull toward its personal best
///     (30% of the time, coefficient 1.5) or toward the global best
///     (coefficient 2.0), the two references being mutually exclusive;
///   * vigilance (r1 < 0.92) - the sign of (fit_i - mean_fit) decides the
///     direction: birds better than the average drift gently toward the
///     centroid with amplitude A1 = exp(-mean_fit/fit_i) 0.5, while worse
///     birds are pushed away from it with amplitude A2 = exp(-gap/mean_fit)
///     plus a small jitter, which keeps sub-populations from collapsing;
///   * flight (r1 >= 0.92) - a large leap toward the champion plus a random
///     excursion, mimicking the sudden relocation of a flock.
/// A swarm factor bsf = e^{-t/2500} decays every amplitude as evaluations
/// are spent, so the rules tighten into local refinement.
pub fn bird_swarm_algorithm(
  config : BswaConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 20 {
    20
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let pbest : Array[Array[Double]] = []
  let pbest_fit : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  pbest.push(best_pos)
  pbest_fit.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    pbest.push(p)
    pbest_fit.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // Swarm centroid and mean fitness drive the vigilance rule.
    let mean_pos : Array[Double] = []
    for _d in 0..<dim {
      mean_pos.push(0.0)
    }
    let mut fit_sum = 0.0
    for i in 0..<pop_size {
      fit_sum = fit_sum + fits[i]
      for d in 0..<dim {
        mean_pos[d] = mean_pos[d] + pop[i][d]
      }
    }
    let inv_pop = 1.0 / pop_size.to_double()
    let mean_fit = fit_sum * inv_pop
    for d in 0..<dim {
      mean_pos[d] = mean_pos[d] * inv_pop
    }

    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let bsf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < 0.92 {
        if r1 < 0.88 {
          // foraging: personal best or global best, mutually exclusive
          if r2 < 0.30 {
            for d in 0..<dim {
              let nv = xi[d] + 1.5 * r3 * (pbest[i][d] - xi[d]) * bsf
              cand.push(max(bounds.min, min(bounds.max, nv)))
            }
          } else {
            for d in 0..<dim {
              let nv = xi[d] + 2.0 * r3 * (best_pos[d] - xi[d]) * bsf
              cand.push(max(bounds.min, min(bounds.max, nv)))
            }
          }
        } else {
          // vigilance: sign of the fitness gap picks approach or retreat
          let gap = fits[i] - mean_fit
          if gap < 0.0 {
            let a1 = @math.exp(-mean_fit / (fits[i] + 1.0e-9)) * 0.5
            for d in 0..<dim {
              let nv = xi[d] + a1 * r2 * (mean_pos[d] - xi[d]) * bsf
              cand.push(max(bounds.min, min(bounds.max, nv)))
            }
          } else {
            let a2 = @math.exp(-gap / (mean_fit + 1.0e-9)) * 0.5
            for d in 0..<dim {
              let nv = xi[d] -
                a2 * r2 * (mean_pos[d] - xi[d]) * bsf +
                r3 * (2.0 * RNG::next(current_rng) - 1.0) * span * 0.02 * bsf
              cand.push(max(bounds.min, min(bounds.max, nv)))
            }
          }
        }
      } else {
        // flight: sudden relocation toward the champion plus a random draw
        for d in 0..<dim {
          let nv = xi[d] +
            (best_pos[d] - xi[d]) * (0.5 + r2) * bsf +
            r3 * (2.0 * RNG::next(current_rng) - 1.0) * span * 0.15 * bsf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < pbest_fit[i] {
        pbest_fit[i] = fv
        pbest[i] = cand
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Tree Growth Algorithm (2018) =====

///|
pub struct TGAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_tga_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TGAConfig {
  TGAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Tree Growth Algorithm (2018).  Every round the forest is ranked by fitness
/// and split into roles, with the growth factor tgf = e^{-t/2500} shrinking
/// all strides as evaluations are spent:
///   * the best third (well-lit trees) grows toward the light, i.e. toward
///     the global best, with a differential kick from two random trees;
///   * the middle trees absorb water from the soil around the better of two
///     neighbouring trees plus a weak pull toward the best, a local
///     differential exploration at tgf^2 scale;
///   * the worst tenth (crowded-out trees) are replaced by a random seedling
///     grafted halfway onto the best-lit tree, so competition keeps the
///     forest diverse instead of letting it collapse.
pub fn tree_growth_algorithm(
  config : TGAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 20 {
    20
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let n_lit = if pop_size / 3 < 3 { 3 } else { pop_size / 3 }
  let n_crowded = if pop_size / 10 < 2 { 2 } else { pop_size / 10 }

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // rank the forest by fitness: order[0] is the best-lit tree
    let order = sort_indices(fits)
    let rank : Array[Int32] = []
    for _k in 0..<pop_size {
      rank.push(0)
    }
    for k in 0..<pop_size {
      rank[order[k]] = k
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let tgf = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if rank[i] < n_lit {
        // growth toward the light (global best) with a differential kick
        let j = (r2 * pop_size.to_double()).to_int() % pop_size
        let k = (r3 * pop_size.to_double()).to_int() % pop_size
        let xj = pop[j]
        let xk = pop[k]
        for d in 0..<dim {
          let nv = xi[d] +
            2.0 * r1 * tgf * (best_pos[d] - xi[d]) +
            0.5 * r3 * (xj[d] - xk[d]) * tgf +
            (2.0 * RNG::next(current_rng) - 1.0) * span * 0.04 * tgf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if rank[i] < pop_size - n_crowded {
        // water absorption: local exploration around the better of two trees
        let j = (r2 * pop_size.to_double()).to_int() % pop_size
        let k = (r3 * pop_size.to_double()).to_int() % pop_size
        let xj = if fits[j] <= fits[k] { pop[j] } else { pop[k] }
        for d in 0..<dim {
          let nv = xi[d] +
            0.9 * r1 * tgf * (xj[d] - xi[d]) +
            0.3 * tgf * (best_pos[d] - xi[d]) +
            (2.0 * r3 - 1.0) * span * 0.08 * tgf * tgf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // competition: a crowded-out tree is replaced by a seedling that is a
        // random draw grafted halfway onto the best-lit tree
        let seedling = random_solution(dim, bounds, current_rng)
        for d in 0..<dim {
          let nv = 0.5 * seedling[d] +
            0.5 * best_pos[d] +
            (2.0 * r2 - 1.0) * span * 0.03 * tgf
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Plant Propagation Algorithm (2014) =====

///|
pub struct PPAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_ppa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> PPAConfig {
  PPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Plant Propagation Algorithm (2014).  Each round the plants are ranked by
/// fitness and reproduce through runners, with the propagation factor
/// ppf = e^{-t/2500} shrinking all distances as evaluations are spent:
///   * a well-adapted plant (low rank ratio rn) sends out many runners over a
///     short distance, so it refines its own neighbourhood;
///   * a poorly adapted plant sends out fewer runners that travel far, which
///     spreads the species into unexplored ground;
///   * 35% of the runners head toward the best plant instead of a random
///     direction, giving the population a pull toward the global optimum, and
///     a runner replaces its parent only when it is better (competition).
pub fn plant_propagation_algorithm(
  config : PPAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  let denom = if pop_size > 1 { (pop_size - 1).to_double() } else { 1.0 }

  while evaluations < max_evals {
    // rank the plants: order[0] is the best adapted
    let order = sort_indices(fits)
    let rank : Array[Int32] = []
    for _k in 0..<pop_size {
      rank.push(0)
    }
    for k in 0..<pop_size {
      rank[order[k]] = k
    }
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let ppf = @math.exp(-evaluations.to_double() / 2500.0)
      let rn = rank[i].to_double() / denom
      // well-adapted plants send many short runners, poor ones few long ones
      let runner_count = 1 + ((1.0 - rn) * 2.0).to_int()
      let reach = span * (0.02 + 0.30 * rn) * ppf
      let mut child_pos : Array[Double] = xi
      let mut child_fit = fits[i]
      for _c in 0..<runner_count {
        if evaluations >= max_evals {
          break
        }
        let cand : Array[Double] = []
        if RNG::next(current_rng) < 0.35 {
          // runner growing toward the best plant
          let w = RNG::next(current_rng) * 0.5 + 0.3
          for d in 0..<dim {
            let nv = xi[d] +
              w * (best_pos[d] - xi[d]) +
              (2.0 * RNG::next(current_rng) - 1.0) * reach * 0.25
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        } else {
          // runner in a random direction
          for d in 0..<dim {
            let nv = xi[d] + (2.0 * RNG::next(current_rng) - 1.0) * reach
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        }
        let fv = evaluate(f, cand)
        evaluations = evaluations + 1
        if fv < child_fit {
          child_fit = fv
          child_pos = cand
        }
        if fv < best_val {
          best_val = fv
          best_pos = cand
        }
      }
      pop[i] = child_pos
      fits[i] = child_fit
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Honey Bee Mating Optimization (2001) =====

///|
pub struct HBMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hbmo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HBMOConfig {
  HBMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Honey Bee Mating Optimization (2001).  The best individual is the queen
/// and every other bee is a drone that may mate with her.  The annealing
/// speed speed = e^{-t/2500} controls both the mating acceptance and the
/// brood care:
///   * a drone with normalized fitness gap g = |fit_i - fit_best| /
///     (|fit_i| + |fit_best| + eps) is accepted with probability
///     p = e^{-g / speed}, so early (fast) rounds accept diverse drones while
///     late (slow) rounds only accept drones almost as good as the queen;
///   * an accepted drone crosses with the queen,
///     brood = queen + w (drone - queen) + jitter, inheriting part of the
///     drone genome;
///   * a rejected drone stays as a worker and patrols its neighbourhood with
///     a weak pull toward the queen plus a fine jitter, so the hive keeps
///     converging even when mating is rare.
/// The brood replaces its drone only when it is better, and the queen is
/// replaced whenever a brood beats her.
pub fn honey_bee_mating_optimization(
  config : HBMOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // annealing mating speed, floored so late rounds still allow mating
    let speed = 0.2 + 0.8 * @math.exp(-evaluations.to_double() / 2500.0)
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let gap = abs_val(fits[i] - best_val) /
        (abs_val(fits[i]) + abs_val(best_val) + 1.0e-9)
      let p_mate = @math.exp(-gap / (speed + 1.0e-9))
      let cand : Array[Double] = []
      if r1 < p_mate {
        // mating: cross the drone genome into the queen
        let w = 0.3 + 0.6 * r2
        for d in 0..<dim {
          let nv = best_pos[d] +
            w * (xi[d] - best_pos[d]) +
            (2.0 * RNG::next(current_rng) - 1.0) * span * 0.03 * speed
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // worker patrol: weak pull toward the queen plus a fine jitter
        for d in 0..<dim {
          let nv = xi[d] +
            0.35 * r2 * (best_pos[d] - xi[d]) +
            (2.0 * RNG::next(current_rng) - 1.0) * span * 0.02 * speed
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      }
      let fv = evaluate(f, cand)
      evaluations = evaluations + 1
      if fv < fits[i] {
        pop[i] = cand
        fits[i] = fv
      }
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Crayfish Optimization Algorithm (2023) =====

///|
pub struct CrayfishConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_crayfish_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CrayfishConfig {
  CrayfishConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
/// Crayfish Optimization Algorithm (2023).  The pond temperature rises from
/// 20 to 30 degrees as evaluations are spent, which raises the chance of a
/// summer resort, while the crayfish factor cof = e^{-t/2500} shrinks all
/// strides:
///   * foraging - move toward the better of two randomly picked crayfish,
///     plus a differential kick and a jitter;
///   * summer resort - contract toward the cave, i.e. the global best, with a
///     fine jitter that tightens the search;
///   * competition - the winner keeps refining near the cave while the loser
///     is pushed away from its rival to stay diverse.
pub fn crayfish_optimization(
  config : CrayfishConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let max_evals = config.max_evaluations
  let current_rng = rng
  let pop_size = if config.population_size < 10 {
    10
  } else {
    config.population_size
  }
  let span = bounds.max - bounds.min
  let max_evals_d = max_evals.to_double()

  let pop : Array[Array[Double]] = []
  let fits : Array[Double] = []
  let mut best_pos : Array[Double] = random_solution(dim, bounds, current_rng)
  let mut best_val = evaluate(f, best_pos)
  let mut evaluations = 1
  pop.push(best_pos)
  fits.push(best_val)
  for _i in 1..<pop_size {
    let p = random_solution(dim, bounds, current_rng)
    let v = evaluate(f, p)
    evaluations = evaluations + 1
    pop.push(p)
    fits.push(v)
    if v < best_val {
      best_val = v
      best_pos = p
    }
  }

  while evaluations < max_evals {
    // pond temperature rises 20 -> 30 degrees
    let temp = 20.0 + 10.0 * evaluations.to_double() / max_evals_d
    let p_resort = 0.3 + 0.2 * (temp - 20.0) / 10.0
    for i in 0..<pop_size {
      if evaluations >= max_evals {
        break
      }
      let xi = pop[i]
      let cof = @math.exp(-evaluations.to_double() / 2500.0)
      let r1 = RNG::next(current_rng)
      let r2 = RNG::next(current_rng)
      let r3 = RNG::next(current_rng)
      let cand : Array[Double] = []
      if r1 < 0.45 {
        // foraging: follow the better of two random crayfish
        let j = (r2 * pop_size.to_double()).to_int() % pop_size
        let k = (r3 * pop_size.to_double()).to_int() % pop_size
        let xj = pop[j]
        let xk = pop[k]
        let xfood = if fits[j] <= fits[k] { xj } else { xk }
        for d in 0..<dim {
          let nv = xi[d] +
            1.8 * cof * (xfood[d] - xi[d]) +
            0.4 * r3 * cof * (xj[d] - xk[d]) +
            (2.0 * RNG::next(current_rng) - 1.0) * span * 0.012 * cof
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else if r1 < 0.45 + p_resort {
        // summer resort: contract toward the cave (global best)
        for d in 0..<dim {
          let nv = xi[d] +
            1.8 * cof * (best_pos[d] - xi[d]) +
            (2.0 * r2 - 1.0) * span * 0.010 * cof
          cand.push(max(bounds.min, min(bounds.max, nv)))
        }
      } else {
        // competition: winner refines, loser retreats from its rival
        let j = (r2 * pop_size.to_double()).to_int() % pop_size
        let xj = pop[j]
        if fits[i] <= fits[j] {
          for d in 0..<dim {
            let nv = xi[d] +
              0.8 * cof * (best_pos[d] - xi[d]) +
              (2.0 * r3 - 1.0) * span * 0.010 * cof
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        } else {
          for d in 0..<dim {
            let nv = xi[d] +
              0.4 * cof * (best_pos[d] - xi[d]) -
              0.3 * r2 * cof * (xj[d] - xi[d]) +
              (2.0 * r3 - 1.0) * span * 0.010 * cof
            cand.push(max(bounds.min, min(bounds.max, nv)))
          }
        }
      }
      pop[i] = cand
      let fv = evaluate(f, cand)
      fits[i] = fv
      evaluations = evaluations + 1
      if fv < best_val {
        best_val = fv
        best_pos = cand
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}
