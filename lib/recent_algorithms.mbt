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
    let _r1 = RNG::next(current_rng)
    let r2 = RNG::next(current_rng)
    let r3 = RNG::next(current_rng)
    let _r4 = RNG::next(current_rng)

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
