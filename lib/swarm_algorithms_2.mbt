// ===== Harris Hawks Optimization (HHO) =====

///|
pub struct HHOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_hho_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> HHOConfig {
  HHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn harris_hawks_optimization(
  config : HHOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let hawks : Array[Array[Double]] = []
  let hawk_values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    hawks.push(sol)
    hawk_values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let _r_escape = RNG::next(current_rng)

      if abs_val(e0) >= 1.0 {
        // Exploration phase
        let r_q = RNG::next(current_rng)
        let r_rand = RNG::next(current_rng)

        if r_q >= 0.5 {
          let r_leader = RNG::next(current_rng)
          let new_pos : Array[Double] = []
          for _d in 0..<dim {
            let pos = bounds.min + r_rand * (bounds.max - bounds.min)
            new_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          let _r_sub = RNG::next(current_rng)
          let final_pos : Array[Double] = []
          for d in 0..<dim {
            let pos = new_pos[d] + r_leader * abs_val(new_pos[d] - hawks[i][d])
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        } else {
          let avg_pos : Array[Double] = []
          for d in 0..<dim {
            let mut sum = 0.0
            for j in 0..<n {
              sum = sum + hawks[j][d]
            }
            avg_pos.push(sum / n.to_double())
          }

          let r_ld1 = RNG::next(current_rng)
          let r_ld2 = RNG::next(current_rng)

          let final_pos : Array[Double] = []
          for d in 0..<dim {
            let pos = r_ld1 * best_pos[d] - r_ld2 * avg_pos[d]
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        }
      } else {
        // Exploitation phase
        let r_escape2 = RNG::next(current_rng)

        if r_escape2 >= 0.5 {
          let r_cf = RNG::next(current_rng)

          let final_pos : Array[Double] = []
          for d in 0..<dim {
            let pos = best_pos[d] +
              abs_val(best_pos[d] - hawks[i][d]) * e0 * r_cf
            final_pos.push(max(bounds.min, min(bounds.max, pos)))
          }
          hawks[i] = final_pos
        } else {
          let r_cf1 = RNG::next(current_rng)
          let r_cf2 = RNG::next(current_rng)

          let final_pos : Array[Double] = []
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
    evaluations,
  }
}

// ===== Flower Pollination Algorithm (FPA) =====

///|
pub struct FPAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  switch_probability : Double
}

///|
pub fn make_fpa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> FPAConfig {
  FPAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    switch_probability: 0.8,
  }
}

///|
pub fn flower_pollination_algorithm(
  config : FPAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let flowers : Array[Array[Double]] = []
  let flower_values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    flowers.push(sol)
    flower_values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if flower_values[i] < best_val {
      best_val = flower_values[i]
      best_pos = flowers[i]
    }
  }

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_switch = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      if r_switch < config.switch_probability {
        // Global pollination
        let r_levy1 = RNG::next(current_rng)
        let r_levy2 = RNG::next(current_rng)

        for d in 0..<dim {
          let levy_step = abs_val(r_levy1 * (best_pos[d] - flowers[i][d]))
          let pos = flowers[i][d] + levy_step * (r_levy2 - 0.5) * 2.0
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Local pollination
        let r_local = RNG::next(current_rng)
        let r_idx1 = RNG::next(current_rng)
        let r_idx2 = RNG::next(current_rng)

        let j = (r_idx1 * n.to_double()).floor().to_int()
        let k = (r_idx2 * n.to_double()).floor().to_int()
        let j_idx = max_int(0, min_int(n - 1, j))
        let k_idx = max_int(0, min_int(n - 1, k))

        for d in 0..<dim {
          let pos = flowers[i][d] +
            r_local * (flowers[j_idx][d] - flowers[k_idx][d])
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
    evaluations,
  }
}

// ===== Dragonfly Algorithm (DA) =====

///|
pub struct DAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  separation_weight : Double
  alignment_weight : Double
  cohesion_weight : Double
}

///|
pub fn make_da_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> DAConfig {
  DAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    separation_weight: 1.0,
    alignment_weight: 1.0,
    cohesion_weight: 1.0,
  }
}

///|
pub fn dragonfly_algorithm(
  config : DAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let range = bounds.max - bounds.min

  let current_rng = rng
  let dragonflies : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let df_values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    dragonflies.push(sol)
    df_values.push(evaluate(f, sol))

    let init_vel : Array[Double] = []
    for _d in 0..<dim {
      let v = RNG::next(current_rng)
      init_vel.push((v * 2.0 - 1.0) * range * 0.01)
    }
    velocities.push(init_vel)
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
    let e_weight = 1.0 - iter

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_sep = RNG::next(current_rng)
      let r_ali = RNG::next(current_rng)
      let r_coh = RNG::next(current_rng)
      let r_food = RNG::next(current_rng)

      // Find actual neighbors and compute metrics
      let sep : Array[Double] = []
      let ali : Array[Double] = []
      let coh : Array[Double] = []
      for d in 0..<dim {
        let mut sep_sum = 0.0
        let mut ali_sum = 0.0
        let mut coh_sum = 0.0
        let mut count = 0

        for j in 0..<n {
          if j == i {
            continue
          }
          sep_sum = sep_sum + (dragonflies[i][d] - dragonflies[j][d])
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
      let new_pos : Array[Double] = []
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
    evaluations,
  }
}

// ===== Crow Search Algorithm (CSA) =====

///|
pub struct CSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  flight_length : Double
  awareness_probability : Double
}

///|
pub fn make_csa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> CSAConfig {
  CSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    flight_length: 2.0,
    awareness_probability: 0.1,
  }
}

///|
pub fn crow_search_algorithm(
  config : CSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let crows : Array[Array[Double]] = []
  let crow_values : Array[Double] = []
  let memory : Array[Array[Double]] = []
  let memory_values : Array[Double] = []

  for i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    crows.push(sol)
    crow_values.push(evaluate(f, sol))
    memory.push(sol)
    memory_values.push(crow_values[i])
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if crow_values[i] < best_val {
      best_val = crow_values[i]
      best_pos = crows[i]
    }
  }

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_awareness = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      if r_awareness >= config.awareness_probability {
        // Follow another crow
        let r_target = RNG::next(current_rng)
        let target = (r_target * n.to_double()).floor().to_int()
        let target_idx = max_int(0, min_int(n - 1, target))

        let r_flight = RNG::next(current_rng)

        for d in 0..<dim {
          let pos = crows[i][d] +
            config.flight_length *
            r_flight *
            (memory[target_idx][d] - crows[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Random search
        let r_rand = RNG::next(current_rng)
        for _d in 0..<dim {
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
    evaluations,
  }
}

// ===== Salp Swarm Algorithm (SalpSA) =====

///|
pub struct SalpSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_salp_sa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SalpSAConfig {
  SalpSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn salp_swarm_algorithm(
  config : SalpSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let salps : Array[Array[Double]] = []
  let salp_values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    salps.push(sol)
    salp_values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  // Find initial best
  let mut _best_idx = 0
  for i in 0..<n {
    if salp_values[i] < best_val {
      best_val = salp_values[i]
      best_pos = salps[i]
      _best_idx = i
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let c1 = 2.0 * exp_approx(-(4.0 * iter) * (4.0 * iter))

    let new_salps : Array[Array[Double]] = []

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_c2 = RNG::next(current_rng)
      let r_c3 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }
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
    evaluations,
  }
}

// ===== Monarch Butterfly Optimization (MBO) =====

///|
pub struct MBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  migration_ratio : Double
  perturbation : Double
}

///|
pub fn make_mbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MBOConfig {
  MBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    migration_ratio: 0.5,
    perturbation: 0.3,
  }
}

///|
pub fn monarch_butterfly_optimization(
  config : MBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let n1 = (n.to_double() * config.migration_ratio).to_int()
  let n2 = n - n1

  let current_rng = rng
  let butterflies1 : Array[Array[Double]] = []
  let values1 : Array[Double] = []
  let butterflies2 : Array[Array[Double]] = []
  let values2 : Array[Double] = []

  for _i in 0..<n1 {
    let sol = random_solution(dim, bounds, current_rng)
    butterflies1.push(sol)
    values1.push(evaluate(f, sol))
  }
  for _i in 0..<n2 {
    let sol = random_solution(dim, bounds, current_rng)
    butterflies2.push(sol)
    values2.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let r_mig = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      if r_mig < config.migration_ratio {
        for d in 0..<dim {
          let r_levy = RNG::next(current_rng)
          let pos = best_pos[d] + r_levy * (butterflies1[i][d] - best_pos[d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        let r_rand = RNG::next(current_rng)
        for _d in 0..<dim {
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
      if evaluations >= config.max_evaluations {
        break
      }

      let r_adj = RNG::next(current_rng)

      let new_pos : Array[Double] = []
      if r_adj < config.perturbation {
        for d in 0..<dim {
          let r_mut = RNG::next(current_rng)
          let pos = butterflies2[i][d] +
            (r_mut - 0.5) * (bounds.max - bounds.min) * 0.1
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        let r_swap = RNG::next(current_rng)
        let partner = (r_swap * n1.to_double()).floor().to_int()
        let p_idx = max_int(0, min_int(n1 - 1, partner))
        for d in 0..<dim {
          let pos = butterflies2[i][d] +
            iter * (butterflies1[p_idx][d] - butterflies2[i][d])
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
    evaluations,
  }
}

// ===== Artificial Ecosystem-based Optimization (AEO) =====

///|
pub struct AEOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_aeo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> AEOConfig {
  AEOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn artificial_ecosystem_optimization(
  config : AEOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let individuals : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    individuals.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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

    // Update positions
    for rank in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let idx = sorted_indices[rank]
      let r_role = RNG::next(current_rng)

      let new_pos : Array[Double] = []

      if rank == 0 {
        // Producer (best) - exploit
        for d in 0..<dim {
          let r_prod = RNG::next(current_rng)
          let pos = individuals[idx][d] +
            r_prod * (best_pos[d] - individuals[idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else if r_role < 0.5 {
        // Consumer - explore toward better solution
        let better_idx = sorted_indices[rank - 1]
        for d in 0..<dim {
          let r_cons = RNG::next(current_rng)
          let pos = individuals[idx][d] +
            r_cons * (individuals[better_idx][d] - individuals[idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Decomposer - random search
        for d in 0..<dim {
          let r_decomp = RNG::next(current_rng)
          let pos = individuals[idx][d] +
            (r_decomp - 0.5) * (bounds.max - bounds.min) * (1.0 - iter)
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
    evaluations,
  }
}

// ===== Invasive Weed Optimization (IWO) =====

///|
pub struct IWOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  initial_population : Int32
  max_population : Int32
  min_seeds : Int32
  max_seeds : Int32
  sigma_initial : Double
  sigma_final : Double
  nonlinear : Double
}

///|
pub fn make_iwo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  init_pop : Int32,
  max_pop : Int32,
  min_seed : Int32,
  max_seed : Int32,
  sig_init : Double,
  sig_final : Double,
  n_nonlinear : Double,
) -> IWOConfig {
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
    nonlinear: n_nonlinear,
  }
}

///|
pub fn invasive_weed_optimization(
  config : IWOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds

  let current_rng = rng
  let population : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<config.initial_population {
    let sol = random_solution(dim, bounds, current_rng)
    population.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = config.initial_population
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<population.length() {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = population[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter_frac = evaluations.to_double() / config.max_evaluations.to_double()
    let sigma = config.sigma_final +
      (config.sigma_initial - config.sigma_final) *
      pow_approx(1.0 - iter_frac, config.nonlinear)

    let new_pop : Array[Array[Double]] = []
    let new_vals : Array[Double] = []

    let sorted_idx = sort_indices(values)
    let n = population.length()

    for rank in 0..<n {
      if evaluations + new_pop.length() >= config.max_evaluations {
        break
      }

      let idx = sorted_idx[rank]
      let norm_rank = (rank.to_double() + 1.0) / n.to_double()
      let seeds_d : Double = config.min_seeds.to_double() +
        norm_rank *
        (config.max_seeds.to_double() - config.min_seeds.to_double())
      let seeds = seeds_d.round().to_int()

      for _s in 0..<seeds {
        if evaluations + new_pop.length() >= config.max_evaluations {
          break
        }

        let child : Array[Double] = []
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)
          let gauss = box_muller(r1, r2)
          let pos = population[idx][d] + sigma * gauss
          child.push(max(bounds.min, min(bounds.max, pos)))
        }

        new_pop.push(child)
        new_vals.push(evaluate(f, child))
      }
    }

    // Combine and keep best max_population
    let combined : Array[Array[Double]] = []
    let combined_vals : Array[Double] = []
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
    evaluations,
  }
}

// ===== Bacteria Foraging Optimization (BFO) =====

///|
pub struct BFOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  chemotactic_steps : Int32
  swim_length : Int32
  tumble_step : Double
}

///|
pub fn make_bfo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  chem_steps : Int32,
  swim_len : Int32,
  tumble : Double,
) -> BFOConfig {
  BFOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    chemotactic_steps: chem_steps,
    swim_length: swim_len,
    tumble_step: tumble,
  }
}

///|
pub fn bacteria_foraging_optimization(
  config : BFOConfig,
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
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = bacteria[i]
    }
  }

  while evaluations < config.max_evaluations {
    for b_idx in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // Tumble: generate random direction
      let delta : Array[Double] = []
      for _d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
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
      for _s in 0..<config.swim_length {
        if evaluations >= config.max_evaluations {
          break
        }

        let new_pos : Array[Double] = []
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
        if evaluations >= config.max_evaluations {
          break
        }
        let worst_idx = idx[n - 1 - i]
        let sol = random_solution(dim, bounds, current_rng)
        bacteria[worst_idx] = sol
        values[worst_idx] = evaluate(f, sol)
        evaluations = evaluations + 1
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Shuffled Frog Leaping Algorithm (SFLA) =====

///|
pub struct SFLAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  memeplexes : Int32
  sub_iterations : Int32
}

///|
pub fn make_sfla_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  memes : Int32,
  sub_iter : Int32,
) -> SFLAConfig {
  SFLAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    memeplexes: memes,
    sub_iterations: sub_iter,
  }
}

///|
pub fn shuffled_frog_leaping_algorithm(
  config : SFLAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let memes = config.memeplexes
  let memeplex_sz = n / memes

  let current_rng = rng
  let frogs : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    frogs.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let memeplex : Array[Int32] = []
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
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
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
        let new_pos2 : Array[Double] = []
        for d in 0..<dim {
          let r = RNG::next(current_rng)
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
          let sol = random_solution(dim, bounds, current_rng)
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
    evaluations,
  }
}

// ===== Glowworm Swarm Optimization (GSO) =====

///|
pub struct GSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  luciferin_decay : Double
  luciferin_enhance : Double
  neighbor_threshold : Int32
  initial_luciferin : Double
  initial_radius : Double
}

///|
pub fn make_gso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  luc_decay : Double,
  luc_enhance : Double,
  neigh_thresh : Int32,
  init_luc : Double,
  init_rad : Double,
) -> GSOConfig {
  GSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    luciferin_decay: luc_decay,
    luciferin_enhance: luc_enhance,
    neighbor_threshold: neigh_thresh,
    initial_luciferin: init_luc,
    initial_radius: init_rad,
  }
}

///|
pub fn glowworm_swarm_optimization(
  config : GSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let worms : Array[Array[Double]] = []
  let values : Array[Double] = []
  let luciferin : Array[Double] = []
  let radius : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    worms.push(sol)
    values.push(evaluate(f, sol))
    luciferin.push(config.initial_luciferin)
    radius.push(config.initial_radius)
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      luciferin[i] = (1.0 - config.luciferin_decay) * luciferin[i] +
        config.luciferin_enhance * (1.0 / (values[i] + 1.0))
    }

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // Find neighbors with higher luciferin within radius
      let neighbors : Array[Int32] = []
      for j in 0..<n {
        if j == i {
          continue
        }
        let dist = euclidean_distance(worms[i], worms[j])
        if dist <= radius[i] && luciferin[j] > luciferin[i] {
          neighbors.push(j)
        }
      }

      if neighbors.length() == 0 {
        // Random movement
        let _r = RNG::next(current_rng)
        let step = radius[i] * 0.3
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let rd = RNG::next(current_rng)
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
        let weights : Array[Double] = []
        let mut total = 0.0
        for j in 0..<neighbors.length() {
          let w = luciferin[neighbors[j]] - luciferin[i]
          weights.push(w)
          total = total + w
        }

        let sel = RNG::next(current_rng)
        let mut cum = 0.0
        let mut chosen = neighbors[0]
        for j in 0..<neighbors.length() {
          cum = cum + weights[j] / total
          if sel <= cum {
            chosen = neighbors[j]
            break
          }
        }

        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let rd = RNG::next(current_rng)
          let step = min(1.0, radius[i] * rd)
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
      }

      // Adaptive sensor range applies to both branches: an isolated glowworm
      // (no neighbours) must widen its reach to find any.
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

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Gravitational Search Algorithm (GSA) =====

///|
pub struct GSAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  gravity_constant : Double
  decay_alpha : Double
}

///|
pub fn make_gsa_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  g0 : Double,
  alpha : Double,
) -> GSAConfig {
  GSAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    gravity_constant: g0,
    decay_alpha: alpha,
  }
}

///|
pub fn gravitational_search_algorithm(
  config : GSAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let agents : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    agents.push(sol)
    values.push(evaluate(f, sol))

    let vel : Array[Double] = []
    for _d in 0..<dim {
      let _r = RNG::next(current_rng)
      vel.push(0.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = agents[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // G decays exponentially over time (standard GSA); decay_alpha scales the rate
    let g = config.gravity_constant *
      exp_approx(-20.0 * config.decay_alpha * iter)

    // Refreshed best/worst of the current population (minimization: best = min fitness)
    let mut cur_best = values[0]
    let mut cur_worst = values[0]
    for i in 1..<n {
      if values[i] < cur_best {
        cur_best = values[i]
      }
      if values[i] > cur_worst {
        cur_worst = values[i]
      }
    }

    // Gravitational masses: best agent gets mass ~1, worst ~0
    let masses : Array[Double] = []
    let mut total_mass = 0.0
    for i in 0..<n {
      let m = (cur_worst - values[i]) / (cur_worst - cur_best + 1.0e-30)
      masses.push(m)
      total_mass = total_mass + m
    }
    let mass_denom = total_mass + 1.0e-30

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // Gravitational acceleration on agent i from all other agents
      let force : Array[Double] = []
      for _d in 0..<dim {
        force.push(0.0)
      }
      for d in 0..<dim {
        let mut f_d = 0.0
        for j in 0..<n {
          if j == i {
            continue
          }
          let dist = euclidean_distance(agents[i], agents[j])
          if dist < 1.0e-10 {
            continue
          }
          let r_force = RNG::next(current_rng)
          f_d = f_d +
            r_force *
            g *
            (masses[j] / mass_denom) *
            (agents[j][d] - agents[i][d]) /
            (dist * dist + 1.0e-10)
        }
        force[d] = f_d
      }

      // Update velocity and position
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r_vel = RNG::next(current_rng)
        let vel = r_vel * velocities[i][d] + force[d]
        velocities[i][d] = vel
        let pos = agents[i][d] + vel
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
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Teaching-Learning-Based Optimization (TLBO) =====

///|
pub struct TLBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_tlbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> TLBOConfig {
  TLBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn teaching_learning_based_optimization(
  config : TLBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let students : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    students.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      if values[i] < values[teacher_idx] {
        teacher_idx = i
      }
      if values[i] > values[worst_idx] {
        worst_idx = i
      }
    }
    let teacher = students[teacher_idx]

    // Calculate mean of all students
    let mean : Array[Double] = []
    for d in 0..<dim {
      let mut sum = 0.0
      for i in 0..<n {
        sum = sum + students[i][d]
      }
      mean.push(sum / n.to_double())
    }

    // Teacher phase
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }
      if i == teacher_idx {
        continue
      }

      let tf_val = RNG::next(current_rng)
      let tf = 1.0 + (if tf_val < 0.5 { 0.0 } else { 1.0 })
      let r1 = RNG::next(current_rng)

      let new_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let partner_r = RNG::next(current_rng)
      let partner = (partner_r * n.to_double()).floor().to_int() % n
      let partner_idx = if partner < 0 { n + partner } else { partner }

      let r1 = RNG::next(current_rng)
      let new_pos : Array[Double] = []

      if values[i] < values[partner_idx] {
        for d in 0..<dim {
          let pos = students[i][d] +
            r1 * (students[i][d] - students[partner_idx][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        for d in 0..<dim {
          let pos = students[i][d] +
            r1 * (students[partner_idx][d] - students[i][d])
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
    evaluations,
  }
}

// ===== Animal Migration Optimization (AMO) =====

///|
pub struct AMOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  migration_interval : Int32
}

///|
pub fn make_amo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  mig_interval : Int32,
) -> AMOConfig {
  AMOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    migration_interval: mig_interval,
  }
}

///|
pub fn animal_migration_optimization(
  config : AMOConfig,
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
  let mut best_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let r_follow = RNG::next(current_rng)
      let new_pos : Array[Double] = []

      if r_follow < 0.7 {
        // Follow best
        let r1 = RNG::next(current_rng)
        for d in 0..<dim {
          let pos = animals[i][d] + r1 * (best_pos[d] - animals[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
      } else {
        // Explore randomly toward a random neighbor
        let r_idx = RNG::next(current_rng)
        let j = (r_idx * n.to_double()).floor().to_int() % n
        let j_idx = if j < 0 { n + j } else { j }
        let r2 = RNG::next(current_rng)
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
    if migration_count % config.migration_interval == 0 &&
      evaluations < config.max_evaluations {
      let idx = sort_indices(values)
      let replace_count = n / 5
      for k in 0..<replace_count {
        if evaluations >= config.max_evaluations {
          break
        }
        let worst_idx = idx[n - 1 - k]
        let sol = random_solution(dim, bounds, current_rng)
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
    evaluations,
  }
}

// ===== Multi-Verse Optimizer (MVO) =====

///|
pub struct MVOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  wormhole_exploitation : Double
  travelling_distance_rate : Double
}

///|
pub fn make_mvo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> MVOConfig {
  MVOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    wormhole_exploitation: 1.0,
    travelling_distance_rate: 0.6,
  }
}

///|
pub fn multi_verse_optimizer(
  config : MVOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let range = bounds.max - bounds.min

  let current_rng = rng
  let universes : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    universes.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = universes[i]
    }
  }

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()

    // Wormhole Existence Probability: decreases from high to low
    let wep = 2.0 * iter * iter * iter + 0.2
    // Travelling Distance Rate: decreases from 1 to 0
    let tdr = 1.0 - iter * iter

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r_wep = RNG::next(current_rng)
        let r_tdr = RNG::next(current_rng)
        let r_idx = RNG::next(current_rng)

        let mut pos : Double = universes[i][d]
        if r_wep < wep {
          // White/black hole tunnel: the worse universe adopts the better one
          let j = (r_idx * n.to_double()).floor().to_int()
          let j_idx = if j >= n { n - 1 } else if j < 0 { 0 } else { j }
          if values[j_idx] < values[i] {
            pos = universes[j_idx][d]
          }
          // Wormhole: exploit best solution
        } else if r_tdr < 0.5 {
          pos = best_pos[d] + tdr * range * r1
        } else {
          pos = best_pos[d] - tdr * range * r1
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
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Krill Herd (KH) =====

///|
pub struct KHConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  max_speed : Double
  foraging_speed : Double
}

///|
pub fn make_kh_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> KHConfig {
  KHConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    max_speed: 0.01,
    foraging_speed: 0.02,
  }
}

///|
pub fn krill_herd(
  config : KHConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let krills : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    krills.push(sol)
    values.push(evaluate(f, sol))

    let vel : Array[Double] = []
    for _d in 0..<dim {
      let _r = RNG::next(current_rng)
      vel.push(0.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      // Find neighbors
      let neighbors : Array[Int32] = []
      let max_neighbors = 2
      for j in 0..<n {
        if j == i {
          continue
        }
        if neighbors.length() >= max_neighbors {
          break
        }
        let dist = euclidean_distance(krills[i], krills[j])
        if dist < 1.0 {
          neighbors.push(j)
        }
      }

      // Direction of motion (induced + foraging)
      let motion : Array[Double] = []
      for d in 0..<dim {
        // Induced motion: follow neighbors
        let mut induced = 0.0
        for ni in 0..<neighbors.length() {
          let j = neighbors[ni]
          let dist = euclidean_distance(krills[i], krills[j])
          if dist > 1.0e-10 {
            induced = induced + (krills[j][d] - krills[i][d]) / dist
          }
        }
        // Foraging motion: toward food (best)
        let food = best_pos[d] - krills[i][d]
        // Random diffusion
        let r_diff = RNG::next(current_rng)

        let dir = induced * config.foraging_speed +
          food * 0.5 +
          r_diff * step_size
        motion.push(dir)
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
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
    evaluations,
  }
}
