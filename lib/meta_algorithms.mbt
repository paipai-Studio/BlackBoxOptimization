// ===== Imperialist Competitive Algorithm (ICA) =====

///|
pub struct ICAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  imperialist_count : Int32
  assimilation_rate : Double
  revolution_rate : Double
}

///|
pub fn make_ica_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> ICAConfig {
  ICAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    imperialist_count: 3,
    assimilation_rate: 0.3,
    revolution_rate: 0.3,
  }
}

///|
pub fn imperialist_competitive_algorithm(
  config : ICAConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let imp_count = config.imperialist_count

  let current_rng = rng
  let countries : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    countries.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = countries[i]
    }
  }

  // Sort and assign imperialists and colonies
  let mut sorted_idx = sort_indices(values)

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Re-sort
    sorted_idx = sort_indices(values)

    // Imperialists: best imp_count countries
    let imperialists : Array[Int32] = []
    for i in 0..<imp_count {
      imperialists.push(sorted_idx[i])
    }

    // Assign colonies to imperialists based on cost
    let colonies : Array[Array[Int32]] = []
    for _i in 0..<imp_count {
      let col : Array[Int32] = []
      colonies.push(col)
    }

    for i in imp_count..<n {
      let idx = sorted_idx[i]
      // Assign to weakest imperialist (simplified: round-robin by probability)
      let r_assign = RNG::next(current_rng)
      let imp_idx = (r_assign * imp_count.to_double()).floor().to_int()
      let target = max_int(0, min_int(imp_count - 1, imp_idx))
      colonies[target].push(idx)
    }

    // Assimilation: colonies move toward their imperialist
    for k in 0..<imp_count {
      if evaluations >= config.max_evaluations {
        break
      }
      let imp_idx = imperialists[k]
      let imp_pos = countries[imp_idx]

      for c in 0..<colonies[k].length() {
        if evaluations >= config.max_evaluations {
          break
        }
        let col_idx = colonies[k][c]

        // Assimilation: move each dimension toward the imperialist by an
        // independent random fraction, plus a decaying exploration offset so
        // colonies keep searching instead of freezing onto the imperialist.
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r_assim = RNG::next(current_rng)
          let r_offset = RNG::next(current_rng)
          let offset = (r_offset - 0.5) *
            2.0 *
            (bounds.max - bounds.min) *
            0.1 *
            (1.0 - iter)
          let pos = countries[col_idx][d] +
            config.assimilation_rate *
            r_assim *
            (imp_pos[d] - countries[col_idx][d]) +
            offset
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
      if evaluations >= config.max_evaluations {
        break
      }
      let r_rev = RNG::next(current_rng)
      if r_rev < config.revolution_rate {
        let dims_to_change = max_int(1, dim / 3)
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r_p = RNG::next(current_rng)
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
    evaluations,
  }
}

// ===== Stochastic Fractal Search (SFS) =====

///|
pub struct SFSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  discovery_iterations : Int32
  learning_iterations : Int32
}

///|
pub fn make_sfs_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SFSConfig {
  SFSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    discovery_iterations: 3,
    learning_iterations: 2,
  }
}

///|
pub fn stochastic_fractal_search(
  config : SFSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let range = bounds.max - bounds.min

  let current_rng = rng
  let points : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    points.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      for d in 0..<dim {
        if evaluations >= config.max_evaluations {
          break
        }

        let sigma = range * 0.1
        let r_gauss = RNG::next(current_rng)
        let r_gauss2 = RNG::next(current_rng)
        let gauss_val = box_muller(r_gauss, r_gauss2)

        // Explore nearby
        let pos1 = points[i][d] + sigma * gauss_val
        let new_pos : Array[Double] = []
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

        if evaluations >= config.max_evaluations {
          break
        }

        // Explore in opposite direction
        let pos2 = points[i][d] - sigma * gauss_val
        let new_pos2 : Array[Double] = []
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
      if evaluations >= config.max_evaluations {
        break
      }

      let r_learn = RNG::next(current_rng)

      // Find a better point
      let mut better_idx = -1
      let mut min_dist = 1.0e308
      for j in 0..<n {
        if j == i {
          continue
        }
        if values[j] < values[i] {
          let dist = euclidean_distance(points[i], points[j])
          if dist < min_dist {
            min_dist = dist
            better_idx = j
          }
        }
      }

      if better_idx >= 0 && min_dist < 1.0e308 {
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let pos = points[i][d] +
            r_learn * (points[better_idx][d] - points[i][d])
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
    evaluations,
  }
}

// ===== Symbiotic Organisms Search (SOS) =====

///|
pub struct SOSConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_sos_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SOSConfig {
  SOSConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn symbiotic_organisms_search(
  config : SOSConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let organisms : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    organisms.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = organisms[i]
    }
  }

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_phase = RNG::next(current_rng)
      let r_j = RNG::next(current_rng)
      let j = ((r_j * (n - 1).to_double()).floor().to_int() + 1) % n

      if r_phase < 0.33 {
        // Mutualism: both organisms benefit
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)
          let mutual_d = (organisms[i][d] + organisms[j][d]) / 2.0
          let pos = organisms[i][d] +
            r1 * (best_pos[d] - mutual_d) +
            r2 * (mutual_d - organisms[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        if new_val < values[i] {
          organisms[i] = new_pos
          values[i] = new_val
          if new_val < best_val {
            best_val = new_val
            best_pos = new_pos
          }
        }
      } else if r_phase < 0.66 {
        // Commensalism: one benefits, other unaffected
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let pos = organisms[i][d] + r1 * (best_pos[d] - organisms[j][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        if new_val < values[i] {
          organisms[i] = new_pos
          values[i] = new_val
          if new_val < best_val {
            best_val = new_val
            best_pos = new_pos
          }
        }
        // Parasitism: one benefits at expense of other
      } else if values[i] < values[j] {
        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let pos = organisms[j][d] + r1 * (organisms[j][d] - best_pos[d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
        organisms[j] = new_pos
        values[j] = evaluate(f, new_pos)
        evaluations = evaluations + 1
        if values[j] < best_val {
          best_val = values[j]
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

// ===== Biogeography-Based Optimization (BBO) =====

///|
pub struct BBOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  mutation_rate : Double
}

///|
pub fn make_bbo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  mut_rate : Double,
) -> BBOConfig {
  BBOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    mutation_rate: mut_rate,
  }
}

///|
pub fn biogeography_based_optimization(
  config : BBOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let mut_rate = config.mutation_rate

  let current_rng = rng
  let habitats : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    habitats.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = habitats[i]
    }
  }

  while evaluations < config.max_evaluations {
    // Sort habitats by fitness and derive rates from rank
    let sorted_idx = sort_indices(values)
    let rank : Array[Int32] = []
    for _r in 0..<n {
      rank.push(0)
    }
    for r in 0..<n {
      rank[sorted_idx[r]] = r
    }

    // Emigration rates: best habitats (rank 0) emigrate most, worst least
    let emigrate_rates : Array[Double] = []
    let immigrate_rates : Array[Double] = []
    let rank_denom = if n > 1 { (n - 1).to_double() } else { 1.0 }
    for i in 0..<n {
      let frac = rank[i].to_double() / rank_denom
      emigrate_rates.push(1.0 - frac)
      immigrate_rates.push(frac)
    }

    // Migration
    let new_habitats : Array[Array[Double]] = []
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_hab : Array[Double] = []
      for d in 0..<dim {
        let r = RNG::next(current_rng)
        if r < immigrate_rates[i] {
          // Roulette-wheel emigrant selection weighted by emigration rate
          // (the previous first-match scan almost always picked rank 0)
          let mut total_rate = 0.0
          for j in 0..<n {
            if j == i {
              continue
            }
            total_rate = total_rate + emigrate_rates[j]
          }
          let mut found_emigrant = -1
          if total_rate > 1.0e-12 {
            let pick = RNG::next(current_rng) * total_rate
            let mut cum = 0.0
            for j in 0..<n {
              if j == i {
                continue
              }
              cum = cum + emigrate_rates[j]
              if pick <= cum {
                found_emigrant = j
                break
              }
            }
          }
          if found_emigrant < 0 {
            for j in 0..<n {
              if j != i {
                found_emigrant = j
                break
              }
            }
          }
          if found_emigrant >= 0 {
            new_hab.push(habitats[found_emigrant][d])
          } else {
            new_hab.push(habitats[i][d])
          }
        } else {
          new_hab.push(habitats[i][d])
        }
      }

      // Mutation
      let final_hab : Array[Double] = []
      for d in 0..<dim {
        let r_mut = RNG::next(current_rng)
        if r_mut < mut_rate {
          let r_rand = RNG::next(current_rng)
          final_hab.push(bounds.min + r_rand * (bounds.max - bounds.min))
        } else {
          final_hab.push(max(bounds.min, min(bounds.max, new_hab[d])))
        }
      }

      let old_val = values[i]
      let new_val = evaluate(f, final_hab)
      evaluations = evaluations + 1

      if new_val < old_val {
        new_habitats.push(final_hab)
        values[i] = new_val

        if new_val < best_val {
          best_val = new_val
          best_pos = final_hab
        }
      } else {
        new_habitats.push(habitats[i])
      }
    }

    if new_habitats.length() > 0 {
      for i in 0..<new_habitats.length() {
        habitats[i] = new_habitats[i]
      }
    }
  }

  OptimizationResult::{
    best_solution: best_pos,
    best_value: best_val,
    evaluations,
  }
}

// ===== Elephant Herding Optimization (EHO) =====

///|
pub struct EHOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
  clans : Int32
  alpha : Double
}

///|
pub fn make_eho_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
  n_clans : Int32,
  a : Double,
) -> EHOConfig {
  EHOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
    clans: n_clans,
    alpha: a,
  }
}

///|
pub fn elephant_herding_optimization(
  config : EHOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let n_clans = config.clans
  let alpha = config.alpha

  let current_rng = rng
  let elephants : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    elephants.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = elephants[i]
    }
  }

  while evaluations < config.max_evaluations {
    // Matriarch update: each clan's best leads
    let clan_size = n / n_clans
    for c in 0..<n_clans {
      if evaluations >= config.max_evaluations {
        break
      }
      let start = c * clan_size
      let end = if start + clan_size < n { start + clan_size } else { n }

      // Find clan best
      let mut clan_best_idx = start
      for i in (start + 1)..<end {
        if values[i] < values[clan_best_idx] {
          clan_best_idx = i
        }
      }

      // Update clan members toward clan best
      for i in start..<end {
        if evaluations >= config.max_evaluations {
          break
        }
        if i == clan_best_idx {
          continue
        }

        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let pos = elephants[i][d] +
            alpha * r1 * (elephants[clan_best_idx][d] - elephants[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }
        elephants[i] = new_pos
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        values[i] = new_val

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      }
    }

    // Group leader update: global best leads clan leaders
    // Replace worst individual in each clan toward best
    for c in 0..<n_clans {
      if evaluations >= config.max_evaluations {
        break
      }
      let start = c * clan_size
      let end = if start + clan_size < n { start + clan_size } else { n }

      // Find clan worst
      let mut clan_worst_idx = start
      for i in (start + 1)..<end {
        if values[i] > values[clan_worst_idx] {
          clan_worst_idx = i
        }
      }

      let r2 = RNG::next(current_rng)
      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r3 = RNG::next(current_rng)
        let pos = best_pos[d] +
          r3 * (bounds.max - bounds.min) * 0.1 * (r2 - 0.5) * 2.0
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }
      elephants[clan_worst_idx] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[clan_worst_idx] = new_val

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

// ===== Equilibrium Optimizer (EO) =====

///|
pub struct EOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_eo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> EOConfig {
  EOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn equilibrium_optimizer(
  config : EOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let particles : Array[Array[Double]] = []
  let values : Array[Double] = []
  let velocities : Array[Array[Double]] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    particles.push(sol)
    values.push(evaluate(f, sol))

    let vel : Array[Double] = []
    for _d in 0..<dim {
      let _r = RNG::next(current_rng)
      vel.push(0.0)
    }
    velocities.push(vel)
  }

  let mut evaluations = n

  // Equilibrium pool: store best 4 particles
  let eq_sorted_idx = sort_indices(values)
  let eq_pool : Array[Int32] = []
  let pool_size = if 4 < eq_sorted_idx.length() {
    4
  } else {
    eq_sorted_idx.length()
  }
  for i in 0..<pool_size {
    eq_pool.push(eq_sorted_idx[i])
  }
  let _eq_pool_len = eq_pool.length()

  let mut best_pos : Array[Double] = particles[eq_sorted_idx[0]]
  let mut best_val = values[eq_sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    // Exponential decay factor
    let a1 = 2.0
    let a2 = 1.0
    let t = (1.0 - iter) * (a2 * iter)
    let exp_decay = exp_approx(-a1 * t)

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      // Random selection from equilibrium pool
      let r_pool = RNG::next(current_rng)
      let pool_idx = (r_pool * eq_pool.length().to_double()).floor().to_int() %
        eq_pool.length()
      let eq_idx = eq_pool[pool_idx]

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)
        let r_gauss = RNG::next(current_rng)
        let gauss = box_muller(r_gauss, r2)

        let candidate = particles[eq_idx][d]
        let diff = abs_val(candidate - particles[i][d])
        let sign = if r1 < 0.5 { -1.0 } else { 1.0 }

        let new_val_d = candidate +
          sign * (r1 * diff * exp_decay) +
          gauss * 0.1 * (r2 * (particles[i][d] - candidate))
        new_pos.push(max(bounds.min, min(bounds.max, new_val_d)))
      }

      particles[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

      // Update equilibrium pool
      if new_val < values[eq_pool[3]] {
        eq_pool[3] = i
        // Sort pool by values
        for p in 0..<4 {
          for q in (p + 1)..<4 {
            if values[eq_pool[q]] < values[eq_pool[p]] {
              let temp = eq_pool[p]
              eq_pool[p] = eq_pool[q]
              eq_pool[q] = temp
            }
          }
        }
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

// ===== Wind Driven Optimization (WDO) =====

///|
pub struct WDOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_wdo_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> WDOConfig {
  WDOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn wind_driven_optimization(
  config : WDOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let friction = 0.1
  let gravity = 0.2
  let wind_factor = 0.4

  let current_rng = rng
  let particles : Array[Array[Double]] = []
  let velocities : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    particles.push(sol)

    let vel : Array[Double] = []
    for _d in 0..<dim {
      let v = RNG::next(current_rng)
      vel.push((v - 0.5) * 2.0)
    }
    velocities.push(vel)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = particles[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    let iter = evaluations.to_double() / config.max_evaluations.to_double()
    let eta = 1.0 - iter

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_gauss = RNG::next(current_rng)
      let r_gauss2 = RNG::next(current_rng)
      let gauss = box_muller(r_gauss, r_gauss2)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)

        let vel = (1.0 - friction) * velocities[i][d] +
          gravity * (best_pos[d] - particles[i][d]) +
          wind_factor * r1 * (best_pos[d] - particles[i][d]) +
          eta * gauss * 0.1 * r2

        velocities[i][d] = vel
        let pos = particles[i][d] + vel
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      particles[i] = new_pos
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

// ===== Jaya Algorithm =====

///|
pub struct JAYAConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_jaya_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> JAYAConfig {
  JAYAConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn jaya_optimization(
  config : JAYAConfig,
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
  let mut sorted_idx = sort_indices(values)
  let mut best_pos : Array[Double] = population[sorted_idx[0]]
  let mut best_val = values[sorted_idx[0]]
  let mut global_best_pos : Array[Double] = best_pos
  let mut global_best_val : Double = best_val
  let mut worst_pos : Array[Double] = population[sorted_idx[n - 1]]
  let mut _worst_val = values[sorted_idx[n - 1]]

  while evaluations < config.max_evaluations {
    sorted_idx = sort_indices(values)
    best_pos = population[sorted_idx[0]]
    best_val = values[sorted_idx[0]]
    worst_pos = population[sorted_idx[n - 1]]
    _worst_val = values[sorted_idx[n - 1]]

    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)
        let r2 = RNG::next(current_rng)

        let cur_pos = population[i][d]
        let pos = cur_pos +
          r1 * (best_pos[d] - cur_pos) -
          r2 * (worst_pos[d] - cur_pos)
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
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
        if new_val < global_best_val {
          global_best_val = new_val
          global_best_pos = new_pos
        }
      }
    }
  }

  OptimizationResult::{
    best_solution: global_best_pos,
    best_value: global_best_val,
    evaluations,
  }
}

// ===== Soccer League Competition (SLC) =====

///|
pub struct SLCConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_slc_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> SLCConfig {
  SLCConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn soccer_league_competition(
  config : SLCConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size
  let n_teams = 4
  let team_size = n / n_teams

  let current_rng = rng
  let players : Array[Array[Double]] = []
  let values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    players.push(sol)
    values.push(evaluate(f, sol))
  }

  let mut evaluations = n
  let mut best_pos : Array[Double] = []
  let mut best_val = 1.0e308

  for i in 0..<n {
    if values[i] < best_val {
      best_val = values[i]
      best_pos = players[i]
    }
  }

  while evaluations < config.max_evaluations {
    for t in 0..<n_teams {
      if evaluations >= config.max_evaluations {
        break
      }
      let start = t * team_size
      let end = if start + team_size < n { start + team_size } else { n }

      let mut team_best_idx = start
      for i in (start + 1)..<end {
        if values[i] < values[team_best_idx] {
          team_best_idx = i
        }
      }

      for i in start..<end {
        if evaluations >= config.max_evaluations {
          break
        }
        if i == team_best_idx {
          continue
        }

        let new_pos : Array[Double] = []
        for d in 0..<dim {
          let r1 = RNG::next(current_rng)
          let r2 = RNG::next(current_rng)

          let pos = players[i][d] +
            r1 * (players[team_best_idx][d] - players[i][d]) +
            r2 * (best_pos[d] - players[i][d])
          new_pos.push(max(bounds.min, min(bounds.max, pos)))
        }

        players[i] = new_pos
        let new_val = evaluate(f, new_pos)
        evaluations = evaluations + 1
        values[i] = new_val

        if new_val < best_val {
          best_val = new_val
          best_pos = new_pos
        }
      }
    }

    let r_exchange = RNG::next(current_rng)
    if r_exchange < 0.3 {
      let r_team1 = RNG::next(current_rng)
      let r_team2 = RNG::next(current_rng)
      let t1 = (r_team1 * n_teams.to_double()).floor().to_int() % n_teams
      let t2 = (r_team2 * n_teams.to_double()).floor().to_int() % n_teams
      if t1 != t2 {
        let s1 = t1 * team_size
        let s2 = t2 * team_size
        let r_p1 = RNG::next(current_rng)
        let r_p2 = RNG::next(current_rng)
        let p1 = s1 +
          (r_p1 * team_size.to_double()).floor().to_int() % team_size
        let p2 = s2 +
          (r_p2 * team_size.to_double()).floor().to_int() % team_size
        if p1 < n && p2 < n {
          let temp_pos = players[p1]
          let temp_val = values[p1]
          players[p1] = players[p2]
          values[p1] = values[p2]
          players[p2] = temp_pos
          values[p2] = temp_val
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

// ===== Gaussian Bare-bones PSO (GBPSO) =====

///|
pub struct GBPSOConfig {
  max_evaluations : Int32
  dimension : Int32
  bounds : Bounds
  population_size : Int32
}

///|
pub fn make_gbpso_config(
  max_evals : Int32,
  dim : Int32,
  b : Bounds,
  pop_sz : Int32,
) -> GBPSOConfig {
  GBPSOConfig::{
    max_evaluations: max_evals,
    dimension: dim,
    bounds: b,
    population_size: pop_sz,
  }
}

///|
pub fn gaussian_bare_bones_pso(
  config : GBPSOConfig,
  f : (Array[Double]) -> Double,
  rng : RNG,
) -> OptimizationResult {
  let dim = config.dimension
  let bounds = config.bounds
  let n = config.population_size

  let current_rng = rng
  let particles : Array[Array[Double]] = []
  let pbest : Array[Array[Double]] = []
  let values : Array[Double] = []
  let pbest_values : Array[Double] = []

  for _i in 0..<n {
    let sol = random_solution(dim, bounds, current_rng)
    particles.push(sol)
    let val = evaluate(f, sol)
    values.push(val)
    pbest.push(sol)
    pbest_values.push(val)
  }

  let mut evaluations = n
  let sorted_idx = sort_indices(values)
  let mut gbest_pos : Array[Double] = particles[sorted_idx[0]]
  let mut gbest_val = values[sorted_idx[0]]

  while evaluations < config.max_evaluations {
    for i in 0..<n {
      if evaluations >= config.max_evaluations {
        break
      }

      let r_gauss = RNG::next(current_rng)
      let r_gauss2 = RNG::next(current_rng)
      let gauss = box_muller(r_gauss, r_gauss2)

      let new_pos : Array[Double] = []
      for d in 0..<dim {
        let r1 = RNG::next(current_rng)

        let pos = gbest_pos[d] + gauss * (gbest_pos[d] - pbest[i][d]) * r1
        new_pos.push(max(bounds.min, min(bounds.max, pos)))
      }

      particles[i] = new_pos
      let new_val = evaluate(f, new_pos)
      evaluations = evaluations + 1
      values[i] = new_val

      if new_val < pbest_values[i] {
        pbest_values[i] = new_val
        pbest[i] = new_pos
      }

      if new_val < gbest_val {
        gbest_val = new_val
        gbest_pos = new_pos
      }
    }
  }

  OptimizationResult::{
    best_solution: gbest_pos,
    best_value: gbest_val,
    evaluations,
  }
}
