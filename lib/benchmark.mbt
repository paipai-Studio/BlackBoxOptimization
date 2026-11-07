// ===== 测试函数 =====

///|
pub fn sphere_function(x : Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<x.length() {
    sum = sum + x[i] * x[i]
  }
  sum
}

///|
pub fn rastrigin_function(x : Array[Double]) -> Double {
  let n : Double = x.length().to_double()
  let mut sum = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum = sum + xi * xi - 10.0 * cos_approx(2.0 * PI * xi)
  }
  sum + 10.0 * n
}

///|
pub fn rosenbrock_function(x : Array[Double]) -> Double {
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

///|
pub fn ackley_function(x : Array[Double]) -> Double {
  let n : Double = x.length().to_double()
  let mut sum_sq = 0.0
  let mut sum_cos = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum_sq = sum_sq + xi * xi
    sum_cos = sum_cos + cos_approx(2.0 * PI * xi)
  }
  -20.0 * exp_approx(-0.2 * (sum_sq / n).sqrt()) -
  exp_approx(sum_cos / n) +
  20.0 +
  2.718281828459045
}

///|
pub fn griewank_function(x : Array[Double]) -> Double {
  let mut sum_sq = 0.0
  let mut product = 1.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum_sq = sum_sq + xi * xi
    product = product * cos_approx(xi / (i + 1).to_double().sqrt())
  }
  sum_sq / 4000.0 - product + 1.0
}

///|
pub fn schwefel_function(x : Array[Double]) -> Double {
  let n : Double = x.length().to_double()
  let mut sum = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum = sum + xi * sin_approx(abs_val(xi).sqrt())
  }
  418.9829 * n - sum
}

///|
pub fn levy_function(x : Array[Double]) -> Double {
  let n = x.length()
  let w : Array[Double] = []
  for i in 0..<n {
    w.push(1.0 + (x[i] - 1.0) / 4.0)
  }
  let head = sin_approx(PI * w[0]) * sin_approx(PI * w[0])
  let mut sum = head
  for i in 0..<(n - 1) {
    let wi = w[i]
    let s = sin_approx(PI * wi + 1.0)
    sum = sum + (wi - 1.0) * (wi - 1.0) * (1.0 + 10.0 * s * s)
  }
  let wn = w[n - 1]
  let tail_s = sin_approx(2.0 * PI * wn)
  sum + (wn - 1.0) * (wn - 1.0) * (1.0 + tail_s * tail_s)
}

///|
pub fn michalewicz_function(x : Array[Double]) -> Double {
  let m = 10
  let mut sum = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    let base = sin_approx((i + 1).to_double() * xi * xi / PI)
    let mut power = 1.0
    for _k in 0..<m {
      power = power * base * base
    }
    sum = sum + sin_approx(xi) * power
  }
  -sum
}

///|
pub fn zakharov_function(x : Array[Double]) -> Double {
  let mut sum_sq = 0.0
  let mut sum_ix = 0.0
  for i in 0..<x.length() {
    let xi = x[i]
    sum_sq = sum_sq + xi * xi
    sum_ix = sum_ix + 0.5 * (i + 1).to_double() * xi
  }
  sum_sq + sum_ix * sum_ix + sum_ix * sum_ix * sum_ix * sum_ix
}

///|
pub fn booth_function(x : Array[Double]) -> Double {
  let x1 = x[0]
  let x2 = x[1]
  (x1 + 2.0 * x2 - 7.0) * (x1 + 2.0 * x2 - 7.0) +
  (2.0 * x1 + x2 - 5.0) * (2.0 * x1 + x2 - 5.0)
}

///|
pub fn matyas_function(x : Array[Double]) -> Double {
  let x1 = x[0]
  let x2 = x[1]
  0.26 * (x1 * x1 + x2 * x2) - 0.48 * x1 * x2
}

// ===== 基准测试框架 =====

///|
pub fn calc_relative_error(
  best_value : Double,
  known_optimal : Double,
) -> Double {
  let diff = abs_val(best_value - known_optimal)
  if known_optimal == 0.0 {
    if best_value == 0.0 {
      0.0
    } else {
      diff
    }
  } else {
    diff / abs_val(known_optimal)
  }
}

///|
pub fn array_sum(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<arr.length() {
    sum = sum + arr[i]
  }
  sum
}

///|
pub fn array_min(arr : Array[Double]) -> Double {
  let mut min_val = arr[0]
  for i in 1..<arr.length() {
    if arr[i] < min_val {
      min_val = arr[i]
    }
  }
  min_val
}

///|
pub fn array_max(arr : Array[Double]) -> Double {
  let mut max_val = arr[0]
  for i in 1..<arr.length() {
    if arr[i] > max_val {
      max_val = arr[i]
    }
  }
  max_val
}

///|
pub fn array_avg(arr : Array[Double]) -> Double {
  array_sum(arr) / arr.length().to_double()
}

///|
pub fn fmt_fixed(x : Double, decimals : Int32) -> String {
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

///|
fn make_factor(acc : Double, n : Int32) -> Double {
  if n == 0 {
    acc
  } else {
    make_factor(acc * 10.0, n - 1)
  }
}

///|
fn pad_zeros(s : String, count : Int32) -> String {
  if count == 0 {
    s
  } else {
    pad_zeros("0" + s, count - 1)
  }
}

///|

///|
/// 算法注册表：单一默认超参来源，新增算法只需在此追加一个条目。
pub struct AlgoEntry {
  name : String
  run : (Int32, Int32, Bounds, RNG, (Array[Double]) -> Double) -> OptimizationResult
}

///|
pub fn default_algorithms() -> Array[AlgoEntry] {
  let entries : Array[AlgoEntry] = []
  entries.push(AlgoEntry::{
    name: "Random Search",
    run: (e, d, b, r, f) => {
      random_search(make_random_search_config(e, d, b), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Hill Climbing",
    run: (e, d, b, r, f) => {
      hill_climbing(make_hill_climbing_config(e, d, b, 0.5, 10), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Simulated Annealing",
    run: (e, d, b, r, f) => {
      simulated_annealing(
        make_simulated_annealing_config(e, d, b, 100.0, 0.995, 0.5),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "PSO",
    run: (e, d, b, r, f) => {
      pso(make_pso_config(e, d, b, 30, 0.7, 1.5, 1.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Differential Evolution",
    run: (e, d, b, r, f) => {
      differential_evolution(make_de_config(e, d, b, 30, 0.8, 0.9), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Genetic Algorithm",
    run: (e, d, b, r, f) => {
      genetic_algorithm(make_ga_config(e, d, b, 50, 0.01, 0.8), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ES",
    run: (e, d, b, r, f) => {
      evolution_strategy(make_es_config(e, d, b, 30, 100, 0.3), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PBIL",
    run: (e, d, b, r, f) => {
      population_based_incremental_learning(
        make_pbil_config(e, d, b, 100, 10, 0.3),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "Nelder-Mead",
    run: (e, d, b, r, f) => {
      nelder_mead(make_nm_config(e, d, b, 1.0, 2.0, 0.5, 0.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Firefly",
    run: (e, d, b, r, f) => {
      firefly_algorithm(make_firefly_config(e, d, b, 25, 0.5, 1.0, 1.0), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ACO",
    run: (e, d, b, r, f) => {
      aco(make_aco_config(e, d, b, 20, 0.1, 1.0, 10.0), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Bee Colony",
    run: (e, d, b, r, f) => {
      bee_colony(make_bee_colony_config(e, d, b, 25, 100), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Tabu Search",
    run: (e, d, b, r, f) => {
      tabu_search(make_tabu_search_config(e, d, b, 50, 30, 0.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Local Search",
    run: (e, d, b, r, f) => {
      local_search(make_local_search_config(e, d, b, 30, 0.5, 0.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CMA-ES",
    run: (e, d, b, r, f) => {
      cma_es(make_cmaes_config(e, d, b, 0.5, 20, 10), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CEM",
    run: (e, d, b, r, f) => {
      cross_entropy_method(make_cem_config(e, d, b, 100, 0.1, 0.7), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Cuckoo Search",
    run: (e, d, b, r, f) => {
      cuckoo_search(make_cuckoo_config(e, d, b, 25, 0.25, 1.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GWO",
    run: (e, d, b, r, f) => {
      grey_wolf_optimizer(make_gwo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Hooke-Jeeves",
    run: (e, d, b, r, f) => {
      hooke_jeeves(make_hooke_jeeves_config(e, d, b, 0.5, 0.5, 1.0e-6), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WOA",
    run: (e, d, b, r, f) => {
      whale_optimization(make_woa_config(e, d, b, 25, 1.0), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Bat Algorithm",
    run: (e, d, b, r, f) => {
      bat_algorithm(make_bat_config(e, d, b, 25, 0.0, 1.0, 0.5, 1.0), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Harmony Search",
    run: (e, d, b, r, f) => {
      harmony_search(make_hs_config(e, d, b, 20, 0.9, 0.3, 0.05), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SCA",
    run: (e, d, b, r, f) => {
      sine_cosine_algorithm(make_sca_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MFO",
    run: (e, d, b, r, f) => {
      moth_flame_optimization(make_mfo_config(e, d, b, 25, 1.0), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ALO",
    run: (e, d, b, r, f) => {
      antlion_optimizer(make_alo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SSA",
    run: (e, d, b, r, f) => {
      sparrow_search_algorithm(make_ssa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HHO",
    run: (e, d, b, r, f) => {
      harris_hawks_optimization(make_hho_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FPA",
    run: (e, d, b, r, f) => {
      flower_pollination_algorithm(make_fpa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DA",
    run: (e, d, b, r, f) => {
      dragonfly_algorithm(make_da_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CSA",
    run: (e, d, b, r, f) => {
      crow_search_algorithm(make_csa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SalpSA",
    run: (e, d, b, r, f) => {
      salp_swarm_algorithm(make_salp_sa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MBO",
    run: (e, d, b, r, f) => {
      monarch_butterfly_optimization(make_mbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AEO",
    run: (e, d, b, r, f) => {
      artificial_ecosystem_optimization(make_aeo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "IWO",
    run: (e, d, b, r, f) => {
      invasive_weed_optimization(
        make_iwo_config(e, d, b, 10, 40, 1, 5, 2.0, 0.1, 3.0),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "BFO",
    run: (e, d, b, r, f) => {
      bacteria_foraging_optimization(
        make_bfo_config(e, d, b, 25, 5, 3, 0.2),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "SFLA",
    run: (e, d, b, r, f) => {
      shuffled_frog_leaping_algorithm(
        make_sfla_config(e, d, b, 30, 5, 10),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "GSO",
    run: (e, d, b, r, f) => {
      glowworm_swarm_optimization(
        make_gso_config(e, d, b, 25, 0.4, 0.6, 5, 5.0, 1.0),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "GSA",
    run: (e, d, b, r, f) => {
      gravitational_search_algorithm(
        make_gsa_config(e, d, b, 25, 100.0, 0.99),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "TLBO",
    run: (e, d, b, r, f) => {
      teaching_learning_based_optimization(make_tlbo_config(e, d, b, 30), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AMO",
    run: (e, d, b, r, f) => {
      animal_migration_optimization(make_amo_config(e, d, b, 25, 5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MVO",
    run: (e, d, b, r, f) => {
      multi_verse_optimizer(make_mvo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "KH",
    run: (e, d, b, r, f) => krill_herd(make_kh_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "ICA",
    run: (e, d, b, r, f) => {
      imperialist_competitive_algorithm(make_ica_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SFS",
    run: (e, d, b, r, f) => {
      stochastic_fractal_search(make_sfs_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SOS",
    run: (e, d, b, r, f) => {
      symbiotic_organisms_search(make_sos_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BBO",
    run: (e, d, b, r, f) => {
      biogeography_based_optimization(make_bbo_config(e, d, b, 25, 0.01), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EHO",
    run: (e, d, b, r, f) => {
      elephant_herding_optimization(make_eho_config(e, d, b, 25, 5, 0.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EO",
    run: (e, d, b, r, f) => {
      equilibrium_optimizer(make_eo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WDO",
    run: (e, d, b, r, f) => {
      wind_driven_optimization(make_wdo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Jaya",
    run: (e, d, b, r, f) => {
      jaya_optimization(make_jaya_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SLC",
    run: (e, d, b, r, f) => {
      soccer_league_competition(make_slc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GBPSO",
    run: (e, d, b, r, f) => {
      gaussian_bare_bones_pso(make_gbpso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MSA",
    run: (e, d, b, r, f) => {
      moth_search_algorithm(make_msa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PFA",
    run: (e, d, b, r, f) => {
      pathfinder_algorithm(make_pfa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CPA",
    run: (e, d, b, r, f) => {
      colony_predation_algorithm(make_cpa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BCO",
    run: (e, d, b, r, f) => {
      bacteria_colony_optimization(make_bco_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FFOA",
    run: (e, d, b, r, f) => {
      flying_fox_optimization(make_ffoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AFSA",
    run: (e, d, b, r, f) => {
      artificial_fish_swarm(make_afsa_config(e, d, b, 25, 2.0, 0.5, 0.5), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FWA",
    run: (e, d, b, r, f) => {
      fireworks_algorithm(make_fwa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BSO",
    run: (e, d, b, r, f) => {
      brain_storm_optimization(make_bso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BSA",
    run: (e, d, b, r, f) => {
      backtracking_search_algorithm(make_bsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PIO",
    run: (e, d, b, r, f) => {
      pigeon_inspired_optimization(make_pio_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AOA",
    run: (e, d, b, r, f) => {
      arithmetic_optimization_algorithm(make_aoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HBA",
    run: (e, d, b, r, f) => {
      honey_badger_algorithm(make_hba_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "POA",
    run: (e, d, b, r, f) => {
      pelican_optimization_algorithm(make_poa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TSA",
    run: (e, d, b, r, f) => {
      tuna_swarm_optimization(make_tsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GJO",
    run: (e, d, b, r, f) => {
      golden_jackal_optimization(make_gjo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GOA",
    run: (e, d, b, r, f) => {
      grasshopper_optimization_algorithm(make_goa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SOA",
    run: (e, d, b, r, f) => {
      seagull_optimization_algorithm(make_soa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CSO",
    run: (e, d, b, r, f) => {
      cat_swarm_optimization(make_cso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GTO",
    run: (e, d, b, r, f) => {
      gorilla_troops_optimizer(make_gto_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WCA",
    run: (e, d, b, r, f) => {
      water_cycle_algorithm(make_wca_config(e, d, b, 25, 4), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SMA",
    run: (e, d, b, r, f) => {
      slime_mould_algorithm(make_sma_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AO",
    run: (e, d, b, r, f) => aquila_optimizer(make_ao_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "ChOA",
    run: (e, d, b, r, f) => {
      chimp_optimization_algorithm(make_choa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DO",
    run: (e, d, b, r, f) => {
      dandelion_optimizer(make_do_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ZOA",
    run: (e, d, b, r, f) => {
      zebra_optimization_algorithm(make_zoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MRFO",
    run: (e, d, b, r, f) => {
      manta_ray_foraging_optimization(make_mrfo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HGS",
    run: (e, d, b, r, f) => {
      hunger_games_search(make_hgs_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RSA",
    run: (e, d, b, r, f) => {
      reptile_search_algorithm(make_rsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MPA",
    run: (e, d, b, r, f) => {
      marine_predators_algorithm(make_mpa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GEO",
    run: (e, d, b, r, f) => {
      golden_eagle_optimizer(make_geo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RUN",
    run: (e, d, b, r, f) => {
      runge_kutta_optimizer(make_run_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AVOA",
    run: (e, d, b, r, f) => {
      african_vultures_optimization(make_avoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WSO",
    run: (e, d, b, r, f) => {
      white_shark_optimizer(make_wso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "COA",
    run: (e, d, b, r, f) => {
      coati_optimization_algorithm(make_coa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "jDE",
    run: (e, d, b, r, f) => jde_algorithm(make_jde_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "SHADE",
    run: (e, d, b, r, f) => {
      shade_algorithm(make_shade_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EP",
    run: (e, d, b, r, f) => {
      evolutionary_programming(make_ep_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "OOA",
    run: (e, d, b, r, f) => {
      osprey_optimization_algorithm(make_ooa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TDO",
    run: (e, d, b, r, f) => {
      tasmanian_devil_optimization(make_tdo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PDO",
    run: (e, d, b, r, f) => {
      prairie_dog_optimization(make_pdo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NGO",
    run: (e, d, b, r, f) => {
      northern_goshawk_optimization(make_ngo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BWO",
    run: (e, d, b, r, f) => {
      beluga_whale_optimization(make_bwo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GBO",
    run: (e, d, b, r, f) => {
      gradient_based_optimizer(make_gbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GO",
    run: (e, d, b, r, f) => growth_optimizer(make_go_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "WaOA",
    run: (e, d, b, r, f) => {
      walrus_optimizer(make_waoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CLONALG",
    run: (e, d, b, r, f) => {
      clonal_selection_algorithm(make_clonalg_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ABC",
    run: (e, d, b, r, f) => {
      artificial_bee_colony(make_abc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EM",
    run: (e, d, b, r, f) => {
      electromagnetism_like_mechanism(make_em_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DBO",
    run: (e, d, b, r, f) => {
      dung_beetle_optimizer(make_dbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TSO",
    run: (e, d, b, r, f) => {
      tunicate_swarm_algorithm(make_tso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "KOA",
    run: (e, d, b, r, f) => {
      kookaburra_optimization_algorithm(make_koa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FHO",
    run: (e, d, b, r, f) => {
      fire_hawk_optimizer(make_fho_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RIME",
    run: (e, d, b, r, f) => {
      rime_optimization_algorithm(make_rime_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "INFO",
    run: (e, d, b, r, f) => {
      info_optimization_algorithm(make_info_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SO",
    run: (e, d, b, r, f) => snake_optimizer(make_so_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "MGO",
    run: (e, d, b, r, f) => {
      mountain_gazelle_optimizer(make_mgo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CHIO",
    run: (e, d, b, r, f) => {
      coronavirus_herd_immunity_optimization(
        make_chio_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "PO",
    run: (e, d, b, r, f) => {
      political_optimizer(make_po_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SCSO",
    run: (e, d, b, r, f) => {
      sand_cat_swarm_optimization(make_scso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SABO",
    run: (e, d, b, r, f) => {
      subtraction_average_based_optimizer(make_sabo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HPO",
    run: (e, d, b, r, f) => {
      hunter_prey_optimizer(make_hpo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AHA",
    run: (e, d, b, r, f) => {
      artificial_hummingbird_algorithm(make_aha_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BKA",
    run: (e, d, b, r, f) => {
      black_winged_kite_algorithm(make_bka_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ARO",
    run: (e, d, b, r, f) => {
      artificial_rabbits_optimization(make_aro_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FOX",
    run: (e, d, b, r, f) => fox_optimizer(make_fox_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "DMO",
    run: (e, d, b, r, f) => {
      dwarf_mongoose_optimization(make_dmo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RSO",
    run: (e, d, b, r, f) => {
      rat_swarm_optimizer(make_rso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EVO",
    run: (e, d, b, r, f) => {
      egyptian_vulture_optimization(make_evo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EEFO",
    run: (e, d, b, r, f) => {
      electric_eel_foraging_optimization(make_eefo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Golden-SA",
    run: (e, d, b, r, f) => {
      golden_sine_algorithm(make_golden_sa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "YSGA",
    run: (e, d, b, r, f) => {
      yellow_saddle_goatfish_algorithm(make_ysga_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FBI",
    run: (e, d, b, r, f) => {
      forensic_based_investigation(make_fbi_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CapSA",
    run: (e, d, b, r, f) => {
      capuchin_search_algorithm(make_capsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SWO",
    run: (e, d, b, r, f) => {
      spider_wasp_optimizer(make_swo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GAO",
    run: (e, d, b, r, f) => {
      green_anaconda_optimization(make_gao_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MAO",
    run: (e, d, b, r, f) => {
      mexican_axolotl_optimization(make_mao_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "UMDA",
    run: (e, d, b, r, f) => {
      univariate_marginal_distribution_algorithm(
        make_umda_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "cGA",
    run: (e, d, b, r, f) => {
      compact_genetic_algorithm(make_cga_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CA",
    run: (e, d, b, r, f) => {
      cultural_algorithm(make_ca_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "COY",
    run: (e, d, b, r, f) => {
      coyote_optimization_algorithm(make_coy_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BES",
    run: (e, d, b, r, f) => {
      bald_eagle_search(make_bes_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GOOSE",
    run: (e, d, b, r, f) => {
      goose_optimization_algorithm(make_goose_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Powell",
    run: (e, d, b, r, f) => {
      powell_conjugate_direction_method(make_powell_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TEO",
    run: (e, d, b, r, f) => {
      thermal_exchange_optimization(make_teo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SHO",
    run: (e, d, b, r, f) => {
      spotted_hyena_optimizer(make_sho_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ROA",
    run: (e, d, b, r, f) => {
      remora_optimization_algorithm(make_roa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HOA",
    run: (e, d, b, r, f) => {
      horse_herd_optimization(make_hoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BWA",
    run: (e, d, b, r, f) => {
      black_widow_optimization_algorithm(make_bwa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CPSO",
    run: (e, d, b, r, f) => {
      constriction_pso(make_cpso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LPSO",
    run: (e, d, b, r, f) => local_best_pso(make_lpso_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "WWO",
    run: (e, d, b, r, f) => {
      water_wave_optimization(make_wwo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WHO",
    run: (e, d, b, r, f) => {
      wild_horse_optimizer(make_who_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GD",
    run: (e, d, b, r, f) => {
      great_deluge_algorithm(make_gd_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "VNS",
    run: (e, d, b, r, f) => {
      variable_neighborhood_search(make_vns_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Rosenbrock",
    run: (e, d, b, r, f) => {
      rosenbrock_method(make_rosenbrock_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BBBC",
    run: (e, d, b, r, f) => {
      big_bang_big_crunch(make_bbbc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BOA",
    run: (e, d, b, r, f) => {
      butterfly_optimization(make_boa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MFA",
    run: (e, d, b, r, f) => mayfly_algorithm(make_mfa_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "SBO",
    run: (e, d, b, r, f) => {
      secretary_bird_optimization(make_sbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TA",
    run: (e, d, b, r, f) => {
      threshold_accepting(make_ta_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ACOR",
    run: (e, d, b, r, f) => {
      ant_colony_continuous(make_acor_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "QPSO",
    run: (e, d, b, r, f) => {
      quantum_behaved_pso(make_qpso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SaDE",
    run: (e, d, b, r, f) => {
      self_adaptive_de(make_sade_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BH",
    run: (e, d, b, r, f) => {
      black_hole_algorithm(make_bh_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LOA",
    run: (e, d, b, r, f) => {
      lion_optimization(make_loa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BHop",
    run: (e, d, b, r, f) => basin_hopping(make_bhop_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "CoDE",
    run: (e, d, b, r, f) => composite_de(make_code_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "ODE",
    run: (e, d, b, r, f) => opposition_de(make_ode_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "CBO",
    run: (e, d, b, r, f) => {
      colliding_bodies_optimization(make_cbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TWO",
    run: (e, d, b, r, f) => {
      teamwork_optimization(make_two_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SPO",
    run: (e, d, b, r, f) => {
      spiral_optimization(make_spo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RTH",
    run: (e, d, b, r, f) => red_tailed_hawk(make_rth_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "CRS",
    run: (e, d, b, r, f) => {
      controlled_random_search(make_crs_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "UPSO",
    run: (e, d, b, r, f) => unified_pso(make_upso_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "RGA",
    run: (e, d, b, r, f) => real_coded_ga(make_rga_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "InS",
    run: (e, d, b, r, f) => interior_search(make_ins_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "BMO",
    run: (e, d, b, r, f) => {
      barnacles_mating_optimizer(make_bmo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NRBO",
    run: (e, d, b, r, f) => {
      newton_raphson_based_optimization(make_nrbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SPSA",
    run: (e, d, b, r, f) => {
      simultaneous_perturbation_spa(make_spsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DIRECT",
    run: (e, d, b, r, f) => {
      dividing_rectangles(make_direct_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Memetic",
    run: (e, d, b, r, f) => {
      memetic_algorithm(make_memetic_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "JADE",
    run: (e, d, b, r, f) => jade(make_jade_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "LCA",
    run: (e, d, b, r, f) => {
      league_championship_algorithm(make_lca_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HGSO",
    run: (e, d, b, r, f) => {
      henry_gas_solubility_optimization(make_hgso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "COOT",
    run: (e, d, b, r, f) => {
      coot_optimization_algorithm(make_coot_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DOA",
    run: (e, d, b, r, f) => {
      dingo_optimization_algorithm(make_doa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "STOA",
    run: (e, d, b, r, f) => {
      sooty_tern_optimization_algorithm(make_stoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GKSO",
    run: (e, d, b, r, f) => {
      genghis_khan_shark_optimizer(make_gkso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FLA",
    run: (e, d, b, r, f) => {
      ficks_law_algorithm(make_fla_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HO",
    run: (e, d, b, r, f) => {
      hippopotamus_optimization_algorithm(make_ho_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NMRA",
    run: (e, d, b, r, f) => {
      naked_mole_rat_algorithm(make_nmra_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FFO",
    run: (e, d, b, r, f) => {
      fennec_fox_optimization_algorithm(make_ffo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EPO",
    run: (e, d, b, r, f) => {
      emperor_penguin_optimization_algorithm(make_epo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CPO",
    run: (e, d, b, r, f) => {
      crested_porcupine_optimizer(make_cpo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MTBO",
    run: (e, d, b, r, f) => {
      mountaineering_team_based_optimization(
        make_mtbo_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "CryStAl",
    run: (e, d, b, r, f) => {
      crystal_structure_algorithm(
        make_crystal_structure_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "NOA",
    run: (e, d, b, r, f) => {
      nutcracker_optimization_algorithm(make_noa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ESOA",
    run: (e, d, b, r, f) => {
      egret_swarm_optimization_algorithm(make_esoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MOA",
    run: (e, d, b, r, f) => {
      meerkat_optimization_algorithm(make_moa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AOS",
    run: (e, d, b, r, f) => {
      atomic_orbital_search(make_aos_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MA",
    run: (e, d, b, r, f) => {
      mayfly_optimization_algorithm(make_ma_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SSO",
    run: (e, d, b, r, f) => {
      social_spider_optimization(make_sso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SLnO",
    run: (e, d, b, r, f) => {
      sea_lion_optimization_algorithm(make_sea_lion_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FOA",
    run: (e, d, b, r, f) => {
      fruit_fly_optimization_algorithm(make_foa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SFO",
    run: (e, d, b, r, f) => {
      sunflower_optimization_algorithm(make_sfo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NNA",
    run: (e, d, b, r, f) => {
      neural_network_algorithm(make_nna_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PBO",
    run: (e, d, b, r, f) => {
      polar_bear_optimization(make_pbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FSS",
    run: (e, d, b, r, f) => {
      fish_school_search(make_fss_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CFO",
    run: (e, d, b, r, f) => {
      central_force_optimization(make_cfo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CGO",
    run: (e, d, b, r, f) => {
      chaos_game_optimization(make_cgo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CDO",
    run: (e, d, b, r, f) => {
      chernobyl_disaster_optimizer(make_cdo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ILS",
    run: (e, d, b, r, f) => {
      iterated_local_search(make_ils_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SCE",
    run: (e, d, b, r, f) => {
      shuffled_complex_evolution(make_sce_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GRASP",
    run: (e, d, b, r, f) => {
      greedy_randomized_adaptive_search(make_grasp_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BO",
    run: (e, d, b, r, f) => bonobo_optimizer(make_bo_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "FSA",
    run: (e, d, b, r, f) => {
      flamingo_search_algorithm(make_fsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SMO",
    run: (e, d, b, r, f) => {
      spider_monkey_optimization(make_smo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HLO",
    run: (e, d, b, r, f) => {
      human_learning_optimization(make_hlo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FDA",
    run: (e, d, b, r, f) => {
      flow_direction_algorithm(make_fda_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CRO",
    run: (e, d, b, r, f) => {
      coral_reefs_optimization(make_cro_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AAA",
    run: (e, d, b, r, f) => {
      artificial_algae_algorithm(make_aaa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "VSA",
    run: (e, d, b, r, f) => {
      vortex_search_algorithm(make_vsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "VPL",
    run: (e, d, b, r, f) => {
      volleyball_premier_league(make_vpl_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WEO",
    run: (e, d, b, r, f) => {
      water_evaporation_optimization(make_weo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SNS",
    run: (e, d, b, r, f) => {
      social_network_search(make_sns_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TCO",
    run: (e, d, b, r, f) => {
      termite_colony_optimization(make_tco_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RO",
    run: (e, d, b, r, f) => ray_optimization(make_ro_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "GbSA",
    run: (e, d, b, r, f) => {
      galaxy_based_search_algorithm(make_gbsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BA",
    run: (e, d, b, r, f) => bees_algorithm(make_ba_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "CO",
    run: (e, d, b, r, f) => cheetah_optimizer(make_co_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "CSS",
    run: (e, d, b, r, f) => {
      charged_system_search(make_css_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MBA",
    run: (e, d, b, r, f) => {
      mine_blast_algorithm(make_mba_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PCA",
    run: (e, d, b, r, f) => {
      particle_collision_algorithm(make_pca_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "IMO",
    run: (e, d, b, r, f) => {
      ions_motion_optimization(make_imo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EMA",
    run: (e, d, b, r, f) => {
      exchange_market_algorithm(make_ema_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WCO",
    run: (e, d, b, r, f) => {
      world_cup_optimization(make_wco_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "AEFA",
    run: (e, d, b, r, f) => {
      artificial_electric_field_algorithm(make_aefa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "OIO",
    run: (e, d, b, r, f) => {
      optics_inspired_optimization(make_oio_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LAPO",
    run: (e, d, b, r, f) => {
      lightning_attachment_procedure_optimization(
        make_lapo_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "FGBO",
    run: (e, d, b, r, f) => {
      football_game_based_optimization(make_fgbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EA",
    run: (e, d, b, r, f) => {
      earthquake_algorithm(make_ea_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LAHC",
    run: (e, d, b, r, f) => {
      late_acceptance_hill_climbing(make_lahc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HuS",
    run: (e, d, b, r, f) => {
      hunting_search_algorithm(make_hus_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RRT",
    run: (e, d, b, r, f) => {
      record_to_record_travel(make_rrt_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NRO",
    run: (e, d, b, r, f) => {
      nuclear_reaction_optimization(make_nro_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HTS",
    run: (e, d, b, r, f) => {
      heat_transfer_search(make_hts_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PGO",
    run: (e, d, b, r, f) => {
      plasma_generation_optimization(make_pgo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "IWD",
    run: (e, d, b, r, f) => {
      intelligent_water_drops(make_iwd_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RFD",
    run: (e, d, b, r, f) => {
      river_formation_dynamics(make_rfd_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LSA",
    run: (e, d, b, r, f) => {
      lightning_search_algorithm(make_lsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SS",
    run: (e, d, b, r, f) => scatter_search(make_ss_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "OBA",
    run: (e, d, b, r, f) => {
      old_bachelor_acceptance(make_oba_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "APO",
    run: (e, d, b, r, f) => {
      artificial_physics_optimization(make_apo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HBO",
    run: (e, d, b, r, f) => {
      heap_based_optimizer(make_hbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SAO",
    run: (e, d, b, r, f) => {
      snow_ablation_optimizer(make_sao_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BRO",
    run: (e, d, b, r, f) => {
      battle_royale_optimization(make_bro_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EFO",
    run: (e, d, b, r, f) => {
      electromagnetic_field_optimization(make_efo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GBMO",
    run: (e, d, b, r, f) => {
      gases_brownian_motion_optimization(make_gbmo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DSA",
    run: (e, d, b, r, f) => {
      differential_search_algorithm(make_dsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GSK",
    run: (e, d, b, r, f) => {
      gaining_sharing_knowledge(make_gsk_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ALNS",
    run: (e, d, b, r, f) => {
      adaptive_large_neighborhood_search(make_alns_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SGO",
    run: (e, d, b, r, f) => {
      social_group_optimization(make_sgo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ASO",
    run: (e, d, b, r, f) => {
      atom_search_optimization(make_aso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BCMO",
    run: (e, d, b, r, f) => {
      balancing_composite_motion_optimization(
        make_bcmo_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "MKE",
    run: (e, d, b, r, f) => {
      monkey_king_evolution(make_mke_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RFO",
    run: (e, d, b, r, f) => {
      red_fox_optimization(make_rfo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SAR",
    run: (e, d, b, r, f) => {
      search_and_rescue_optimization(make_sar_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ChSA",
    run: (e, d, b, r, f) => {
      chameleon_swarm_algorithm(make_chsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RDA",
    run: (e, d, b, r, f) => {
      red_deer_algorithm(make_rda_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SqSA",
    run: (e, d, b, r, f) => {
      squirrel_search_algorithm(make_sqsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Puma",
    run: (e, d, b, r, f) => puma_optimizer(make_puma_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "StMO",
    run: (e, d, b, r, f) => {
      starling_murmuration_optimizer(make_stmo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SFOA",
    run: (e, d, b, r, f) => {
      swordfish_optimization(make_sfoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ABO",
    run: (e, d, b, r, f) => {
      african_buffalo_optimization(make_abo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Gannet",
    run: (e, d, b, r, f) => {
      gannet_optimization(make_gannet_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Narwhal",
    run: (e, d, b, r, f) => {
      narwhal_optimizer(make_narwhal_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Serval",
    run: (e, d, b, r, f) => {
      serval_optimization(make_serval_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ChSO",
    run: (e, d, b, r, f) => {
      chicken_swarm_optimization(make_chso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PR",
    run: (e, d, b, r, f) => path_relinking(make_pr_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "CrSO",
    run: (e, d, b, r, f) => {
      cockroach_swarm_optimization(make_crso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "JS",
    run: (e, d, b, r, f) => {
      jellyfish_search_optimizer(make_js_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WSA",
    run: (e, d, b, r, f) => {
      water_strider_algorithm(make_wsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SeHO",
    run: (e, d, b, r, f) => {
      sea_horse_optimizer(make_seho_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WBO",
    run: (e, d, b, r, f) => {
      wildebeest_herd_optimization(make_wbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "OWL",
    run: (e, d, b, r, f) => owl_optimization(make_owl_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "RRO",
    run: (e, d, b, r, f) => {
      raven_roosting_optimization(make_rro_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SaFO",
    run: (e, d, b, r, f) => {
      sailfish_optimizer(make_safo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GPC",
    run: (e, d, b, r, f) => {
      giza_pyramids_construction(make_gpc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SSD",
    run: (e, d, b, r, f) => {
      social_ski_driver(make_ssd_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FDO",
    run: (e, d, b, r, f) => {
      fitness_dependent_optimizer(make_fdo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CFA",
    run: (e, d, b, r, f) => {
      cuttlefish_optimization(make_cfa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "OPA",
    run: (e, d, b, r, f) => {
      orca_predation_algorithm(make_opa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "APuf",
    run: (e, d, b, r, f) => {
      arctic_puffin_optimization(make_apuf_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WPA",
    run: (e, d, b, r, f) => {
      wolf_pack_algorithm(make_wpa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Racc",
    run: (e, d, b, r, f) => {
      raccoon_optimization(make_racc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Monkey",
    run: (e, d, b, r, f) => {
      monkey_search_algorithm(make_monkey_search_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BSwA",
    run: (e, d, b, r, f) => {
      bird_swarm_algorithm(make_bswa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TGA",
    run: (e, d, b, r, f) => {
      tree_growth_algorithm(make_tga_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PPA",
    run: (e, d, b, r, f) => {
      plant_propagation_algorithm(make_ppa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "HBMO",
    run: (e, d, b, r, f) => {
      honey_bee_mating_optimization(make_hbmo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Crayfish",
    run: (e, d, b, r, f) => {
      crayfish_optimization(make_crayfish_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "QSA",
    run: (e, d, b, r, f) => {
      queuing_search_algorithm(make_qsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Starfish",
    run: (e, d, b, r, f) => {
      starfish_optimization(make_starfish_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TBO",
    run: (e, d, b, r, f) => {
      tiger_beetle_optimization(make_tbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Lyrebird",
    run: (e, d, b, r, f) => {
      lyrebird_optimization(make_lyrebird_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SEOA",
    run: (e, d, b, r, f) => {
      social_emotional_optimization(make_seoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GFO",
    run: (e, d, b, r, f) => {
      guineafowl_optimization(make_guineafowl_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DEO",
    run: (e, d, b, r, f) => {
      dolphin_optimization(make_dolphin_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BAS",
    run: (e, d, b, r, f) => {
      beetle_antennae_search(make_bas_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ManSA",
    run: (e, d, b, r, f) => {
      mantis_search_algorithm(make_mantis_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EWA",
    run: (e, d, b, r, f) => {
      earthworm_optimization(make_ewa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GrSA",
    run: (e, d, b, r, f) => {
      group_search_optimization(make_grsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FrOA",
    run: (e, d, b, r, f) => {
      forest_optimization(make_froa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ChRO",
    run: (e, d, b, r, f) => {
      chemical_reaction_optimization(make_chro_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SPBO",
    run: (e, d, b, r, f) => {
      student_psychological_optimization(make_spbo_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PTO",
    run: (e, d, b, r, f) => {
      parallel_tempering_optimization(make_pto_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CLPSO",
    run: (e, d, b, r, f) => {
      comprehensive_learning_pso(make_clpso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ExtO",
    run: (e, d, b, r, f) => {
      extremal_optimization(make_exto_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BFGS",
    run: (e, d, b, r, f) => bfgs_optimizer(make_bfgs_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "SOMA",
    run: (e, d, b, r, f) => {
      self_organizing_migration(make_soma_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SDS",
    run: (e, d, b, r, f) => {
      stochastic_diffusion_search(make_sds_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DAnl",
    run: (e, d, b, r, f) => {
      dual_annealing(make_dual_annealing_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SQA",
    run: (e, d, b, r, f) => {
      simulated_quantum_annealing(make_sqa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "COBYLA",
    run: (e, d, b, r, f) => {
      cobyla_optimizer(make_cobyla_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NS",
    run: (e, d, b, r, f) => {
      nested_sampling(make_nested_sampling_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MADS",
    run: (e, d, b, r, f) => {
      mesh_adaptive_direct_search(make_mads_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "NEWUOA",
    run: (e, d, b, r, f) => {
      newuoa_optimizer(make_newuoa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PA",
    run: (e, d, b, r, f) => {
      population_annealing(make_population_annealing_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "OAIN",
    run: (e, d, b, r, f) => opt_ainet(make_oain_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "xNES",
    run: (e, d, b, r, f) => {
      exponential_natural_evolution_strategy(
        make_xnes_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "Demon",
    run: (e, d, b, r, f) => {
      demon_algorithm(make_demon_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MCS",
    run: (e, d, b, r, f) => {
      multilevel_coordinate_search(make_mcs_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WL",
    run: (e, d, b, r, f) => {
      wang_landau_sampling(make_wl_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LSHADE",
    run: (e, d, b, r, f) => {
      lshade_algorithm(make_lshade_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "1+1ES",
    run: (e, d, b, r, f) => {
      one_plus_one_evolution_strategy(
        make_one_plus_one_es_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "PRAXIS",
    run: (e, d, b, r, f) => {
      praxis_optimizer(make_praxis_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BOBYQA",
    run: (e, d, b, r, f) => {
      bobyqa_optimizer(make_bobyqa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "ASA",
    run: (e, d, b, r, f) => {
      adaptive_simulated_annealing(make_asa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GSS",
    run: (e, d, b, r, f) => {
      golden_section_search(make_gss_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Subplex",
    run: (e, d, b, r, f) => {
      subplex_optimizer(make_subplex_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MLSL",
    run: (e, d, b, r, f) => {
      multi_level_single_linkage(make_mlsl_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LBB",
    run: (e, d, b, r, f) => {
      lipschitz_branch_and_bound(make_lipschitz_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MTS",
    run: (e, d, b, r, f) => {
      multiple_trajectory_search(make_mts_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "TPE",
    run: (e, d, b, r, f) => {
      tree_parzen_estimator(make_tpe_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "EGO",
    run: (e, d, b, r, f) => {
      efficient_global_optimization(make_ego_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GR",
    run: (e, d, b, r, f) => {
      gutmann_ramm_optimization(make_gr_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "sepCMA",
    run: (e, d, b, r, f) => {
      separable_cma_es(make_sepcma_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "IFilt",
    run: (e, d, b, r, f) => {
      implicit_filtering(make_ifilt_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CG",
    run: (e, d, b, r, f) => {
      conjugate_gradient_method(make_cg_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DHO",
    run: (e, d, b, r, f) => {
      deer_hunting_optimization(make_dho_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CHC",
    run: (e, d, b, r, f) => {
      cross_generation_elitist_algorithm(make_chc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MDS",
    run: (e, d, b, r, f) => {
      multi_directional_search(make_mds_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CC",
    run: (e, d, b, r, f) => {
      cooperative_coevolution(make_cc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "MGA",
    run: (e, d, b, r, f) => {
      micro_genetic_algorithm(make_mga_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "CellGA",
    run: (e, d, b, r, f) => {
      cellular_genetic_algorithm(make_cellga_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GOMEA",
    run: (e, d, b, r, f) => {
      gene_pool_optimal_mixing(make_gomea_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "SPSO",
    run: (e, d, b, r, f) => species_pso(make_spso_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "OSA",
    run: (e, d, b, r, f) => {
      ostrich_optimization(make_osa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BB",
    run: (e, d, b, r, f) => barzilai_borwein(make_bb_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "Anderson",
    run: (e, d, b, r, f) => {
      anderson_acceleration(make_anderson_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "GWSO",
    run: (e, d, b, r, f) => {
      great_white_shark_optimizer(make_gwso_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RSM",
    run: (e, d, b, r, f) => {
      response_surface_methodology(make_rsm_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PC",
    run: (e, d, b, r, f) => {
      probability_collectives(make_pc_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "WMA",
    run: (e, d, b, r, f) => {
      woodpecker_mating_algorithm(make_wma_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "BMDA",
    run: (e, d, b, r, f) => {
      bivariate_marginal_distribution_algorithm(
        make_bmda_config(e, d, b, 25),
        f,
        r,
      )
    },
  })
  entries.push(AlgoEntry::{
    name: "PGSA",
    run: (e, d, b, r, f) => {
      plant_growth_simulation(make_pgsa_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "Crowd",
    run: (e, d, b, r, f) => {
      crowding_genetic_algorithm(make_crowd_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "PBT",
    run: (e, d, b, r, f) => {
      population_based_training(make_pbt_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "FS",
    run: (e, d, b, r, f) => fitness_sharing(make_fs_config(e, d, b, 25), f, r),
  })
  entries.push(AlgoEntry::{
    name: "Bandit",
    run: (e, d, b, r, f) => {
      bandit_allocation_algorithm(make_bandit_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "LA",
    run: (e, d, b, r, f) => {
      learning_automata_algorithm(make_la_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "RBR",
    run: (e, d, b, r, f) => {
      rank_based_replacement_algorithm(make_rbr_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "OGA",
    run: (e, d, b, r, f) => {
      orthogonal_genetic_algorithm(make_oga_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "G3",
    run: (e, d, b, r, f) => {
      generalized_generation_gap_pcx(make_g3_config(e, d, b, 25), f, r)
    },
  })
  entries.push(AlgoEntry::{
    name: "DG2",
    run: (e, d, b, r, f) => {
      differential_grouping_decomposition(make_dg2_config(e, d, b, 25), f, r)
    },
  })
  entries
}

///|
fn pad_right(s : String, width : Int32) -> String {
  let mut out = s
  while out.length() < width {
    out = out + " "
  }
  out
}

///|
pub fn run_benchmark_at_evals(
  f : (Array[Double]) -> Double,
  _func_name : String,
  dim : Int32,
  bounds : Bounds,
  known_optimal : Double,
  max_evals : Int32,
  rng : RNG,
) -> Array[BenchmarkResult] {
  let results : Array[BenchmarkResult] = []
  let algos = default_algorithms()
  for i in 0..<algos.length() {
    let algo = algos[i]
    let result = (algo.run)(max_evals, dim, bounds, rng, f)
    let error = calc_relative_error(result.best_value, known_optimal)
    results.push(BenchmarkResult::{
      algorithm_name: algo.name,
      best_value: result.best_value,
      evaluations: result.evaluations,
      relative_error: error,
    })
  }
  results
}

///|
pub fn print_benchmark_table(
  results : Array[BenchmarkResult],
  func_name : String,
  max_evals : Int32,
) -> Unit {
  println("  \{func_name} (evals = \{max_evals})")
  println(
    "  +--------------------+-------------------+------------+-------------------+",
  )
  println(
    "  | Algorithm          | Best Value        | Evals      | Relative Error    |",
  )
  println(
    "  +--------------------+-------------------+------------+-------------------+",
  )
  for i in 0..<results.length() {
    let r = results[i]
    let err_str = fmt_fixed(r.relative_error, 6)
    println(
      "  | \{r.algorithm_name} | \{r.best_value} | \{r.evaluations} | \{err_str} |",
    )
  }
  println(
    "  +--------------------+-------------------+------------+-------------------+",
  )
  println("")
}

///|
pub fn run_multiple_trials(
  f : (Array[Double]) -> Double,
  func_name : String,
  dim : Int32,
  bounds : Bounds,
  known_optimal : Double,
  max_evals : Int32,
  num_trials : Int32,
  base_seed : Double,
) -> Unit {
  println("=== 多次运行统计 (num_trials = \{num_trials}) ===")
  println("  函数: \{func_name}")
  println("")

  let algos = default_algorithms()
  let all_values : Array[Array[Double]] = []
  for _i in 0..<algos.length() {
    all_values.push([])
  }

  for trial in 0..<num_trials {
    let rng = make_rng(base_seed + trial.to_double() * 1000.0)
    for i in 0..<algos.length() {
      let result = (algos[i].run)(max_evals, dim, bounds, rng, f)
      all_values[i].push(result.best_value)
    }
  }

  println(
    "  Algorithm          | Min Value           | Max Value           | Avg Value           | Avg Error",
  )
  println(
    "  +-------------------+---------------------+---------------------+---------------------+------------------+",
  )

  for i in 0..<algos.length() {
    let vals = all_values[i]
    let v_min = array_min(vals)
    let v_max = array_max(vals)
    let v_avg = array_avg(vals)
    let v_err = calc_relative_error(v_avg, known_optimal)
    println(
      "  | \{pad_right(algos[i].name, 18)} | \{v_min} | \{v_max} | \{v_avg} | \{v_err} |",
    )
  }
  println("")
}

///|
pub fn run_evals_comparison(
  f : (Array[Double]) -> Double,
  func_name : String,
  dim : Int32,
  bounds : Bounds,
  known_optimal : Double,
  rng : RNG,
) -> Unit {
  println("=== 评估次数对比 ===")
  println("  函数: \{func_name}")
  println("  维度: \{dim}")
  println("  已知最优值: \{known_optimal}")
  println("")

  let eval_budgets : Array[Int32] = [100, 500, 1000, 5000, 10000]

  for i in 0..<eval_budgets.length() {
    let max_evals = eval_budgets[i]
    let results = run_benchmark_at_evals(
      f, func_name, dim, bounds, known_optimal, max_evals, rng,
    )
    print_benchmark_table(results, func_name, max_evals)
  }
}
