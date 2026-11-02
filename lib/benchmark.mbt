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
