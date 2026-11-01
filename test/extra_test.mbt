// 早期算法与基础设施的补充测试

///|
test "PSO produces valid result" {
  let config = @lib.make_pso_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.7,
    1.5,
    1.5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.pso(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PSO improves with more evaluations" {
  let config1 = @lib.make_pso_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.7,
    1.5,
    1.5,
  )
  let result1 = @lib.pso(config1, @lib.sphere_function, @lib.make_rng(42.0))
  let config2 = @lib.make_pso_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.7,
    1.5,
    1.5,
  )
  let result2 = @lib.pso(config2, @lib.sphere_function, @lib.make_rng(42.0))
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DE produces valid result" {
  let config = @lib.make_de_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.8,
    0.9,
  )
  let result = @lib.differential_evolution(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DE improves with more evaluations" {
  let config1 = @lib.make_de_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.8,
    0.9,
  )
  let result1 = @lib.differential_evolution(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_de_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.8,
    0.9,
  )
  let result2 = @lib.differential_evolution(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GA produces valid result" {
  let config = @lib.make_ga_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    0.01,
    0.8,
  )
  let result = @lib.genetic_algorithm(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GA improves with more evaluations" {
  let config1 = @lib.make_ga_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    0.01,
    0.8,
  )
  let result1 = @lib.genetic_algorithm(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_ga_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    0.01,
    0.8,
  )
  let result2 = @lib.genetic_algorithm(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Nelder-Mead produces valid result" {
  let config = @lib.make_nm_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    1.0,
    2.0,
    0.5,
    0.5,
  )
  let result = @lib.nelder_mead(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Nelder-Mead improves with more evaluations" {
  let config1 = @lib.make_nm_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    1.0,
    2.0,
    0.5,
    0.5,
  )
  let result1 = @lib.nelder_mead(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_nm_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    1.0,
    2.0,
    0.5,
    0.5,
  )
  let result2 = @lib.nelder_mead(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Firefly produces valid result" {
  let config = @lib.make_firefly_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    0.5,
    1.0,
    1.0,
  )
  let result = @lib.firefly_algorithm(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Firefly improves with more evaluations" {
  let config1 = @lib.make_firefly_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    0.5,
    1.0,
    1.0,
  )
  let result1 = @lib.firefly_algorithm(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_firefly_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    0.5,
    1.0,
    1.0,
  )
  let result2 = @lib.firefly_algorithm(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ACO produces valid result" {
  let config = @lib.make_aco_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    0.1,
    1.0,
    10.0,
  )
  let result = @lib.aco(config, @lib.sphere_function, @lib.make_rng(42.0))
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ACO improves with more evaluations" {
  let config1 = @lib.make_aco_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    0.1,
    1.0,
    10.0,
  )
  let result1 = @lib.aco(config1, @lib.sphere_function, @lib.make_rng(42.0))
  let config2 = @lib.make_aco_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    0.1,
    1.0,
    10.0,
  )
  let result2 = @lib.aco(config2, @lib.sphere_function, @lib.make_rng(42.0))
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Bee Colony produces valid result" {
  let config = @lib.make_bee_colony_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    100,
  )
  let result = @lib.bee_colony(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Bee Colony improves with more evaluations" {
  let config1 = @lib.make_bee_colony_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    100,
  )
  let result1 = @lib.bee_colony(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_bee_colony_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    100,
  )
  let result2 = @lib.bee_colony(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Tabu Search produces valid result" {
  let config = @lib.make_tabu_search_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    30,
    0.5,
  )
  let result = @lib.tabu_search(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Tabu Search improves with more evaluations" {
  let config1 = @lib.make_tabu_search_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    30,
    0.5,
  )
  let result1 = @lib.tabu_search(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_tabu_search_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    30,
    0.5,
  )
  let result2 = @lib.tabu_search(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Local Search produces valid result" {
  let config = @lib.make_local_search_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.5,
    0.5,
  )
  let result = @lib.local_search(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Local Search improves with more evaluations" {
  let config1 = @lib.make_local_search_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.5,
    0.5,
  )
  let result1 = @lib.local_search(
    config1,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  let config2 = @lib.make_local_search_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.5,
    0.5,
  )
  let result2 = @lib.local_search(
    config2,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CMA-ES produces valid result" {
  let config = @lib.make_cmaes_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    0.5,
    20,
    10,
  )
  let result = @lib.cma_es(config, @lib.sphere_function, @lib.make_rng(42.0))
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CMA-ES improves with more evaluations" {
  let config1 = @lib.make_cmaes_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    0.5,
    20,
    10,
  )
  let result1 = @lib.cma_es(config1, @lib.sphere_function, @lib.make_rng(42.0))
  let config2 = @lib.make_cmaes_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    0.5,
    20,
    10,
  )
  let result2 = @lib.cma_es(config2, @lib.sphere_function, @lib.make_rng(42.0))
  assert_true(result2.best_value <= result1.best_value)
}

// ===== 数值精度断言（修复后的收敛性保障）=====

///|
test "DE converges on sphere" {
  let config = @lib.make_de_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.8,
    0.9,
  )
  let result = @lib.differential_evolution(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.best_value < 0.01)
}

///|
test "PSO converges on sphere" {
  let config = @lib.make_pso_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    0.7,
    1.5,
    1.5,
  )
  let result = @lib.pso(config, @lib.sphere_function, @lib.make_rng(42.0))
  assert_true(result.best_value < 0.01)
}

///|
test "DA converges on sphere" {
  let config = @lib.make_da_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let result = @lib.dragonfly_algorithm(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.best_value < 0.05)
}

///|
test "FFOA converges on sphere" {
  let config = @lib.make_ffoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let result = @lib.flying_fox_optimization(
    config,
    @lib.sphere_function,
    @lib.make_rng(42.0),
  )
  assert_true(result.best_value < 0.01)
}

///|
test "math approximations now use native precision" {
  let cos_pi = @lib.cos_approx(@lib.PI)
  assert_true(@lib.abs_val(cos_pi + 1.0) < 1.0e-12)
  let sin_half = @lib.sin_approx(@lib.PI / 2.0)
  assert_true(@lib.abs_val(sin_half - 1.0) < 1.0e-12)
  let e1 = @lib.exp_approx(1.0)
  assert_true(@lib.abs_val(e1 - 2.718281828459045) < 1.0e-12)
  let e10 = @lib.exp_approx(10.0)
  assert_true(@lib.abs_val(e10 - 22026.465794806718) < 1.0e-6)
  let l2 = @lib.pow_approx(2.0, 10.0)
  assert_true(@lib.abs_val(l2 - 1024.0) < 1.0e-9)
}

// ===== RNG 质量测试 =====

///|
test "RNG different seeds produce different sequences" {
  let rng1 = @lib.make_rng(1.0)
  let v1 = @lib.RNG::next(rng1)
  let v2 = @lib.RNG::next(rng1)
  let rng3 = @lib.make_rng(2.0)
  let v3 = @lib.RNG::next(rng3)
  let v4 = @lib.RNG::next(rng3)
  assert_true(v1 != v3 || v2 != v4)
}

///|
test "RNG zero seed does not deadlock" {
  let rng1 = @lib.make_rng(0.0)
  let v1 = @lib.RNG::next(rng1)
  let v2 = @lib.RNG::next(rng1)
  assert_true(v1 != v2)
}

///|
test "RNG long sequence stays uniform" {
  let rng = @lib.make_rng(7.0)
  let mut sum = 0.0
  let n = 10000
  for _i in 0..<n {
    let v = @lib.RNG::next(rng)
    assert_true(v >= 0.0)
    assert_true(v < 1.0)
    sum = sum + v
  }
  let avg = sum / n.to_double()
  assert_true(avg > 0.48)
  assert_true(avg < 0.52)
}

// ===== 注册表与基准框架测试 =====

///|
test "default_algorithms has complete registry" {
  let algos = @lib.default_algorithms()
  assert_true(algos.length() >= 75)
  for i in 0..<algos.length() {
    assert_true(algos[i].name.length() > 0)
  }
}

///|
test "registry runs every algorithm on sphere" {
  let algos = @lib.default_algorithms()
  let bounds = @lib.make_bounds(-5.0, 5.0)
  for i in 0..<algos.length() {
    let result = (algos[i].run)(
      1000,
      3,
      bounds,
      @lib.make_rng(42.0),
      @lib.sphere_function,
    )
    assert_true(result.best_value >= 0.0)
    assert_true(result.evaluations > 0)
    assert_true(result.evaluations <= 2000)
  }
}

///|
test "fmt_fixed produces fixed decimal strings" {
  let s = @lib.fmt_fixed(0.5, 3)
  assert_true(s.length() >= 5)
}

///|
test "run_benchmark_at_evals returns one result per algorithm" {
  let results = @lib.run_benchmark_at_evals(
    @lib.sphere_function,
    "Sphere",
    3,
    @lib.make_bounds(-5.0, 5.0),
    0.0,
    200,
    @lib.make_rng(42.0),
  )
  let algos = @lib.default_algorithms()
  assert_eq(results.length(), algos.length())
}
