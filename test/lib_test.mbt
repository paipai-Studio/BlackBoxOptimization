// BlackBox Optimization 测试

///|
test "cos_approx basic values" {
  assert_eq(@lib.cos_approx(0.0), 1.0)
  let cos_pi = @lib.cos_approx(@lib.PI)
  assert_true(cos_pi > -1.05)
  assert_true(cos_pi < -0.95)
  let cos_pi_half = @lib.cos_approx(@lib.PI / 2.0)
  assert_true(cos_pi_half > -0.05)
  assert_true(cos_pi_half < 0.05)
}

///|
test "cos_approx periodicity" {
  let result = @lib.cos_approx(3.0 * @lib.PI)
  assert_true(result > -1.05)
  assert_true(result < -0.95)
}

///|
test "exp_approx basic values" {
  assert_eq(@lib.exp_approx(0.0), 1.0)
  let result = @lib.exp_approx(1.0)
  assert_true(result > 2.5)
  assert_true(result < 3.0)
}

///|
test "exp_approx negative values" {
  let result = @lib.exp_approx(-1.0)
  assert_true(result > 0.0)
  assert_true(result < 1.0)
}

///|
test "max and min functions" {
  assert_eq(@lib.max(3.0, 5.0), 5.0)
  assert_eq(@lib.max(-1.0, -5.0), -1.0)
  assert_eq(@lib.min(3.0, 5.0), 3.0)
  assert_eq(@lib.min(-1.0, -5.0), -5.0)
}

///|
test "abs_val function" {
  assert_eq(@lib.abs_val(5.0), 5.0)
  assert_eq(@lib.abs_val(-5.0), 5.0)
  assert_eq(@lib.abs_val(0.0), 0.0)
}

///|
test "reduce_to_pi" {
  let result = @lib.reduce_to_pi(0.0)
  assert_eq(result, 0.0)
}

///|
test "RNG generates values in [0,1)" {
  let rng = @lib.make_rng(42.0)
  let val = @lib.RNG::next(rng)
  assert_true(val >= 0.0)
  assert_true(val < 1.0)
}

///|
test "RNG deterministic with same seed" {
  let rng1 = @lib.make_rng(123.0)
  let val1 = @lib.RNG::next(rng1)
  let val2 = @lib.RNG::next(rng1)

  let rng3 = @lib.make_rng(123.0)
  let val3 = @lib.RNG::next(rng3)
  let val4 = @lib.RNG::next(rng3)

  assert_eq(val1, val3)
  assert_eq(val2, val4)
}

///|
test "sphere function at origin" {
  let zeros : Array[Double] = [0.0, 0.0, 0.0]
  let result = @lib.sphere_function(zeros)
  assert_eq(result, 0.0)
}

///|
test "sphere function at non-origin" {
  let input : Array[Double] = [1.0, 2.0, 3.0]
  let result = @lib.sphere_function(input)
  assert_eq(result, 14.0)
}

///|
test "rastrigin function at origin" {
  let zeros : Array[Double] = [0.0, 0.0, 0.0]
  let result = @lib.rastrigin_function(zeros)
  assert_eq(result, 0.0)
}

///|
test "rastrigin function non-negative" {
  let input : Array[Double] = [1.0, 2.0, 3.0]
  let result = @lib.rastrigin_function(input)
  assert_true(result >= 0.0)
}

///|
test "rosenbrock function at origin" {
  let input : Array[Double] = [1.0, 1.0, 1.0]
  let result = @lib.rosenbrock_function(input)
  assert_eq(result, 0.0)
}

///|
test "rosenbrock function at non-optimal" {
  let input : Array[Double] = [0.0, 0.0, 0.0]
  let result = @lib.rosenbrock_function(input)
  assert_true(result > 0.0)
}

///|
test "array_sum function" {
  let arr : Array[Double] = [1.0, 2.0, 3.0, 4.0]
  let result = @lib.array_sum(arr)
  assert_eq(result, 10.0)
}

///|
test "array_min function" {
  let arr : Array[Double] = [3.0, 1.0, 4.0, 1.0, 5.0]
  let result = @lib.array_min(arr)
  assert_eq(result, 1.0)
}

///|
test "array_max function" {
  let arr : Array[Double] = [3.0, 1.0, 4.0, 1.0, 5.0]
  let result = @lib.array_max(arr)
  assert_eq(result, 5.0)
}

///|
test "array_avg function" {
  let arr : Array[Double] = [2.0, 4.0, 6.0]
  let result = @lib.array_avg(arr)
  assert_eq(result, 4.0)
}

///|
test "calc_relative_error with zero optimal" {
  let error1 = @lib.calc_relative_error(0.0, 0.0)
  assert_eq(error1, 0.0)

  let error2 = @lib.calc_relative_error(5.0, 0.0)
  assert_eq(error2, 5.0)
}

///|
test "calc_relative_error with non-zero optimal" {
  let error = @lib.calc_relative_error(110.0, 100.0)
  assert_eq(error, 0.1)
}

///|
test "random_solution generates within bounds" {
  let bounds = @lib.make_bounds(-5.0, 5.0)
  let rng = @lib.make_rng(42.0)
  let solution = @lib.random_solution(5, bounds, rng)

  for i in 0..<solution.length() {
    assert_true(solution[i] >= -5.0)
    assert_true(solution[i] <= 5.0)
  }
}

///|
test "random_search improves with more evaluations" {
  let config1 = @lib.make_random_search_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.random_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_random_search_config(
    10000,
    3,
    @lib.make_bounds(-5.0, 5.0),
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.random_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "hill_climbing produces valid result" {
  let config = @lib.make_hill_climbing_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    0.5,
    5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.hill_climbing(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "simulated_annealing produces valid result" {
  let config = @lib.make_simulated_annealing_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    100.0,
    0.995,
    0.5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.simulated_annealing(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "evaluate function works" {
  let f = @lib.sphere_function
  let solution : Array[Double] = [1.0, 2.0]
  let result = @lib.evaluate(f, solution)
  assert_eq(result, 5.0)
}

///|
test "CEM produces valid result" {
  let config = @lib.make_cem_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    0.1,
    0.7,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.cross_entropy_method(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "cuckoo_search produces valid result" {
  let config = @lib.make_cuckoo_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    0.25,
    1.5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.cuckoo_search(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GWO produces valid result" {
  let config = @lib.make_gwo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.grey_wolf_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Hooke-Jeeves produces valid result" {
  let config = @lib.make_hooke_jeeves_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    0.5,
    0.5,
    1.0e-6,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.hooke_jeeves(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CEM improves with more evaluations" {
  let config1 = @lib.make_cem_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    0.1,
    0.7,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cross_entropy_method(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_cem_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    50,
    0.1,
    0.7,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cross_entropy_method(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GWO improves with more evaluations" {
  let config1 = @lib.make_gwo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.grey_wolf_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_gwo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.grey_wolf_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WOA produces valid result" {
  let config = @lib.make_woa_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    1.0,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.whale_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Bat Algorithm produces valid result" {
  let config = @lib.make_bat_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    0.0,
    1.0,
    0.5,
    1.0,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.bat_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Harmony Search produces valid result" {
  let config = @lib.make_hs_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    0.9,
    0.3,
    0.05,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.harmony_search(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Harmony Search improves with more evaluations" {
  let config1 = @lib.make_hs_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    0.9,
    0.3,
    0.05,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.harmony_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_hs_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    0.9,
    0.3,
    0.05,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.harmony_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "pow_approx works correctly" {
  let result1 = @lib.pow_approx(2.0, 3.0)
  assert_true(result1 > 7.0 && result1 < 9.0)

  let result2 = @lib.pow_approx(4.0, 0.5)
  assert_true(result2 > 1.5 && result2 < 2.5)
}

///|
test "sin_approx basic values" {
  assert_true(@lib.abs_val(@lib.sin_approx(0.0)) < 0.001)
  let sin_half = @lib.sin_approx(@lib.PI / 2.0)
  assert_true(sin_half > 0.99 && sin_half < 1.01)
  let sin_pi = @lib.sin_approx(@lib.PI)
  assert_true(@lib.abs_val(sin_pi) < 0.01)
}

///|
test "SCA produces valid result" {
  let config = @lib.make_sca_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sine_cosine_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SCA improves with more evaluations" {
  let config1 = @lib.make_sca_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sine_cosine_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_sca_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sine_cosine_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MFO produces valid result" {
  let config = @lib.make_mfo_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    1.0,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.moth_flame_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MFO improves with more evaluations" {
  let config1 = @lib.make_mfo_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    1.0,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.moth_flame_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mfo_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    15,
    1.0,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.moth_flame_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ALO produces valid result" {
  let config = @lib.make_alo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.antlion_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ALO improves with more evaluations" {
  let config1 = @lib.make_alo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.antlion_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_alo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.antlion_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SSA produces valid result" {
  let config = @lib.make_ssa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sparrow_search_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SSA improves with more evaluations" {
  let config1 = @lib.make_ssa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sparrow_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ssa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sparrow_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HHO produces valid result" {
  let config = @lib.make_hho_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.harris_hawks_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HHO improves with more evaluations" {
  let config1 = @lib.make_hho_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.harris_hawks_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_hho_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.harris_hawks_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FPA produces valid result" {
  let config = @lib.make_fpa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.flower_pollination_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FPA improves with more evaluations" {
  let config1 = @lib.make_fpa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.flower_pollination_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_fpa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.flower_pollination_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DA produces valid result" {
  let config = @lib.make_da_config(500, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dragonfly_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DA improves with more evaluations" {
  let config1 = @lib.make_da_config(100, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dragonfly_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_da_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 15)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dragonfly_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "IWO produces valid result" {
  let config = @lib.make_iwo_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    10,
    40,
    1,
    5,
    2.0,
    0.1,
    3.0,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.invasive_weed_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "IWO improves with more evaluations" {
  let config1 = @lib.make_iwo_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    10,
    40,
    1,
    5,
    2.0,
    0.1,
    3.0,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.invasive_weed_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_iwo_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    10,
    40,
    1,
    5,
    2.0,
    0.1,
    3.0,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.invasive_weed_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BFO produces valid result" {
  let config = @lib.make_bfo_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    3,
    2,
    0.2,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.bacteria_foraging_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BFO improves with more evaluations" {
  let config1 = @lib.make_bfo_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    3,
    2,
    0.2,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bacteria_foraging_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bfo_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    3,
    2,
    0.2,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bacteria_foraging_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SFLA produces valid result" {
  let config = @lib.make_sfla_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    4,
    5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.shuffled_frog_leaping_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SFLA improves with more evaluations" {
  let config1 = @lib.make_sfla_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    4,
    5,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.shuffled_frog_leaping_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sfla_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    4,
    5,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.shuffled_frog_leaping_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GSO produces valid result" {
  let config = @lib.make_gso_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    0.4,
    0.6,
    5,
    5.0,
    1.0,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.glowworm_swarm_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GSO improves with more evaluations" {
  let config1 = @lib.make_gso_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    0.4,
    0.6,
    5,
    5.0,
    1.0,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.glowworm_swarm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gso_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    0.4,
    0.6,
    5,
    5.0,
    1.0,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.glowworm_swarm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GSA produces valid result" {
  let config = @lib.make_gsa_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    100.0,
    0.99,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.gravitational_search_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GSA improves with more evaluations" {
  let config1 = @lib.make_gsa_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    100.0,
    0.99,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gravitational_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gsa_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    100.0,
    0.99,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gravitational_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TLBO produces valid result" {
  let config = @lib.make_tlbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.teaching_learning_based_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TLBO improves with more evaluations" {
  let config1 = @lib.make_tlbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.teaching_learning_based_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_tlbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.teaching_learning_based_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AMO produces valid result" {
  let config = @lib.make_amo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 20, 5)
  let rng = @lib.make_rng(42.0)
  let result = @lib.animal_migration_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AMO improves with more evaluations" {
  let config1 = @lib.make_amo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 20, 5)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.animal_migration_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_amo_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    20,
    5,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.animal_migration_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MVO produces valid result" {
  let config = @lib.make_mvo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.multi_verse_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MVO improves with more evaluations" {
  let config1 = @lib.make_mvo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.multi_verse_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_mvo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.multi_verse_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "KH produces valid result" {
  let config = @lib.make_kh_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.krill_herd(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "KH improves with more evaluations" {
  let config1 = @lib.make_kh_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.krill_herd(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_kh_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.krill_herd(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ICA produces valid result" {
  let config = @lib.make_ica_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.imperialist_competitive_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ICA improves with more evaluations" {
  let config1 = @lib.make_ica_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.imperialist_competitive_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ica_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.imperialist_competitive_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SFS produces valid result" {
  let config = @lib.make_sfs_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.stochastic_fractal_search(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SFS improves with more evaluations" {
  let config1 = @lib.make_sfs_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.stochastic_fractal_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sfs_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.stochastic_fractal_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SOS produces valid result" {
  let config = @lib.make_sos_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.symbiotic_organisms_search(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SOS improves with more evaluations" {
  let config1 = @lib.make_sos_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.symbiotic_organisms_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sos_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.symbiotic_organisms_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BBO produces valid result" {
  let config = @lib.make_bbo_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    0.01,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.biogeography_based_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BBO improves with more evaluations" {
  let config1 = @lib.make_bbo_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    0.01,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.biogeography_based_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bbo_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    0.01,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.biogeography_based_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EHO produces valid result" {
  let config = @lib.make_eho_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    5,
    0.5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.elephant_herding_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EHO improves with more evaluations" {
  let config1 = @lib.make_eho_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    5,
    0.5,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.elephant_herding_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_eho_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    5,
    0.5,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.elephant_herding_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EO produces valid result" {
  let config = @lib.make_eo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.equilibrium_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EO improves with more evaluations" {
  let config1 = @lib.make_eo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.equilibrium_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_eo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.equilibrium_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WDO produces valid result" {
  let config = @lib.make_wdo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.wind_driven_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WDO improves with more evaluations" {
  let config1 = @lib.make_wdo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.wind_driven_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_wdo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.wind_driven_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MSA produces valid result" {
  let config = @lib.make_msa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.moth_search_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MSA improves with more evaluations" {
  let config1 = @lib.make_msa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.moth_search_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_msa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.moth_search_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PFA produces valid result" {
  let config = @lib.make_pfa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.pathfinder_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PFA improves with more evaluations" {
  let config1 = @lib.make_pfa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.pathfinder_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_pfa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.pathfinder_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CPA produces valid result" {
  let config = @lib.make_cpa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.colony_predation_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CPA improves with more evaluations" {
  let config1 = @lib.make_cpa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.colony_predation_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cpa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.colony_predation_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BCO produces valid result" {
  let config = @lib.make_bco_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bacteria_colony_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BCO improves with more evaluations" {
  let config1 = @lib.make_bco_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bacteria_colony_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bco_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bacteria_colony_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FFOA produces valid result" {
  let config = @lib.make_ffoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.flying_fox_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FFOA improves with more evaluations" {
  let config1 = @lib.make_ffoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.flying_fox_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ffoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.flying_fox_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AFSA produces valid result" {
  let config = @lib.make_afsa_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    2.0,
    0.5,
    0.5,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_fish_swarm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AFSA improves with more evaluations" {
  let config1 = @lib.make_afsa_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    2.0,
    0.5,
    0.5,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_fish_swarm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_afsa_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    2.0,
    0.5,
    0.5,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_fish_swarm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FWA produces valid result" {
  let config = @lib.make_fwa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fireworks_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FWA improves with more evaluations" {
  let config1 = @lib.make_fwa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fireworks_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_fwa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fireworks_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BSO produces valid result" {
  let config = @lib.make_bso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.brain_storm_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BSO improves with more evaluations" {
  let config1 = @lib.make_bso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.brain_storm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.brain_storm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BSA produces valid result" {
  let config = @lib.make_bsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.backtracking_search_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BSA improves with more evaluations" {
  let config1 = @lib.make_bsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.backtracking_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.backtracking_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PIO produces valid result" {
  let config = @lib.make_pio_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.pigeon_inspired_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PIO improves with more evaluations" {
  let config1 = @lib.make_pio_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.pigeon_inspired_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pio_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.pigeon_inspired_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AOA produces valid result" {
  let config = @lib.make_aoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.arithmetic_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AOA improves with more evaluations" {
  let config1 = @lib.make_aoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.arithmetic_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_aoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.arithmetic_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HBA produces valid result" {
  let config = @lib.make_hba_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.honey_badger_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HBA improves with more evaluations" {
  let config1 = @lib.make_hba_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.honey_badger_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_hba_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.honey_badger_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "POA produces valid result" {
  let config = @lib.make_poa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.pelican_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "POA improves with more evaluations" {
  let config1 = @lib.make_poa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.pelican_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_poa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.pelican_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TSA produces valid result" {
  let config = @lib.make_tsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.tuna_swarm_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TSA improves with more evaluations" {
  let config1 = @lib.make_tsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.tuna_swarm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_tsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.tuna_swarm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GJO produces valid result" {
  let config = @lib.make_gjo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.golden_jackal_optimization(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GJO improves with more evaluations" {
  let config1 = @lib.make_gjo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.golden_jackal_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gjo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.golden_jackal_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Jaya produces valid result" {
  let config = @lib.make_jaya_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.jaya_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Jaya improves with more evaluations" {
  let config1 = @lib.make_jaya_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.jaya_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_jaya_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.jaya_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SLC produces valid result" {
  let config = @lib.make_slc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.soccer_league_competition(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SLC improves with more evaluations" {
  let config1 = @lib.make_slc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.soccer_league_competition(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_slc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.soccer_league_competition(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GBPSO produces valid result" {
  let config = @lib.make_gbpso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gaussian_bare_bones_pso(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GBPSO improves with more evaluations" {
  let config1 = @lib.make_gbpso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gaussian_bare_bones_pso(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gbpso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gaussian_bare_bones_pso(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GOA produces valid result" {
  let config = @lib.make_goa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.grasshopper_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GOA improves with more evaluations" {
  let config1 = @lib.make_goa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.grasshopper_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_goa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.grasshopper_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SOA produces valid result" {
  let config = @lib.make_soa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.seagull_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SOA improves with more evaluations" {
  let config1 = @lib.make_soa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.seagull_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_soa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.seagull_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CSO produces valid result" {
  let config = @lib.make_cso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cat_swarm_optimization(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CSO improves with more evaluations" {
  let config1 = @lib.make_cso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cat_swarm_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_cso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cat_swarm_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GTO produces valid result" {
  let config = @lib.make_gto_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gorilla_troops_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GTO improves with more evaluations" {
  let config1 = @lib.make_gto_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gorilla_troops_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gto_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gorilla_troops_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WCA produces valid result" {
  let config = @lib.make_wca_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25, 4)
  let rng = @lib.make_rng(42.0)
  let result = @lib.water_cycle_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WCA improves with more evaluations" {
  let config1 = @lib.make_wca_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25, 4)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.water_cycle_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_wca_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
    4,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.water_cycle_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SMA produces valid result" {
  let config = @lib.make_sma_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.slime_mould_algorithm(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SMA improves with more evaluations" {
  let config1 = @lib.make_sma_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.slime_mould_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_sma_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.slime_mould_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AO produces valid result" {
  let config = @lib.make_ao_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.aquila_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AO improves with more evaluations" {
  let config1 = @lib.make_ao_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.aquila_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ao_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.aquila_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ChOA produces valid result" {
  let config = @lib.make_choa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.chimp_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ChOA improves with more evaluations" {
  let config1 = @lib.make_choa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.chimp_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_choa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.chimp_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DO produces valid result" {
  let config = @lib.make_do_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dandelion_optimizer(config, @lib.sphere_function, rng)

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DO improves with more evaluations" {
  let config1 = @lib.make_do_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dandelion_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_do_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dandelion_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ZOA produces valid result" {
  let config = @lib.make_zoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.zebra_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )

  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ZOA improves with more evaluations" {
  let config1 = @lib.make_zoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.zebra_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_zoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.zebra_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ES produces valid result" {
  let config = @lib.make_es_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    100,
    0.3,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.evolution_strategy(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ES improves with more evaluations" {
  let config1 = @lib.make_es_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    100,
    0.3,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.evolution_strategy(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_es_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    30,
    100,
    0.3,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.evolution_strategy(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PBIL produces valid result" {
  let config = @lib.make_pbil_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    100,
    10,
    0.3,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.population_based_incremental_learning(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PBIL improves with more evaluations" {
  let config1 = @lib.make_pbil_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    100,
    10,
    0.3,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.population_based_incremental_learning(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pbil_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    100,
    10,
    0.3,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.population_based_incremental_learning(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MRFO produces valid result" {
  let config = @lib.make_mrfo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.manta_ray_foraging_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MRFO improves with more evaluations" {
  let config1 = @lib.make_mrfo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.manta_ray_foraging_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mrfo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.manta_ray_foraging_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HGS produces valid result" {
  let config = @lib.make_hgs_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.hunger_games_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HGS improves with more evaluations" {
  let config1 = @lib.make_hgs_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.hunger_games_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_hgs_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.hunger_games_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RSA produces valid result" {
  let config = @lib.make_rsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.reptile_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RSA improves with more evaluations" {
  let config1 = @lib.make_rsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.reptile_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_rsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.reptile_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MPA produces valid result" {
  let config = @lib.make_mpa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.marine_predators_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MPA improves with more evaluations" {
  let config1 = @lib.make_mpa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.marine_predators_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mpa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.marine_predators_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GEO produces valid result" {
  let config = @lib.make_geo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.golden_eagle_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GEO improves with more evaluations" {
  let config1 = @lib.make_geo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.golden_eagle_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_geo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.golden_eagle_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RUN produces valid result" {
  let config = @lib.make_run_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.runge_kutta_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RUN improves with more evaluations" {
  let config1 = @lib.make_run_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.runge_kutta_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_run_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.runge_kutta_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AVOA produces valid result" {
  let config = @lib.make_avoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.african_vultures_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AVOA improves with more evaluations" {
  let config1 = @lib.make_avoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.african_vultures_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_avoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.african_vultures_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WSO produces valid result" {
  let config = @lib.make_wso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.white_shark_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WSO improves with more evaluations" {
  let config1 = @lib.make_wso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.white_shark_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_wso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.white_shark_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "COA produces valid result" {
  let config = @lib.make_coa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.coati_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "COA improves with more evaluations" {
  let config1 = @lib.make_coa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.coati_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_coa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.coati_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "jDE produces valid result" {
  let config = @lib.make_jde_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.jde_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "jDE improves with more evaluations" {
  let config1 = @lib.make_jde_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.jde_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_jde_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.jde_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SHADE produces valid result" {
  let config = @lib.make_shade_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.shade_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SHADE improves with more evaluations" {
  let config1 = @lib.make_shade_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.shade_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_shade_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.shade_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EP produces valid result" {
  let config = @lib.make_ep_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.evolutionary_programming(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EP improves with more evaluations" {
  let config1 = @lib.make_ep_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.evolutionary_programming(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ep_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.evolutionary_programming(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "OOA produces valid result" {
  let config = @lib.make_ooa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.osprey_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "OOA improves with more evaluations" {
  let config1 = @lib.make_ooa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.osprey_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ooa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.osprey_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TDO produces valid result" {
  let config = @lib.make_tdo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.tasmanian_devil_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TDO improves with more evaluations" {
  let config1 = @lib.make_tdo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.tasmanian_devil_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_tdo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.tasmanian_devil_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PDO produces valid result" {
  let config = @lib.make_pdo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.prairie_dog_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PDO improves with more evaluations" {
  let config1 = @lib.make_pdo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.prairie_dog_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pdo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.prairie_dog_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NGO produces valid result" {
  let config = @lib.make_ngo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.northern_goshawk_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NGO improves with more evaluations" {
  let config1 = @lib.make_ngo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.northern_goshawk_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ngo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.northern_goshawk_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BWO produces valid result" {
  let config = @lib.make_bwo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.beluga_whale_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BWO improves with more evaluations" {
  let config1 = @lib.make_bwo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.beluga_whale_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bwo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.beluga_whale_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GBO produces valid result" {
  let config = @lib.make_gbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gradient_based_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GBO improves with more evaluations" {
  let config1 = @lib.make_gbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gradient_based_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gradient_based_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GO produces valid result" {
  let config = @lib.make_go_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.growth_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GO improves with more evaluations" {
  let config1 = @lib.make_go_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.growth_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_go_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.growth_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WaOA produces valid result" {
  let config = @lib.make_waoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.walrus_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WaOA improves with more evaluations" {
  let config1 = @lib.make_waoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.walrus_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_waoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.walrus_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CLONALG produces valid result" {
  let config = @lib.make_clonalg_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.clonal_selection_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CLONALG improves with more evaluations" {
  let config1 = @lib.make_clonalg_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.clonal_selection_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_clonalg_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.clonal_selection_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ABC produces valid result" {
  let config = @lib.make_abc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_bee_colony(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ABC improves with more evaluations" {
  let config1 = @lib.make_abc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_bee_colony(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_abc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_bee_colony(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EM produces valid result" {
  let config = @lib.make_em_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.electromagnetism_like_mechanism(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EM improves with more evaluations" {
  let config1 = @lib.make_em_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.electromagnetism_like_mechanism(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_em_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.electromagnetism_like_mechanism(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DBO produces valid result" {
  let config = @lib.make_dbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dung_beetle_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DBO improves with more evaluations" {
  let config1 = @lib.make_dbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dung_beetle_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_dbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dung_beetle_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TSO produces valid result" {
  let config = @lib.make_tso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.tunicate_swarm_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TSO improves with more evaluations" {
  let config1 = @lib.make_tso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.tunicate_swarm_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_tso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.tunicate_swarm_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "KOA produces valid result" {
  let config = @lib.make_koa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.kookaburra_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "KOA improves with more evaluations" {
  let config1 = @lib.make_koa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.kookaburra_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_koa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.kookaburra_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FHO produces valid result" {
  let config = @lib.make_fho_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fire_hawk_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FHO improves with more evaluations" {
  let config1 = @lib.make_fho_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fire_hawk_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_fho_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fire_hawk_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RIME produces valid result" {
  let config = @lib.make_rime_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.rime_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RIME improves with more evaluations" {
  let config1 = @lib.make_rime_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.rime_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_rime_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.rime_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "INFO produces valid result" {
  let config = @lib.make_info_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.info_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "INFO improves with more evaluations" {
  let config1 = @lib.make_info_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.info_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_info_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.info_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SO produces valid result" {
  let config = @lib.make_so_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.snake_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SO improves with more evaluations" {
  let config1 = @lib.make_so_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.snake_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_so_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.snake_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MGO produces valid result" {
  let config = @lib.make_mgo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mountain_gazelle_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MGO improves with more evaluations" {
  let config1 = @lib.make_mgo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mountain_gazelle_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mgo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mountain_gazelle_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CHIO produces valid result" {
  let config = @lib.make_chio_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.coronavirus_herd_immunity_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CHIO improves with more evaluations" {
  let config1 = @lib.make_chio_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.coronavirus_herd_immunity_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_chio_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.coronavirus_herd_immunity_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PO produces valid result" {
  let config = @lib.make_po_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.political_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PO improves with more evaluations" {
  let config1 = @lib.make_po_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.political_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_po_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.political_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SCSO produces valid result" {
  let config = @lib.make_scso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sand_cat_swarm_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SCSO improves with more evaluations" {
  let config1 = @lib.make_scso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sand_cat_swarm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_scso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sand_cat_swarm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SABO produces valid result" {
  let config = @lib.make_sabo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.subtraction_average_based_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SABO improves with more evaluations" {
  let config1 = @lib.make_sabo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.subtraction_average_based_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sabo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.subtraction_average_based_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HPO produces valid result" {
  let config = @lib.make_hpo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.hunter_prey_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HPO improves with more evaluations" {
  let config1 = @lib.make_hpo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.hunter_prey_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_hpo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.hunter_prey_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AHA produces valid result" {
  let config = @lib.make_aha_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_hummingbird_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AHA improves with more evaluations" {
  let config1 = @lib.make_aha_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_hummingbird_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_aha_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_hummingbird_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BKA produces valid result" {
  let config = @lib.make_bka_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.black_winged_kite_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BKA improves with more evaluations" {
  let config1 = @lib.make_bka_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.black_winged_kite_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bka_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.black_winged_kite_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ARO produces valid result" {
  let config = @lib.make_aro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_rabbits_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ARO improves with more evaluations" {
  let config1 = @lib.make_aro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_rabbits_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_aro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_rabbits_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FOX produces valid result" {
  let config = @lib.make_fox_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fox_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FOX improves with more evaluations" {
  let config1 = @lib.make_fox_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fox_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_fox_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fox_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DMO produces valid result" {
  let config = @lib.make_dmo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dwarf_mongoose_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DMO improves with more evaluations" {
  let config1 = @lib.make_dmo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dwarf_mongoose_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_dmo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dwarf_mongoose_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RSO produces valid result" {
  let config = @lib.make_rso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.rat_swarm_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RSO improves with more evaluations" {
  let config1 = @lib.make_rso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.rat_swarm_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_rso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.rat_swarm_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EVO produces valid result" {
  let config = @lib.make_evo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.egyptian_vulture_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EVO improves with more evaluations" {
  let config1 = @lib.make_evo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.egyptian_vulture_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_evo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.egyptian_vulture_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EEFO produces valid result" {
  let config = @lib.make_eefo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.electric_eel_foraging_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EEFO improves with more evaluations" {
  let config1 = @lib.make_eefo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.electric_eel_foraging_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_eefo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.electric_eel_foraging_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Golden-SA produces valid result" {
  let config = @lib.make_golden_sa_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.golden_sine_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Golden-SA improves with more evaluations" {
  let config1 = @lib.make_golden_sa_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.golden_sine_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_golden_sa_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.golden_sine_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "YSGA produces valid result" {
  let config = @lib.make_ysga_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.yellow_saddle_goatfish_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "YSGA improves with more evaluations" {
  let config1 = @lib.make_ysga_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.yellow_saddle_goatfish_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ysga_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.yellow_saddle_goatfish_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FBI produces valid result" {
  let config = @lib.make_fbi_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.forensic_based_investigation(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FBI improves with more evaluations" {
  let config1 = @lib.make_fbi_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.forensic_based_investigation(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_fbi_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.forensic_based_investigation(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CapSA produces valid result" {
  let config = @lib.make_capsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.capuchin_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CapSA improves with more evaluations" {
  let config1 = @lib.make_capsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.capuchin_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_capsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.capuchin_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SWO produces valid result" {
  let config = @lib.make_swo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.spider_wasp_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SWO improves with more evaluations" {
  let config1 = @lib.make_swo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.spider_wasp_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_swo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.spider_wasp_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GAO produces valid result" {
  let config = @lib.make_gao_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.green_anaconda_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GAO improves with more evaluations" {
  let config1 = @lib.make_gao_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.green_anaconda_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gao_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.green_anaconda_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MAO produces valid result" {
  let config = @lib.make_mao_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mexican_axolotl_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MAO improves with more evaluations" {
  let config1 = @lib.make_mao_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mexican_axolotl_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mao_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mexican_axolotl_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "UMDA produces valid result" {
  let config = @lib.make_umda_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.univariate_marginal_distribution_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "UMDA improves with more evaluations" {
  let config1 = @lib.make_umda_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.univariate_marginal_distribution_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_umda_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.univariate_marginal_distribution_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "cGA produces valid result" {
  let config = @lib.make_cga_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.compact_genetic_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "cGA improves with more evaluations" {
  let config1 = @lib.make_cga_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.compact_genetic_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cga_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.compact_genetic_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CA produces valid result" {
  let config = @lib.make_ca_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cultural_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CA improves with more evaluations" {
  let config1 = @lib.make_ca_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cultural_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ca_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cultural_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}
