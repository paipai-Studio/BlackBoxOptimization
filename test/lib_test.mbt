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

///|
test "COY produces valid result" {
  let config = @lib.make_coy_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.coyote_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "COY improves with more evaluations" {
  let config1 = @lib.make_coy_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.coyote_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_coy_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.coyote_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BES produces valid result" {
  let config = @lib.make_bes_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bald_eagle_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BES improves with more evaluations" {
  let config1 = @lib.make_bes_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bald_eagle_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bes_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bald_eagle_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GOOSE produces valid result" {
  let config = @lib.make_goose_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.goose_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GOOSE improves with more evaluations" {
  let config1 = @lib.make_goose_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.goose_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_goose_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.goose_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Powell produces valid result" {
  let config = @lib.make_powell_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.powell_conjugate_direction_method(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Powell improves with more evaluations" {
  let config1 = @lib.make_powell_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.powell_conjugate_direction_method(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_powell_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.powell_conjugate_direction_method(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TEO produces valid result" {
  let config = @lib.make_teo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.thermal_exchange_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TEO improves with more evaluations" {
  let config1 = @lib.make_teo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.thermal_exchange_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_teo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.thermal_exchange_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SHO produces valid result" {
  let config = @lib.make_sho_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.spotted_hyena_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SHO improves with more evaluations" {
  let config1 = @lib.make_sho_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.spotted_hyena_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sho_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.spotted_hyena_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ROA produces valid result" {
  let config = @lib.make_roa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.remora_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ROA improves with more evaluations" {
  let config1 = @lib.make_roa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.remora_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_roa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.remora_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HOA produces valid result" {
  let config = @lib.make_hoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.horse_herd_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HOA improves with more evaluations" {
  let config1 = @lib.make_hoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.horse_herd_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_hoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.horse_herd_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BWA produces valid result" {
  let config = @lib.make_bwa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.black_widow_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BWA improves with more evaluations" {
  let config1 = @lib.make_bwa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.black_widow_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bwa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.black_widow_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CPSO produces valid result" {
  let config = @lib.make_cpso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.constriction_pso(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CPSO improves with more evaluations" {
  let config1 = @lib.make_cpso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.constriction_pso(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_cpso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.constriction_pso(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LPSO produces valid result" {
  let config = @lib.make_lpso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.local_best_pso(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LPSO improves with more evaluations" {
  let config1 = @lib.make_lpso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.local_best_pso(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_lpso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.local_best_pso(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WWO produces valid result" {
  let config = @lib.make_wwo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.water_wave_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WWO improves with more evaluations" {
  let config1 = @lib.make_wwo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.water_wave_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_wwo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.water_wave_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WHO produces valid result" {
  let config = @lib.make_who_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.wild_horse_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WHO improves with more evaluations" {
  let config1 = @lib.make_who_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.wild_horse_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_who_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.wild_horse_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GD produces valid result" {
  let config = @lib.make_gd_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.great_deluge_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GD improves with more evaluations" {
  let config1 = @lib.make_gd_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.great_deluge_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_gd_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.great_deluge_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "VNS produces valid result" {
  let config = @lib.make_vns_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.variable_neighborhood_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "VNS improves with more evaluations" {
  let config1 = @lib.make_vns_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.variable_neighborhood_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_vns_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.variable_neighborhood_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Rosenbrock method produces valid result" {
  let config = @lib.make_rosenbrock_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.rosenbrock_method(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Rosenbrock method improves with more evaluations" {
  let config1 = @lib.make_rosenbrock_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.rosenbrock_method(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_rosenbrock_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.rosenbrock_method(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BBBC produces valid result" {
  let config = @lib.make_bbbc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.big_bang_big_crunch(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BBBC improves with more evaluations" {
  let config1 = @lib.make_bbbc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.big_bang_big_crunch(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bbbc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.big_bang_big_crunch(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BOA produces valid result" {
  let config = @lib.make_boa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.butterfly_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BOA improves with more evaluations" {
  let config1 = @lib.make_boa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.butterfly_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_boa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.butterfly_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MFA produces valid result" {
  let config = @lib.make_mfa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mayfly_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MFA improves with more evaluations" {
  let config1 = @lib.make_mfa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mayfly_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_mfa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mayfly_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SBO produces valid result" {
  let config = @lib.make_sbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.secretary_bird_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SBO improves with more evaluations" {
  let config1 = @lib.make_sbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.secretary_bird_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.secretary_bird_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TA produces valid result" {
  let config = @lib.make_ta_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.threshold_accepting(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TA improves with more evaluations" {
  let config1 = @lib.make_ta_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.threshold_accepting(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ta_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.threshold_accepting(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ACOR produces valid result" {
  let config = @lib.make_acor_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.ant_colony_continuous(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ACOR improves with more evaluations" {
  let config1 = @lib.make_acor_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.ant_colony_continuous(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_acor_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.ant_colony_continuous(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "QPSO produces valid result" {
  let config = @lib.make_qpso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.quantum_behaved_pso(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "QPSO improves with more evaluations" {
  let config1 = @lib.make_qpso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.quantum_behaved_pso(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_qpso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.quantum_behaved_pso(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SaDE produces valid result" {
  let config = @lib.make_sade_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.self_adaptive_de(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SaDE improves with more evaluations" {
  let config1 = @lib.make_sade_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.self_adaptive_de(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_sade_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.self_adaptive_de(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BH produces valid result" {
  let config = @lib.make_bh_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.black_hole_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BH improves with more evaluations" {
  let config1 = @lib.make_bh_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.black_hole_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bh_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.black_hole_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LOA produces valid result" {
  let config = @lib.make_loa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.lion_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LOA improves with more evaluations" {
  let config1 = @lib.make_loa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.lion_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_loa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.lion_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BHop produces valid result" {
  let config = @lib.make_bhop_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.basin_hopping(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BHop improves with more evaluations" {
  let config1 = @lib.make_bhop_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.basin_hopping(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bhop_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.basin_hopping(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CoDE produces valid result" {
  let config = @lib.make_code_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.composite_de(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CoDE improves with more evaluations" {
  let config1 = @lib.make_code_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.composite_de(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_code_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.composite_de(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ODE produces valid result" {
  let config = @lib.make_ode_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.opposition_de(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ODE improves with more evaluations" {
  let config1 = @lib.make_ode_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.opposition_de(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ode_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.opposition_de(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CBO produces valid result" {
  let config = @lib.make_cbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.colliding_bodies_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CBO improves with more evaluations" {
  let config1 = @lib.make_cbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.colliding_bodies_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.colliding_bodies_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TWO produces valid result" {
  let config = @lib.make_two_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.teamwork_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TWO improves with more evaluations" {
  let config1 = @lib.make_two_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.teamwork_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_two_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.teamwork_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SPO produces valid result" {
  let config = @lib.make_spo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.spiral_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SPO improves with more evaluations" {
  let config1 = @lib.make_spo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.spiral_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_spo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.spiral_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RTH produces valid result" {
  let config = @lib.make_rth_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.red_tailed_hawk(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RTH improves with more evaluations" {
  let config1 = @lib.make_rth_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.red_tailed_hawk(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_rth_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.red_tailed_hawk(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CRS produces valid result" {
  let config = @lib.make_crs_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.controlled_random_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CRS improves with more evaluations" {
  let config1 = @lib.make_crs_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.controlled_random_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_crs_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.controlled_random_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "UPSO produces valid result" {
  let config = @lib.make_upso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.unified_pso(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "UPSO improves with more evaluations" {
  let config1 = @lib.make_upso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.unified_pso(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_upso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.unified_pso(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RGA produces valid result" {
  let config = @lib.make_rga_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.real_coded_ga(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RGA improves with more evaluations" {
  let config1 = @lib.make_rga_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.real_coded_ga(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_rga_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.real_coded_ga(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "InS produces valid result" {
  let config = @lib.make_ins_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.interior_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "InS improves with more evaluations" {
  let config1 = @lib.make_ins_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.interior_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ins_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.interior_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BMO produces valid result" {
  let config = @lib.make_bmo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.barnacles_mating_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BMO improves with more evaluations" {
  let config1 = @lib.make_bmo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.barnacles_mating_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bmo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.barnacles_mating_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NRBO produces valid result" {
  let config = @lib.make_nrbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.newton_raphson_based_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NRBO improves with more evaluations" {
  let config1 = @lib.make_nrbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.newton_raphson_based_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_nrbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.newton_raphson_based_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SPSA produces valid result" {
  let config = @lib.make_spsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.simultaneous_perturbation_spa(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SPSA improves with more evaluations" {
  let config1 = @lib.make_spsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.simultaneous_perturbation_spa(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_spsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.simultaneous_perturbation_spa(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BAS produces valid result" {
  let config = @lib.make_bas_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.beetle_antennae_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BAS improves with more evaluations" {
  let config1 = @lib.make_bas_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.beetle_antennae_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bas_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.beetle_antennae_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ManSA produces valid result" {
  let config = @lib.make_mantis_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mantis_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ManSA improves with more evaluations" {
  let config1 = @lib.make_mantis_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mantis_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mantis_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mantis_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EWA produces valid result" {
  let config = @lib.make_ewa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.earthworm_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EWA improves with more evaluations" {
  let config1 = @lib.make_ewa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.earthworm_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ewa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.earthworm_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GrSA produces valid result" {
  let config = @lib.make_grsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.group_search_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GrSA improves with more evaluations" {
  let config1 = @lib.make_grsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.group_search_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_grsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.group_search_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FrOA produces valid result" {
  let config = @lib.make_froa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.forest_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FrOA improves with more evaluations" {
  let config1 = @lib.make_froa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.forest_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_froa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.forest_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ChRO produces valid result" {
  let config = @lib.make_chro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.chemical_reaction_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ChRO improves with more evaluations" {
  let config1 = @lib.make_chro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.chemical_reaction_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_chro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.chemical_reaction_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SPBO produces valid result" {
  let config = @lib.make_spbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.student_psychological_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SPBO improves with more evaluations" {
  let config1 = @lib.make_spbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.student_psychological_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_spbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.student_psychological_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PTO produces valid result" {
  let config = @lib.make_pto_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.parallel_tempering_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PTO improves with more evaluations" {
  let config1 = @lib.make_pto_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.parallel_tempering_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pto_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.parallel_tempering_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CLPSO produces valid result" {
  let config = @lib.make_clpso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.comprehensive_learning_pso(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CLPSO improves with more evaluations" {
  let config1 = @lib.make_clpso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.comprehensive_learning_pso(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_clpso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.comprehensive_learning_pso(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BFGS produces valid result" {
  let config = @lib.make_bfgs_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bfgs_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BFGS improves with more evaluations" {
  let config1 = @lib.make_bfgs_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bfgs_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bfgs_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bfgs_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ExtO produces valid result" {
  let config = @lib.make_exto_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.extremal_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ExtO improves with more evaluations" {
  let config1 = @lib.make_exto_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.extremal_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_exto_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.extremal_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DIRECT produces valid result" {
  let config = @lib.make_direct_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dividing_rectangles(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DIRECT improves with more evaluations" {
  let config1 = @lib.make_direct_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dividing_rectangles(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_direct_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dividing_rectangles(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Memetic produces valid result" {
  let config = @lib.make_memetic_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.memetic_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Memetic improves with more evaluations" {
  let config1 = @lib.make_memetic_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.memetic_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_memetic_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.memetic_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "JADE produces valid result" {
  let config = @lib.make_jade_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.jade(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "JADE improves with more evaluations" {
  let config1 = @lib.make_jade_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.jade(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_jade_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.jade(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LCA produces valid result" {
  let config = @lib.make_lca_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.league_championship_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LCA improves with more evaluations" {
  let config1 = @lib.make_lca_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.league_championship_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_lca_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.league_championship_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HGSO produces valid result" {
  let config = @lib.make_hgso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.henry_gas_solubility_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HGSO improves with more evaluations" {
  let config1 = @lib.make_hgso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.henry_gas_solubility_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_hgso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.henry_gas_solubility_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "COOT produces valid result" {
  let config = @lib.make_coot_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.coot_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "COOT improves with more evaluations" {
  let config1 = @lib.make_coot_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.coot_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_coot_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.coot_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DOA produces valid result" {
  let config = @lib.make_doa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dingo_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DOA improves with more evaluations" {
  let config1 = @lib.make_doa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dingo_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_doa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dingo_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "STOA produces valid result" {
  let config = @lib.make_stoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sooty_tern_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "STOA improves with more evaluations" {
  let config1 = @lib.make_stoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sooty_tern_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_stoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sooty_tern_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GKSO produces valid result" {
  let config = @lib.make_gkso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.genghis_khan_shark_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GKSO improves with more evaluations" {
  let config1 = @lib.make_gkso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.genghis_khan_shark_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gkso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.genghis_khan_shark_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FLA produces valid result" {
  let config = @lib.make_fla_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.ficks_law_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FLA improves with more evaluations" {
  let config1 = @lib.make_fla_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.ficks_law_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_fla_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.ficks_law_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HO produces valid result" {
  let config = @lib.make_ho_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.hippopotamus_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HO improves with more evaluations" {
  let config1 = @lib.make_ho_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.hippopotamus_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ho_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.hippopotamus_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NMRA produces valid result" {
  let config = @lib.make_nmra_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.naked_mole_rat_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NMRA improves with more evaluations" {
  let config1 = @lib.make_nmra_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.naked_mole_rat_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_nmra_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.naked_mole_rat_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FFO produces valid result" {
  let config = @lib.make_ffo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fennec_fox_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FFO improves with more evaluations" {
  let config1 = @lib.make_ffo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fennec_fox_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ffo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fennec_fox_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EPO produces valid result" {
  let config = @lib.make_epo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.emperor_penguin_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EPO improves with more evaluations" {
  let config1 = @lib.make_epo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.emperor_penguin_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_epo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.emperor_penguin_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CPO produces valid result" {
  let config = @lib.make_cpo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.crested_porcupine_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CPO improves with more evaluations" {
  let config1 = @lib.make_cpo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.crested_porcupine_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cpo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.crested_porcupine_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MTBO produces valid result" {
  let config = @lib.make_mtbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mountaineering_team_based_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MTBO improves with more evaluations" {
  let config1 = @lib.make_mtbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mountaineering_team_based_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mtbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mountaineering_team_based_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CryStAl produces valid result" {
  let config = @lib.make_crystal_structure_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.crystal_structure_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CryStAl improves with more evaluations" {
  let config1 = @lib.make_crystal_structure_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.crystal_structure_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_crystal_structure_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.crystal_structure_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NOA produces valid result" {
  let config = @lib.make_noa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.nutcracker_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NOA improves with more evaluations" {
  let config1 = @lib.make_noa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.nutcracker_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_noa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.nutcracker_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ESOA produces valid result" {
  let config = @lib.make_esoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.egret_swarm_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ESOA improves with more evaluations" {
  let config1 = @lib.make_esoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.egret_swarm_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_esoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.egret_swarm_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MOA produces valid result" {
  let config = @lib.make_moa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.meerkat_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MOA improves with more evaluations" {
  let config1 = @lib.make_moa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.meerkat_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_moa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.meerkat_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AOS produces valid result" {
  let config = @lib.make_aos_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.atomic_orbital_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AOS improves with more evaluations" {
  let config1 = @lib.make_aos_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.atomic_orbital_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_aos_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.atomic_orbital_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MA produces valid result" {
  let config = @lib.make_ma_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mayfly_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MA improves with more evaluations" {
  let config1 = @lib.make_ma_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mayfly_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ma_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mayfly_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SSO produces valid result" {
  let config = @lib.make_sso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.social_spider_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SSO improves with more evaluations" {
  let config1 = @lib.make_sso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.social_spider_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.social_spider_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SLnO produces valid result" {
  let config = @lib.make_sea_lion_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.sea_lion_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SLnO improves with more evaluations" {
  let config1 = @lib.make_sea_lion_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sea_lion_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sea_lion_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sea_lion_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FOA produces valid result" {
  let config = @lib.make_foa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fruit_fly_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FOA improves with more evaluations" {
  let config1 = @lib.make_foa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fruit_fly_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_foa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fruit_fly_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SFO produces valid result" {
  let config = @lib.make_sfo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sunflower_optimization_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SFO improves with more evaluations" {
  let config1 = @lib.make_sfo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sunflower_optimization_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sfo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sunflower_optimization_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NNA produces valid result" {
  let config = @lib.make_nna_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.neural_network_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NNA improves with more evaluations" {
  let config1 = @lib.make_nna_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.neural_network_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_nna_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.neural_network_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PBO produces valid result" {
  let config = @lib.make_pbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.polar_bear_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PBO improves with more evaluations" {
  let config1 = @lib.make_pbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.polar_bear_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.polar_bear_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FSS produces valid result" {
  let config = @lib.make_fss_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fish_school_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FSS improves with more evaluations" {
  let config1 = @lib.make_fss_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fish_school_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_fss_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fish_school_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CFO produces valid result" {
  let config = @lib.make_cfo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.central_force_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CFO improves with more evaluations" {
  let config1 = @lib.make_cfo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.central_force_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cfo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.central_force_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CGO produces valid result" {
  let config = @lib.make_cgo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.chaos_game_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CGO improves with more evaluations" {
  let config1 = @lib.make_cgo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.chaos_game_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cgo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.chaos_game_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CDO produces valid result" {
  let config = @lib.make_cdo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.chernobyl_disaster_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CDO improves with more evaluations" {
  let config1 = @lib.make_cdo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.chernobyl_disaster_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cdo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.chernobyl_disaster_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ILS produces valid result" {
  let config = @lib.make_ils_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.iterated_local_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ILS improves with more evaluations" {
  let config1 = @lib.make_ils_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.iterated_local_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ils_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.iterated_local_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SCE produces valid result" {
  let config = @lib.make_sce_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.shuffled_complex_evolution(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SCE improves with more evaluations" {
  let config1 = @lib.make_sce_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.shuffled_complex_evolution(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sce_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.shuffled_complex_evolution(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GRASP produces valid result" {
  let config = @lib.make_grasp_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.greedy_randomized_adaptive_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GRASP improves with more evaluations" {
  let config1 = @lib.make_grasp_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.greedy_randomized_adaptive_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_grasp_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.greedy_randomized_adaptive_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BO produces valid result" {
  let config = @lib.make_bo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bonobo_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BO improves with more evaluations" {
  let config1 = @lib.make_bo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bonobo_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bonobo_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FSA produces valid result" {
  let config = @lib.make_fsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.flamingo_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FSA improves with more evaluations" {
  let config1 = @lib.make_fsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.flamingo_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_fsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.flamingo_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SMO produces valid result" {
  let config = @lib.make_smo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.spider_monkey_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SMO improves with more evaluations" {
  let config1 = @lib.make_smo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.spider_monkey_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_smo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.spider_monkey_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HLO produces valid result" {
  let config = @lib.make_hlo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.human_learning_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HLO improves with more evaluations" {
  let config1 = @lib.make_hlo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.human_learning_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_hlo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.human_learning_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FDA produces valid result" {
  let config = @lib.make_fda_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.flow_direction_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FDA improves with more evaluations" {
  let config1 = @lib.make_fda_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.flow_direction_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_fda_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.flow_direction_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CRO produces valid result" {
  let config = @lib.make_cro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.coral_reefs_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CRO improves with more evaluations" {
  let config1 = @lib.make_cro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.coral_reefs_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.coral_reefs_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AAA produces valid result" {
  let config = @lib.make_aaa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_algae_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AAA improves with more evaluations" {
  let config1 = @lib.make_aaa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_algae_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_aaa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_algae_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "VSA produces valid result" {
  let config = @lib.make_vsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.vortex_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "VSA improves with more evaluations" {
  let config1 = @lib.make_vsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.vortex_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_vsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.vortex_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "VPL produces valid result" {
  let config = @lib.make_vpl_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.volleyball_premier_league(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "VPL improves with more evaluations" {
  let config1 = @lib.make_vpl_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.volleyball_premier_league(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_vpl_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.volleyball_premier_league(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WEO produces valid result" {
  let config = @lib.make_weo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.water_evaporation_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WEO improves with more evaluations" {
  let config1 = @lib.make_weo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.water_evaporation_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_weo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.water_evaporation_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SNS produces valid result" {
  let config = @lib.make_sns_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.social_network_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SNS improves with more evaluations" {
  let config1 = @lib.make_sns_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.social_network_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_sns_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.social_network_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TCO produces valid result" {
  let config = @lib.make_tco_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.termite_colony_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TCO improves with more evaluations" {
  let config1 = @lib.make_tco_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.termite_colony_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_tco_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.termite_colony_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RO produces valid result" {
  let config = @lib.make_ro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.ray_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RO improves with more evaluations" {
  let config1 = @lib.make_ro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.ray_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.ray_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GbSA produces valid result" {
  let config = @lib.make_gbsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.galaxy_based_search_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GbSA improves with more evaluations" {
  let config1 = @lib.make_gbsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.galaxy_based_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gbsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.galaxy_based_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BA produces valid result" {
  let config = @lib.make_ba_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bees_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BA improves with more evaluations" {
  let config1 = @lib.make_ba_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bees_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ba_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bees_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CO produces valid result" {
  let config = @lib.make_co_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cheetah_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CO improves with more evaluations" {
  let config1 = @lib.make_co_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cheetah_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_co_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cheetah_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CSS produces valid result" {
  let config = @lib.make_css_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.charged_system_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CSS improves with more evaluations" {
  let config1 = @lib.make_css_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.charged_system_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_css_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.charged_system_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MBA produces valid result" {
  let config = @lib.make_mba_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mine_blast_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MBA improves with more evaluations" {
  let config1 = @lib.make_mba_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mine_blast_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_mba_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mine_blast_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PCA produces valid result" {
  let config = @lib.make_pca_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.particle_collision_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PCA improves with more evaluations" {
  let config1 = @lib.make_pca_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.particle_collision_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pca_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.particle_collision_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "IMO produces valid result" {
  let config = @lib.make_imo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.ions_motion_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "IMO improves with more evaluations" {
  let config1 = @lib.make_imo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.ions_motion_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_imo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.ions_motion_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EMA produces valid result" {
  let config = @lib.make_ema_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.exchange_market_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EMA improves with more evaluations" {
  let config1 = @lib.make_ema_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.exchange_market_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ema_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.exchange_market_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WCO produces valid result" {
  let config = @lib.make_wco_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.world_cup_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WCO improves with more evaluations" {
  let config1 = @lib.make_wco_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.world_cup_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_wco_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.world_cup_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "AEFA produces valid result" {
  let config = @lib.make_aefa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_electric_field_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "AEFA improves with more evaluations" {
  let config1 = @lib.make_aefa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_electric_field_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_aefa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_electric_field_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "OIO produces valid result" {
  let config = @lib.make_oio_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.optics_inspired_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "OIO improves with more evaluations" {
  let config1 = @lib.make_oio_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.optics_inspired_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_oio_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.optics_inspired_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LAPO produces valid result" {
  let config = @lib.make_lapo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.lightning_attachment_procedure_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LAPO improves with more evaluations" {
  let config1 = @lib.make_lapo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.lightning_attachment_procedure_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_lapo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.lightning_attachment_procedure_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FGBO produces valid result" {
  let config = @lib.make_fgbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.football_game_based_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FGBO improves with more evaluations" {
  let config1 = @lib.make_fgbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.football_game_based_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_fgbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.football_game_based_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EA produces valid result" {
  let config = @lib.make_ea_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.earthquake_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EA improves with more evaluations" {
  let config1 = @lib.make_ea_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.earthquake_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ea_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.earthquake_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LAHC produces valid result" {
  let config = @lib.make_lahc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.late_acceptance_hill_climbing(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LAHC improves with more evaluations" {
  let config1 = @lib.make_lahc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.late_acceptance_hill_climbing(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_lahc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.late_acceptance_hill_climbing(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HuS produces valid result" {
  let config = @lib.make_hus_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.hunting_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HuS improves with more evaluations" {
  let config1 = @lib.make_hus_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.hunting_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_hus_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.hunting_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RRT produces valid result" {
  let config = @lib.make_rrt_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.record_to_record_travel(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RRT improves with more evaluations" {
  let config1 = @lib.make_rrt_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.record_to_record_travel(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_rrt_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.record_to_record_travel(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NRO produces valid result" {
  let config = @lib.make_nro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.nuclear_reaction_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NRO improves with more evaluations" {
  let config1 = @lib.make_nro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.nuclear_reaction_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_nro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.nuclear_reaction_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HTS produces valid result" {
  let config = @lib.make_hts_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.heat_transfer_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HTS improves with more evaluations" {
  let config1 = @lib.make_hts_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.heat_transfer_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_hts_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.heat_transfer_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PGO produces valid result" {
  let config = @lib.make_pgo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.plasma_generation_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PGO improves with more evaluations" {
  let config1 = @lib.make_pgo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.plasma_generation_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_pgo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.plasma_generation_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "IWD produces valid result" {
  let config = @lib.make_iwd_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.intelligent_water_drops(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "IWD improves with more evaluations" {
  let config1 = @lib.make_iwd_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.intelligent_water_drops(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_iwd_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.intelligent_water_drops(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RFD produces valid result" {
  let config = @lib.make_rfd_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.river_formation_dynamics(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RFD improves with more evaluations" {
  let config1 = @lib.make_rfd_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.river_formation_dynamics(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_rfd_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.river_formation_dynamics(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LSA produces valid result" {
  let config = @lib.make_lsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.lightning_search_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LSA improves with more evaluations" {
  let config1 = @lib.make_lsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.lightning_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_lsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.lightning_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SS produces valid result" {
  let config = @lib.make_ss_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.scatter_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SS improves with more evaluations" {
  let config1 = @lib.make_ss_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.scatter_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ss_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.scatter_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "OBA produces valid result" {
  let config = @lib.make_oba_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.old_bachelor_acceptance(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "OBA improves with more evaluations" {
  let config1 = @lib.make_oba_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.old_bachelor_acceptance(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_oba_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.old_bachelor_acceptance(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "APO produces valid result" {
  let config = @lib.make_apo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.artificial_physics_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "APO improves with more evaluations" {
  let config1 = @lib.make_apo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.artificial_physics_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_apo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.artificial_physics_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HBO produces valid result" {
  let config = @lib.make_hbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.heap_based_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HBO improves with more evaluations" {
  let config1 = @lib.make_hbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.heap_based_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_hbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.heap_based_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SAO produces valid result" {
  let config = @lib.make_sao_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.snow_ablation_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SAO improves with more evaluations" {
  let config1 = @lib.make_sao_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.snow_ablation_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sao_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.snow_ablation_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BRO produces valid result" {
  let config = @lib.make_bro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.battle_royale_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BRO improves with more evaluations" {
  let config1 = @lib.make_bro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.battle_royale_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.battle_royale_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EFO produces valid result" {
  let config = @lib.make_efo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.electromagnetic_field_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EFO improves with more evaluations" {
  let config1 = @lib.make_efo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.electromagnetic_field_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_efo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.electromagnetic_field_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GBMO produces valid result" {
  let config = @lib.make_gbmo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gases_brownian_motion_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GBMO improves with more evaluations" {
  let config1 = @lib.make_gbmo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gases_brownian_motion_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gbmo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gases_brownian_motion_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DSA produces valid result" {
  let config = @lib.make_dsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.differential_search_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DSA improves with more evaluations" {
  let config1 = @lib.make_dsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.differential_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_dsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.differential_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GSK produces valid result" {
  let config = @lib.make_gsk_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gaining_sharing_knowledge(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GSK improves with more evaluations" {
  let config1 = @lib.make_gsk_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gaining_sharing_knowledge(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gsk_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gaining_sharing_knowledge(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ALNS produces valid result" {
  let config = @lib.make_alns_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.adaptive_large_neighborhood_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ALNS improves with more evaluations" {
  let config1 = @lib.make_alns_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.adaptive_large_neighborhood_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_alns_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.adaptive_large_neighborhood_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SGO produces valid result" {
  let config = @lib.make_sgo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.social_group_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SGO improves with more evaluations" {
  let config1 = @lib.make_sgo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.social_group_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sgo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.social_group_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ASO produces valid result" {
  let config = @lib.make_aso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.atom_search_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ASO improves with more evaluations" {
  let config1 = @lib.make_aso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.atom_search_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_aso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.atom_search_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BCMO produces valid result" {
  let config = @lib.make_bcmo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.balancing_composite_motion_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BCMO improves with more evaluations" {
  let config1 = @lib.make_bcmo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.balancing_composite_motion_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_bcmo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.balancing_composite_motion_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MKE produces valid result" {
  let config = @lib.make_mke_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.monkey_king_evolution(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MKE improves with more evaluations" {
  let config1 = @lib.make_mke_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.monkey_king_evolution(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_mke_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.monkey_king_evolution(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RFO produces valid result" {
  let config = @lib.make_rfo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.red_fox_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RFO improves with more evaluations" {
  let config1 = @lib.make_rfo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.red_fox_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_rfo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.red_fox_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SAR produces valid result" {
  let config = @lib.make_sar_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.search_and_rescue_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SAR improves with more evaluations" {
  let config1 = @lib.make_sar_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.search_and_rescue_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sar_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.search_and_rescue_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ChSA produces valid result" {
  let config = @lib.make_chsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.chameleon_swarm_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ChSA improves with more evaluations" {
  let config1 = @lib.make_chsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.chameleon_swarm_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_chsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.chameleon_swarm_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RDA produces valid result" {
  let config = @lib.make_rda_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.red_deer_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RDA improves with more evaluations" {
  let config1 = @lib.make_rda_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.red_deer_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_rda_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.red_deer_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SqSA produces valid result" {
  let config = @lib.make_sqsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.squirrel_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SqSA improves with more evaluations" {
  let config1 = @lib.make_sqsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.squirrel_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sqsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.squirrel_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Puma produces valid result" {
  let config = @lib.make_puma_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.puma_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Puma improves with more evaluations" {
  let config1 = @lib.make_puma_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.puma_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_puma_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.puma_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "StMO produces valid result" {
  let config = @lib.make_stmo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.starling_murmuration_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "StMO improves with more evaluations" {
  let config1 = @lib.make_stmo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.starling_murmuration_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_stmo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.starling_murmuration_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SFOA produces valid result" {
  let config = @lib.make_sfoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.swordfish_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SFOA improves with more evaluations" {
  let config1 = @lib.make_sfoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.swordfish_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_sfoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.swordfish_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ABO produces valid result" {
  let config = @lib.make_abo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.african_buffalo_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ABO improves with more evaluations" {
  let config1 = @lib.make_abo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.african_buffalo_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_abo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.african_buffalo_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Gannet produces valid result" {
  let config = @lib.make_gannet_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gannet_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Gannet improves with more evaluations" {
  let config1 = @lib.make_gannet_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gannet_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_gannet_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gannet_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Narwhal produces valid result" {
  let config = @lib.make_narwhal_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.narwhal_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Narwhal improves with more evaluations" {
  let config1 = @lib.make_narwhal_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.narwhal_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_narwhal_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.narwhal_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Serval produces valid result" {
  let config = @lib.make_serval_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.serval_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Serval improves with more evaluations" {
  let config1 = @lib.make_serval_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.serval_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_serval_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.serval_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ChSO produces valid result" {
  let config = @lib.make_chso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.chicken_swarm_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ChSO improves with more evaluations" {
  let config1 = @lib.make_chso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.chicken_swarm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_chso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.chicken_swarm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PR produces valid result" {
  let config = @lib.make_pr_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.path_relinking(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PR improves with more evaluations" {
  let config1 = @lib.make_pr_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.path_relinking(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_pr_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.path_relinking(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CrSO produces valid result" {
  let config = @lib.make_crso_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cockroach_swarm_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CrSO improves with more evaluations" {
  let config1 = @lib.make_crso_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cockroach_swarm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_crso_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cockroach_swarm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "JS produces valid result" {
  let config = @lib.make_js_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.jellyfish_search_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "JS improves with more evaluations" {
  let config1 = @lib.make_js_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.jellyfish_search_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_js_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.jellyfish_search_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WSA produces valid result" {
  let config = @lib.make_wsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.water_strider_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WSA improves with more evaluations" {
  let config1 = @lib.make_wsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.water_strider_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_wsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.water_strider_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SeHO produces valid result" {
  let config = @lib.make_seho_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sea_horse_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SeHO improves with more evaluations" {
  let config1 = @lib.make_seho_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sea_horse_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_seho_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sea_horse_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WBO produces valid result" {
  let config = @lib.make_wbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.wildebeest_herd_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WBO improves with more evaluations" {
  let config1 = @lib.make_wbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.wildebeest_herd_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_wbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.wildebeest_herd_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "OWL produces valid result" {
  let config = @lib.make_owl_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.owl_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "OWL improves with more evaluations" {
  let config1 = @lib.make_owl_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.owl_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_owl_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.owl_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "RRO produces valid result" {
  let config = @lib.make_rro_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.raven_roosting_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "RRO improves with more evaluations" {
  let config1 = @lib.make_rro_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.raven_roosting_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_rro_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.raven_roosting_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SaFO produces valid result" {
  let config = @lib.make_safo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.sailfish_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SaFO improves with more evaluations" {
  let config1 = @lib.make_safo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.sailfish_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_safo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.sailfish_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GPC produces valid result" {
  let config = @lib.make_gpc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.giza_pyramids_construction(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GPC improves with more evaluations" {
  let config1 = @lib.make_gpc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.giza_pyramids_construction(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gpc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.giza_pyramids_construction(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SSD produces valid result" {
  let config = @lib.make_ssd_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.social_ski_driver(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SSD improves with more evaluations" {
  let config1 = @lib.make_ssd_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.social_ski_driver(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ssd_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.social_ski_driver(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "FDO produces valid result" {
  let config = @lib.make_fdo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.fitness_dependent_optimizer(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "FDO improves with more evaluations" {
  let config1 = @lib.make_fdo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.fitness_dependent_optimizer(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_fdo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.fitness_dependent_optimizer(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CFA produces valid result" {
  let config = @lib.make_cfa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cuttlefish_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CFA improves with more evaluations" {
  let config1 = @lib.make_cfa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cuttlefish_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cfa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cuttlefish_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "OPA produces valid result" {
  let config = @lib.make_opa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.orca_predation_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "OPA improves with more evaluations" {
  let config1 = @lib.make_opa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.orca_predation_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_opa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.orca_predation_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "APuf produces valid result" {
  let config = @lib.make_apuf_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.arctic_puffin_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "APuf improves with more evaluations" {
  let config1 = @lib.make_apuf_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.arctic_puffin_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_apuf_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.arctic_puffin_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WPA produces valid result" {
  let config = @lib.make_wpa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.wolf_pack_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WPA improves with more evaluations" {
  let config1 = @lib.make_wpa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.wolf_pack_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_wpa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.wolf_pack_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Racc produces valid result" {
  let config = @lib.make_racc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.raccoon_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Racc improves with more evaluations" {
  let config1 = @lib.make_racc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.raccoon_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_racc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.raccoon_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Monkey produces valid result" {
  let config = @lib.make_monkey_search_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.monkey_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Monkey improves with more evaluations" {
  let config1 = @lib.make_monkey_search_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.monkey_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_monkey_search_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.monkey_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BSwA produces valid result" {
  let config = @lib.make_bswa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bird_swarm_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BSwA improves with more evaluations" {
  let config1 = @lib.make_bswa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bird_swarm_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bswa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bird_swarm_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TGA produces valid result" {
  let config = @lib.make_tga_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.tree_growth_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TGA improves with more evaluations" {
  let config1 = @lib.make_tga_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.tree_growth_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_tga_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.tree_growth_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PPA produces valid result" {
  let config = @lib.make_ppa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.plant_propagation_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PPA improves with more evaluations" {
  let config1 = @lib.make_ppa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.plant_propagation_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ppa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.plant_propagation_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "HBMO produces valid result" {
  let config = @lib.make_hbmo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.honey_bee_mating_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "HBMO improves with more evaluations" {
  let config1 = @lib.make_hbmo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.honey_bee_mating_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_hbmo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.honey_bee_mating_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Crayfish produces valid result" {
  let config = @lib.make_crayfish_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.crayfish_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Crayfish improves with more evaluations" {
  let config1 = @lib.make_crayfish_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.crayfish_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_crayfish_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.crayfish_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "QSA produces valid result" {
  let config = @lib.make_qsa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.queuing_search_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "QSA improves with more evaluations" {
  let config1 = @lib.make_qsa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.queuing_search_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_qsa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.queuing_search_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Starfish produces valid result" {
  let config = @lib.make_starfish_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.starfish_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Starfish improves with more evaluations" {
  let config1 = @lib.make_starfish_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.starfish_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_starfish_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.starfish_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TBO produces valid result" {
  let config = @lib.make_tbo_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.tiger_beetle_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TBO improves with more evaluations" {
  let config1 = @lib.make_tbo_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.tiger_beetle_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_tbo_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.tiger_beetle_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Lyrebird produces valid result" {
  let config = @lib.make_lyrebird_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.lyrebird_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Lyrebird improves with more evaluations" {
  let config1 = @lib.make_lyrebird_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.lyrebird_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_lyrebird_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.lyrebird_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SEOA produces valid result" {
  let config = @lib.make_seoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.social_emotional_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SEOA improves with more evaluations" {
  let config1 = @lib.make_seoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.social_emotional_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_seoa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.social_emotional_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GFO produces valid result" {
  let config = @lib.make_guineafowl_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.guineafowl_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GFO improves with more evaluations" {
  let config1 = @lib.make_guineafowl_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.guineafowl_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_guineafowl_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.guineafowl_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DEO produces valid result" {
  let config = @lib.make_dolphin_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.dolphin_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DEO improves with more evaluations" {
  let config1 = @lib.make_dolphin_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dolphin_optimization(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_dolphin_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dolphin_optimization(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SOMA produces valid result" {
  let config = @lib.make_soma_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.self_organizing_migration(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SOMA improves with more evaluations" {
  let config1 = @lib.make_soma_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.self_organizing_migration(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_soma_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.self_organizing_migration(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "SDS produces valid result" {
  let config = @lib.make_sds_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.stochastic_diffusion_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SDS improves with more evaluations" {
  let config1 = @lib.make_sds_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.stochastic_diffusion_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sds_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.stochastic_diffusion_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DAnl produces valid result" {
  let config = @lib.make_dual_annealing_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.dual_annealing(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SQA produces valid result" {
  let config = @lib.make_sqa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.simulated_quantum_annealing(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "SQA improves with more evaluations" {
  let config1 = @lib.make_sqa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.simulated_quantum_annealing(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_sqa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.simulated_quantum_annealing(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DAnl improves with more evaluations" {
  let config1 = @lib.make_dual_annealing_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.dual_annealing(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_dual_annealing_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.dual_annealing(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "COBYLA produces valid result" {
  let config = @lib.make_cobyla_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cobyla_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "COBYLA improves with more evaluations" {
  let config1 = @lib.make_cobyla_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cobyla_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_cobyla_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cobyla_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NS produces valid result" {
  let config = @lib.make_nested_sampling_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.nested_sampling(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NS improves with more evaluations" {
  let config1 = @lib.make_nested_sampling_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.nested_sampling(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_nested_sampling_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.nested_sampling(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MADS produces valid result" {
  let config = @lib.make_mads_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.mesh_adaptive_direct_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MADS improves with more evaluations" {
  let config1 = @lib.make_mads_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.mesh_adaptive_direct_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mads_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.mesh_adaptive_direct_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "NEWUOA produces valid result" {
  let config = @lib.make_newuoa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.newuoa_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "NEWUOA improves with more evaluations" {
  let config1 = @lib.make_newuoa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.newuoa_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_newuoa_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.newuoa_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PA produces valid result" {
  let config = @lib.make_population_annealing_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.population_annealing(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PA improves with more evaluations" {
  let config1 = @lib.make_population_annealing_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.population_annealing(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_population_annealing_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.population_annealing(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "OAIN produces valid result" {
  let config = @lib.make_oain_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.opt_ainet(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "OAIN improves with more evaluations" {
  let config1 = @lib.make_oain_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.opt_ainet(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_oain_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.opt_ainet(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "xNES produces valid result" {
  let config = @lib.make_xnes_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.exponential_natural_evolution_strategy(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "xNES improves with more evaluations" {
  let config1 = @lib.make_xnes_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.exponential_natural_evolution_strategy(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_xnes_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.exponential_natural_evolution_strategy(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Demon produces valid result" {
  let config = @lib.make_demon_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.demon_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Demon improves with more evaluations" {
  let config1 = @lib.make_demon_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.demon_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_demon_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.demon_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MCS produces valid result" {
  let config = @lib.make_mcs_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.multilevel_coordinate_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MCS improves with more evaluations" {
  let config1 = @lib.make_mcs_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.multilevel_coordinate_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mcs_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.multilevel_coordinate_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "WL produces valid result" {
  let config = @lib.make_wl_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.wang_landau_sampling(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "WL improves with more evaluations" {
  let config1 = @lib.make_wl_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.wang_landau_sampling(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_wl_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.wang_landau_sampling(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LSHADE produces valid result" {
  let config = @lib.make_lshade_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.lshade_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LSHADE improves with more evaluations" {
  let config1 = @lib.make_lshade_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.lshade_algorithm(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_lshade_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.lshade_algorithm(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "1+1ES produces valid result" {
  let config = @lib.make_one_plus_one_es_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.one_plus_one_evolution_strategy(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "1+1ES improves with more evaluations" {
  let config1 = @lib.make_one_plus_one_es_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.one_plus_one_evolution_strategy(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_one_plus_one_es_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.one_plus_one_evolution_strategy(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "PRAXIS produces valid result" {
  let config = @lib.make_praxis_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.praxis_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "PRAXIS improves with more evaluations" {
  let config1 = @lib.make_praxis_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.praxis_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_praxis_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.praxis_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "BOBYQA produces valid result" {
  let config = @lib.make_bobyqa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.bobyqa_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "BOBYQA improves with more evaluations" {
  let config1 = @lib.make_bobyqa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.bobyqa_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_bobyqa_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.bobyqa_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "ASA produces valid result" {
  let config = @lib.make_asa_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.adaptive_simulated_annealing(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "ASA improves with more evaluations" {
  let config1 = @lib.make_asa_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.adaptive_simulated_annealing(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_asa_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.adaptive_simulated_annealing(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GSS produces valid result" {
  let config = @lib.make_gss_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.golden_section_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GSS improves with more evaluations" {
  let config1 = @lib.make_gss_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.golden_section_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_gss_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.golden_section_search(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "Subplex produces valid result" {
  let config = @lib.make_subplex_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.subplex_optimizer(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "Subplex improves with more evaluations" {
  let config1 = @lib.make_subplex_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.subplex_optimizer(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_subplex_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.subplex_optimizer(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MLSL produces valid result" {
  let config = @lib.make_mlsl_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.multi_level_single_linkage(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MLSL improves with more evaluations" {
  let config1 = @lib.make_mlsl_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.multi_level_single_linkage(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mlsl_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.multi_level_single_linkage(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "LBB produces valid result" {
  let config = @lib.make_lipschitz_config(
    500,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng = @lib.make_rng(42.0)
  let result = @lib.lipschitz_branch_and_bound(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "LBB improves with more evaluations" {
  let config1 = @lib.make_lipschitz_config(
    100,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.lipschitz_branch_and_bound(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_lipschitz_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.lipschitz_branch_and_bound(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MTS finds good solution on sphere" {
  let config = @lib.make_mts_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.multiple_trajectory_search(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MTS improves with more evaluations" {
  let config1 = @lib.make_mts_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.multiple_trajectory_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mts_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.multiple_trajectory_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "TPE finds good solution on sphere" {
  let config = @lib.make_tpe_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.tree_parzen_estimator(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "TPE improves with more evaluations" {
  let config1 = @lib.make_tpe_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.tree_parzen_estimator(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_tpe_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.tree_parzen_estimator(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "EGO finds good solution on sphere" {
  let config = @lib.make_ego_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.efficient_global_optimization(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "EGO improves with more evaluations" {
  let config1 = @lib.make_ego_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.efficient_global_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_ego_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.efficient_global_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "GR finds good solution on sphere" {
  let config = @lib.make_gr_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.gutmann_ramm_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "GR improves with more evaluations" {
  let config1 = @lib.make_gr_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.gutmann_ramm_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_gr_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.gutmann_ramm_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "sepCMA finds good solution on sphere" {
  let config = @lib.make_sepcma_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.separable_cma_es(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "sepCMA improves with more evaluations" {
  let config1 = @lib.make_sepcma_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.separable_cma_es(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_sepcma_config(
    5000,
    3,
    @lib.make_bounds(-5.0, 5.0),
    25,
  )
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.separable_cma_es(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "IFilt finds good solution on sphere" {
  let config = @lib.make_ifilt_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.implicit_filtering(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "IFilt improves with more evaluations" {
  let config1 = @lib.make_ifilt_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.implicit_filtering(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_ifilt_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.implicit_filtering(config2, @lib.sphere_function, rng2)

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CG finds good solution on sphere" {
  let config = @lib.make_cg_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.conjugate_gradient_method(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CG improves with more evaluations" {
  let config1 = @lib.make_cg_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.conjugate_gradient_method(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cg_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.conjugate_gradient_method(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "DHO finds good solution on sphere" {
  let config = @lib.make_dho_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.deer_hunting_optimization(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "DHO improves with more evaluations" {
  let config1 = @lib.make_dho_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.deer_hunting_optimization(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_dho_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.deer_hunting_optimization(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CHC finds good solution on sphere" {
  let config = @lib.make_chc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cross_generation_elitist_algorithm(
    config, @lib.sphere_function, rng,
  )
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CHC improves with more evaluations" {
  let config1 = @lib.make_chc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cross_generation_elitist_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_chc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cross_generation_elitist_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MDS finds good solution on sphere" {
  let config = @lib.make_mds_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.multi_directional_search(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MDS improves with more evaluations" {
  let config1 = @lib.make_mds_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.multi_directional_search(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mds_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.multi_directional_search(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "CC finds good solution on sphere" {
  let config = @lib.make_cc_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.cooperative_coevolution(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "CC improves with more evaluations" {
  let config1 = @lib.make_cc_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.cooperative_coevolution(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_cc_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.cooperative_coevolution(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}

///|
test "MGA finds good solution on sphere" {
  let config = @lib.make_mga_config(500, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng = @lib.make_rng(42.0)
  let result = @lib.micro_genetic_algorithm(config, @lib.sphere_function, rng)
  assert_true(result.evaluations > 0)
  assert_true(result.best_value >= 0.0)
}

///|
test "MGA improves with more evaluations" {
  let config1 = @lib.make_mga_config(100, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.micro_genetic_algorithm(
    config1, @lib.sphere_function, rng1,
  )

  let config2 = @lib.make_mga_config(5000, 3, @lib.make_bounds(-5.0, 5.0), 25)
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.micro_genetic_algorithm(
    config2, @lib.sphere_function, rng2,
  )

  assert_true(result2.best_value <= result1.best_value)
}
