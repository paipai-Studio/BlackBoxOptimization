// BlackBox Optimization 测试


fn approx_eq(a: Double, b: Double, tol: Double) -> Bool {
  let diff = a - b
  let abs_diff = if diff < 0.0 { -diff } else { diff }
  abs_diff < tol
}

test "cos_approx basic values" {
  assert_eq!(@lib.cos_approx(0.0), 1.0)
  assert_eq!(approx_eq(@lib.cos_approx(@lib.PI), -1.0, 0.05), true)
  assert_eq!(@lib.cos_approx(@lib.PI / 2.0).round(), 0.0)
}

test "cos_approx periodicity" {
  let result = @lib.cos_approx(3.0 * @lib.PI)
  assert_eq!(approx_eq(result, -1.0, 0.05), true)
}

test "exp_approx basic values" {
  assert_eq!(@lib.exp_approx(0.0), 1.0)
  let result = @lib.exp_approx(1.0)
  assert_eq!(result.round(), 3.0)
}

test "exp_approx negative values" {
  let result = @lib.exp_approx(-1.0)
  assert_eq!(result > 0.0, true)
  assert_eq!(result < 1.0, true)
}

test "max and min functions" {
  assert_eq!(@lib.max(3.0, 5.0), 5.0)
  assert_eq!(@lib.max(-1.0, -5.0), -1.0)
  assert_eq!(@lib.min(3.0, 5.0), 3.0)
  assert_eq!(@lib.min(-1.0, -5.0), -5.0)
}

test "abs_val function" {
  assert_eq!(@lib.abs_val(5.0), 5.0)
  assert_eq!(@lib.abs_val(-5.0), 5.0)
  assert_eq!(@lib.abs_val(0.0), 0.0)
}

test "reduce_to_pi" {
  let result = @lib.reduce_to_pi(0.0)
  assert_eq!(result, 0.0)
}

test "RNG generates values in [0,1)" {
  let rng = @lib.make_rng(42.0)
  let (_, val) = @lib.RNG::next(rng)
  assert_eq!(val >= 0.0, true)
  assert_eq!(val < 1.0, true)
}

test "RNG deterministic with same seed" {
  let rng1 = @lib.make_rng(123.0)
  let (rng2, val1) = @lib.RNG::next(rng1)
  let (_, val2) = @lib.RNG::next(rng2)

  let rng3 = @lib.make_rng(123.0)
  let (rng4, val3) = @lib.RNG::next(rng3)
  let (_, val4) = @lib.RNG::next(rng4)

  assert_eq!(val1, val3)
  assert_eq!(val2, val4)
}

test "sphere function at origin" {
  let zeros: Array[Double] = [0.0, 0.0, 0.0]
  let result = @lib.sphere_function(zeros)
  assert_eq!(result, 0.0)
}

test "sphere function at non-origin" {
  let input: Array[Double] = [1.0, 2.0, 3.0]
  let result = @lib.sphere_function(input)
  assert_eq!(result, 14.0)
}

test "rastrigin function at origin" {
  let zeros: Array[Double] = [0.0, 0.0, 0.0]
  let result = @lib.rastrigin_function(zeros)
  assert_eq!(result, 0.0)
}

test "rastrigin function non-negative" {
  let input: Array[Double] = [1.0, 2.0, 3.0]
  let result = @lib.rastrigin_function(input)
  assert_eq!(result >= 0.0, true)
}

test "rosenbrock function at optimum" {
  let input: Array[Double] = [1.0, 1.0, 1.0]
  let result = @lib.rosenbrock_function(input)
  assert_eq!(result, 0.0)
}

test "rosenbrock function at non-optimal" {
  let input: Array[Double] = [0.0, 0.0, 0.0]
  let result = @lib.rosenbrock_function(input)
  assert_eq!(result > 0.0, true)
}

test "array_sum function" {
  let arr: Array[Double] = [1.0, 2.0, 3.0, 4.0]
  let result = @lib.array_sum(arr)
  assert_eq!(result, 10.0)
}

test "array_min function" {
  let arr: Array[Double] = [3.0, 1.0, 4.0, 1.0, 5.0]
  let result = @lib.array_min(arr)
  assert_eq!(result, 1.0)
}

test "array_max function" {
  let arr: Array[Double] = [3.0, 1.0, 4.0, 1.0, 5.0]
  let result = @lib.array_max(arr)
  assert_eq!(result, 5.0)
}

test "array_avg function" {
  let arr: Array[Double] = [2.0, 4.0, 6.0]
  let result = @lib.array_avg(arr)
  assert_eq!(result, 4.0)
}

test "calc_relative_error with zero optimal" {
  let error1 = @lib.calc_relative_error(0.0, 0.0)
  assert_eq!(error1, 0.0)

  let error2 = @lib.calc_relative_error(5.0, 0.0)
  assert_eq!(error2, 5.0)
}

test "calc_relative_error with non-zero optimal" {
  let error = @lib.calc_relative_error(110.0, 100.0)
  assert_eq!(error, 0.1)
}

test "random_solution generates within bounds" {
  let bounds = @lib.make_bounds(-5.0, 5.0)
  let rng = @lib.make_rng(42.0)
  let (_, solution) = @lib.random_solution(5, bounds, rng)

  for i in 0..<solution.length() {
    assert_eq!(solution[i] >= -5.0, true)
    assert_eq!(solution[i] <= 5.0, true)
  }
}

test "random_search improves with more evaluations" {
  let config1 = @lib.make_random_search_config(100, 3, @lib.make_bounds(-5.0, 5.0))
  let rng1 = @lib.make_rng(42.0)
  let result1 = @lib.random_search(config1, @lib.sphere_function, rng1)

  let config2 = @lib.make_random_search_config(10000, 3, @lib.make_bounds(-5.0, 5.0))
  let rng2 = @lib.make_rng(42.0)
  let result2 = @lib.random_search(config2, @lib.sphere_function, rng2)

  assert_eq!(result2.best_value <= result1.best_value, true)
}

test "hill_climbing produces valid result" {
  let config = @lib.make_hill_climbing_config(500, 3, @lib.make_bounds(-5.0, 5.0), 0.5, 5)
  let rng = @lib.make_rng(42.0)
  let result = @lib.hill_climbing(config, @lib.sphere_function, rng)

  assert_eq!(result.evaluations > 0, true)
  assert_eq!(result.best_value >= 0.0, true)
}

test "simulated_annealing produces valid result" {
  let config = @lib.make_simulated_annealing_config(500, 3, @lib.make_bounds(-5.0, 5.0), 100.0, 0.995, 0.5)
  let rng = @lib.make_rng(42.0)
  let result = @lib.simulated_annealing(config, @lib.sphere_function, rng)

  assert_eq!(result.evaluations > 0, true)
  assert_eq!(result.best_value >= 0.0, true)
}

test "evaluate function works" {
  let f = @lib.sphere_function
  let solution: Array[Double] = [1.0, 2.0]
  let result = @lib.evaluate(f, solution)
  assert_eq!(result, 5.0)
}
