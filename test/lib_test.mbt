// BlackBox Optimization 测试

import lib

test "cos_approx basic values" {
  assert_eq!(lib.cos_approx(0.0), 1.0, "cos(0) should be 1")
  assert_eq!(lib.cos_approx(lib.PI), -1.0, "cos(pi) should be -1")
  assert_eq!(lib.cos_approx(lib.PI / 2.0).round(), 0.0, "cos(pi/2) should be ~0")
}

test "cos_approx periodicity" {
  let result = lib.cos_approx(3.0 * lib.PI)
  assert_eq!(result, -1.0, "cos(3*pi) should be -1")
}

test "exp_approx basic values" {
  assert_eq!(lib.exp_approx(0.0), 1.0, "exp(0) should be 1")
  let result = lib.exp_approx(1.0)
  assert_eq!(result.round(), 3.0, "exp(1) should be ~2.718")
}

test "exp_approx negative values" {
  let result = lib.exp_approx(-1.0)
  assert!(result > 0.0, "exp(-1) should be positive")
  assert!(result < 1.0, "exp(-1) should be less than 1")
}

test "max and min functions" {
  assert_eq!(lib.max(3.0, 5.0), 5.0, "max(3, 5) should be 5")
  assert_eq!(lib.max(-1.0, -5.0), -1.0, "max(-1, -5) should be -1")
  assert_eq!(lib.min(3.0, 5.0), 3.0, "min(3, 5) should be 3")
  assert_eq!(lib.min(-1.0, -5.0), -5.0, "min(-1, -5) should be -5")
}

test "abs_val function" {
  assert_eq!(lib.abs_val(5.0), 5.0, "abs(5) should be 5")
  assert_eq!(lib.abs_val(-5.0), 5.0, "abs(-5) should be 5")
  assert_eq!(lib.abs_val(0.0), 0.0, "abs(0) should be 0")
}

test "reduce_to_pi" {
  let result = lib.reduce_to_pi(0.0)
  assert_eq!(result, 0.0, "reduce_to_pi(0) should be 0")
}

test "RNG generates values in [0,1)" {
  let rng = lib.RNG::{ seed: 42.0 }
  let (_, val) = lib.RNG::next(rng)
  assert!(val >= 0.0, "RNG value should be >= 0")
  assert!(val < 1.0, "RNG value should be < 1")
}

test "RNG deterministic with same seed" {
  let rng1 = lib.RNG::{ seed: 123.0 }
  let (rng2, val1) = lib.RNG::next(rng1)
  let (_, val2) = lib.RNG::next(rng2)

  let rng3 = lib.RNG::{ seed: 123.0 }
  let (rng4, val3) = lib.RNG::next(rng3)
  let (_, val4) = lib.RNG::next(rng4)

  assert_eq!(val1, val3, "Same seed should produce same first value")
  assert_eq!(val2, val4, "Same seed should produce same second value")
}

test "sphere function at origin" {
  let zeros: Array[Double] = [0.0, 0.0, 0.0]
  let result = lib.sphere_function(zeros)
  assert_eq!(result, 0.0, "Sphere at origin should be 0")
}

test "sphere function at non-origin" {
  let input: Array[Double] = [1.0, 2.0, 3.0]
  let result = lib.sphere_function(input)
  assert_eq!(result, 14.0, "Sphere(1,2,3) = 1+4+9 = 14")
}

test "rastrigin function at origin" {
  let zeros: Array[Double] = [0.0, 0.0, 0.0]
  let result = lib.rastrigin_function(zeros)
  assert_eq!(result, 0.0, "Rastrigin at origin should be 0")
}

test "rastrigin function non-negative" {
  let input: Array[Double] = [1.0, 2.0, 3.0]
  let result = lib.rastrigin_function(input)
  assert!(result >= 0.0, "Rastrigin function should be non-negative")
}

test "rosenbrock function at origin" {
  let input: Array[Double] = [1.0, 1.0, 1.0]
  let result = lib.rosenbrock_function(input)
  assert_eq!(result, 0.0, "Rosenbrock at (1,1,1) should be 0")
}

test "rosenbrock function at non-optimal" {
  let input: Array[Double] = [0.0, 0.0, 0.0]
  let result = lib.rosenbrock_function(input)
  assert!(result > 0.0, "Rosenbrock at (0,0,0) should be > 0")
}

test "array_sum function" {
  let arr: Array[Double] = [1.0, 2.0, 3.0, 4.0]
  let result = lib.array_sum(arr)
  assert_eq!(result, 10.0, "Sum of [1,2,3,4] should be 10")
}

test "array_min function" {
  let arr: Array[Double] = [3.0, 1.0, 4.0, 1.0, 5.0]
  let result = lib.array_min(arr)
  assert_eq!(result, 1.0, "Min of [3,1,4,1,5] should be 1")
}

test "array_max function" {
  let arr: Array[Double] = [3.0, 1.0, 4.0, 1.0, 5.0]
  let result = lib.array_max(arr)
  assert_eq!(result, 5.0, "Max of [3,1,4,1,5] should be 5")
}

test "array_avg function" {
  let arr: Array[Double] = [2.0, 4.0, 6.0]
  let result = lib.array_avg(arr)
  assert_eq!(result, 4.0, "Avg of [2,4,6] should be 4")
}

test "calc_relative_error with zero optimal" {
  let error1 = lib.calc_relative_error(0.0, 0.0)
  assert_eq!(error1, 0.0, "Error should be 0 when both are 0")

  let error2 = lib.calc_relative_error(5.0, 0.0)
  assert_eq!(error2, 5.0, "Error should be absolute diff when optimal is 0")
}

test "calc_relative_error with non-zero optimal" {
  let error = lib.calc_relative_error(110.0, 100.0)
  assert_eq!(error, 0.1, "Error should be 0.1 for 110 vs 100")
}

test "random_solution generates within bounds" {
  let bounds = lib.Bounds::{ min: -5.0, max: 5.0 }
  let rng = lib.RNG::{ seed: 42.0 }
  let (_, solution) = lib.random_solution(5, bounds, rng)

  for i in 0..<solution.length() {
    assert!(solution[i] >= -5.0, "Value should be >= min bound")
    assert!(solution[i] <= 5.0, "Value should be <= max bound")
  }
}

test "random_search improves with more evaluations" {
  let config1 = lib.RandomSearchConfig::{
    max_evaluations: 100,
    dimension: 3,
    bounds: lib.Bounds::{ min: -5.0, max: 5.0 }
  }
  let rng1 = lib.RNG::{ seed: 42.0 }
  let result1 = lib.random_search(config1, lib.sphere_function, rng1)

  let config2 = lib.RandomSearchConfig::{
    max_evaluations: 10000,
    dimension: 3,
    bounds: lib.Bounds::{ min: -5.0, max: 5.0 }
  }
  let rng2 = lib.RNG::{ seed: 42.0 }
  let result2 = lib.random_search(config2, lib.sphere_function, rng2)

  assert!(result2.best_value <= result1.best_value, "More evaluations should not worsen result")
}

test "hill_climbing produces valid result" {
  let config = lib.HillClimbingConfig::{
    max_evaluations: 500,
    dimension: 3,
    bounds: lib.Bounds::{ min: -5.0, max: 5.0 },
    step_size: 0.5,
    restarts: 5
  }
  let rng = lib.RNG::{ seed: 42.0 }
  let result = lib.hill_climbing(config, lib.sphere_function, rng)

  assert!(result.evaluations > 0, "Should have some evaluations")
  assert!(result.best_value >= 0.0, "Sphere function value should be non-negative")
}

test "simulated_annealing produces valid result" {
  let config = lib.SimulatedAnnealingConfig::{
    max_evaluations: 500,
    dimension: 3,
    bounds: lib.Bounds::{ min: -5.0, max: 5.0 },
    initial_temperature: 100.0,
    cooling_rate: 0.995,
    step_size: 0.5
  }
  let rng = lib.RNG::{ seed: 42.0 }
  let result = lib.simulated_annealing(config, lib.sphere_function, rng)

  assert!(result.evaluations > 0, "Should have some evaluations")
  assert!(result.best_value >= 0.0, "Sphere function value should be non-negative")
}

test "evaluate function works" {
  let f = lib.sphere_function
  let solution: Array[Double] = [1.0, 2.0]
  let result = lib.evaluate(f, solution)
  assert_eq!(result, 5.0, "evaluate(sphere, [1,2]) should be 5")
}