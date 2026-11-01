// Black-Box Optimization Framework - 核心库

// ===== 数学常量 =====

///|
pub const PI : Double = 3.14159265358979323846

///|
pub const TWO_PI : Double = 6.28318530717958647692

// ===== 工具函数 =====

///|
pub struct RNG {
  mut seed : UInt64
}

///|
pub fn RNG::next(self : RNG) -> Double {
  let mut s = self.seed
  s = s ^ (s >> 12)
  s = s ^ (s << 25)
  s = s ^ (s >> 27)
  let result = s * 0x2545F4914F6CDD1DUL
  self.seed = s
  (result >> 11).to_double() / 9007199254740992.0
}

///|
pub fn make_rng(seed : Double) -> RNG {
  let mut z = seed.to_uint64() + 0x9E3779B97F4A7C15UL
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL
  RNG::{ seed: (z ^ (z >> 31)) | 1UL, }
}

///|
pub fn reduce_to_pi(x : Double) -> Double {
  let k = (x / TWO_PI).round()
  x - k * TWO_PI
}

///|
pub fn cos_approx(x : Double) -> Double {
  @math.cos(x)
}

///|
pub fn sin_approx(x : Double) -> Double {
  @math.sin(x)
}

///|
pub fn exp_approx(x : Double) -> Double {
  @math.exp(x)
}

///|
pub fn pow_approx(base : Double, exp : Double) -> Double {
  @math.pow(if base <= 0.0 { 0.0001 } else { base }, exp)
}

///|
pub fn max(a : Double, b : Double) -> Double {
  if a > b {
    a
  } else {
    b
  }
}

///|
pub fn min(a : Double, b : Double) -> Double {
  if a < b {
    a
  } else {
    b
  }
}

///|
fn max_int(a : Int32, b : Int32) -> Int32 {
  if a > b {
    a
  } else {
    b
  }
}

///|
fn min_int(a : Int32, b : Int32) -> Int32 {
  if a < b {
    a
  } else {
    b
  }
}

///|
pub fn abs_val(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

// ===== 核心类型定义 =====

///|
pub struct Bounds {
  min : Double
  max : Double
}

///|
pub fn make_bounds(min : Double, max : Double) -> Bounds {
  Bounds::{ min, max, }
}

///|
pub struct OptimizationResult {
  best_solution : Array[Double]
  best_value : Double
  evaluations : Int32
}

///|
pub struct BenchmarkResult {
  algorithm_name : String
  best_value : Double
  evaluations : Int32
  relative_error : Double
}

// ===== 辅助函数 =====

///|
pub fn random_solution(
  dim : Int32,
  bounds : Bounds,
  rng : RNG,
) -> Array[Double] {
  let solution : Array[Double] = []
  for _i in 0..<dim {
    let r = RNG::next(rng)
    let rand_val : Double = bounds.min + r * (bounds.max - bounds.min)
    solution.push(rand_val)
  }
  solution
}

///|
pub fn evaluate(
  f : (Array[Double]) -> Double,
  solution : Array[Double],
) -> Double {
  f(solution)
}
