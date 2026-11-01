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
/// Per-dimension bounds for problems with mixed-scale variables
/// (e.g. x_0 in [-100,100] while x_1 in [-1,1]).
pub struct MultiBounds {
  mins : Array[Double]
  maxs : Array[Double]
}

///|
pub fn make_multi_bounds(
  mins : Array[Double],
  maxs : Array[Double],
) -> MultiBounds {
  MultiBounds::{ mins, maxs, }
}

///|
/// Build uniform per-dimension bounds from a single Bounds object.
pub fn MultiBounds::from_bounds(dim : Int32, bounds : Bounds) -> MultiBounds {
  let mins : Array[Double] = []
  let maxs : Array[Double] = []
  for _i in 0..<dim {
    mins.push(bounds.min)
    maxs.push(bounds.max)
  }
  MultiBounds::{ mins, maxs, }
}

///|
pub fn MultiBounds::dimension(self : MultiBounds) -> Int32 {
  min_int(self.mins.length(), self.maxs.length())
}

///|
pub fn MultiBounds::get_min(self : MultiBounds, d : Int32) -> Double {
  self.mins[d]
}

///|
pub fn MultiBounds::get_max(self : MultiBounds, d : Int32) -> Double {
  self.maxs[d]
}

///|
/// Clamp a vector against the per-dimension bounds.
pub fn MultiBounds::clamp(
  self : MultiBounds,
  x : Array[Double],
) -> Array[Double] {
  let n = min_int(self.mins.length(), x.length())
  let out : Array[Double] = []
  for d in 0..<n {
    out.push(max(self.mins[d], min(self.maxs[d], x[d])))
  }
  out
}

///|
/// Random solution within per-dimension bounds.
pub fn random_solution_mb(
  dim : Int32,
  mb : MultiBounds,
  rng : RNG,
) -> Array[Double] {
  let n = min_int(dim, mb.mins.length())
  let solution : Array[Double] = []
  for i in 0..<n {
    let r = RNG::next(rng)
    solution.push(mb.mins[i] + r * (mb.maxs[i] - mb.mins[i]))
  }
  solution
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
