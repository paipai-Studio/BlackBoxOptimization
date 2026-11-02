///|
fn main {
  let bounds = @lib.make_bounds(-5.0, 5.0)
  let dim : Int32 = 3
  let funcs : Array[(String, (Array[Double]) -> Double)] = [
    ("Sphere", @lib.sphere_function),
    ("Ackley", @lib.ackley_function),
    ("Rosenbrock", @lib.rosenbrock_function),
    ("Levy", @lib.levy_function),
    ("Zakharov", @lib.zakharov_function),
    ("Booth", @lib.booth_function),
    ("Matyas", @lib.matyas_function),
    ("Rastrigin", @lib.rastrigin_function),
  ]
  let seeds = 30
  let algos : Array[
    (
      String,
      (Int32, Int32, @lib.Bounds, Int32, (Array[Double]) -> Double, @lib.RNG) -> @lib.OptimizationResult,
    ),
  ] = [
    (
      "Rosenbrock",
      fn(e, d, b, p, f, r) {
        @lib.rosenbrock_method(@lib.make_rosenbrock_config(e, d, b, p), f, r)
      },
    ),
  ]
  for algo in algos {
    let (aname, run) = algo
    for func in funcs {
      let (fname, f) = func
      let mut sum1 = 0.0
      let mut sum5 = 0.0
      let mut flags = 0
      let mut worst1 = 0.0
      let mut worst5 = 0.0
      for s in 0..<seeds {
        let r1 = run(
          1000,
          dim,
          bounds,
          25,
          f,
          @lib.make_rng((100 + s).to_double()),
        )
        let r5 = run(
          5000,
          dim,
          bounds,
          25,
          f,
          @lib.make_rng((100 + s).to_double()),
        )
        sum1 = sum1 + r1.best_value
        sum5 = sum5 + r5.best_value
        if r5.best_value > r1.best_value {
          flags = flags + 1
        }
        if r1.best_value > worst1 {
          worst1 = r1.best_value
        }
        if r5.best_value > worst5 {
          worst5 = r5.best_value
        }
      }
      let m1 = sum1 / seeds.to_double()
      let m5 = sum5 / seeds.to_double()
      println(
        "\{aname}|\{fname}|mean1k=\{m1}|mean5k=\{m5}|worst1k=\{worst1}|worst5k=\{worst5}|regress_seeds=\{flags}/\{seeds}",
      )
    }
  }
  println("SCAN_DONE")
}
