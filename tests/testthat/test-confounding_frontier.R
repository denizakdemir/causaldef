# =============================================================================
# Tests for confounding_frontier()
# =============================================================================

test_that("confounding_frontier returns correct class", {
  frontier <- confounding_frontier(
    alpha_range = c(-1, 1),
    gamma_range = c(-1, 1),
    grid_size = 10
  )
  
  expect_s3_class(frontier, "confounding_frontier")
  expect_true("grid" %in% names(frontier))
  expect_true("params" %in% names(frontier))
})

test_that("confounding_frontier creates correct grid", {
  frontier <- confounding_frontier(
    alpha_range = c(-1, 1),
    gamma_range = c(-1, 1),
    grid_size = 10
  )
  
  expect_equal(nrow(frontier$grid), 100)  # 10 x 10
  expect_true("alpha" %in% names(frontier$grid))
  expect_true("gamma" %in% names(frontier$grid))
  expect_true("delta" %in% names(frontier$grid))
})

test_that("confounding_frontier follows confounding lower bound behavior", {
  # Confounding lower bound: delta increases with |alpha*gamma| and is ~0 near alpha=0 or gamma=0
  # Use odd grid size to include alpha=0 and gamma=0 exactly
  frontier <- confounding_frontier(
    alpha_range = c(-2, 2),
    gamma_range = c(-2, 2),
    grid_size = 21  # Odd size ensures 0 is included
  )
  
  # At alpha=0 or gamma=0, delta should be ~0
  # Use wider tolerance since grid may not hit exact zero
  zero_alpha <- frontier$grid[abs(frontier$grid$alpha) < 0.2, ]
  zero_gamma <- frontier$grid[abs(frontier$grid$gamma) < 0.2, ]
  
  # At exact zero, delta should be zero
  expect_true(mean(zero_alpha$delta) < 0.1)
  expect_true(mean(zero_gamma$delta) < 0.1)
  
  # At high alpha AND gamma, delta should be higher
  high_both <- frontier$grid[
    abs(frontier$grid$alpha) > 1.5 & abs(frontier$grid$gamma) > 1.5, 
  ]
  
  expect_true(mean(high_both$delta) > 0.1)
})

test_that("confounding_frontier accepts spec for parameter estimation", {
  df <- .simulate_confounded(200, seed = 1)
  spec <- causal_spec(df, "A", "Y")
  
  frontier <- confounding_frontier(
    spec = spec,
    alpha_range = c(-1, 1),
    gamma_range = c(-1, 1),
    grid_size = 10
  )
  
  expect_s3_class(frontier, "confounding_frontier")
  expect_true(!is.null(frontier$params$sigma_A))
  expect_true(!is.null(frontier$params$sigma_Y))
})

test_that("confounding_frontier delta is symmetric in sign", {
  frontier <- confounding_frontier(
    alpha_range = c(-2, 2),
    gamma_range = c(-2, 2),
    grid_size = 21
  )
  
  # delta(alpha, gamma) should equal delta(-alpha, -gamma)
  pos <- frontier$grid[frontier$grid$alpha > 0 & frontier$grid$gamma > 0, ]
  neg <- frontier$grid[frontier$grid$alpha < 0 & frontier$grid$gamma < 0, ]
  
  # Compare mean deltas (should be similar)
  expect_equal(mean(pos$delta), mean(neg$delta), tolerance = 0.01)
})

test_that("print.confounding_frontier works", {
  frontier <- confounding_frontier(grid_size = 10)
  
  expect_output(print(frontier), "Confounding Frontier")
})
