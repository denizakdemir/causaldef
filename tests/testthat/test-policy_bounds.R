# =============================================================================
# Tests for policy_regret_bound()
# =============================================================================

test_that("policy_regret_bound returns correct class", {
  result <- policy_regret_bound(
    deficiency = 0.1,
    utility_range = c(0, 1)
  )
  
  expect_s3_class(result, "policy_bound")
  expect_true("safety_floor" %in% names(result))
  expect_true("delta" %in% names(result))
})

test_that("policy_regret_bound accepts deficiency object", {
  spec <- .make_test_spec(n = 100)
  def_result <- estimate_deficiency(spec, methods = "iptw", n_boot = 10)
  
  bound <- policy_regret_bound(
    deficiency = def_result,
    utility_range = c(0, 1)
  )
  
  expect_s3_class(bound, "policy_bound")
  expect_equal(bound$delta, min(def_result$estimates))
})

test_that("policy_regret_bound can use upper CI for delta", {
  spec <- .make_test_spec(n = 100)
  def_result <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 20, verbose = FALSE)
  
  bound_upper <- policy_regret_bound(
    deficiency = def_result,
    utility_range = c(0, 1),
    delta_mode = "upper"
  )
  
  expect_s3_class(bound_upper, "policy_bound")
  expect_true(bound_upper$delta >= min(def_result$estimates))
  expect_equal(bound_upper$delta, min(def_result$ci[, 2]))
})

test_that("policy regret bounds match manuscript", {
  # Transfer penalty (upper bound additive term): M * delta
  # Minimax floor (lower bound): (M/2) * delta
  delta <- 0.1
  M <- 10  # utility range diff
  
  bound <- policy_regret_bound(
    deficiency = delta,
    utility_range = c(0, M)
  )
  
  expected_floor <- M * delta
  expect_equal(bound$safety_floor, expected_floor)
  expect_equal(bound$transfer_penalty, expected_floor)
  expect_equal(bound$minimax_floor, 0.5 * M * delta)
})

test_that("regret bound includes obs_regret", {
  delta <- 0.1
  obs_regret <- 0.05
  M <- 1
  
  bound <- policy_regret_bound(
    deficiency = delta,
    utility_range = c(0, M),
    obs_regret = obs_regret
  )
  
  # regret_bound = obs_regret + M*delta
  expected_bound <- obs_regret + M * delta
  expect_equal(bound$regret_bound, expected_bound)
})

test_that("policy_regret_bound validates inputs", {
  # Invalid delta
  expect_error(
    policy_regret_bound(deficiency = -0.1, utility_range = c(0, 1))
  )
  
  # Invalid utility range
  expect_error(
    policy_regret_bound(deficiency = 0.1, utility_range = c(1, 0))
  )
})

test_that("print.policy_bound works", {
  bound <- policy_regret_bound(
    deficiency = 0.1,
    utility_range = c(0, 1),
    obs_regret = 0.05
  )
  
  expect_output(print(bound), "Policy Regret")
  expect_output(print(bound), "Transfer penalty")
})
