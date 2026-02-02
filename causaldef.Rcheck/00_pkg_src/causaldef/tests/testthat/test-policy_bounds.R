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

test_that("safety floor follows Theorem 3.2", {
  # Theorem 3.2: safety_floor = 2 * M * delta
  delta <- 0.1
  M <- 10  # utility range diff
  
  bound <- policy_regret_bound(
    deficiency = delta,
    utility_range = c(0, M)
  )
  
  expected_floor <- 2 * M * delta
  expect_equal(bound$safety_floor, expected_floor)
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
  
  # Theorem 3.2: regret_bound = obs_regret + 2*M*delta
  expected_bound <- obs_regret + 2 * M * delta
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
  expect_output(print(bound), "Safety")
})
