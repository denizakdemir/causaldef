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

test_that("policy_regret_bound warns when selecting the minimum across methods", {
  spec <- .make_test_spec(n = 100)
  def_result <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 10, verbose = FALSE)

  expect_warning(
    bound <- policy_regret_bound(
      deficiency = def_result,
      utility_range = c(0, 1)
    ),
    "optimistic"
  )

  expect_equal(bound$delta_selection, "optimistic_minimum")
})

test_that("policy_regret_bound can use upper CI for delta", {
  spec <- .make_test_spec(n = 100)
  def_result <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 20, verbose = FALSE)
  
  bound_upper <- policy_regret_bound(
    deficiency = def_result,
    utility_range = c(0, 1),
    delta_mode = "upper",
    method = "iptw"
  )
  
  expect_s3_class(bound_upper, "policy_bound")
  expect_true(bound_upper$delta >= def_result$estimates["iptw"])
  expect_equal(bound_upper$delta, unname(def_result$ci["iptw", 2]))
  expect_equal(bound_upper$delta_selection, "explicit_method")
})

test_that("policy regret bounds match manuscript", {
  # Policy Regret Transfer: value error |A| M delta, regret penalty 2 |A| M delta
  # Minimax Safety Floor: (M/4) delta
  delta <- 0.1
  M <- 10  # utility range diff

  bound <- policy_regret_bound(
    deficiency = delta,
    utility_range = c(0, M)
  )

  expect_equal(bound$n_actions, 2L)
  expect_equal(bound$value_error, 2 * M * delta)
  expect_equal(bound$transfer_penalty, 2 * 2 * M * delta)
  expect_equal(bound$safety_floor, bound$transfer_penalty)
  expect_equal(bound$minimax_floor, 0.25 * M * delta)
})

test_that("constant policies drop the action-count factor", {
  bound <- policy_regret_bound(deficiency = 0.1, utility_range = c(0, 1),
                               policy_type = "constant", n_actions = 5)
  expect_equal(bound$value_error, 0.1)
  expect_equal(bound$transfer_penalty, 0.2)
})

test_that("covariate policies scale with the number of actions", {
  b3 <- policy_regret_bound(deficiency = 0.1, utility_range = c(0, 1), n_actions = 3)
  expect_equal(b3$transfer_penalty, 2 * 3 * 0.1)
  expect_error(policy_regret_bound(deficiency = 0.1, n_actions = 1))
})

test_that("factor 2 example from the manuscript is respected", {
  # Two actions, trivial covariates: worst-case regret 2 M d while simulated regret 0.
  d <- 0.1
  b <- policy_regret_bound(deficiency = d, utility_range = c(0, 1),
                           policy_type = "constant", obs_regret = 0)
  expect_gte(b$regret_bound, 2 * d)
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
  
  # regret_bound = obs_regret + 2 * n_actions * M * delta (default: 2 actions)
  expected_bound <- obs_regret + 4 * M * delta
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
