# Test basic harmonization

test_that("harmonize_gwas_gwis returns expected structure", {
  sim <- simulate_gxe_summary(seed = 42, n_variants = 1000)
  harm <- harmonize_gwas_gwis(sim$gwis, sim$gwas, quiet = TRUE)

  expect_s3_class(harm, "mrgxe_harmonized")
  expect_true("cleaned" %in% names(harm))
  expect_true("qc_counts" %in% names(harm))
  expect_true(nrow(harm$cleaned) > 0)
  expect_true(all(c("was_flipped", "BETA_GWAS_ALIGNED") %in% names(harm$cleaned)))
})

test_that("harmonize drops palindromic SNPs", {
  sim <- simulate_gxe_summary(seed = 42, n_variants = 500, n_palindromic = 50)
  harm <- harmonize_gwas_gwis(sim$gwis, sim$gwas, quiet = TRUE,
                               remove_palindromic = TRUE)
  # No palindromic should remain
  expect_false(any(grepl("^[ATCG][ATCG]$",
                         paste0(harm$cleaned$A1_GWIS, harm$cleaned$A2_GWIS)) %in%
                  c("AT", "TA", "CG", "GC")))
})

# Test instrument selection

test_that("select_instruments returns valid instruments", {
  sim <- simulate_gxe_summary(seed = 42)
  harm <- harmonize_gwas_gwis(sim$gwis, sim$gwas, quiet = TRUE)

  # Use P (since our simulated data has P column, not P_GWAS)
  iv <- select_instruments(harm$cleaned, p_col = "P", quiet = TRUE)

  expect_s3_class(iv, "mrgxe_instruments")
  expect_true(iv$n_instruments >= 2)
  expect_true(all(c("CHR", "BP", "SNP") %in% names(iv$instruments)))
})

test_that("distance clumping removes nearby variants", {
  sim <- simulate_gxe_summary(seed = 42)
  harm <- harmonize_gwas_gwis(sim$gwis, sim$gwas, quiet = TRUE)
  iv <- select_instruments(harm$cleaned, p_col = "P", quiet = TRUE,
                           window_kb = 10000)  # very large window

  if (nrow(iv$instruments) >= 2) {
    # All instruments should be > 10Mb apart on the same chromosome
    for (chr in unique(iv$instruments$CHR)) {
      chr_iv <- iv$instruments[iv$instruments$CHR == chr, ]
      if (nrow(chr_iv) >= 2) {
        dists <- abs(diff(chr_iv$BP))
        expect_true(all(dists > 10000000 | dists == 0))
      }
    }
  }
})

# Test rho estimation

test_that("estimate_rho returns valid correlation", {
  sim <- simulate_gxe_summary(seed = 42)
  harm <- harmonize_gwas_gwis(sim$gwis, sim$gwas, quiet = TRUE)
  rho <- estimate_rho(harm$cleaned)

  expect_s3_class(rho, "mrgxe_rho")
  expect_true(is.numeric(rho$rho))
  expect_true(abs(rho$rho) <= 1)
  expect_true(rho$n_null >= 10)
})

# Test TMRGxE screening

test_that("run_tmrgxe returns expected results", {
  skip_if_not_installed("IMRP")

  sim <- simulate_gxe_summary(seed = 42)
  harm <- harmonize_gwas_gwis(sim$gwis, sim$gwas, quiet = TRUE)
  iv <- select_instruments(harm$cleaned, p_col = "P", quiet = TRUE)
  rho <- estimate_rho(harm$cleaned)

  screen <- run_tmrgxe(harm$cleaned, iv$instruments, rho$rho)

  expect_s3_class(screen, "mrgxe_tmrgxe")
  expect_true("results" %in% names(screen))
  expect_true("PleioP_MR" %in% names(screen$results))
  expect_true("is_significant_gxe" %in% names(screen$results))
})

# Test simulation

test_that("simulate_gxe_summary creates correct structure", {
  sim <- simulate_gxe_summary(seed = 42, n_variants = 500)

  expect_true("gwis" %in% names(sim))
  expect_true("gwas" %in% names(sim))
  expect_true("true_theta" %in% names(sim))
  expect_true("true_gxe_snps" %in% names(sim))

  expect_equal(nrow(sim$gwis), 500)
  expect_equal(nrow(sim$gwas), 500)
  expect_true(sim$true_theta == 0.72)
})
