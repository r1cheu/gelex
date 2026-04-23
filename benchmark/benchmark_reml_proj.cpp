/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */

// Benchmark: REML projection P = V^{-1} - ViX * C_inv * ViX'
//
// Three strategies compared:
//   A. materialized  — form P explicitly, use P * X, (P*K).trace()
//   B. lazy          — keep V, ViX, C_inv factored; apply P on demand
//   C. hybrid (what gelex uses now) — lazy during REML loop, materialize
//      once into V's buffer after convergence, then dense GEMM for scan
//
// Reports both FLOP-dominated timing and matrix memory footprint.

#include <cstdio>
#include <string>
#include <vector>

#include <nanobench.h>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>

namespace
{

constexpr Eigen::Index N = 4000;     // individuals
constexpr Eigen::Index P = 5;        // fixed-effect columns
constexpr Eigen::Index NCOMP = 3;    // variance components (residual + 2 K)
constexpr Eigen::Index CHUNK = 512;  // GWAS chunk size

auto make_spd(Eigen::Index n) -> Eigen::MatrixXd
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
    A = A * A.transpose();
    A.diagonal().array() += static_cast<double>(n);
    return A;
}

void report_memory()
{
    const auto bytes = [](Eigen::Index r, Eigen::Index c)
    { return static_cast<double>(r) * static_cast<double>(c) * 8.0 / 1e6; };

    std::printf(
        "\n=== Memory footprint (MB) — n=%ld, p=%ld ===\n",
        static_cast<long>(N),
        static_cast<long>(P));
    std::printf("  V^{-1}       (n x n): %7.2f MB\n", bytes(N, N));
    std::printf(
        "  proj = P     (n x n): %7.2f MB  [materialized only]\n", bytes(N, N));
    std::printf("  ViX          (n x p): %7.2f MB\n", bytes(N, P));
    std::printf("  C_inv        (p x p): %7.2f MB\n", bytes(P, P));
    std::printf(
        "\n  Strategy A (materialized):   V + P          = %.2f MB\n",
        2.0 * bytes(N, N));
    std::printf(
        "  Strategy B (pure lazy):      V + ViX + C_inv = %.2f MB\n",
        bytes(N, N) + bytes(N, P) + bytes(P, P));
    std::printf(
        "  Strategy C (hybrid, in-place after fit): V (becomes P) = %.2f MB\n",
        bytes(N, N));
    std::printf("\n");
}

}  // namespace

TEST_CASE("REML projection strategies", "[!benchmark][math][reml]")
{
    // Build representative state.
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(N, P);
    Eigen::MatrixXd V_inv = make_spd(N).inverse();
    Eigen::VectorXd y = Eigen::VectorXd::Random(N);
    std::vector<Eigen::MatrixXd> Ks;
    Ks.push_back(make_spd(N));  // K1
    Ks.push_back(make_spd(N));  // K2

    // Factored pieces (shared by all strategies).
    Eigen::MatrixXd ViX = V_inv * X;
    Eigen::MatrixXd XtViX = X.transpose() * ViX;
    Eigen::MatrixXd C_inv = XtViX.llt().solve(Eigen::MatrixXd::Identity(P, P));
    Eigen::VectorXd Py = V_inv * y - ViX * (C_inv * (ViX.transpose() * y));

    // Materialized P (for strategy A / C scan).
    Eigen::MatrixXd P_mat = V_inv;
    P_mat.noalias() -= ViX * C_inv * ViX.transpose();

    // dvpy for AI Hessian (one column per variance component).
    Eigen::MatrixXd dvpy(N, NCOMP);
    dvpy.col(0) = Py;
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(Ks.size()); ++i)
    {
        dvpy.col(i + 1).noalias() = Ks[i] * Py;
    }

    // SNP chunk.
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(N, CHUNK);

    report_memory();

    ankerl::nanobench::Bench loop;
    loop.title("REML loop hot path (one AI iteration)")
        .unit("iter")
        .warmup(3)
        .minEpochIterations(10);

    // --- Strategy A: materialized P ---
    loop.run(
        "A materialized: form P + tr(P) + tr(P*K)*2 + P*dvpy",
        [&]()
        {
            Eigen::MatrixXd Pm = V_inv;
            Pm.noalias() -= ViX * C_inv * ViX.transpose();
            double tr = Pm.trace();
            for (const auto& K : Ks)
            {
                tr += (Pm * K).trace();
            }
            Eigen::MatrixXd pd = Pm * dvpy;
            ankerl::nanobench::doNotOptimizeAway(tr);
            ankerl::nanobench::doNotOptimizeAway(pd);
        });

    // --- Strategy B: lazy P (Frobenius identities) ---
    loop.run(
        "B lazy:         trace_proj() + trace_proj_k()*2 + P*dvpy (factored)",
        [&]()
        {
            // tr(P) = tr(V^{-1}) - tr(C_inv * ViX' * ViX)
            Eigen::MatrixXd wtw = ViX.transpose() * ViX;
            double tr = V_inv.trace() - C_inv.cwiseProduct(wtw).sum();

            // tr(P * K_i) = tr(V^{-1} * K_i) - tr(C_inv * ViX' * K_i * ViX)
            for (const auto& K : Ks)
            {
                Eigen::MatrixXd wtkw = ViX.transpose() * K * ViX;
                tr += V_inv.cwiseProduct(K).sum()
                      - C_inv.cwiseProduct(wtkw).sum();
            }

            // P * dvpy = V^{-1} * dvpy - ViX * C_inv * (ViX' * dvpy)
            Eigen::MatrixXd pd = V_inv * dvpy;
            Eigen::MatrixXd w_dvpy = ViX.transpose() * dvpy;
            pd.noalias() -= ViX * (C_inv * w_dvpy);

            ankerl::nanobench::doNotOptimizeAway(tr);
            ankerl::nanobench::doNotOptimizeAway(pd);
        });

    ankerl::nanobench::Bench scan;
    scan.title("GWAS scan per chunk (m=" + std::to_string(CHUNK) + ")")
        .unit("chunk")
        .warmup(3)
        .minEpochIterations(20);

    // --- Scan with materialized P (current gelex hybrid post-fit) ---
    scan.run(
        "C materialized scan: W = P * Z (single GEMM)",
        [&]()
        {
            Eigen::MatrixXd W(N, CHUNK);
            W.noalias() = P_mat * Z;
            ankerl::nanobench::doNotOptimizeAway(W);
        });

    // --- Scan with lazy apply ---
    scan.run(
        "B lazy scan:         W = V*Z - ViX*C_inv*(ViX'*Z)",
        [&]()
        {
            Eigen::MatrixXd W = V_inv * Z;
            W.noalias() -= ViX * (C_inv * (ViX.transpose() * Z));
            ankerl::nanobench::doNotOptimizeAway(W);
        });
}
