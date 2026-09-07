// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_INFRA_DETAIL_EIGEN_THREAD_GUARD_H_
#define GELEX_INFRA_DETAIL_EIGEN_THREAD_GUARD_H_

#include <Eigen/Core>

namespace gelex::detail
{

class EigenThreadGuard
{
   public:
    EigenThreadGuard() : old_thread_count_(Eigen::nbThreads())
    {
        Eigen::setNbThreads(1);
    }

    ~EigenThreadGuard() { Eigen::setNbThreads(old_thread_count_); }

    EigenThreadGuard(const EigenThreadGuard&) = delete;
    EigenThreadGuard& operator=(const EigenThreadGuard&) = delete;
    EigenThreadGuard(EigenThreadGuard&&) = delete;
    EigenThreadGuard& operator=(EigenThreadGuard&&) = delete;

   private:
    int old_thread_count_;
};

}  // namespace gelex::detail

#endif  // GELEX_INFRA_DETAIL_EIGEN_THREAD_GUARD_H_
