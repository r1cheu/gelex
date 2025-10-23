# Gelex Bayes 模块重构计划

## 项目概述

**目标**: 简化 BayesModel, BayesState, MCMCSamples, MCMCResult 的设计，提高代码可读性和可维护性

**时间估计**: 2-3 天（包括测试）

**风险等级**: 中等（核心模块重构，但接口变化可控）

---

## 一、当前问题分析

### 1.1 核心问题
- ❌ BayesModel/State 过度使用模板元编程（tuple + optional + concept）
- ❌ MCMCSamples 使用 bool 标志位，与其他类结构不一致
- ❌ MCMCResult 滥用继承（DominantSummary : AdditiveSummary : RandomSummary）
- ❌ 访问接口复杂：`get<T>()` 返回 `expected<reference_wrapper<T>>`
- ❌ 大量重复的条件判断代码

### 1.2 影响范围
**核心文件**（6 个）:
- `include/gelex/model/bayes/model.h` (194 行)
- `src/model/bayes/model.cpp` (176 行)
- `include/gelex/estimator/bayes/samples.h` (145 行)
- `src/estimator/bayes/samples.cpp` (108 行)
- `include/gelex/estimator/bayes/result.h` (144 行)
- `src/estimator/bayes/result.cpp` (115 行)

**受影响文件**（28 个）:
- 采样器：11 个文件 (common, fixed, random, additive/*, dominant, pi)
- MCMC 核心：mcmc.h, params.h
- 输出/日志：result_writer.*, bayes_logger.*
- 其他：prior_manager.*, predictor.*, fit_command.cpp

---

## 二、设计目标

### 2.1 核心原则
1. ✅ **KISS**: 简单直接，易于理解
2. ✅ **一致性**: 所有层次使用相同的设计模式
3. ✅ **组合优于继承**: 避免不必要的继承关系
4. ✅ **明确的接口**: 清晰的命名和返回类型
5. ✅ **零成本抽象**: 避免不必要的性能开销

### 2.2 设计决策

#### ✅ 统一使用直接成员 + optional
```cpp
// 所有类都使用这个模式
class XxxModel {
private:
    std::optional<FixedEffect> fixed_;
    std::vector<RandomEffect> random_;
    std::optional<AdditiveEffect> additive_;
    std::optional<DominantEffect> dominant_;
    ResidualEffect residual_;  // always exists
};
```

#### ✅ 简单明确的访问方法
```cpp
// 返回指针，nullptr 表示不存在
const FixedEffect* fixed_effect() const {
    return fixed_ ? &*fixed_ : nullptr;
}

// bool 检查
bool has_fixed_effect() const {
    return fixed_.has_value();
}
```

#### ✅ 组合而非继承
```cpp
// ❌ 错误
struct DominantSummary : AdditiveSummary : RandomSummary { ... };

// ✅ 正确
struct DominantSummary {
    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary ratios;
};
```

---

## 三、详细重构步骤

### Phase 1: 准备工作（0.5 天）

#### 1.1 备份和分支
```bash
git checkout -b refactor/simplify-bayes-design
git add -A && git commit -m "checkpoint: before refactor"
```

#### 1.2 建立测试基线
```bash
# 运行现有测试，记录输出
pixi run test > /tmp/test_baseline.txt 2>&1

# 确保编译通过
pixi run build
```

#### 1.3 创建临时兼容层（可选）
- 如果需要分阶段迁移，可以先保留旧接口作为 deprecated

---

### Phase 2: 重构 BayesModel 和 BayesState（1 天）

#### 2.1 修改 BayesModel 头文件
**文件**: `include/gelex/model/bayes/model.h`

**变更**:
```cpp
// ❌ 删除
namespace detail {
    template <typename T>
    concept BayesEffectType = ...;
}

// ❌ 删除
using EffectsTuple = std::tuple<...>;

// ✅ 新增简单成员
private:
    std::optional<bayes::FixedEffect> fixed_;
    std::vector<bayes::RandomEffect> random_;
    std::optional<bayes::AdditiveEffect> additive_;
    std::optional<bayes::DominantEffect> dominant_;
    bayes::Residual residual_;

// ✅ 新增简单访问方法
public:
    const bayes::FixedEffect* fixed_effect() const;
    bayes::FixedEffect* fixed_effect();
    bool has_fixed_effect() const;

    const std::vector<bayes::RandomEffect>& random_effects() const;
    bool has_random_effects() const;

    const bayes::AdditiveEffect* additive_effect() const;
    bool has_additive_effect() const;

    const bayes::DominantEffect* dominant_effect() const;
    bool has_dominant_effect() const;

    const bayes::Residual& residual() const;
```

#### 2.2 修改 BayesModel 实现
**文件**: `src/model/bayes/model.cpp`

**变更**:
- 实现新的访问方法
- 修改 `create()` 方法中的初始化逻辑
- 移除 `get<T>()` 模板方法实现

#### 2.3 修改 BayesState
**文件**: `include/gelex/model/bayes/model.h`, `src/model/bayes/model.cpp`

**变更**: 完全对应 BayesModel 的结构

```cpp
class BayesState {
public:
    static auto create(const BayesModel& model) -> std::expected<BayesState, Error>;

    bayes::FixedState* fixed_state();
    const bayes::FixedState* fixed_state() const;

    std::vector<bayes::RandomState>& random_states();
    const std::vector<bayes::RandomState>& random_states() const;

    bayes::AdditiveState* additive_state();
    const bayes::AdditiveState* additive_state() const;

    bayes::DominantState* dominant_state();
    const bayes::DominantState* dominant_state() const;

    bayes::ResidualState& residual_state();
    const bayes::ResidualState& residual_state() const;

    void compute_heritability();

private:
    std::optional<bayes::FixedState> fixed_;
    std::vector<bayes::RandomState> random_;
    std::optional<bayes::AdditiveState> additive_;
    std::optional<bayes::DominantState> dominant_;
    bayes::ResidualState residual_;
};
```

#### 2.4 更新所有采样器
**文件**:
- `src/model/bayes/samplers/common.cpp` (Fixed, Random, Residual)
- `src/model/bayes/samplers/additive/*.cpp` (A, B, C, RR, RRD)
- `src/model/bayes/samplers/dominant.cpp`
- `src/model/bayes/samplers/pi.cpp`

**变更**: 将所有的 `get<T>()` 调用替换为新的访问方法

```cpp
// ❌ 旧代码
auto effect = model.get<bayes::FixedEffect>();
if (!effect) return;
const auto& e = effect->get();

// ✅ 新代码
const auto* effect = model.fixed_effect();
if (!effect) return;
```

#### 2.5 测试点
```bash
# 每个子步骤后编译测试
pixi run build
pixi run test
```

---

### Phase 3: 重构 MCMCSamples（0.5 天）

#### 3.1 修改 Samples 结构体
**文件**: `include/gelex/estimator/bayes/samples.h`

**变更**:
```cpp
// ✅ 简化结构，移除继承
struct EffectSamples {
    Samples coeffs;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct VarianceEffectSamples {
    Samples coeffs;
    Samples variance;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct AdditiveSamples {
    Samples coeffs;
    Samples variance;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct DominantSamples {
    Samples coeffs;
    Samples variance;
    Samples ratios;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct PiSamples {
    Samples prop;
    explicit operator bool() const { return !prop.empty(); }
};
```

#### 3.2 修改 MCMCSamples
**文件**: `include/gelex/estimator/bayes/samples.h`, `src/estimator/bayes/samples.cpp`

**变更**:
```cpp
class MCMCSamples {
public:
    explicit MCMCSamples(const MCMCParams& params, const BayesState& state);

    void store(const BayesState& state, Index record_idx, Index chain_idx);

    // ✅ 返回指针，nullptr 表示不存在
    const EffectSamples* fixed() const { return fixed_.get(); }
    const std::vector<VarianceEffectSamples>& random() const { return random_; }
    const AdditiveSamples* additive() const { return additive_.get(); }
    const DominantSamples* dominant() const { return dominant_.get(); }
    const PiSamples* pi() const { return pi_.get(); }
    const VarianceEffectSamples& residual() const { return residual_; }

private:
    // ❌ 移除所有 store_* bool 标志位

    // ✅ 使用 unique_ptr 表示可选
    std::unique_ptr<EffectSamples> fixed_;
    std::vector<VarianceEffectSamples> random_;
    std::unique_ptr<AdditiveSamples> additive_;
    std::unique_ptr<DominantSamples> dominant_;
    std::unique_ptr<PiSamples> pi_;
    VarianceEffectSamples residual_;  // always exists
};
```

#### 3.3 简化 store() 方法
```cpp
void MCMCSamples::store(
    const BayesState& state,
    Index record_idx,
    Index chain_idx)
{
    // ✅ 直接检查指针，无需额外标志位
    if (fixed_) {
        fixed_->coeffs[chain_idx].col(record_idx)
            = state.fixed_state()->coeffs;
    }

    if (additive_) {
        const auto* s = state.additive_state();
        additive_->coeffs[chain_idx].col(record_idx) = s->coeffs;
        additive_->variance[chain_idx](0, record_idx) = s->variance;
    }

    // ... 其他类似
}
```

---

### Phase 4: 重构 MCMCResult（0.5 天）

#### 4.1 修改 Summary 结构体
**文件**: `include/gelex/estimator/bayes/result.h`

**变更**: 移除继承，使用组合

```cpp
// ❌ 删除继承链
// struct DominantSummary : AdditiveSummary : RandomSummary : FixedSummary

// ✅ 使用组合
struct EffectSummary {
    PosteriorSummary coeffs;
};

struct VarianceEffectSummary {
    PosteriorSummary coeffs;
    PosteriorSummary variance;
};

struct AdditiveSummary {
    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary pve;
};

struct DominantSummary {
    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary ratios;
    PosteriorSummary pve;
};
```

#### 4.2 修改 MCMCResult
**文件**: `include/gelex/estimator/bayes/result.h`, `src/estimator/bayes/result.cpp`

**变更**:
```cpp
class MCMCResult {
public:
    void compute(std::optional<double> prob = std::nullopt);

    // ✅ 返回指针
    const EffectSummary* fixed() const { return fixed_.get(); }
    const std::vector<VarianceEffectSummary>& random() const { return random_; }
    const AdditiveSummary* additive() const { return additive_.get(); }
    const DominantSummary* dominant() const { return dominant_.get(); }
    const PosteriorSummary& residual() const { return residual_; }

private:
    std::unique_ptr<EffectSummary> fixed_;
    std::vector<VarianceEffectSummary> random_;
    std::unique_ptr<AdditiveSummary> additive_;
    std::unique_ptr<DominantSummary> dominant_;
    PosteriorSummary residual_;

    // ... 其他成员保持不变
};
```

---

### Phase 5: 更新依赖文件（0.5 天）

#### 5.1 更新 MCMC 主循环
**文件**: `include/gelex/estimator/bayes/mcmc.h`

**变更**: 更新 `update_indicators()` 中的访问方式

```cpp
// ❌ 旧代码
if (model.has_additive_effect()) {
    auto effect = model.get<bayes::AdditiveEffect>();
    if (effect) {
        const auto& e = effect->get();
        // ...
    }
}

// ✅ 新代码
if (const auto* effect = model.additive_effect()) {
    indicator.update(chain, sigma_squared("_a"), state.additive_state()->variance);
}
```

#### 5.2 更新结果输出
**文件**: `src/estimator/bayes/result_writer.cpp`

#### 5.3 更新日志输出
**文件**: `src/logger/bayes_logger.cpp`

#### 5.4 更新先验管理
**文件**: `src/model/bayes/prior_manager.cpp`

#### 5.5 更新预测器
**文件**: `src/predictor/bayes/predictor.cpp`

#### 5.6 更新命令行入口
**文件**: `src/cli/fit_command.cpp`

---

### Phase 6: 测试和验证（0.5 天）

#### 6.1 编译测试
```bash
pixi run build
```

#### 6.2 运行单元测试
```bash
pixi run test
```

#### 6.3 运行集成测试
```bash
# 使用示例数据运行完整流程
./gelex fit --bfile data/test --pheno data/pheno.tsv \
    --method RR --out test_output

# 对比结果与重构前是否一致
diff test_output.* baseline_output.*
```

#### 6.4 性能测试
```bash
# 记录重构前后的运行时间
time ./gelex fit --bfile large_data --method RR --out perf_test
```

---

## 四、风险管理

### 4.1 潜在风险
1. **接口变更破坏现有代码**
   - 缓解措施：先在分支上完成，充分测试后再合并

2. **性能回退**
   - 缓解措施：使用 benchmark 对比重构前后性能

3. **引入新 bug**
   - 缓解措施：保留原始实现作为参考，对比输出结果

### 4.2 回滚策略
```bash
# 如果出现严重问题，可以回滚
git checkout main
git branch -D refactor/simplify-bayes-design
```

---

## 五、代码审查清单

### 5.1 设计审查
- [ ] 所有类使用一致的设计模式
- [ ] 访问方法命名清晰明确
- [ ] 消除了不必要的模板元编程
- [ ] 移除了继承滥用
- [ ] 没有引入新的复杂性

### 5.2 功能审查
- [ ] 所有采样器正常工作
- [ ] MCMC 采样结果正确
- [ ] 后验统计量计算正确
- [ ] 文件输出格式不变

### 5.3 性能审查
- [ ] 编译时间没有显著增加
- [ ] 运行时性能没有回退
- [ ] 内存使用合理

### 5.4 测试审查
- [ ] 所有现有测试通过
- [ ] 新增针对重构的测试（如有必要）

---

## 六、预期成果

### 6.1 代码改进
- 减少约 200-300 行模板元编程代码
- 访问方法调用简化 50%+
- 消除所有 `store_*` bool 标志位
- 继承层级减少 3 层

### 6.2 可读性提升
```cpp
// ❌ 重构前
auto effect = model.get<bayes::AdditiveEffect>();
if (!effect) return std::unexpected(effect.error());
const auto& e = effect->get();
auto state = states.get<bayes::AdditiveState>();
if (!state) return;
const auto& s = state->get();

// ✅ 重构后
const auto* effect = model.additive_effect();
if (!effect) return;
auto* state = states.additive_state();
if (!state) return;
```

### 6.3 维护性提升
- 新增 effect 类型只需添加成员和访问方法
- 错误信息清晰易懂
- IDE 代码补全友好

---

## 七、时间表

| 阶段 | 任务 | 时间 | 负责人 |
|-----|------|------|--------|
| Phase 1 | 准备工作 | 0.5 天 | - |
| Phase 2 | BayesModel/State | 1 天 | - |
| Phase 3 | MCMCSamples | 0.5 天 | - |
| Phase 4 | MCMCResult | 0.5 天 | - |
| Phase 5 | 依赖文件 | 0.5 天 | - |
| Phase 6 | 测试验证 | 0.5 天 | - |
| **总计** | | **3.5 天** | |

---

## 八、后续优化（可选）

### 8.1 进一步简化
- 考虑将 `BayesModel`, `BayesState`, `MCMCSamples` 合并为一个类族
- 使用 CRTP 减少代码重复

### 8.2 性能优化
- 如果访问频繁，考虑缓存指针
- 使用 `std::span` 代替部分 vector 传递

### 8.3 文档更新
- 更新 API 文档
- 添加设计模式说明
- 提供迁移指南（如果有外部用户）

---

## 九、附录

### A. 关键文件列表

**核心文件（必须修改）**:
```
include/gelex/model/bayes/model.h
src/model/bayes/model.cpp
include/gelex/estimator/bayes/samples.h
src/estimator/bayes/samples.cpp
include/gelex/estimator/bayes/result.h
src/estimator/bayes/result.cpp
```

**采样器文件（必须修改）**:
```
src/model/bayes/samplers/common.{h,cpp}
src/model/bayes/samplers/additive/{a,b,c,rr,rrd}.{h,cpp}
src/model/bayes/samplers/dominant.{h,cpp}
src/model/bayes/samplers/pi.{h,cpp}
```

**其他文件（必须修改）**:
```
include/gelex/estimator/bayes/mcmc.h
src/estimator/bayes/result_writer.cpp
src/logger/bayes_logger.{h,cpp}
src/model/bayes/prior_manager.cpp
src/predictor/bayes/predictor.cpp
src/cli/fit_command.cpp
```

### B. 参考资料
- C++ Core Guidelines: https://isocpp.github.io/CppCoreGuidelines/
- "Effective Modern C++" by Scott Meyers
- "C++ Software Design" by Klaus Iglberger

---

**文档版本**: v1.0
**创建日期**: 2025-10-23
**最后更新**: 2025-10-23
