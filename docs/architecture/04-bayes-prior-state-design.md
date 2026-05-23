# Bayes PriorState Architecture

> `GeneticPriorState` 是 `GeneticPrior` 的可变运行态对应物。
> 它拥有 prior-local chain state，也拥有 sample / checkpoint 的 record schema。

- Module: `gelex::bayes`
- Status: Beta
- 当前真源：`prior_state.h`、`state_capabilities.h`、`gaussian_prior_state.h`

---

## 目标

本文记录 Bayes prior architecture 的 runtime-state 侧。它补充
`02-bayes-prior-design.md`：02 停在 immutable prior contract，本文从
`GeneticPrior` 创建 mutable state 的时刻开始。

State 层按数据职责组织，不按 method name 组织。一个 BayesR 风格的 state 并不是因为
recipe 名字叫 BayesR 才特殊，而是因为它同时拥有 variance、proportion、component
state。

---

## 数据流

```mermaid
flowchart LR
    Prior["GeneticPrior<br/>immutable block"] -->|make_state(num_markers, num_individuals)| State["GeneticPriorState"]

    Prior --> PCaps["immutable capabilities<br/>variance specs / proportion specs / multipliers"]
    PCaps --> State

    State --> SCaps["mutable capabilities<br/>variance / proportion / component"]
    SCaps --> Kernel["sampler binding"]
    SCaps --> SampleRecords["sample records"]
    SCaps --> CheckpointRecords["checkpoint records"]

    Kernel --> HotLoop["marker / iteration hot path"]
```

动态发现应该发生在 setup 阶段：

1. prior 创建 shape 正确的 concrete state object。
2. kernel binder 查询一次所需 state capabilities。
3. sampler 在 hot loop 中使用 typed references 或 spans。
4. writer / checkpoint code 消费 state 暴露的 records。

marker-level hot path 不应该反复发现一个 state 是否有 variance、proportion、
component slot。

---

## State Root

`GeneticPriorState` 是 mutable identity root。它有三个职责：

- 提供 dynamic capability query 与 requirement checks；
- 暴露 sample records；
- 暴露 checkpoint records，包括 read-only visitor 与 mutable visitor。

record visitation 放在 state 上是有意设计。state 拥有可变 shape，因此它最清楚哪些
records 存在，以及哪些 records 是 sample-only、哪些必须进入 checkpoint。把这套逻辑
移到外部，会产生另一套必须同步每个 concrete state shape 的多态映射。

因此 record visitor function 属于 state contract，而不是独立 writer concern。

---

## Mutable Capabilities

State capabilities 对应下游实际需要的 mutable data。

| Capability | Mutable data | 说明 |
|---|---|---|
| `VarianceStateCap` | 一个或多个 `Eigen::VectorXd` variance slots | slot 数量跟随 prior modes/specs |
| `ProportionStateCap` | assignment、count、value、update policy | assignment 与 value 必须同属一个状态组件 |
| `ComponentStateCap` | component GEBV vectors 与 component variance | mixture-style sampler 消费的 component slot |

没有 mutable multiplier capability。Multiplier 是 immutable prior data，属于
`MultiplierSpecCap`。

`Component` 是下游 sampler 消费的领域词。state 暴露的是 generic component slot，
不把某个 concrete prior family 的名字带进 capability。

---

## State Components

`ProportionState` 组合必须一致更新的字段：

- `assignment`：marker-to-class allocation；
- `count`：由 assignment 派生的 class counts；
- `value`：当前 class probabilities；
- `update`：`value` 是 fixed 还是 sampled。

这些字段是一个状态组件。把 assignment 与 value 拆成不同 capability，会允许 consumer
观察或更新 half state。

`ComponentState` 组合 per-component runtime work buffers：

- `gebv`：每个 active component 一个 vector；
- `gebv_var`：每个 component 的 variance summary。

component buffer 数量在 state construction 时由 prior shape 推导。Scaled mixture prior
中 zero multiplier class 不需要 active component buffer。Joint mixture prior 当前使用
一个 component slot 表达 active joint contribution。

---

## Concrete State Shapes

Concrete state 是 mutable capabilities 的小组合：

| State | Capabilities | Runtime meaning |
|---|---|---|
| `GaussianState` | variance | Gaussian marker variance only |
| `SpikeSlabGaussianState` | variance + proportion | 2-class allocation 与 optional sampled proportion |
| `ScaledMixtureGaussianState` | variance + proportion + component | multi-class scaled component state |
| `JointMixtureGaussianState` | variance + proportion + component | joint A/D allocation state |

Concrete class 拥有 storage。Capability interface 只暴露 spans 与 record helper 行为。
这样数据所有权保持局部，不需要一个塞满 optional fields 的中心 state object。

---

## Record Schema

Record exposure 按 capability 组合：

| Capability | Sample records | Checkpoint records |
|---|---|---|
| variance | variance value | variance value |
| proportion | assignment；只有 sampled update 才输出 value | assignment、count、value |
| component | component variance summary | component variance summary 与 GEBV buffers |

这个区分是有意的：

- sample records 面向 posterior output；
- checkpoint records 必须足够恢复 chain；
- fixed proportion value 不需要作为 sampled posterior value 输出；
- checkpoint data 不从 immutable prior 反推缺失的 mutable pieces。

Concrete state 组合自己实现的 capability record visitors。这样 record presence 与
sampler 绑定看到的 capability surface 保持一致。

---

## Prior 到 State 的构造

State construction 是 immutable spec 到 mutable runtime shape 的桥。

prior 提供：

- marker variance specs；
- proportion specs 与 update policy；
- relevant fixed multipliers；
- covered genetic modes。

state constructor 接收：

- `num_markers`，用于 marker-length vectors 与 assignments；
- `num_individuals`，用于 component GEBV buffers；
- 初始化 mutable values 所需的 prior specs。

初始化规则：

- variance specs 按 `MarkerVarianceScope` 展开成 concrete variance vectors；
- proportion state 初始时所有 marker 分配到 class 0，`value` 等于 prior initial value；
- component buffers zero-initialized，并按 active component count 定长。

这些都是 runtime facts，不属于 `BayesRecipeConfig` 或 `BayesPrior`。

---

## 下游消费

Sampler 应绑定最小必要 surface：

- Gaussian variance update 需要 `VarianceStateCap`；
- spike-slab allocation 需要 `VarianceStateCap` 与 `ProportionStateCap`；
- scaled mixture update 还需要 `ComponentStateCap`，以及 prior 侧 immutable multiplier；
- joint mixture update 需要 joint prior layout 与同样的 mutable state capabilities。

Writer / checkpoint code 应消费 record visitors，而不是重新实现 capability-specific
schema logic。这样文件 schema 决策与 state shape 保持在同一个 owner 下。

---

## 边界规则

State 层遵循这些规则：

- state object 内不按 recipe name 分支。
- 除非出现真正的跨格式需求，否则不把 state shape 复制到外部 writer registry。
- assignment / count / value 不拆到不同 owner。
- multiplier 不变成 mutable state。
- concrete state 能组合 capability 时，不引入 monolithic optional-field state。
- marker-level hot loop 中不做 capability lookup。

这套设计接受 state root 上的局部多态，因为 state 本来就是 mutable shape 的 owner。
它避免增加第二套只用于重新发现相同 shape 的多态层。

---

## 迁移说明

当前代码仍有 legacy `InferenceState` 与 engine paths。目标 runtime 方向是显式区分
model state 与 prior state：

- model/effect state 拥有 coefficients 以及 residual/random runtime values；
- `GeneticPriorState` 拥有 prior-local mutable values；
- sampler setup 通过 capability 绑定 `BayesModel`、`BayesPrior`、model state、
  prior state；
- legacy `BayesMethod` 不继续承载新的 runtime state。

重构旧 `BayesState` / `*Effect` 路径时保持这个分离：effect design data 留在
`BayesModel`，coefficient state 留在 model runtime state，prior-local mutable state
留在 `GeneticPriorState` 后面。
