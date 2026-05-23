# Bayes Recipe / Prior Architecture

> `BayesRecipe` 是编译层。
> `BayesPrior` 是下游推断消费的不可变概率合同。

- Module: `gelex::bayes`
- Status: Beta
- 当前真源：`include/gelex/model/bayes/` 与 `src/model/bayes/`

---

## 目标

本文记录 `BayesRecipe`、`BayesPrior`、`GeneticPrior` 的架构边界。它不是
header API 的重复；精确接口以代码为准，本文只说明对象为什么这样拆、数据如何流动、
以及后续扩展应该放在哪一层。

目标路径是：

1. 用户侧 method preset 与 override 进入 `BayesRecipe`。
2. `BayesRecipe::make_prior(model)` 将输入编译成 `BayesPrior`。
3. 推断代码读取 `BayesPrior` 与 `GeneticPrior` capability，而不是读取 recipe 名字。
4. 可变链状态由 `GeneticPrior::make_state(...)` 创建，详见
   `04-bayes-prior-state-design.md`。

当前 legacy MCMC wiring 仍然存在。这部分是迁移路径，不是新架构的真源。

---

## 职责划分

| 层 | 拥有 | 不拥有 |
|---|---|---|
| `BayesModel` | phenotype、fixed/random/genetic design data、模型方差事实 | method default、prior topology、sampler state |
| `BayesRecipeConfig` | method-agnostic user override shape | preset topology、final prior semantics、CLI syntax |
| `BayesRecipe` | preset selection、config validation、到 `BayesPrior` 的 compile boundary | runtime state、kernel selection、checkpoint schema |
| `BayesRecipeImpl` | preset-specific topology、default resolution | public ownership、mutable chain state |
| `BayesPrior` | model-level immutable prior contract：random、genetic blocks、residual | method name、作为逻辑真源的 recipe provenance、mutable state |
| `GeneticPrior` | 单个 genetic prior block 及其 immutable capabilities | coefficients、assignments、sampled values、writer paths |

单一真源规则：

- method topology 属于 concrete recipe implementation；
- final probability facts 属于 `BayesPrior` / `GeneticPrior`；
- mutable inference facts 属于 `GeneticPriorState`；
- CLI spelling 与 parser 细节留在 model layer 之外。

---

## 数据流

```mermaid
flowchart LR
    CLI["CLI / caller"] --> Config["BayesRecipeConfig<br/>override bag"]
    Config --> Recipe["BayesRecipe<br/>preset + config"]
    Model["BayesModel<br/>data + design"] --> Recipe
    Recipe -->|make_prior(model)| Prior["BayesPrior<br/>random + genetics + residual"]
    Prior --> Genetic["GeneticPrior blocks<br/>immutable capabilities"]
    Genetic -->|make_state(shapes)| State["GeneticPriorState<br/>mutable capabilities"]

    Genetic --> Binder["kernel binder / writer binder"]
    State --> Binder
    Binder --> Runtime["MCMC / VI runtime"]
```

关键边界是 `Recipe -> Prior` 的 compile step。跨过这个边界后，下游不应该再问
“原始 method 是 BayesA / BayesB / BayesC / BayesR / BayesCD 吗”，而应该问
“这个 prior 暴露了哪些 immutable capabilities”。

---

## Recipe 层

`BayesRecipe` 组合 `BayesRecipePreset` 与 `BayesRecipeConfig`。

`BayesRecipeConfig` 是 override bag：

- `modes` 声明请求的 genetic modes，并限制为 `{A}`、`{D}` 或 `{A, D}`；
- `additive` 与 `dominance` 保存 per-effect optional overrides，例如
  heritability、proportion、multiplier、proportion update；
- `joint_proportion` 与 `joint_proportion_update` 用于 joint method；
- `random_variance_proportion` 控制 random variance 的默认翻译。

它不表达 method topology。例如一个 preset 使用 per-marker variance、
per-effect variance、independent proportion，还是 joint proportion，这些事实属于
concrete recipe implementation。

Concrete recipe implementation 将用户意图翻译成 final prior facts：

| Preset | Genetic prior shape | Compile-time facts |
|---|---|---|
| `BayesRR` | Gaussian | per-effect marker variance |
| `BayesA` | Gaussian | per-marker marker variance |
| `BayesB` | Spike-slab Gaussian | per-marker variance、2-class proportion |
| `BayesC` | Spike-slab Gaussian | per-effect variance、2-class proportion |
| `BayesR` | Scaled mixture Gaussian | per-effect variance、paired proportion and multiplier |
| `BayesCD` | Joint mixture Gaussian | `{A, D}` joint block、4-class joint proportion |

Default resolution 也属于 recipe 层：

- effect heritability 缺失时使用该 genetic mode 的 recipe default；
- BayesB/C 缺省为 two-class simplex；
- BayesR 的 proportion 与 multiplier 作为一组默认值处理；
- BayesCD 使用 `[both-off, A-only, D-only, both-on]` 的 joint layout；
- marker variance target 在构造 prior 前由 model variance facts 与 active marker
  weight 推导。

这样 `BayesPrior` 不需要知道 method name，但 final prior 仍然是完整显式的概率合同。

---

## Prior 层

`BayesPrior` 是 immutable model-level probability contract，包含：

- 一个 random variance prior；
- 一个有序的 owned `GeneticPrior` block 列表；
- 一个 residual variance prior。

Genetic block 的顺序有语义。它是 state creation、kernel binding、sample records、
checkpoint binding 的下游遍历顺序。consumer 应该保留该顺序，而不是按 method name
或 mode 自行排序。

`GeneticPrior` 是单个 genetic block 的 identity root。它回答两个问题：

- 该 block 覆盖哪些 `GeneticMode`；
- 在给定 marker / individual shape 时如何创建匹配的 mutable prior state。

Concrete prior class 拥有 immutable data 并组合 capability。Capability 是正交的：

| Capability | 暴露的 immutable data | 典型 consumer |
|---|---|---|
| `VarianceSpecCap` | marker variance specs | state creation、variance kernels |
| `ProportionSpecCap` | initial proportion 与 update policy | state creation、mixture kernels |
| `MultiplierSpecCap` | fixed scale multipliers | BayesR-style kernels |

Capability model 避免中心 enum / variant。新增 prior 时，concrete prior 只组合自己
实际拥有的能力；下游在 setup 阶段按 capability 绑定。

---

## 边界规则

这套拆分让 SOLID 压力停留在局部：

- SRP：recipe 编译用户意图；prior 保存不可变概率事实；prior state 保存可变链事实。
- OCP：新增 prior 行为通常只增加 concrete `GeneticPrior`，必要时增加 capability，
  不需要改已有 prior class。
- LSP：下游依赖 `GeneticPrior` + capabilities，而不是把某个 concrete class 当成
  method 的替身。
- ISP：variance、proportion、multiplier 分开，consumer 只 require 自己需要的能力。
- DIP：kernel / writer 依赖 prior/state capability contract，而不是 recipe preset。

实践规则：

- `make_prior(model)` 之后不要再按 `BayesRecipePreset` 分支。
- `BayesPrior` 不保存 assignment、sampled proportion、component GEBV、
  marker variance value 这类 mutable facts。
- CLI ingress check 与 parser syntax 不进入 recipe/prior 类型。
- capability query 已能回答的问题，不再增加一套 method-name catalog。
- concrete prior/state 已自然拥有的数据，不为了统一外形强行搬到中心 wrapper。

---

## 扩展路径

新增 recipe 时按以下顺序推进：

1. 复用或新增能表达 final probability contract 的 `GeneticPrior` shape。
2. 只有当下游需要一种新的 immutable data 时，才新增 capability。
3. 在 recipe compile logic 中把 config 与 model facts 翻译成 final prior。
4. 在匹配的 prior-state 层补 state construction 与 record exposure。
5. kernel 在 setup 阶段绑定 capability，不在 marker-level hot loop 中查询。

新增 prior 而不新增 recipe 也是合法的。`BayesPrior` 是下游合同；`BayesRecipe` 只是
产生它的一种方式。

---

## 迁移说明

当前代码仍包含 legacy `BayesMethod`、legacy prior constants、以及
`InferenceState`-based MCMC/VI state。这些类型不应继续塑造新架构。

迁移期间：

- 新的 recipe tests 应通过检查生成的 `BayesPrior` 来验证 public compile contract；
- 旧 engine call site 可以暂时保留，但新代码不应继续扩大 legacy method surface；
- 重构 `BayesModel`、effect objects、runtime state 时，最终方向应是
  `BayesModel + BayesPrior -> BayesState/PriorState`，而不是
  `BayesMethod -> legacy InferenceState`。
