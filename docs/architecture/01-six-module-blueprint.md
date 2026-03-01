# Gelex 七模块蓝图（最终版）

## 1. 目标

仓库固定为七模块架构，并由 CMake target 强制边界：

- `infra`
- `types`
- `io`
- `data`
- `model`
- `algo`
- `pipeline`

`include/gelex/**` 顶层目录固定为：
`infra`、`types`、`io`、`data`、`model`、`algo`、`pipeline`。

允许保留单文件入口：`include/gelex/exception.h`。

## 2. 依赖方向与白名单

总体方向（单向、无环）：

`infra <- types <- io <- data <- model <- algo <- pipeline`

允许跨层向下依赖，但必须命中白名单：

- `infra`: 无
- `types`: `infra`
- `io`: `infra`, `types`
- `data`: `infra`, `types`, `io`
- `model`: `infra`, `types`, `data`
- `algo`: `infra`, `types`, `data`, `model`
- `pipeline`: `infra`, `types`, `io`, `data`, `model`, `algo`

上层（`apps`/`tests`/`benchmark`）只允许通过 `pipeline` 或明确公开 API 访问。

## 3. 模块职责与归属

### 3.1 infra

职责：日志、计时、通用统计、格式化等基础设施能力，不包含领域业务。

允许依赖：STL、Eigen、`gelex/exception.h`。

禁止依赖：`types/io/data/model/algo/pipeline`。

目录归属：

- `include/gelex/infra/**`
- `src/infra/**`

### 3.2 types

职责：值对象、配置结构、结果 DTO、异常、枚举、轻量策略声明。

允许依赖：STL、Eigen、`gelex/exception.h`。

禁止依赖：`io/data/model/algo/pipeline`。

目录归属：

- `include/gelex/types/**`
- `src/types/**`
- `include/gelex/exception.h`

### 3.3 io

职责：文本/二进制读写、格式解析、文件序列化，不包含领域业务判断。

允许依赖：`infra`、`types`。

禁止依赖：`data/model/algo/pipeline`。

目录归属：

- `include/gelex/io/**`
- `src/io/**`

### 3.4 data

职责：DataFrame、基因型矩阵、GRM、样本/位点对齐与数据准备。

允许依赖：`infra`、`types`、`io`。

禁止依赖：`model/algo/pipeline`。

目录归属：

- `include/gelex/data/**`
- `src/data/**`

### 3.5 model

职责：统计模型状态与参数表达（频率学/贝叶斯），不做 I/O 与流程编排。

允许依赖：`infra`、`types`、`data`。

禁止依赖：`pipeline`。

目录归属：

- `include/gelex/model/**`
- `src/model/**`

### 3.6 algo

职责：估计器、采样器、优化器、统计推断（REML/MCMC/GWAS 计算）。

允许依赖：`infra`、`types`、`data`、`model`。

禁止依赖：`pipeline`。

目录归属：

- `include/gelex/algo/**`
- `src/algo/**`

### 3.7 pipeline

职责：用例编排（`load -> build -> fit -> write`）、预测与模拟流程组织。

允许依赖：`infra`、`types`、`io`、`data`、`model`、`algo`。

目录归属：

- `include/gelex/pipeline/**`
- `src/pipeline/**`

## 4. 结构与 include 硬约束

1. 禁止新增目录：
   `include/gelex/{estimator,optim,gwas,predict,sim,logger,utils,detail}/**`。
2. 禁止新增 `#include "gelex/logger.h"`，统一使用 `gelex/infra/logger.h`。
3. `apps/**`、`tests/**`、`benchmark/**` 禁止 include `src/**` 私有头。
4. 公共头禁止直接 include `include/gelex/internal/**`。
5. 所有模块必须满足白名单依赖，禁止反向依赖与循环依赖。

## 5. CMake Target 定义

每个模块一个 target：

- `gelex_infra`
- `gelex_types`
- `gelex_io`
- `gelex_data`
- `gelex_model`
- `gelex_algo`
- `gelex_pipeline`

`target_link_libraries` 必须按白名单声明；通过 `PUBLIC/PRIVATE` include
范围阻断越层 include。

## 6. 完成标准（Definition of Done）

1. `include/gelex/**` 顶层目录仅有七模块目录与必要单文件入口（`exception.h`）。
2. `apps/tests/benchmark` 不 include `src/**` 私有头。
3. 公共头不 include `include/gelex/internal/**`。
4. 七个 target 依赖图满足单向且无环。
5. 验证命令全部通过：
   - `python tools/arch/check_module_boundaries.py --format text`
   - `pixi run build-debug`
   - `pixi run test`
   - `pre-commit run -a`
