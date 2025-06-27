# gelexy

[![GitHub issues](https://img.shields.io/github/issues/r1cheu/gelexy?color=green)](https://github.com/r1cheu/gelexy/issues/new)
Gelexy is a Python package used for genomic prediction, leveraging C++ to handle the dense computation part, achieving the best balance between ease of use and performance.

> [!NOTE]
> This project is undergoing rapid development.

## Installation

Currently, Only installation from source is supported.

```bash
git clone https://github.com/r1cheu/gelexy.git
mamba create -f environment.yml

```

## Quick Start

### GBLUP

```Python
import gelexy as gx

mk_m = gx.make_gblup("~/project/gelexy/notebooks/data/f1_train.tsv")
a = gx.make_grm("/home/rlchen/project/gelexy/notebooks/data/f1_train.bed")
d = gx.make_grm(
    "/home/rlchen/project/gelexy/notebooks/data/f1_train.bed", method="dom"
)
model = mk_m.make("Yield_per_plant ~ 1 + g[a] + g[d]", {"a": a, "d": d})

est = gx.Estimator(max_iter=20)
est.fit(model, em_init=False)
```
