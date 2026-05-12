# Drafting · Architecture Spec Template

工程蓝图风格的 HTML 设计文档模板。米色纸张、黑边粗框、网格背景、信号橙重音、EB Garamond 斜体作为强调符号。

## 起步

```bash
cp docs/architecture/templates/drafting.html docs/architecture/NN-your-doc.html
```

打开新文件，替换内容即可。不需要 build 步骤，浏览器直接打开。

## 文件结构

```
docs/architecture/
├── 01-six-module-blueprint.md
├── 02-bayes-prior-design.html      ← 已有文档
├── NN-your-doc.html                ← 复制后的新文档
└── templates/
    ├── drafting.html               ← 模板源文件
    └── README.md
```

## 顶部三块

### 1. Title block（标题块）

```html
<div class="titleblock">
  <div class="main">
    <div class="stamp"><span>Architecture Spec · NN · Module</span></div>
    <h1>主标题 <span class="accent">关键词.</span></h1>
    <p class="lede">引言。一两句话，<strong>粗体</strong>用于核心概念。</p>
  </div>
  <dl class="data">
    <div class="row"><dt>Sheet</dt><dd>NN / Name</dd></div>
    <div class="row"><dt>Module</dt><dd>code::path</dd></div>
    <div class="row"><dt>Stability</dt><dd>Beta <span class="tag">// note</span></dd></div>
    ...
  </dl>
</div>
```

- `.stamp` — 顶部小印章（橙色），写文档分类
- `h1 > .accent` — 主标题里需要变色的关键词，外包 `<span class="accent">`
- `.lede` — EB Garamond 斜体引言
- `.data` — 右栏元数据，自由增减 `.row`

### 2. Sheetbar（图签条）

```html
<div class="sheetbar">
  <div class="lhs"><span>Rev. NN</span> · Architecture Spec · 2026</div>
  <div>scale 1 : 1 · sheet N of M</div>
</div>
```

可删可留。蓝图风格的版本/页码条。

## 章节

```html
<section class="section">
  <div class="gutter"><div class="id">§ 01 / NAME</div></div>
  <div class="body">
    <h2>章节标题 <span class="accent">关键词</span></h2>
    <p class="dek">小标题。EB Garamond 斜体，限制宽度。</p>

    <p>正文段落。<code>inline code</code> 嵌入技术词汇。</p>

    <h3>子标题</h3>      <!-- 前面会自动加橙色短线 -->
    <h4>SECTION LABEL</h4> <!-- 全大写 mono 标签 -->
  </div>
</section>
```

- `.gutter > .id` — 左侧竖排的章节编号（可旋转大写）
- `h2 > .accent` — 章节标题里的强调词
- `.dek` — 章节副标题（斜体引子）

## 组件

### Callout

```html
<div class="callout">           <!-- 默认（橙）-->
<div class="callout deny">      <!-- 禁止 / 约束（橙）-->
<div class="callout note">      <!-- 注解（蓝）-->
<div class="callout req">       <!-- 允许 / 要求（绿）-->
  <div class="label">标签文本</div>
  <p>内容</p>
</div>
```

左上角带方块角点，颜色对应种类。

### Modules grid（4 格模块图）

```html
<div class="modules">
  <div class="m"><div class="tag">A · TAG</div><h5>Name</h5><p>说明</p></div>
  <div class="m focus">...</div>   <!-- .focus 加深底色 + 底部三角指示 -->
  <div class="m">...</div>
  <div class="m">...</div>
</div>
```

四列等宽，移动端自动堆叠。`.focus` 标记当前文档的"主角"。

### 代码块

```html
<pre data-lang="C++"><code><span class="k">auto</span> x = <span class="fn">f</span>();</code></pre>
```

`data-lang` 在右上角显示语言。Token 类：

| class | 颜色 | 用途 |
|---|---|---|
| `.k`  | signal | 关键字 (`auto`, `class`, `if`) |
| `.t`  | cobalt | 类型 (`BayesState`) |
| `.s`  | moss   | 字符串字面量 |
| `.c`  | ink-soft 斜体 | 注释 |
| `.fn` | signal 加粗 | 函数名 |

### Pull quote

```html
<div class="pull">
  <span class="mark">⌐ epigraph</span>
  引文内容，EB Garamond 斜体，跨两条横线。
</div>
```

### 结尾签名

```html
<div class="end">
  <div>Drawing No. <span class="sig">GLX-NN-NAME</span></div>
  <div>Signed · gelex::module</div>
</div>
```

## 配色（CSS 变量）

集中在 `:root`，调色只需改这几个：

```css
--paper:    #ebe6dc;   /* 主背景 */
--ink:      #0a1a26;   /* 主文字、粗边框 */
--signal:   #d65420;   /* 主强调色（橙）★ 改这里换主色 */
--signal-2: #f08a4a;   /* 主强调色浅版 */
--cobalt:   #1a4f7a;   /* note callout + 类型 token */
--moss:     #4a6b3a;   /* req callout + 字符串 token */
```

换主色时若新色相和 `--cobalt` 或 `--moss` 冲撞，要同步改副色，保证 callout 三色可区分。

## 字体

通过 Google Fonts 加载：

- **Familjen Grotesk** — 标题、正文（瑞典字体，几何感强、字距紧凑）
- **EB Garamond** — 斜体强调、`.lede`、`.dek`、`.pull`、`spec h5`
- **Fragment Mono** — 所有 mono 标签、代码、stamps（带打字机抖动）

需要离线时，把 `<link rel="stylesheet">` 替换成自托管字体文件。

## 约定

- 章节编号用 `§ NN`，与文档系列编号（`02-bayes-prior-design.html` 的 `02`）独立
- `<title>` 写"文档主题 · Architecture Spec"
- `data-lang` 写规范语言名：`C++`、`Bash`、`Python`、`TOML`
- 引用代码符号用 `<code>` 而非 `<em>`，斜体保留给概念性强调
- 不要在正文里嵌入 emoji（与现有 Gelex 文档风格一致）
