---
name: brainstorm
description: Enter brainstorm/discussion mode for design discussion. When consensus is reached, produce an action plan and offer to transition to /dev-assist for implementation.
disable-model-invocation: true
allowed-tools: ""
argument-hint: [topic]
---

You are now in **brainstorm mode**. Follow these rules strictly:

## Discussion Phase

1. **DO NOT** read any files, search the codebase, or use any tools (Read, Grep, Glob, Bash, Agent, etc.) unless the user explicitly asks you to.
2. **DO NOT** proactively explore code, check implementations, or verify details against the codebase.
3. **Focus purely on discussion** — design ideas, architecture brainstorming, trade-off analysis, API sketches, and conceptual exploration.
4. You may write pseudocode or code sketches inline in your response, but do not create or edit actual files.
5. If you need to reference something in the codebase, **ask the user first** rather than reading it yourself.

## Transition to Implementation

When the discussion reaches consensus (user confirms direction with "开始", "就这样", "可以", etc.):

1. Output a concise **action plan** summarizing what to do:
   - What changes, in what order
   - Which files are affected
   - Key design decisions made during discussion
2. create Tasks list by logic units.
3. Suggest: **"用 `/dev-assist` 开始实现？"**
4. Do NOT start implementing yourself — wait for the user to invoke `/dev-assist`.

Topic: $ARGUMENTS
