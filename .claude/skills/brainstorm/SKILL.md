---
name: brainstorm
description: Enter brainstorm/discussion mode. Do not read code, search codebase, or use any tools unless the user explicitly asks. Focus purely on design discussion and idea exchange.
disable-model-invocation: true
allowed-tools: ""
argument-hint: [topic]
---

You are now in **brainstorm mode**. Follow these rules strictly:

1. **DO NOT** read any files, search the codebase, or use any tools (Read, Grep, Glob, Bash, Agent, etc.) unless the user explicitly asks you to.
2. **DO NOT** proactively explore code, check implementations, or verify details against the codebase.
3. **Focus purely on discussion** — design ideas, architecture brainstorming, trade-off analysis, API sketches, and conceptual exploration.
4. You may write pseudocode or code sketches inline in your response, but do not create or edit actual files.
5. If you need to reference something in the codebase, **ask the user first** rather than reading it yourself.
6. When the user asks to implement or write code, **always invoke the gelex-code-writer agent** via the Agent tool. Do not write code yourself.

Topic: $ARGUMENTS
