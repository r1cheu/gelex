---
name: dev-assist
description: Enter incremental development mode. Work through fixes/changes one at a time, verify each before moving on. Optimized for the user's concise, fix-by-fix workflow.
disable-model-invocation: true
argument-hint: [task description or "next" to continue]
---

You are now in **dev-assist mode** for incremental development. Follow these rules:

## Workflow

1. **One change at a time** — never batch multiple fixes. Complete one, verify build, then ask "next?".
2. **Build after every change** — run `pixi r build-debug` after each edit. Do not proceed until green.
3. **Wait before LSP queries** — after editing, wait a few seconds before checking clangd diagnostics.
4. **Short status updates** — after each fix, report only: what changed, build result. No summaries of what was already known.

## Communication Protocol

- The user communicates in very concise Chinese (sometimes 1-2 characters like "开始", "要的", "下一个").
- Mirror this brevity. Lead with action, not explanation.
- When the user says "下一个" or "next", move to the next pending issue without re-explaining context.
- Do not restate what the user said or what you already know.

## When the User Gives "Why" but Not "How"

The user tends to explain motivation rather than implementation steps. When you receive a directive that is clear on *why* but ambiguous on *how*:

1. **State your proposed approach in 1-2 sentences** before implementing.
2. Wait for confirmation ("要的", "开始", etc.) or correction.
3. This prevents wasted work from misinterpreting intent.

Example:
> User: "这一段明显可以接着化简"
> You: "我打算把 X 提取成 Y，用 Z 替换原来的循环。可以吗？"

## Code Changes

- **Always delegate to gelex-code-writer agent** for non-trivial implementations (new classes, multi-file refactors).
- For small, targeted edits (rename, delete a line, fix a type), edit directly.
- After editing, if the change touches headers, grep for affected call sites and fix them too.
- Never leave broken call sites — fix test call sites after refactoring.

## Revert Protocol

- If the user says "取消" or "回退", revert immediately without asking for confirmation.
- Revert cleanly — remove all traces, not just the main change.

## Progress Tracking

- Use tasks to track multi-step work within the session.
- Mark each task done immediately after completion, not in batches.
- When starting a session with pending work, list remaining items briefly.

Topic: $ARGUMENTS
