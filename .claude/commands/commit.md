---
description: Generates a commit message from staged (or unstaged) changes.
argument-hint: [optional context...]
allowed-tools: Bash(git diff:*)
---

# Task: Generate a Conventional Commit Message

Based on the following code changes, please generate a concise and descriptive commit message that follows the Conventional Commits specification.

The template is as follows:

<type>[optional scope]: <description>
[optional body]

The commit types are:
feat: Commits, which adds a new feature
fix: Commits, that fixes a bug
refactor: refactored code that neither fixes a bug nor adds a feature but rewrites/restructures your code.
chore : Changes that do not relate to a fix or feature and don’t modify src or test files basically miscellaneous commits (for example, updating dependencies or modifying .gitignore file)
perf : Commits are special refactor commits, geared towards improving performance.
ci : Continuous integration related.
ops : Commits, that affect operational components like infrastructure, deployment, backup , recovery …
build : Changes that affect the build system build tool, ci pipeline, dependencies, project version, …
docs : Commits, that affect documentation, such as the README.
style : changes that do not affect the meaning of the code, likely related to code formatting such as white-space, missing semi-colons, etc.
revert: reverts a previous commit.
test:commits that add missing tests or correct existing tests

rules you must follow:
Limit the subject line to 50 characters
Capitalize the subject/description line
Do not end the subject line with a period
Separate the subject from the body with a blank line
Wrap the body at 72 characters
Use the body to explain what and why
Use the imperative mood in the subject line let it seem like you’re giving a command eg “feat: Add unit tests for user authentication”. Using the imperative mood in commit messages makes them more consistent and commands-like, which is helpful in understanding the actions taken.
Do not add info of Claude or AI in the commit message.

## Staged Changes:

!git diff --cached

## Unstaged Changes:

!git diff

## Optional User Context:

$ARGUMENTS
