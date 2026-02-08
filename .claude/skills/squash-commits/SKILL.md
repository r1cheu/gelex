---
name: squash-commits
description: Reorganize unpushed commits into clean, logical commits. Use when the user wants to squash commits, combine fixup commits, clean up commit history before pushing, or reorganize commits into logical units. Triggers on "squash commits", "clean up commits", or "reorganize history".
---

# Squash/Reorganize Unpushed Commits

Combine multiple small unpushed commits into clean, logical commits. Can create a single commit or multiple well-organized commits.

## Workflow

### Step 1: Check Repository State

```bash
git status
git log --oneline @{u}..HEAD
```

Ensure working directory is clean. If there are uncommitted changes, ask user to commit or stash first.

### Step 2: Analyze Unpushed Commits

Display commits and analyze what they do:

```bash
# Show unpushed commits with details
git log @{u}..HEAD --oneline

# Show full commit messages
git log @{u}..HEAD --format="%h %s%n%b%n---"

# Count commits
git rev-list --count @{u}..HEAD
```

Group commits by logical purpose (e.g., feature additions, bug fixes, refactoring, tests).

### Step 3: Ask User for Reorganization Plan

Present the commits and ask how to reorganize:

**Option A: Squash all into one commit**
- Simplest option
- Good when all commits serve one logical purpose

**Option B: Squash into multiple logical commits**
- Example: One commit for feature, one for tests, one for docs
- Better for complex changes with distinct aspects

**Option C: Interactive rebase**
- User wants full control over squashing, reordering, editing

### Step 4: Execute Reorganization

**For Option A (Single Commit):**

```bash
COMMIT_COUNT=$(git rev-list --count @{u}..HEAD)
git reset --soft HEAD~$COMMIT_COUNT
git commit -m "type: descriptive message"
```

**For Option B (Multiple Clean Commits):**

Use interactive rebase:

```bash
git rebase -i @{u}
```

In the editor:
- Group related commits together
- Mark first commit of each group as `pick`
- Mark others in group as `squash` or `fixup`
- Optionally `reword` to improve messages

Example:
```
pick abc1234 Add feature X
squash def5678 Fix typo in feature X
squash ghi9012 Add feature X tests
pick jkl3456 Update documentation
squash mno7890 Fix doc formatting
```

**For Option C (Interactive Rebase):**

```bash
git rebase -i @{u}
```

Let user handle the interactive rebase directly.

### Step 5: Craft Commit Messages

For each new/combined commit:
- Follow project conventions (conventional commits: `feat:`, `fix:`, `refactor:`, etc.)
- Describe the logical change, not implementation details
- Reference related commits if helpful

### Step 6: Verify Result

```bash
# View new commit history
git log --oneline @{u}..HEAD

# Review each commit
git log @{u}..HEAD --format="%h %s%n%b%n---"

# Check what changed overall
git diff @{u}..HEAD --stat
```

Ensure:
- All changes are preserved
- Commits are logical and well-described
- History is clean and readable

## Common Patterns

### Pattern 1: Feature + Fixups → Single Feature Commit

```
Before:
- Add feature X
- Fix typo
- Fix bug in feature X
- Add missing test

After:
- feat: add feature X with tests
```

### Pattern 2: Multiple Logical Units

```
Before:
- Add feature X
- Fix typo in X
- Add tests for X
- Update docs
- Fix doc typo

After:
- feat: add feature X
- test: add tests for feature X
- docs: update documentation for feature X
```

### Pattern 3: Keep Some, Squash Others

```
Before:
- refactor: extract helper function
- Add feature X using helper
- Fix feature X bug
- Add feature X tests

After:
- refactor: extract helper function
- feat: add feature X with tests
```

## Arguments

- `squash-commits` - Interactive reorganization of unpushed commits
- `squash-commits all` - Squash all unpushed into single commit
- `squash-commits <N>` - Squash last N commits

## Safety & Rollback

Always working with unpushed commits only (safe, no force push needed).

If something goes wrong:
```bash
git reflog
git reset --hard HEAD@{1}
```

## Notes

- Default behavior: analyze commits and suggest logical grouping
- Preserve all code changes, only reorganize commits
- Follow project's commit message conventions
- Use conventional commits format when appropriate
