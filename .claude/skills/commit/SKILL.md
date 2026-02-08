---
name: git-commit
description: Generate professional conventional commit messages with emojis. Use when users ask to create, generate, write, or format commit messages from diffs, code changes, git status output, or change descriptions. Also trigger for git-related tasks like "write a commit for this", "what should my commit message be", or "create commit messages for these changes".
---

# Git Commit Message Generator

Generate clear, professional commit messages following conventional commit standards with appropriate emojis.

## Core Principles

1. **Conventional Commits Format**: Every commit follows `<emoji> <type>: <description>`
2. **Present Tense Imperative**: Write as commands ("add feature" not "added feature")
3. **Atomic Commits**: Each message represents a single logical change
4. **Clear and Concise**: Subject line under 72 characters when possible
5. **Context-Aware**: Incorporate project-specific guidelines when available

## Input Analysis

When the user provides changes, analyze them to identify:

- **Type of change**: Feature, fix, refactor, docs, etc.
- **Scope**: What area of the codebase is affected
- **Breaking changes**: Are there API or behavior changes?
- **Multiple concerns**: Should this be split into multiple commits?

### Input Types You'll Receive

- Git diffs (`git diff` output)
- File change summaries
- Descriptions of modifications
- Git status output
- Multiple files with unrelated changes

## Commit Types and Emojis

Use these standard types with their corresponding emojis:

| Emoji | Type            | When to Use                                     |
| ----- | --------------- | ----------------------------------------------- |
| ✨    | `feat`          | New features or functionality                   |
| 🐛    | `fix`           | Bug fixes (not warnings/security)               |
| 🚨    | `fix`           | Fix compiler/linter warnings                    |
| 🔒️   | `fix`           | Fix security vulnerabilities                    |
| 📝    | `docs`          | Documentation only changes                      |
| 💄    | `style`         | Code formatting, whitespace, missing semicolons |
| ♻️    | `refactor`      | Code restructuring without behavior change      |
| ⚡️   | `perf`          | Performance improvements                        |
| ✅    | `test`          | Adding or updating tests                        |
| 🔧    | `chore`         | Build process, dependencies, tooling            |
| 🚀    | `ci`            | CI/CD configuration and scripts                 |
| 🗑️    | `revert`        | Reverting previous commits                      |
| 💥    | `feat` or `fix` | Breaking changes (add ! after type)             |

## Message Structure

### Subject Line (Required)

```
<emoji> <type>[optional scope]: <description>

Examples:
✨ feat: add user authentication system
🐛 fix(api): resolve null pointer in user endpoint
♻️ refactor(auth): simplify token validation logic
📝 docs: update API documentation for v2 endpoints
```

### Body (Optional, for Complex Changes)

Add a body when:

- The change needs additional context
- There's important rationale to explain
- The impact isn't obvious from the subject
- There are breaking changes to describe

Format:

```
<emoji> <type>: <subject>

<blank line>
<detailed explanation of what and why>

<blank line>
BREAKING CHANGE: <description of breaking change if applicable>
```

## Commit Splitting Logic

**When to suggest multiple commits:**

1. **Different concerns**: Changes to unrelated parts of the codebase
   - Example: Database migration + frontend styling = 2 commits

2. **Different types**: Mixing features, fixes, and refactoring
   - Example: New feature + bug fix = 2 commits

3. **File patterns**: Source code vs documentation vs configuration
   - Example: Code changes + README update = 2 commits

4. **Logical grouping**: Changes easier to review separately
   - Example: API endpoint + tests for that endpoint = potentially 1 commit
   - Example: Multiple unrelated API endpoints = multiple commits

5. **Size considerations**: Very large diffs that would be clearer if broken down
   - Use judgment: Is this a cohesive change or multiple changes bundled together?

**When suggesting splits, provide:**

- Clear commit message for each suggested commit
- List of files/changes for each commit
- Explanation of why the split improves clarity

## Workflow

### 1. Analyze Input

Carefully examine what the user provided:

```python
# Check for:
- Is this a diff, description, or file list?
- Are there project guidelines (CLAUDE.md, Agent.md, CONTRIBUTING.md)?
- How many distinct logical changes are present?
- What's the primary purpose of each change?
```

### 2. Identify Project Context

If available, look for project-specific guidelines:

- Check for CLAUDE.md, Agent.md, or CONTRIBUTING.md
- Look for existing commit patterns in git log if provided
- Incorporate project-specific scopes or conventions

### 3. Determine Commit Type

Ask yourself:

- Does this add new functionality? → `feat`
- Does this fix incorrect behavior? → `fix`
- Does this restructure without changing behavior? → `refactor`
- Does this only change documentation? → `docs`
- Does this improve performance? → `perf`
- Does this add/update tests? → `test`
- Is this tooling/build/deps? → `chore`

### 4. Evaluate Breaking Changes

If the change:

- Modifies public API
- Changes expected behavior
- Requires users to update their code
- Removes features

Then add `!` after the type: `feat!:` or `fix!:`
And use the 💥 emoji

### 5. Write the Message

**Subject line:**

- Start with emoji + type
- Add scope in parentheses if helpful: `fix(api):`
- Use imperative mood: "add", "fix", "update"
- Keep under 72 characters
- Don't end with a period

**Body (if needed):**

- Explain what and why (not how - code shows how)
- Wrap at 72 characters
- Separate from subject with blank line

### 6. Self-Verify

Before outputting, check:

- ✓ Emoji matches the type
- ✓ Type is correct for the change
- ✓ Present tense, imperative mood
- ✓ Clear and specific (not vague)
- ✓ Subject under 72 chars
- ✓ Breaking changes properly marked
- ✓ Multiple concerns addressed (split if needed)

## Handling Unclear Input

If the input is vague or insufficient:

**Ask for clarification:**

- "Could you provide the git diff or describe what changed?"
- "What files were modified and how?"
- "Is this a new feature, bug fix, or something else?"
- "Are there any breaking changes?"

**Never guess** - accuracy is more important than speed.

## Examples

### Example 1: Simple Feature

**Input**: "Added login functionality with JWT tokens"

**Output**:

```
✨ feat: add JWT-based user authentication

Implement login system with:
- JWT token generation and validation
- Password hashing with bcrypt
- Session management middleware
```

### Example 2: Bug Fix

**Input**:

```diff
- if (user.age > 18) {
+ if (user.age >= 18) {
```

**Output**:

```
🐛 fix: correct age validation logic

Change age check from > 18 to >= 18 to properly include 18-year-olds
```

### Example 3: Multiple Commits Needed

**Input**: "I updated the README, fixed a null pointer bug in auth.js, and added a new API endpoint for user profiles"

**Output**:

```
I recommend splitting this into 3 separate commits:

1. 📝 docs: update README
   Files: README.md

2. 🐛 fix(auth): resolve null pointer exception
   Files: auth.js

3. ✨ feat(api): add user profile endpoint
   Files: routes/profile.js, controllers/profileController.js

This makes each change easier to review and revert if needed.
```

### Example 4: Breaking Change

**Input**: "Renamed the 'getUserData' function to 'fetchUserProfile' and changed its return type"

**Output**:

```
💥 refactor!: rename getUserData to fetchUserProfile

BREAKING CHANGE: The getUserData function has been renamed to
fetchUserProfile. Additionally, the return type has changed from
object to Promise<UserProfile>. Update all calls accordingly.

Migration:
- getUserData(id) → await fetchUserProfile(id)
```

### Example 5: Refactor

**Input**: "Extracted duplicate code in user and admin services into a shared utility module"

**Output**:

```
♻️ refactor: extract shared validation logic to utils

Move common validation functions from user and admin services
to shared utils/validation.js to reduce code duplication
```

## Special Cases

### Merge Commits

For merge commits, use:

```
🔀 merge: merge feature/user-auth into main
```

### Reverting

When reverting:

```
🗑️ revert: revert "feat: add experimental cache layer"

This reverts commit abc1234 due to performance issues in production.
```

### Initial Commit

For the very first commit:

```
🎉 chore: initial commit
```

### Release/Version Bumps

```
🔖 chore: bump version to 2.1.0
```

## Output Format

Always provide the commit message in a clean, ready-to-use format:

```
<emoji> <type>: <description>
```

Or with body:

```
<emoji> <type>: <description>

<body>
```

**Do not** include:

- Extra commentary about the message
- "Here's your commit message:"
- Markdown code fences (unless specifically requested)
- Git commands like `git commit -m`

Just output the message itself, ready to copy and paste or use directly.

## Advanced: Reading Project Guidelines

If the user's repository has a CLAUDE.md, Agent.md, or CONTRIBUTING.md file, read it to understand:

- Project-specific commit conventions
- Custom scopes they use
- Special prefixes or suffixes
- Team preferences for commit granularity
- Examples from the project

Incorporate these guidelines to make commits that match the project's existing style.

## Final Checklist

Before every output, verify:

1. ✓ Correct emoji for the type
2. ✓ Conventional commit format followed
3. ✓ Present tense, imperative mood
4. ✓ Subject line is clear and concise
5. ✓ Breaking changes properly indicated
6. ✓ Split suggested if multiple concerns present
7. ✓ Ready to use (no extra formatting or commentary)

---

Remember: Your goal is to help developers communicate their changes clearly and consistently. A good commit message tells the story of what changed and why, making the codebase easier to understand and maintain.
