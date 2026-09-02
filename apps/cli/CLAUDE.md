# CLI Subcommand Architecture

Computational subcommands normally live in `apps/cli/<name>/` and separate the
following responsibilities when they exist:

- **Args** `args.h/.cpp` — `setup_<name>_command`: create the subcommand, allocate `auto config = std::make_shared<<Name>Config>()`, bind options directly to `config` members via `add_option(name, config->member, ...)`, and attach the CLI11 callback (the callback captures `config`, tying its lifetime to the `App` since registration and `parse()` live in different scopes)
- **Config** `config.h` — typed `<Name>Config` struct holding the command's options with their default values; bound by Args, consumed by Command
- **Command** `command.h/.cpp` — `<name>_execute(const <Name>Config&)`: takes only the typed config (no `CLI::App`), constructs reporters, runs core operations, and directly reports command-visible state. Command preamble (logging init, banner, `report_command_line`, timing) lives in `cli::execute_cli_command` (the harness), which is the only layer that touches `CLI::App` at run time
- **Reporter** `reporter.h/.cpp` — CLI presentation methods; add `as_observer()` only when a core component emits internal progress/state

These are responsibility boundaries, not a mandatory file template. Specialized
commands such as `completion` may omit layers that have no responsibility, and a
command may add cohesive `compute` or `io` components when needed.

## Reporting

- Use observer events only for state produced inside core components, such as algorithm iteration progress, checkpoint progress, reader chunk progress, or internal state pointers.
- Do not create events for data already available in `command.cpp`: banner, config, dataset intersection, prior summaries, output files, and completion summaries should call reporter/printer methods directly.
- Define reusable core observer events in `gelex/infra/logging/` only when the event is not CLI presentation data.

## Registration

In `main.cpp`: call `setup_<name>_command(program, exit_code)`.
In `apps/cli/<name>/args.cpp`: create `program.add_subcommand("<name>")`,
register options, and attach the CLI11 callback.
`apps/CMakeLists.txt` discovers `apps/cli/**/*.cpp` automatically.
