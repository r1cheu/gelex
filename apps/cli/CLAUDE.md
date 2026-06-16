# CLI Subcommand Architecture

Each subcommand lives in `apps/cli/<name>/` with four layers:

- **Args** `args.h/.cpp` — `setup_<name>_command`: create the subcommand, register flags, and attach the CLI11 callback
- **Config** `config.h/.cpp` — `make_*_config`
- **Command** `command.h/.cpp` — `<name>_execute`: calls `make_config`, constructs reporters, runs core operations, and directly reports command-visible state
- **Reporter** `reporter.h/.cpp` — CLI presentation methods plus `as_observer()` only when a core component emits internal progress/state

## Reporting

- Use observer events only for state produced inside core components, such as algorithm iteration progress, checkpoint progress, reader chunk progress, or internal state pointers.
- Do not create events for data already available in `command.cpp`: banner, config, dataset intersection, prior summaries, output files, and completion summaries should call reporter/printer methods directly.
- Define reusable core observer events in `gelex/infra/logging/` only when the event is not CLI presentation data.

## Reporter Headers

- Forward-declare `spdlog::logger` in reporter headers; `#include <spdlog/logger.h>` only in `.cpp` files.
- Do not forward-declare event structs when their event header is already included.

## Registration

In `main.cpp`: call `setup_<name>_command(program, exit_code)`.
In `apps/cli/<name>/args.cpp`: create `program.add_subcommand("<name>")`,
register options, and attach the CLI11 callback.
In `CMakeLists.txt`: append all `.cpp` files to `CLI_SOURCES`.
