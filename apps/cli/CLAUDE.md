# CLI Subcommand Architecture

Each subcommand lives in `apps/cli/<name>/` with four layers:

- **Args** `args.h/.cpp` — `setup_<name>_command`: create the subcommand, register flags, and attach the CLI11 callback
- **Config** `config.h/.cpp` — `make_*_config`
- **Command** `command.h/.cpp` — `<name>_execute`: calls `make_config`, constructs reporter + engine, runs engine with event visitor
- **Reporter** `reporter.h/.cpp` — one `on_event(const XxxEvent&)` overload per engine event type (optional if no progress/logging needed)

## Events

Engine emits a `std::variant<Event1, Event2, ...>` type. Command layer dispatches
via `std::visit` to reporter overloads. Define event structs in `gelex/infra/logging/`.

## Reporter Headers

- Forward-declare `spdlog::logger` in reporter headers; `#include <spdlog/logger.h>` only in `.cpp` files.
- Do not forward-declare event structs when their event header is already included.

## Registration

In `main.cpp`: call `setup_<name>_command(program, exit_code)`.
In `apps/cli/<name>/args.cpp`: create `program.add_subcommand("<name>")`,
register options, and attach the CLI11 callback.
In `CMakeLists.txt`: append all `.cpp` files to `CLI_SOURCES`.
