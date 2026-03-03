# CLI Subcommand Architecture

Each subcommand lives in `apps/cli/<name>/` with four layers:

- **Args** `<name>_args.h/.cpp` — register flags on `ArgumentParser`
- **Config** `<name>_config.h/.cpp` — satisfies `CliConfig` concept so `make_config<T>()` calls both
- **Command** `<name>_command.h/.cpp` — `int <name>_execute(parser)`: calls `make_config`, constructs reporter + engine, runs engine with event visitor
- **Reporter** `<name>_reporter.h/.cpp` — one `on_event(const XxxEvent&)` overload per engine event type (optional if no progress/logging needed)

## Events

Engine emits a `std::variant<Event1, Event2, ...>` type. Command layer dispatches
via `std::visit` to reporter overloads. Define event structs in `gelex/infra/logging/`.

## Registration

In `main.cpp`: declare `argparse::ArgumentParser <name>("<name>")`, add a
`CommandDescriptor{name, &parser, setup_fn, execute_fn}` entry.
In `CMakeLists.txt`: append all `.cpp` files to `CLI_SOURCES`.
