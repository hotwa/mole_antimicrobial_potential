# Implementation Plans

This directory stores retained implementation plans for substantial repository
changes that were executed in this repo.

Use it for:

- historical execution plans for major refactors or throughput work
- agent-oriented checklists that explain why a change was structured a certain
  way
- plan artifacts that are still useful for reconstructing implementation
  intent

Do not use it for:

- canonical user documentation
- API or CLI reference material
- one-off scratch notes with no lasting architectural value

Retention rule:

- keep plans here only when they describe a meaningful implementation strategy
  that future agents may need to understand
- move purely disposable or abandoned planning notes to `.trash/`

Canonical operational docs remain under `docs/`, especially:

- `docs/repo_layout.md`
- `docs/cli_reference.md`
- `docs/batch_screening_input_format.md`
