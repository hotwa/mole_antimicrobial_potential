# .trash

This directory is the quarantine area for files that should not stay in the
repository root or in the canonical implementation layout, but that were kept
intentionally instead of being deleted.

Use it for:

- one-off patch helpers
- scratch scripts
- temporary sample inputs
- runtime leftovers such as logs
- caches that were moved out of the main tree

Rules:

- Do not treat files here as canonical entrypoints.
- Do not import from `.trash/` in production code.
- If a file becomes useful again, move it back into a proper directory such as
  `scripts/`, `docs/`, or `test/`.
