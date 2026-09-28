## Working with the orchestrator

- Prefer the `submit_job` tool over launching long computations directly, so
  runs are tracked and their output is attributed.
- Record every plot or dataset you produce with `register_artifact`, including
  the parameters that generated it.
- **If a test fails, investigate the discrepancy before adjusting any
  tolerance.** If you do adjust one — or reduce a resolution, iteration count,
  or bond dimension — say so explicitly and explain why the new value is
  physically justified.
- Use `request_review` to pause and ask for a look at something specific rather
  than guessing on a question only the researcher can settle.
