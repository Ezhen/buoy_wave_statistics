# Tests

```bash
pip install pytest --break-system-packages   # not yet in requirements.txt - add it there
pytest tests/ -v
```

`test_stage01.py` covers the two Stage 01 correctness bugs found and
fixed in the 2026-07-19 session (see `CHANGELOG.md`): era-mismatch
fabrication, timestamp-jitter false missingness, and the gap-erasure
regression caught while fixing the first one. These were originally
one-off synthetic scripts built by hand in that session and thrown away
after use — this file exists so the next regression in this stage gets
caught automatically instead of requiring someone to rebuild the same
cases from scratch again.

Add a new `test_stageNN.py` here whenever a stage gets a synthetic
known-answer validation built for it (per this project's own
established discipline — see `PLAN_next_session.md`'s Working
Discipline #2) rather than letting the validation live only in a
throwaway script or a chat transcript.
