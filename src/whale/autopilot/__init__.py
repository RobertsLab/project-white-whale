"""Autopilot: a scheduled controller that moves each study one step per tick.

Step A (this package): deterministic and side-effect free. `whale autopilot
tick --dry-run` reads where every study is, decides the one transition each
would make next, matches failed runs against a playbook of known failures and
prints what it would do. Nothing is launched, pushed, commented on or sent to
an API. See docs/autopilot.md for what steps B-F add.
"""
from whale.autopilot.model import STATES, Decision, RunEvidence, StudyState  # noqa: F401
