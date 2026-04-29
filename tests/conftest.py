"""pytest fixtures shared across the workflow's unit tests.

Adds workflow/scripts to sys.path so tests can import the user-facing
scripts as plain Python modules without a package install.
"""
import os.path as op
import sys

REPO_ROOT = op.dirname(op.dirname(op.abspath(__file__)))
SCRIPTS_DIR = op.join(REPO_ROOT, "workflow", "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)
