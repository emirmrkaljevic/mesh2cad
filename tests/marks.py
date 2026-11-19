"""One representative case from a noise or seed sweep stays in the default run.

The other values are marked ``slow``. ``pytest -m slow`` runs those, and the
weekly workflow runs the whole suite.
"""
import pytest


def sweep(values):
    """Keep the first value. Mark every later value slow."""
    first, *rest = values
    return [first, *(pytest.param(value, marks=pytest.mark.slow) for value in rest)]
