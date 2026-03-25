import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import pyopenms  # noqa: F401

    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")
