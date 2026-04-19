"""Plotting submodule (scanpy-style ``pl`` namespace).

All functions accept a fitted :class:`tmesplit.TMESplit` model and return a
matplotlib :class:`Axes` (or list of Axes). Matplotlib + seaborn are optional
dependencies — install with ``pip install tmesplitpy[plot]``.
"""

from tmesplit.pl._activities import activities
from tmesplit.pl._network import network
from tmesplit.pl._programs import programs
from tmesplit.pl._significance import significance

__all__ = ["programs", "activities", "network", "significance"]
