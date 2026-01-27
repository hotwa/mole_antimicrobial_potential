"""Predictor singleton management."""

from typing import Optional

from src.predictor import AntimicrobialPredictor

_predictor: Optional[AntimicrobialPredictor] = None


def get_predictor() -> AntimicrobialPredictor:
    """Return the singleton predictor instance."""
    global _predictor
    if _predictor is None:
        _predictor = AntimicrobialPredictor()
    return _predictor
