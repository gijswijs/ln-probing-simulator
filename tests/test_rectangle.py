import pytest

from rectangle import Rectangle


def test_box():
    "Cube with width 6, depth 2, height 4 should have a volume of 7 * 3 * 5, since all boundaries are inclusive"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.S() == 7 * 3 * 5
