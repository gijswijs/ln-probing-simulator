import pytest

from rectangle import Rectangle


def test_rectangle():
    "Triangle with width 6, depth 2, height 4 should have a volume of 7 * 3 * 5, since all boundaries are inclusive"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.S() == 7 * 3 * 5


def test_cut_2d():
    "A cut of 9 of this 3-dimensional box should return 33"
    rectangle = Rectangle([0, 0], [5, 5])
    assert rectangle.cut(9) == 33


def test_cut_leq():
    "A cut of 4 of this 2-dimensional box should return 13"
    rectangle = Rectangle([0, 0], [3, 3])
    assert rectangle.cut(4, "<=") == 13


def test_cut_geq():
    "A cut of 4 of this 2-dimensional box should return 6"
    rectangle = Rectangle([0, 0], [3, 3])
    assert rectangle.cut(4, ">=") == 6


def test_cut_greater():
    "A cut of 4 of this 2-dimensional box should return 6"
    rectangle = Rectangle([0, 0], [3, 3])
    assert rectangle.cut(4, ">") == 3


def test_cut_less():
    "A cut of 4 of this 2-dimensional box should return 6"
    rectangle = Rectangle([0, 0], [3, 3])
    assert rectangle.cut(4, "<") == 10


def test_cut_4d():
    "A cut of 12 of this 3-dimensional box should return all 4**4 -1 = 255"
    rectangle = Rectangle([0, 0, 0, 0], [3, 3, 3, 3])
    assert rectangle.cut(12) == 255


def test_cut_a():
    "A cut of 1 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [1, 1, 1])
    assert rectangle.cut(1) == 1


def test_cut_b():
    "A cut of 1 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [1, 1, 1])
    assert rectangle.cut(2) == 4


def test_cut_c():
    "A cut of 1 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [1, 1, 1])
    assert rectangle.cut(4) == 8


def test_cut_d():
    "A cut of 1 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [3, 1, 2])
    assert rectangle.cut(6) == 23


def test_cut_1():
    "A cut of 1 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(1) == 1


def test_cut_2():
    "A cut of 2 of this 3-dimensional box should return 4"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(2) == 4


def test_cut_3():
    "A cut of 3 of this 3-dimensional box should return 9"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 10 because the # of lattice points in a pyramid with
    # sides 2 is 10
    assert rectangle.cut(3) == 10


def test_cut_4():
    "A cut of 4 of this 3-dimensional box should return 9"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 9 because the # of lattice points in a pyramid with
    # sides 3 is 20, but a pyramid of size 0 sticks outside of the
    # rectangle, because it has depth 2. That accounts for 1 lattice
    # points that should be deducted.
    assert rectangle.cut(4) == 19


def test_cut_5():
    "A cut of 5 of this 3-dimensional box should return 24"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 24 because the # of lattice points in a pyramid with
    # size 4 is 35, but we have a pyramids sticking out of the
    # rectangle, with size 1. That accounts for 4 lattice points that
    # should be deducted.
    assert rectangle.cut(5) == 31


def test_cut_7():
    "A cut of 7 of this 3-dimensional box should return 38"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 60 because the # of lattice points in a pyramid with
    # size 6 is 84, but we have 2 pyramids sticking out of the
    # rectangle, with size 3 and 1, whose lattice points sum to 20 +
    # 4 = 24.
    assert rectangle.cut(7) == 60


def test_cut_8():
    "A cut of 8 of this 3-dimensional box should return 38"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 74 because the # of lattice points in a pyramid with
    # size 7 is 120, but we have 3 pyramids sticking out of the
    # rectangle, with size 4, 2 and 0, whose lattice points sum to 35 + 10 +
    # 1 = 46. We also have overlap of 4 lattice points
    assert rectangle.cut(8) == 74


def test_cut_11():
    "A cut of 11 of a 3-dimensional box should return 101"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 101 because the # of lattice points in a pyramid with
    # size 10 is 286, but we have 3 pyramids sticking out of the
    # rectangle, with size 7, 5 and 3, whose lattice points sum to 120 +
    # 56 + 20 = 196. We have two overlapping pyramids of size 0 and 2,
    # so 1 + 10 = 11 overlapping lattice points
    assert rectangle.cut(11) == 101


def test_cut_11_translated():
    "A cut of 11 of a 3-dimensional box should return 101"
    rectangle = Rectangle([2, 5, 6], [8, 7, 10])
    # We expect 101 because the # of lattice points in a pyramid with
    # size 10 is 286, but we have 3 pyramids sticking out of the
    # rectangle, with size 7, 5 and 3, whose lattice points sum to 120 +
    # 56 + 20 = 196. We have two overlapping pyramids of size 0 and 2,
    # so 1 + 10 = 11 overlapping lattice points
    assert rectangle.cut(24) == 101


def test_cut_12():
    "A cut of 12 of a 3-dimensional box should return 104"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 104 because the # of lattice points in a pyramid with
    # size 11 is 364, but we have 3 pyramids sticking out of the
    # rectangle, with size 8, 6 and 4, whose lattice points sum to 156 +
    # 84 + 35 = 284. We have two overlapping pyramids of size 1 and 3,
    # so 4 + 20 = 24 overlapping lattice points -> 364 - 284 + 24
    assert rectangle.cut(12) == 104


def test_cut_13():
    "A cut of 13 of a 3-dimensional box should return 105"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    # We expect 105 because 13 is bigger than the max_value (6+2+4) so
    # all lattice points inside the box should be returned.
    assert rectangle.cut(13) == 105
