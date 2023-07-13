import pytest

from rectangle import Rectangle


def test_box():
    "Cube with width 6, depth 2, height 4 should have a volume of 7 * 3 * 5, since all boundaries are inclusive"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.S() == 7 * 3 * 5


def test_cut_2d():
    "A h_l (strict) of 8 of this 2-dimensional box should return S() - 33 = 36 - 33 = 3"
    rectangle = Rectangle([0, 0], [5, 5])
    assert rectangle.cut(8, 10) == 3


def test_cut_too_big():
    "An overzised cut should raise an ValueError"
    rectangle = Rectangle([0], [2000000])
    with pytest.raises(ValueError):
        rectangle.cut(2249999, 2250000)


@pytest.mark.parametrize(
    "cut, result",
    [(-1, 4), (0, 3), (1, 1)],
)
def test_cut_unit_square(cut, result):
    "h_l with unit square"
    rectangle = Rectangle(
        [0, 0],
        [1, 1],
    )
    assert rectangle.cut(cut, 2) == result


@pytest.mark.parametrize(
    "cut, result",
    [(-1, 8), (0, 7), (1, 4), (2, 1)],
)
def test_cut_unit_cube(cut, result):
    "This cut returned a value bigger than self.S() which is wrong."
    rectangle = Rectangle(
        [0, 0, 0],
        [1, 1, 1],
    )
    assert rectangle.cut(cut, 3) == result


def test_cut_too_big2():
    "This cut returned a value bigger than self.S() which is wrong."
    rectangle = Rectangle(
        [0, 0, 0, 0, 0, 0, 0],
        [
            15_000_000,
            231_922,
            2_000_000,
            4_000_000,
            15_000_000,
            2_000_000,
            1_000_000,
        ],
    )
    assert rectangle.cut(19_615_960, sum(rectangle.u_vertex)) < rectangle.S()


# def test_cut_leq():
#     "A cut of 4 of this 2-dimensional box should return 13"
#     rectangle = Rectangle([0, 0], [3, 3])
#     assert rectangle.cut(4, "<=") == 13


# def test_cut_geq():
#     "A cut of 4 of this 2-dimensional box should return 6"
#     rectangle = Rectangle([0, 0], [3, 3])
#     assert rectangle.cut(4, ">=") == 6


# def test_cut_greater():
#     "A cut of 4 of this 2-dimensional box should return 6"
#     rectangle = Rectangle([0, 0], [3, 3])
#     assert rectangle.cut(4, ">") == 3


# def test_cut_less():
#     "A cut of 4 of this 2-dimensional box should return 6"
#     rectangle = Rectangle([0, 0], [3, 3])
#     assert rectangle.cut(4, "<") == 10


def test_cut_4d():
    "A h_l of 11 (strict) (4 * 3 - 1) of this 4-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0, 0], [3, 3, 3, 3])
    assert rectangle.cut(11, 12) == 1


def test_cut_3d():
    "A h_l of 8 (3 * 3 - 1) of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [3, 3, 3])
    assert rectangle.cut(8, 9) == 1


def test_cut_a():
    "A h_l of 0 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [1, 1, 1])
    assert rectangle.cut(0, 3) == 7


def test_cut_b():
    "A h_l of 1 of this 3-dimensional box should return 4"
    rectangle = Rectangle([0, 0, 0], [1, 1, 1])
    assert rectangle.cut(1, 3) == 4


# def test_cut_c():
#     "A cut of 8 of this 3-dimensional box should return 8"
#     rectangle = Rectangle([0, 0, 0], [1, 1, 1])
#     assert rectangle.cut(4) == 8


def test_cut_d():
    "A h_l of 5 of this 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [3, 1, 2])
    assert rectangle.cut(5, 6) == 1


def test_cut_1():
    "A h_l of 0 of this 3-dimensional box should return 104"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(0, 12) == 104


def test_cut_2():
    "A h_l of 1 of this 3-dimensional box should return 101"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(1, 12) == 101


def test_cut_3():
    "A h_l of 2 of this 3-dimensional box should return 95"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(2, 12) == 95


def test_cut_3_inverse():
    "A h_u of 9 of this 3-dimensional box should return 95"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(-1, 9) == 95


def test_cut_3_double():
    "A h_l of 2 and a h_u of 9 of this 3-dimensional box should return 85"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(2, 9) == 85


def test_cut_4():
    "A h_l of 3 of this 3-dimensional box should return 86"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(3, 12) == 86


def test_cut_5():
    "A h_l of 4 of this 3-dimensional box should return 74"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(4, 12) == 74


def test_cut_7():
    "A h_l of 6 of this 3-dimensional box should return 45"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(6, 12) == 45


def test_cut_8():
    "A h_l of 7 of this 3-dimensional box should return 31"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(7, 12) == 31


def test_cut_11():
    "A h_l of 10 of a 3-dimensional box should return 4"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(10, 12) == 4


def test_cut_11_translated():
    "A h_l of 10 of a 3-dimensional box should return 4"
    rectangle = Rectangle([2, 5, 6], [8, 7, 10])
    assert rectangle.cut(23, 25) == 4


def test_cut_LattE_error_1():
    "A h_l of 10 of a 3-dimensional box should return 4"
    rectangle = Rectangle([0, 425937], [249999, 425939])
    assert rectangle.cut(-1, 249999) == 250000


def test_cut_12():
    "A h_l of 11 of a 3-dimensional box should return 1"
    rectangle = Rectangle([0, 0, 0], [6, 2, 4])
    assert rectangle.cut(11, 12) == 1


# def test_cut_13():
#     "A cut of 13 of a 3-dimensional box should return 105"
#     rectangle = Rectangle([0, 0, 0], [6, 2, 4])
#     # We expect 105 because 13 is bigger than the max_value (6+2+4) so
#     # all lattice points inside the box should be returned.
#     assert rectangle.cut(13) == 105
