import pytest

from hop import Hop, dir0, dir1
from rectangle import Rectangle


@pytest.fixture
def hop_instance():
    hop = Hop([150000], [0], [0], [], [70000])
    return hop


@pytest.fixture
def hop_pss_instance():
    hop = Hop([150000, 80000], [0, 1], [0, 1], [1], [70000, 30000], pss=True)
    return hop


def test_hop_init(hop_instance, hop_pss_instance):
    assert hop_instance.N == 1
    assert hop_pss_instance.N == 2
    assert hop_instance.h_u == 150000
    assert hop_pss_instance.h_u == 230000


def test_can_forward(hop_instance: Hop):
    can_forward_dir0 = hop_instance.can_forward(dir0)
    can_forward_dir1 = hop_instance.can_forward(dir1)

    assert can_forward_dir0
    assert can_forward_dir1


def test_probe_hop_with_single_channel(hop_instance: Hop):
    pss = False
    hop_instance.set_h_and_g(pss)
    hop_instance.probe(dir0, 20000, pss)
    b_l_dir0 = hop_instance.b_l
    h_l_dir0 = hop_instance.h_l
    g_l_dir0 = hop_instance.g_l
    b_u_dir0 = hop_instance.b_u
    h_u_dir0 = hop_instance.h_u
    g_u_dir0 = hop_instance.g_u
    hop_instance.reset_estimates(pss)
    hop_instance.probe(dir1, 20000, pss)
    b_l_dir1 = hop_instance.b_l
    h_l_dir1 = hop_instance.h_l
    g_l_dir1 = hop_instance.g_l
    b_u_dir1 = hop_instance.b_u
    h_u_dir1 = hop_instance.h_u
    g_u_dir1 = hop_instance.g_u
    assert h_u_dir0 == g_u_dir1
    assert g_u_dir0 == h_u_dir1
    assert h_l_dir0 == g_l_dir1
    assert g_l_dir0 == h_l_dir1
    assert b_l_dir0[0] == g_l_dir1
    assert b_u_dir0[0] == g_u_dir1
    assert b_l_dir1[0] == g_l_dir0
    assert b_u_dir1[0] == g_u_dir0


def test_hop_LattE_error():
    pss = True
    hop = Hop([249999, 425939], [0], [1], [], [155000, 425939], pss=pss)
    hop.b_l[1] = 425938
    hop.probe(dir0, 249999, pss)
    assert hop.S_F == 249999


@pytest.mark.parametrize(
    "h_l, result",
    [(-1, 4), (0, 3), (1, 1)],
)
def test_sfgp_unit_square(h_l, result):
    "S_F_generic_pss different h_ls with unit square"
    pss = True
    hop = Hop([1, 1], [0, 1], [0, 1], pss=pss)
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    assert hop.S_F_generic_pss(h_l, 2, -1, 2, R_b) == result


@pytest.mark.parametrize(
    "h_l, result",
    [(-1, 8), (0, 7), (1, 4), (2, 1)],
)
def test_sfgp_unit_cube(h_l, result):
    "S_F_generic_pss different h_ls with unit cube"
    pss = True
    hop = Hop([1, 1, 1], [0, 1, 2], [0, 1, 2], pss=pss)
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    assert hop.S_F_generic_pss(h_l, 3, -1, 3, R_b) == result


@pytest.mark.parametrize(
    "h_l, h_u, g_l, g_u, result",
    [
        (0, 12, -1, 12, 104),
        (1, 12, -1, 12, 101),
        (2, 12, -1, 12, 95),
        (-1, 9, -1, 12, 95),
        (2, 9, -1, 12, 85),
        (-1, 12, 0, 12, 104),
        (-1, 12, 1, 12, 101),
        (-1, 12, 2, 12, 95),
        (-1, 12, -1, 9, 95),
        (-1, 12, 2, 9, 85),
        (3, 12, -1, 12, 86),
        (4, 12, -1, 12, 74),
        (6, 12, -1, 12, 45),
        (7, 12, -1, 12, 31),
        (10, 12, -1, 12, 4),
        (11, 12, -1, 12, 1),
    ],
)
def test_sfgp(h_l, h_u, g_l, g_u, result):
    "S_F_generic_pss different h_ls with unit cube"
    pss = True
    hop = Hop([6, 2, 4], [0, 1, 2], [0, 1, 2], pss=pss)
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    assert hop.S_F_generic_pss(h_l, h_u, g_l, g_u, R_b) == result


@pytest.mark.parametrize(
    "h_l, h_u, g_l, g_u, c, b_l, result",
    [
        (8, 10, -1, 10, [5] * 2, [-1] * 2, 3),
        (11, 12, -1, 12, [3] * 4, [-1] * 4, 1),
        (8, 9, -1, 9, [3] * 3, [-1] * 3, 1),
        (0, 3, -1, 3, [1] * 3, [-1] * 3, 7),
        (1, 3, -1, 3, [1] * 3, [-1] * 3, 4),
        (5, 6, -1, 6, [3, 1, 2], [-1] * 3, 1),
        (23, 25, -1, 25, [8, 7, 10], [1, 4, 5], 4),
    ],
)
def test_sfgp_with_b_l(h_l, h_u, g_l, g_u, c, b_l, result):
    "S_F_generic_pss different h_ls with unit cube"
    pss = True
    hop = Hop(c, range(len(c)), range(len(c)), pss=pss)
    hop.b_l = b_l
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    assert hop.S_F_generic_pss(h_l, h_u, g_l, g_u, R_b) == result


def test_sfgp_too_big():
    "S_F_generic_pss cannot return a value bigger than self.S()"
    pss = True
    c = [
        15_000_000,
        231_922,
        2_000_000,
        4_000_000,
        15_000_000,
        2_000_000,
        1_000_000,
    ]
    hop = Hop(c, range(len(c)), range(len(c)), pss=pss)
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    assert hop.S_F_generic_pss(19_615_960, sum(c), -1, sum(c), R_b) < R_b.S()


def test_sfgp_valueerror():
    "Out of bounds h_l should raise an ValueError"
    pss = True
    hop = Hop([2000000], [0], [0], pss=pss)
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    with pytest.raises(ValueError):
        hop.S_F_generic_pss(2249999, 2250000, -1, 2000000, R_b)


def test_probe_hop_with_multiple_channels_no_pss():
    pss = False
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [70000, 50000])
    hop.set_h_and_g(pss)
    hop.probe(dir0, 65000, pss)
    # The below assert is valid with the original code, but I disagree
    # with it. If an probing amount is so large it can only be forwarded
    # by the largest channel, you can update the b_l of that channel
    # accordingly. b_l[0] should be 64999.
    # FIXME: Fix this error for non-pss probing??
    assert hop.b_l == [-1, -1]
    hop.probe(dir1, 39000, pss)
    # b_l[0] is now adjusted because of the failed probe in dir1. But it
    # should still be 64999, which is bigger
    assert hop.b_l == [61000, 21000]
    assert hop.h_l == 64999
    assert hop.h_u == 100000
    assert hop.g_l == -1
    assert hop.g_u == 38999


def test_probe_hop_with_multiple_channels_pss():
    pss = True
    hop = Hop(
        [40000, 50000, 80000, 100000],
        [0, 1, 2, 3],
        [0, 1, 2, 3],
        [],
        [31000, 32000, 35000, 90000],
    )
    hop.set_h_and_g(pss)
    hop.probe(dir0, 180000, pss)
    assert hop.h_l == 179999
    assert hop.h_u == 270000
    assert hop.h == 188000
    assert hop.g_u == 90000
    assert hop.g_l == -1
    assert hop.uncertainty == 61.037464446274
    # assert hop.S_F == 14_400_988_024_500_260_001 1.3e+19
    assert hop.S_F == 2_366_506_259_887_662_501
    assert hop.b_l == [-1, -1, -1, 9999]
    assert hop.b_u == [40000, 50000, 80000, 100000]


def test_probe_hop_with_error_pss():
    pss = True
    hop = Hop(
        [1000000, 5568000, 500000, 5568000, 500000],
        [0, 3, 4],
        [0, 1, 2],
        [1, 2, 3, 4],
        [355622, 3126295, 469077, 4079402, 130907],
        pss,
    )
    hop.set_h_and_g(pss)
    hop.probe(dir0, 3117005, pss)
    # assert hop.h_l == 179999
    # assert hop.h_u == 270000
    # assert hop.h == 188000
    # assert hop.g_u == 90000
    # assert hop.g_l == -1
    # assert hop.uncertainty == 61.03751921524554
    # # assert hop.S_F == 14_400_988_024_500_260_001 1.3e+19
    # assert hop.S_F == 2_366_596_101_171_202_501
    # assert hop.b_l == [-1, -1, -1, 9999]
    # assert hop.b_u == [40000, 50000, 80000, 100000]


@pytest.mark.parametrize(
    "C, B, pss, direction",
    [
        ([100000, 60000], [26000, 23000], True, dir0),
        ([100000, 60000], [74000, 37000], True, dir1),
        ([100000, 60000], [74000, 37000], True, dir0),
        ([100000, 60000], [26000, 23000], True, dir1),
    ],
)
def test_next_a(C, B, pss, direction):
    hop = Hop(C, [0, 1], [0, 1], [], B, pss=pss)
    hop.set_h_and_g(pss)
    initial_S_F = hop.S_F
    amount = hop.next_a(direction, False, False, pss)
    print(amount)
    hop.probe(direction, amount, pss)
    new_S_F = hop.S_F
    assert round(initial_S_F / new_S_F) == 2


def test_S_F_probe_fail():
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 140_001
    R_b = Rectangle([b_l_i + 1 for b_l_i in hop.b_l], hop.b_u)
    # Amount fails, because 26_000 + 23_000 < 140_001. So amount - 1 is
    # now h_u. In this case (both channels are bidrectional) this
    # directly affects g_l as well, because g_l = sum(c) - h_u - 1 =
    # 160000 - 140000 - 1 = 19999.
    S_F = hop.S_F_generic_pss(
        -1, amount - 1, sum(hop.c) - amount, sum(hop.c), R_b
    )
    hop.probe(dir0, amount, True)
    assert hop.S_F == S_F


def test_S_F_probe_fail_small():
    hop = Hop([100, 60], [0], [0, 1], [], [26, 23], pss=True)
    assert hop.S_F == 101 * 61
    amount = 51
    hop.probe(dir0, amount, True)
    assert hop.S_F == 51 * 61


def test_S_F_probe_fail_2():
    hop = Hop([100000, 60000], [0], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 80001
    hop.probe(dir0, amount, True)
    assert hop.S_F == 80001 * 60001


def test_S_F_probe_success():
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 20000
    hop.probe(dir0, amount, True)
    assert hop.S_F == 100001 * 60001 - 20001 * 10000


def test_S_F_dir1_probe_success():
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 20000
    hop.probe(dir1, amount, True)
    assert hop.S_F == 100001 * 60001 - 20001 * 10000


def test_S_F_dir1_probe_fail():
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 140001
    hop.probe(dir1, amount, True)
    assert hop.S_F == 100001 * 60001 - 20001 * 10000


@pytest.mark.parametrize(
    "C, B, pss, bs, jamming",
    [
        ([100000, 60000], [26000, 23000], True, False, False),
        ([100000, 60000], [74000, 37000], True, False, False),
        ([100000, 60000], [74000, 37000], True, False, False),
        ([100000, 60000], [26000, 23000], True, False, False),
    ],
)
def test_next_dir(C, B, pss, bs, jamming):
    # TODO: I doubt if this has any use in a PSS setting. I guess it
    # works, but the calculations seem meaningless to me
    hop = Hop(C, [0, 1], [0, 1], [], B, pss=pss)
    hop.set_h_and_g(pss)
    direction = hop.next_dir(bs, jamming, pss=pss)
    print("dir0" if direction == dir0 else "dir1")
    assert True
