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
    assert hop.uncertainty == 61.03726544366254
    # assert hop.S_F == 14_400_988_024_500_260_001 1.3e+19
    assert hop.S_F == 2_366_179_851_025_390_001
    assert hop.b_l == [-1, -1, -1, 9999]
    assert hop.b_u == [40000, 50000, 80000, 100000]


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
    amount = 140001
    hop.probe(dir0, amount, True)
    rect = Rectangle([0, 0], [100001, 60001])
    cut = rect.cut(160000 - amount)
    assert hop.S_F == 100001 * 60001 - cut


def test_S_F_probe_fail_small():
    hop = Hop([100, 60], [0], [0, 1], [], [26, 23], pss=True)
    assert hop.S_F == 101 * 61
    amount = 51
    hop.probe(dir0, amount, True)
    rect = Rectangle([0, 0], [101, 61])
    cut = rect.cut(amount)
    assert hop.S_F == 51 * 61 - cut


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
    assert hop.S_F == 100001 * 60001 - 20000 * 10000


def test_S_F_dir1_probe_success():
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 20000
    hop.probe(dir1, amount, True)
    assert hop.S_F == 100001 * 60001 - 20000 * 10000


def test_S_F_dir1_probe_fail():
    hop = Hop([100000, 60000], [0, 1], [0, 1], [], [26000, 23000], pss=True)
    assert hop.S_F == 100001 * 60001
    amount = 140001
    hop.probe(dir1, amount, True)
    assert hop.S_F == 100001 * 60001 - 20000 * 10000


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
    hop = Hop(C, [0, 1], [0, 1], [], B)
    hop.set_h_and_g(pss)
    direction = hop.next_dir(bs, jamming, pss=pss)
    print("dir0" if direction == dir0 else "dir1")
    assert True
