from math import log2

from hop import Hop, dir0, dir1
from synthetic import generate_hop, probe_hop_without_jamming

BITCOIN = 100 * 1000 * 1000
MIN_CAPACITY_SYNTHETIC = 0.01 * BITCOIN
MAX_CAPACITY_SYNTHETIC = 10 * BITCOIN


def test_generate_hop():
    hop: Hop = generate_hop(
        1,
        1,
        0.01 * BITCOIN,
        0.01 * BITCOIN,
        0,
        [
            0.004 * BITCOIN,
        ],
    )
    """
    The uncertainty is the amount of bits required to represent the
    difference between the lower bound and the upper bound. If no probes
    have been done, this is equal to the capacity + 1. This is because
    the lower bound is strict, so innitially it's -1.
    """
    assert hop.uncertainty == log2(0.01 * BITCOIN + 1)


def test_probe_hop_without_jamming_single_channel_pss():
    pss = True
    hop: Hop = Hop([150000], [0], [0], [70000], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 18
    assert hop.b_l[0] == 69999
    assert hop.b_u[0] == 70000
    assert hop.h_l == 69999
    assert hop.h_u == 70000
    assert hop.g_l == 79999
    assert hop.g_u == 80000


def test_probe_hop_without_jamming_single_channel():
    pss = False
    hop: Hop = Hop([150000], [0], [0], [70000], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 18
    assert hop.b_l[0] == 69999
    assert hop.b_u[0] == 70000
    assert hop.h_l == 69999
    assert hop.h_u == 70000
    assert hop.g_l == 79999
    assert hop.g_u == 80000


def test_probe_hop_without_jamming_multi_channel():
    pss = False
    hop: Hop = Hop([100000, 60000], [0, 1], [0, 1], [72345, 23458], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 33
    assert hop.b_l == [72344, 23457]
    assert hop.b_u == [72345, 23458]
    assert hop.h_l == 72344
    assert hop.h_u == 72345
    assert hop.g_l == 36541
    assert hop.g_u == 36542


def test_probe_hop_without_jamming_multi_channel_pss():
    pss = True
    hop: Hop = Hop([100000, 60000], [0, 1], [0, 1], [72345, 23458], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 34
    assert hop.b_l == [35802, -1]
    assert hop.b_u == [95803, 60000]
    assert hop.h_l == 95802
    assert hop.h_u == 95803
    assert hop.g_l == 64196
    assert hop.g_u == 64197


def test_probe_hop_without_jamming_multi_channel_pss_empty():
    pss = True
    hop: Hop = Hop([100000, 60000], [0, 1], [0, 1], [0, 0], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 17
    assert hop.b_l == [-1, -1]
    assert hop.b_u == [0, 0]
    assert hop.h_l == -1
    assert hop.h_u == 0
    assert hop.g_l == 159999
    assert hop.g_u == 160000


def test_probe_hop_without_jamming_multi_channel_pss_full():
    pss = True
    hop: Hop = Hop([100000, 60000], [0, 1], [0, 1], [100000, 60000], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 34
    assert hop.b_l == [99999, 59999]
    assert hop.b_u == [100000, 60000]
    assert hop.h_l == 159999
    assert hop.h_u == 160000
    assert hop.g_l == -1
    assert hop.g_u == 0


def test_probe_hop_without_jamming_multi_channel_pss_near_empty():
    pss = True
    hop: Hop = Hop([100000, 60000], [0, 1], [0, 1], [20000, 10000], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 17
    assert hop.b_l == [-1, -1]
    assert hop.b_u == [30000, 30000]
    assert hop.h_l == 29999
    assert hop.h_u == 30000
    assert hop.g_l == 129999
    assert hop.g_u == 130000


def test_probe_hop_without_jamming_multi_channel_pss_near_full():
    pss = True
    hop: Hop = Hop([100000, 60000], [0, 1], [0, 1], [90000, 50000], pss=pss)
    num_probes = probe_hop_without_jamming(hop, True, pss)
    assert num_probes == 34
    assert hop.b_l == [79999, 39999]
    assert hop.b_u == [100000, 60000]
    assert hop.h_l == 139999
    assert hop.h_u == 140000
    assert hop.g_l == 19999
    assert hop.g_u == 20000
