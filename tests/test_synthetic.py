from math import log2

import pytest

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


@pytest.mark.parametrize(
    "C, B, e_dir0, e_dir1, pss, bs, success",
    [
        (
            [150000],  # C
            [70000],  # B
            [0],  # e_dir0
            [0],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 18,
                "b_l": [69999],
                "b_u": [70000],
                "h_l": 69999,
                "h_u": 70000,
                "g_l": 79999,
                "g_u": 80000,
            },  # success
        ),
        (
            [150000],  # C
            [70000],  # B
            [0],  # e_dir0
            [0],  # e_dir1
            False,  # pss
            True,  # bs
            {
                "num_probes": 18,
                "b_l": [69999],
                "b_u": [70000],
                "h_l": 69999,
                "h_u": 70000,
                "g_l": 79999,
                "g_u": 80000,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            False,  # pss
            True,  # bs
            {
                "num_probes": 33,
                "b_l": [72344, 23457],
                "b_u": [72345, 23458],
                "h_l": 72344,
                "h_u": 72345,
                "g_l": 36541,
                "g_u": 36542,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 18,
                "b_l": [35802, -1],
                "b_u": [95803, 60000],
                "h_l": 95802,
                "h_u": 95803,
                "g_l": 64196,
                "g_u": 64197,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [0, 0],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 17,
                "b_l": [-1, -1],
                "b_u": [0, 0],
                "h_l": -1,
                "h_u": 0,
                "g_l": 159999,
                "g_u": 160000,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [100000, 60000],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 17,
                "b_l": [99999, 59999],
                "b_u": [100000, 60000],
                "h_l": 159999,
                "h_u": 160000,
                "g_l": -1,
                "g_u": 0,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [20000, 10000],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 17,
                "b_l": [-1, -1],
                "b_u": [30000, 30000],
                "h_l": 29999,
                "h_u": 30000,
                "g_l": 129999,
                "g_u": 130000,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [90000, 50000],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 17,
                "b_l": [79999, 39999],
                "b_u": [100000, 60000],
                "h_l": 139999,
                "h_u": 140000,
                "g_l": 19999,
                "g_u": 20000,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [0, 1],  # e_dir0
            [0],  # e_dir1
            False,  # pss
            True,  # bs
            {
                "num_probes": 34,
                "b_l": [72344, -1],
                "b_u": [72345, 60000],
                "h_l": 72344,
                "h_u": 72345,
                "g_l": 27654,
                "g_u": 27655,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [0, 1],  # e_dir0
            [0],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 33,
                "b_l": [72344, 23457],
                "b_u": [72345, 23458],
                "h_l": 95802,
                "h_u": 95803,
                "g_l": 27654,
                "g_u": 27655,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [1],  # e_dir0
            [0, 1],  # e_dir1
            False,  # pss
            True,  # bs
            {
                "num_probes": 33,
                "b_l": [63457, 23457],
                "b_u": [100000, 23458],
                "h_l": 23457,
                "h_u": 23458,
                "g_l": 36541,
                "g_u": 36542,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [1],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 33,
                "b_l": [72344, 23457],
                "b_u": [72345, 23458],
                "h_l": 23457,
                "h_u": 23458,
                "g_l": 64196,
                "g_u": 64197,
            },  # success
        ),
        (
            [100000, 60000, 40000],  # C
            [72345, 23458, 20000],  # B
            [0, 1, 2],  # e_dir0
            [0, 1],  # e_dir1
            False,  # pss
            True,  # bs
            {
                "num_probes": 34,
                "b_l": [63457, 23457, -1],
                "b_u": [72345, 60000, 40000],
                "h_l": 72344,
                "h_u": 72345,
                "g_l": 36541,
                "g_u": 36542,
            },  # success
        ),
        (
            [100000, 60000, 40000],  # C
            [72345, 23458, 20000],  # B
            [0, 1, 2],  # e_dir0
            [0, 1],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 33,
                "b_l": [35802, -1, -1],
                "b_u": [95803, 60000, 40000],
                "h_l": 115802,
                "h_u": 115803,
                "g_l": 64196,
                "g_u": 64197,
            },  # success
        ),
        (
            [2017461, 2017461, 2017461],  # C
            [701210, 717798, 1172118],  # B
            [0, 1, 2],  # e_dir0
            [0, 1, 2],  # e_dir1
            True,  # pss
            True,  # bs
            {
                "num_probes": 42,
                "b_l": [-1, -1, -1],
                "b_u": [2017461, 2017461, 2017461],
                "h_l": 2591125,
                "h_u": 2591126,
                "g_l": 3461256,
                "g_u": 3461257,
            },  # success
        ),
        (
            [100000, 60000],  # C
            [72345, 23458],  # B
            [0, 1],  # e_dir0
            [0, 1],  # e_dir1
            False,  # pss
            False,  # bs
            {
                "num_probes": 32,
                "b_l": [72344, 23457],
                "b_u": [72345, 23458],
                "h_l": 72344,
                "h_u": 72345,
                "g_l": 36541,
                "g_u": 36542,
            },  # success
        ),
        (
            [2017461, 2017461, 2017461],  # C
            [701210, 717798, 1172118],  # B
            [0, 1, 2],  # e_dir0
            [0, 1, 2],  # e_dir1
            True,  # pss
            False,  # bs
            {
                "num_probes": 44,
                "b_l": [-1, -1, -1],
                "b_u": [2017461, 2017461, 2017461],
                "h_l": 2591125,
                "h_u": 2591126,
                "g_l": 3461256,
                "g_u": 3461257,
            },  # success
        ),
    ],
)
def test_probe_hop_without_jamming(C, B, e_dir0, e_dir1, pss, bs, success):
    hop: Hop = Hop(C, e_dir0, e_dir1, [], B)
    hop.set_h_and_g(pss)
    num_probes = probe_hop_without_jamming(hop, bs, pss)
    assert num_probes == success["num_probes"]
    assert hop.b_l == success["b_l"]
    assert hop.b_u == success["b_u"]
    assert hop.h_l == success["h_l"]
    assert hop.h_u == success["h_u"]
    assert hop.g_l == success["g_l"]
    assert hop.g_u == success["g_u"]
