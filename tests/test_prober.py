import pytest

from hop import Hop, dir0, dir1
from prober import Prober

SNAPSHOT_FILENAME = "./snapshots/listchannels-2021-12-09.json"
ENTRY_CHANNEL_CAPACITY = 10 * 100 * 1000 * 1000
# top 10 nodes by degree as per https://1ml.com/node?order=channelcount
ENTRY_NODES = [
    "02ad6fb8d693dc1e4569bcedefadf5f72a931ae027dc0f0c544b34c1c6f3b9a02b",
    "03864ef025fde8fb587d989186ce6a4a186895ee44a926bfc370e2c366597a3f8f",
    "0217890e3aad8d35bc054f43acc00084b25229ecff0ab68debd82883ad65ee8266",
    "0331f80652fb840239df8dc99205792bba2e559a05469915804c08420230e23c7c",
    "0242a4ae0c5bef18048fbecf995094b74bfb0f7391418d71ed394784373f41e4f3",
    "03bb88ccc444534da7b5b64b4f7b15e1eccb18e102db0e400d4b9cfe93763aa26d",
    "03abf6f44c355dec0d5aa155bdbdd6e0c8fefe318eff402de65c6eb2e1be55dc3e",
    "02004c625d622245606a1ea2c1c69cfb4516b703b47945a3647713c05fe4aaeb1c",
    "0395033b252c6f40e3756984162d68174e2bd8060a129c0d3462a9370471c6d28f",
    "0390b5d4492dc2f5318e5233ab2cebf6d48914881a33ef6a9c6bcdbb433ad986d0",
]


@pytest.fixture
def prober():
    prober = Prober(
        SNAPSHOT_FILENAME, "PROBER", ENTRY_NODES, ENTRY_CHANNEL_CAPACITY
    )
    return prober


@pytest.mark.parametrize(
    "target_node_pair, bs, jamming, pss",
    [
        (
            (
                "037386ef907ad3e95682516bc100d971e12211f42f782ff9b0e7fd"
                + "7ac9050920a2",
                "035bddb48a2b891a8e2a704234cb13405d63c62fcd8aa1c6f5ea10"
                + "c4c027654af7",
            ),  # target_node_pair
            True,  # bs
            False,  # Jamming
            False,  # pss
        ),
        (
            (
                "037386ef907ad3e95682516bc100d971e12211f42f782ff9b0e7fd"
                + "7ac9050920a2",
                "035bddb48a2b891a8e2a704234cb13405d63c62fcd8aa1c6f5ea10"
                + "c4c027654af7",
            ),  # target_node_pair
            True,  # bs
            False,  # Jamming
            True,  # pss
        ),
    ],
)
def test_probe_hop(target_node_pair, bs, jamming, pss, prober: Prober):
    # target_hops_node_pair =
    # prober.choose_target_hops_with_n_channels(1, 1)
    prober.probe_hop(target_node_pair, bs=bs, jamming=jamming, pss=pss)
    target_hop: Hop = prober.lnhopgraph[target_node_pair[0]][
        target_node_pair[1]
    ]["hop"]
    assert target_hop.g_l == target_hop.g - 1
    assert target_hop.g_u == target_hop.g
    assert target_hop.h_l == target_hop.h - 1
    assert target_hop.h_u == target_hop.h
    assert target_hop.b_l == list(map(lambda x: x - 1, target_hop.b))
    assert target_hop.b_u == target_hop.b
