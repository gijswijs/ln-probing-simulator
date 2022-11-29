import os

import pytest

from hop import Hop
from prober import Prober

SNAPSHOT_FILENAME = "./tests/data/channels.json"
ENTRY_CHANNEL_CAPACITY = 10 * 100 * 1000 * 1000
# top 10 nodes by degree as per https://1ml.com/node?order=channelcount
ENTRY_NODES = ["Alice", "Bob", "Carol"]


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
                "Carol",
                "Judy",
            ),  # target_node_pair
            True,  # bs
            False,  # Jamming
            False,  # pss
        ),
        (
            (
                "Harold",
                "Alice",
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
    target_hop: Hop = (
        prober.psshopgraph[target_node_pair[0]][target_node_pair[1]]["hop"]
        if pss
        else prober.lnhopgraph[target_node_pair[0]][target_node_pair[1]]["hop"]
    )
    assert target_hop.g_l == target_hop.g - 1
    assert target_hop.g_u == target_hop.g
    assert target_hop.h_l == target_hop.h - 1
    assert target_hop.h_u == target_hop.h
    assert target_hop.b_l == list(map(lambda x: x - 1, target_hop.b))
    assert target_hop.b_u == target_hop.b


def test_channels_json():
    filename = "./tests/data/channels.json"
    entry_nodes = ["Alice", "Bob", "Carol"]
    Prober(filename, "PROBER", entry_nodes, ENTRY_CHANNEL_CAPACITY)

    # paths = list(prober.paths_for_amount(("Dave", "Erin"), 2000))

    assert True


def test_export_graph():
    filename = "./tests/data/channels.json"
    entry_nodes = ["Alice", "Bob", "Carol"]
    prober = Prober(filename, "PROBER", entry_nodes, ENTRY_CHANNEL_CAPACITY)
    filename_lnhopgraph = "./tests/data/hopgraph-2022-11-25.gml"
    filename_psshopgraph = "./tests/data/psshopgraph-2022-11-25.gml"
    prober.export_graph(prober.lnhopgraph, filename_lnhopgraph)
    prober.export_graph(prober.psshopgraph, filename_psshopgraph)

    assert os.path.exists(filename_lnhopgraph)
    assert os.path.exists(filename_psshopgraph)


def test_import_graph():
    filename = [
        "./tests/data/hopgraph-2022-11-25.gml",
        "./tests/data/psshopgraph-2022-11-25.gml",
    ]
    prober = Prober(filename, "PROBER", None, None, gml_file=True)
    assert prober.lnhopgraph["Alice"]["Bob"]["hop"].g_u == 60000
    assert prober.lnhopgraph["Alice"]["Bob"]["hop"].h_u == 100000
    assert prober.psshopgraph["Alice"]["Bob"]["hop"].g_u == 424885
    assert prober.psshopgraph["Alice"]["Bob"]["hop"].h_u == 524885
