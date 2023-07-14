import pytest

from experiments import run_one_instance_of_experiment
from prober import Prober

SNAPSHOT_FILENAME = "./tests/data/channels.json"
ENTRY_CHANNEL_CAPACITY = 10 * 100 * 1000 * 1000
# top 10 nodes by degree as per https://1ml.com/node?order=channelcount
ENTRY_NODES = ["Alice", "Bob", "Carol"]
# very high-dimensional hops are rare in snapshots - too few to run
# experiments on
MAX_MAX_NUM_CHANNELS = 5


@pytest.fixture
def prober():
    prober = Prober(
        SNAPSHOT_FILENAME, "PROBER", ENTRY_NODES, ENTRY_CHANNEL_CAPACITY
    )
    return prober


@pytest.mark.parametrize(
    "remote_probing, bs, pss",
    [
        (False, True, False),
        (False, True, True),
        (True, True, False),
        (True, True, True),
    ],
)
def test_run_one_instance_of_experiment(remote_probing, bs, pss, prober):
    gains_line, speed_line = run_one_instance_of_experiment(
        False, remote_probing, bs, pss, [1], 1, prober, 1
    )
    assert True
