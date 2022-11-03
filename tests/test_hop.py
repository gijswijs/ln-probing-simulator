import pytest

from hop import Hop


def test_can_forward(hop_instance: Hop):
    can_forward_dir0 = hop_instance.can_forward(0)
    can_forward_dir1 = hop_instance.can_forward(1)

    assert can_forward_dir0
    assert can_forward_dir1


@pytest.fixture
def hop_instance():
    hop = Hop([150000], [0], [0])
    return hop
