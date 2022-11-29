import matplotlib.pyplot as plt
import networkx as nx
from networkx import Graph, MultiDiGraph, MultiGraph

from graph import (
    create_multigraph_from_snapshot,
    ln_hopgraph_to_pss_hopgraph,
    ln_multigraph_to_hop_graph,
)
from hop import Hop, dir0, dir1

JSON_FILE = "./tests/data/channels.json"


def draw_graph(G, filename, arrowstyle="-"):
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos)
    labels = {node: node[0:1] for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels)

    ax = plt.gca()

    arcs = dict()

    for e in G.edges:
        nodes = sorted(list(e[0:2]))
        if nodes[0] not in arcs:
            arcs[nodes[0]] = {nodes[1]: 0}
        elif nodes[1] not in arcs[nodes[0]]:
            arcs[nodes[0]][nodes[1]] = 0
        else:
            arcs[nodes[0]][nodes[1]] += 0.3
        ax.annotate(
            "",
            xy=pos[e[1]],
            xycoords="data",
            xytext=pos[e[0]],
            textcoords="data",
            arrowprops=dict(
                arrowstyle=arrowstyle,
                color="0.5",
                shrinkA=10,
                shrinkB=10,
                patchA=None,
                patchB=None,
                connectionstyle="arc3,rad=rrr".replace(
                    "rrr", str(arcs[nodes[0]][nodes[1]])
                ),
            ),
        )

    ax.margins(0.20)
    plt.axis("off")
    plt.savefig(filename)
    plt.cla()


def test_create_multigraph_from_snapshot():
    multigraph: MultiGraph = create_multigraph_from_snapshot(JSON_FILE)
    draw_graph(multigraph, "./tests/data/graph.png")

    multidigraph: MultiDiGraph = multigraph.to_directed()
    draw_graph(multidigraph, "./tests/data/digraph.png", "->")

    def filter_edge(n1, n2, key):
        direction = dir0 if n1 < n2 else dir1
        if direction == dir0:
            return multidigraph[n1][n2][key].get("dir0_enabled", False)
        else:
            return multidigraph[n1][n2][key].get("dir1_enabled", False)

    routinggraph = nx.subgraph_view(multidigraph, filter_edge=filter_edge)
    draw_graph(routinggraph, "./tests/data/routinggraph.png", "->")
    assert True


def test_ln_hopgraph_to_pss_hopgraph():
    multigraph: MultiGraph = create_multigraph_from_snapshot(JSON_FILE)
    lnhopgraph: Graph = ln_multigraph_to_hop_graph(multigraph)
    psshopgraph: Graph = ln_hopgraph_to_pss_hopgraph(lnhopgraph)
    assert lnhopgraph["Alice"]["Bob"]["hop"].h_u == 100000
    assert (
        psshopgraph["Alice"]["Bob"]["hop"].h_u
        == 100000 + 60000 + 20000 + 344885
    )
    assert lnhopgraph["Alice"]["Grace"]["hop"].h_u == 1000099
    assert (
        psshopgraph["Alice"]["Grace"]["hop"].h_u
        == 1000099 + 100000 + 60000 + 20000
    )
    assert lnhopgraph["Bob"]["Grace"]["hop"].h_u == 344885
    assert psshopgraph["Bob"]["Grace"]["hop"].h_u == 344885 + 60000 + 20000
