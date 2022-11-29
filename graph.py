#! ./env/bin/python3

"""
This file is part of Lightning Network Probing Simulator.

Copyright © 2020-2021 University of Luxembourg

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

SPDX-FileType: SOURCE SPDX-FileCopyrightText: 2020-2021 University of
Luxembourg SPDX-License-Identifier: MIT
"""

"""
  Auxiliary operations with the LN graph.
"""

import json

import networkx as nx

from hop import Hop, dir0, dir1


class Channel:
    def __init__(
        self, source, destination, capacity, dir0_enabled, dir1_enabled
    ):
        self.source = source
        self.destination = destination
        self.capacity = capacity
        self.dir0_enabled = dir0_enabled
        self.dir1_enabled = dir1_enabled


def create_multigraph_from_snapshot(snapshot_filename):
    """
    Create a NetworkX multigraph from a clightning's listchannels.json
    snapshot. Multigraph means each edge corresponds to an edge
    (parallel edges allowed).

    Parameters:
    - snapshot_filename: path to the snapshot

    Return:
    - the multigraph (the maximal connected component only).
    """
    print("Creating LN graph from file:", snapshot_filename, "...")
    with open(snapshot_filename, "r") as snapshot_file:
        network = json.load(snapshot_file)
    edges_set, nodes_set = set(), set()
    edges, channels = [], dict()
    # cid -> Channel
    for channel_direction in network["channels"]:
        cid = channel_direction["short_channel_id"]
        direction = (
            channel_direction["source"] < channel_direction["destination"]
        )
        if direction == dir0:
            source = channel_direction["source"]
            destination = channel_direction["destination"]
        else:
            source = channel_direction["destination"]
            destination = channel_direction["source"]
        if cid not in channels:
            # print("creating new channel for", cid)
            dir0_enabled, dir1_enabled = (
                (channel_direction["active"], False)
                if direction == dir0
                else (False, channel_direction["active"])
            )
            channel = Channel(
                source,
                destination,
                channel_direction["satoshis"],
                dir0_enabled,
                dir1_enabled,
            )
            channels[cid] = channel
        else:
            # print("updating existing channels for", cid)
            channel = channels[cid]
            if direction == dir0:
                channel.dir0_enabled = channel_direction["active"]
            else:
                channel.dir1_enabled = channel_direction["active"]
    # count how many uni-directional channels we have
    num_bidirectional = sum(
        [
            1
            for cid in channels
            if channels[cid].dir0_enabled and channels[cid].dir1_enabled
        ]
    )
    print("Total channels:", len(channels))
    print("Bidirectional channels:", num_bidirectional)
    for cid in channels:
        channel = channels[cid]
        edges.append(
            (
                channel.source,
                channel.destination,
                cid,
                {
                    "capacity": channel.capacity,
                    "dir0_enabled": channel.dir0_enabled,
                    "dir1_enabled": channel.dir1_enabled,
                },
            )
        )
        edges_set.add(cid)
        nodes_set.add(source)
        nodes_set.add(destination)
    nodes = list(nodes_set)
    g = nx.MultiGraph()
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    print(
        "LN snapshot contains:",
        g.number_of_nodes(),
        "nodes,",
        g.number_of_edges(),
        "channels.",
    )
    # continue with the largest connected component
    print(
        "Components:",
        nx.number_connected_components(g),
        "\nContinuing with the largest component.",
    )
    largest_cc = max(nx.connected_components(g), key=len)
    g = g.subgraph(largest_cc).copy()
    print(
        "LN graph created with",
        g.number_of_nodes(),
        "nodes,",
        g.number_of_edges(),
        "channels.",
    )
    return g


def ln_multigraph_to_hop_graph(ln_multigraph):
    """
    Generate a hopgraph from an LN multigraph. A hopgraph doesn't allow
    parallel edges. Instead, parallel channels are encoded in edge
    attributes.

    Parameters:
    - ln_multigraph: LN model multigraph

    Return:
    - hop_graph: a non-directed graph where each edge models a hop
    """
    hop_graph = nx.Graph()
    # initialize hop graph with nodes and empty edge attributes
    for n1, n2 in ln_multigraph.edges():
        hop_graph.add_nodes_from([n1, n2])
        hop_graph.add_edge(n1, n2)
    for n1, n2, k, d in ln_multigraph.edges(keys=True, data=True):
        multi_edge = ln_multigraph[n1][n2]
        cids = [cid for cid in multi_edge]
        capacities, e_dir0, e_dir1 = [], [], []
        for i, cid in enumerate(cids):
            capacities.append(multi_edge[cid]["capacity"])
            if multi_edge[cid]["dir0_enabled"]:
                e_dir0.append(i)
            if multi_edge[cid]["dir1_enabled"]:
                e_dir1.append(i)
        hop_graph[n1][n2]["hop"] = Hop(capacities, e_dir0, e_dir1)
    return hop_graph


def ln_hopgraph_to_pss_hopgraph(ln_hopgraph, max_length=2):
    """
    Generate a pss_hopgraph from an ln_hopgraph. A pss_hopgraph is a
    hopgraph that includes all routes equal or shorter than max_length
    from n1 to n2.

    Parameters:
    - ln_hopgraph: LN model ln_hopgraph
    - max_length: maximum length of route from n1 to n2. Defaults to 2,
      which means only one intermediary hop between n1 and n2.

    Return:
    - pss_hopgraph: a non-directed graph where each edge models a pss
      hop
    """
    pss_hopgraph = ln_hopgraph.copy()
    pss_hopdigraph = pss_hopgraph.to_directed()

    for _, _, d in pss_hopgraph.edges(data=True):
        hop = d["hop"]
        hop.set_h_and_g(pss=True)

    def filter_edge(n1, n2):
        """
        Return True if the edge is kept, False if it is excluded.

        Parameters:
        - n1, n2: node IDs of the vertices
        """
        hop = ln_hopgraph[n1][n2]["hop"]
        direction = dir0 if n1 < n2 else dir1
        if direction == dir0:
            return hop.can_forward(dir0)
        else:
            return hop.can_forward(dir1)

    def find_alt_paths(n1, n2):
        paths = nx.shortest_simple_paths(routing_graph, source=n1, target=n2)
        next_path = []
        while True:
            try:
                next_path = next(paths)
                if len(next_path) <= max_length + 1:
                    if len(next_path) > 2:
                        node_pairs = [p for p in zip(next_path, next_path[1:])]
                        capacity = None
                        for u, v in node_pairs:
                            direction = dir0 if u < v else dir1
                            if direction == dir0:
                                new_capacity = pss_hopgraph[u][v]["hop"].h_u
                            else:
                                new_capacity = pss_hopgraph[u][v]["hop"].g_u
                            if capacity is None or new_capacity < capacity:
                                capacity = new_capacity
                        capacities.append(capacity)
                        direction = dir0 if n1 < n2 else dir1
                        if direction == dir0:
                            e_dir0.append(len(capacities) - 1)
                        else:
                            e_dir1.append(len(capacities) - 1)
                        alt.append(len(capacities) - 1)
                else:
                    break
            except (nx.exception.NetworkXNoPath, StopIteration):
                break

    routing_graph = nx.subgraph_view(
        pss_hopdigraph,
        filter_edge=filter_edge,
    )

    for n1, n2, d in pss_hopgraph.edges(data=True):
        hop = d["hop"]
        hop.set_h_and_g(pss=True)
        capacities = hop.c.copy()
        e_dir0 = hop.e[dir0].copy()
        e_dir1 = hop.e[dir1].copy()
        alt = []
        find_alt_paths(n1, n2)
        find_alt_paths(n2, n1)
        pss_hopgraph[n1][n2]["hop_pss"] = Hop(
            capacities, e_dir0, e_dir1, alt, pss=True
        )

    hops_with_alt_routes = 0
    for n1, n2, d in pss_hopgraph.edges(data=True):
        if d["hop_pss"].N > d["hop"].N:
            hops_with_alt_routes += 1
        d["hop"] = d.pop("hop_pss")

    for n1, n2, d in ln_hopgraph.edges(data=True):
        d["hop"].set_h_and_g(pss=False)

    print(
        "PSS hopgraph created with",
        pss_hopgraph.number_of_nodes(),
        "nodes,",
        pss_hopgraph.number_of_edges(),
        "hops.",
    )

    print(hops_with_alt_routes, "hops with at least one alternative route")

    return pss_hopgraph
