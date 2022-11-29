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
  Generation of synthetic hops and their (direct) probing.
"""

from random import randint, random

from hop import Hop, dir0, dir1


def generate_hop(
    min_N,
    max_N,
    min_capacity,
    max_capacity,
    probability_bidirectional,
    balances=None,
):
    """
    Generate a hop.

    Parameters:
    - min_N: minimum number of channels
    - max_N: maximum number of channels
    - min_capacity: minimum capacity of one channel
    - max_capacity: maximum capacity of one channel
    - probability_bidirectional: probability that a channel is enabled
      in both directions
    - balances: channel balances (generated randomly if None)

    Return:
    - a Hop instance
    """
    N = randint(min_N, max_N)
    capacities = [randint(min_capacity, max_capacity) for _ in range(N)]
    # avoid generating hops disabled in both directions (we can't probe
    # them anyway)
    hop_enabled_in_one_direction = False
    while not hop_enabled_in_one_direction:
        enabled_dir0, enabled_dir1 = [], []
        for i in range(N):
            is_bidirectional = random() < probability_bidirectional
            if is_bidirectional:
                enabled_dir0.append(i)
                enabled_dir1.append(i)
            else:
                if random() < 0.5:
                    enabled_dir0.append(i)
                else:
                    enabled_dir1.append(i)
        hop_enabled_in_one_direction = enabled_dir0 or enabled_dir1
    # print("Generating hop: capacities", capacities, "enabled_dir0",
    # enabled_dir0, "enabled_dir1", enabled_dir1)
    return Hop(capacities, enabled_dir0, enabled_dir1, [], balances)


def generate_hops(
    num_target_hops,
    N,
    min_capacity,
    max_capacity,
    probability_bidirectional=1,
):
    """
    Generate num_target_hops random hops.

    Parameters:
    - num_target_hops: the number of hops to generate
    - max_N: maximum number of channel per hop
    - max_capacity: maximum capacity per channel
    - probability_bidirectional: probability that a channel is enabled
      in a given direction

    Return:
    - a list of generated hops
    """
    return [
        generate_hop(
            N, N, min_capacity, max_capacity, probability_bidirectional
        )
        for _ in range(num_target_hops)
    ]


def probe_single_hop(hop: Hop, bs, jamming, pss=False):
    """
    Do a series of (direct) probes until the hop is fully probed.

    Parameters:
    - hop: the target hop
    - bs: amount choice method
    - jamming: do jamming-enhanced probing after h and g are fully
      probed

    Return:
    - gain: achieved information gain
    - num_probes: the total number of probes done
    - num_jams: the total number of jams done
    """
    initial_uncertainty = hop.uncertainty_pss if pss else hop.uncertainty
    num_probes, num_jams = probe_hop_without_jamming(hop, bs, pss=pss), 0
    if jamming:
        for i in range(hop.N):
            # print("\nJamming-enhanced probing channel", i)
            hop.unjam(i, direction=dir0)
            hop.unjam(i, direction=dir1)
            # TODO: can we jam in one direction only? (fewer jams)
            num_jams += hop.jam_all_except_in_direction(i, direction=dir0)
            num_jams += hop.jam_all_except_in_direction(i, direction=dir1)
            num_probes_i, num_jams_i = jam_hop_and_probe_single_channel(
                hop, bs, i
            )
            num_probes += num_probes_i
            num_jams += num_jams_i
    hop.unjam_all()
    final_uncertainty = hop.uncertainty_pss if pss else hop.uncertainty
    gain = initial_uncertainty - final_uncertainty
    # return gain in bits and used number of probes
    return gain, num_probes, num_jams


def probe_hop_without_jamming(hop: Hop, bs, pss=False):
    """
    Probe a hop without jamming.

    Parameters:
    - hop: the target hop
    - bs: amount choice method

    Return:
    - num_probes: the total number of probes done
    """
    num_probes = 0
    while hop.worth_probing_h() or hop.worth_probing_g():
        chosen_dir = hop.next_dir(bs, jamming=False, pss=pss)
        if chosen_dir is None:
            print("Hop is disabled in both directions, cannot probe")
            break
        amount = hop.next_a(chosen_dir, bs, jamming=False, pss=pss)
        hop.probe(chosen_dir, amount, pss)
        num_probes += 1
    return num_probes


def jam_hop_and_probe_single_channel(hop, bs, i):
    """
    Jam all channels in a hop except one and go jamming-enhanced
    probing.

    Parameters:
    - hop: the target hop
    - bs: amount choice method
    - i: the index of the channel to leave unjammed

    Return:
    - num_probes: the total number of probes done
    - num_jams: the total number of jams done
    """
    num_probes, num_jams = 0, 0
    while hop.worth_probing_channel(i):
        chosen_dir = hop.next_dir(bs, jamming=True)
        if chosen_dir is None:
            print("Hop is disabled in both directions, cannot probe")
            break
        amount = hop.next_a(chosen_dir, bs, jamming=True)
        hop.probe(chosen_dir, amount)
        num_probes += 1
    return num_probes, num_jams


def probe_hops_direct(hops, bs, jamming, pss=False):
    """
    Probe each hop from a list of hops.

    Parameters:
    - hops: a list of target hops
    - bs: amount choice method
    - jamming: do jamming-enhanced probing after h and g are fully
      probed

    Return:
    - total_gain: total information gain (total resolved
    uncertainty to initial uncertainty)
    - probing_speed: average bits per probe obtained
    """
    for hop in hops:
        hop.set_h_and_g(pss)
    initial_uncertainty_total = (
        sum([hop.uncertainty_pss for hop in hops])
        if pss
        else sum([hop.uncertainty for hop in hops])
    )
    gains, probes_list = [], []
    for hop in hops:
        gain, probes, jams = probe_single_hop(
            hop, bs=bs, jamming=jamming, pss=pss
        )
        gains.append(gain)
        # count jams as probes (they are payments too!)
        probes_list.append(probes + jams)
    # print("\nProbed with method:", "bs" if bs else "nbs", "with
    # jamming" if jamming else "without jamming")
    final_uncertainty_total = (
        sum([hop.uncertainty_pss for hop in hops])
        if pss
        else sum([hop.uncertainty for hop in hops])
    )
    # print("Final uncertainty:", final_uncertainty_total) print("Total
    # gain:   ", round(sum(gains),2), "after", sum(probes_list),
    # "probes") print("Average per hop:  ",
    # round(sum(gains)/len(gains),2), "after",
    # sum(probes_list)/len(probes_list), "probes")
    total_gain_bits = initial_uncertainty_total - final_uncertainty_total
    probing_speed = total_gain_bits / sum(probes_list)
    total_gain = total_gain_bits / initial_uncertainty_total
    return total_gain, probing_speed
