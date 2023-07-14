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
  A model of a hop with parallel channels.
"""
import subprocess
import tempfile
from itertools import product
from math import log2
from random import randrange

import numpy as np

from rectangle import ProbingRectangle, Rectangle

# We encode channel direction as a boolean. Direction 0 is from the
# alphanumerically lower node ID to the higher, direction 1 is the
# opposite.
dir0 = True
dir1 = False


class Hop:
    def __init__(
        self,
        capacities,
        e_dir0,
        e_dir1,
        alts=[],
        balances=None,
        granularity=1,
        pss=False,
    ):
        """
        Initialize a hop.
        For backwards compatibility we initialize without the assumption
        of pss. Before a probe we need to run hop.set_h_and_g(pss).

        Parameters:
        - capacities: a list of capacities
        - e_dir0: a list of indices of channels enabled in dir0
        - e_dir1: a list of indices of channels enabled in dir1
        - balances: a list of balances (if None, balances are generated
          randomly)
        """
        self.N = len(capacities)
        assert self.N > 0
        # ensure validity of indices
        if len(e_dir0):
            assert max(e_dir0) <= self.N
        if len(e_dir1):
            assert max(e_dir1) <= self.N
        self.c = capacities
        self.e = {dir0: e_dir0, dir1: e_dir1}  # enabled
        self.a = alts  # channels in this hop that represent alternative routes
        self.j = {dir0: [], dir1: []}  # jammed
        self.pss = pss

        if balances:
            # if balances are provided, check their consistency w.r.t.
            # capacities
            assert all(0 <= b <= c for b, c in zip(balances, capacities))
            self.b = balances
        else:
            # for each channel, pick a balance randomly between zero and
            # capacity
            self.b = [randrange(self.c[i]) for i in range(self.N)]

        self.granularity = granularity
        self.uncertainty = None  # will be set later

        self.set_h_and_g(pss=pss)

    def set_h_and_g(self, pss=False):
        # h is how much a hop can forward in dir0, if no channels are
        # jammed
        if self.can_forward(dir0) and not pss:
            self.h = max(
                [b for i, b in enumerate(self.b) if i in self.e[dir0]]
            )
        elif self.can_forward(dir0):
            self.h = sum(
                [b for i, b in enumerate(self.b) if i in self.e[dir0]]
            )
        else:
            self.h = 0
        # g is how much a hop can forward in dir1, if no channels are
        # jammed

        if self.can_forward(dir1) and not pss:
            self.g = max(
                [
                    self.c[i] - b
                    for i, b in enumerate(self.b)
                    if i in self.e[dir1]
                ]
            )
        elif self.can_forward(dir1):
            self.g = sum(
                [
                    self.c[i] - b
                    for i, b in enumerate(self.b)
                    if i in self.e[dir1]
                ]
            )
        else:
            self.g = 0

        self.reset_estimates(pss)

    def can_forward(self, direction):
        # there is at least one channel enabled and not jammed in this
        # direction
        return not all(i in self.j[direction] for i in self.e[direction])

    def jam(self, channel_index, direction):
        num_jams = 0
        if channel_index not in self.j[direction]:
            # print("jamming channel", channel_index, "in direction",
            # "dir0" if direction else "dir1")
            self.j[direction].append(channel_index)
            num_jams += 1
        return num_jams

    def jam_all_except_in_direction(self, channel_index, direction):
        num_jams = 0
        for i in self.e[direction]:
            if i != channel_index:
                num_jams += self.jam(i, direction)
        return num_jams

    def unjam(self, channel_index, direction):
        if channel_index in self.j[direction]:
            # print("unjamming channel", channel_index, "in direction",
            # "dir0" if direction else "dir1")
            self.j[direction].remove(channel_index)

    def unjam_all_in_direction(self, direction):
        for i in self.e[direction]:
            self.unjam(i, direction)

    def unjam_all(self):
        self.unjam_all_in_direction(dir0)
        self.unjam_all_in_direction(dir1)

    def get_corner_points(self):
        """
        Get the corner points of R_u_u that are not yet excluded from F.
        We could similarly get all points, but it's very slow. We use
        this as a shortcut to stop probing when only one point remains.

        Return:
        - points: a list of points (each point is a list of self.N
          coordinates).
        """
        R_b = Rectangle([b_l_i + 1 for b_l_i in self.b_l], self.b_u)
        R_u_u = self.R_h_u.intersect_with(self.R_g_u).intersect_with(R_b)
        R_u_l = self.R_h_u.intersect_with(self.R_g_l).intersect_with(R_b)
        R_l_u = self.R_h_l.intersect_with(self.R_g_u).intersect_with(R_b)
        ranges = [
            [R_u_u.l_vertex[i], R_u_u.u_vertex[i]]
            for i in range(len(R_u_u.l_vertex))
        ]
        points = []
        points_left = self.S_F
        for p in product(*ranges):
            if not R_u_l.contains_point(p) and not R_l_u.contains_point(p):
                points.append(p)
                points_left -= 1
            if points_left == 0:
                break
        return points

    def update_dependent_hop_properties(self, pss=False):
        """
        Hop is defined by the current bounds (h_l, h_u, g_l, g_u).
        Rectangle-related properties are fully determined by these
        bounds. We set these properties here. This function MUST be
        called after every bound update (such as a probe).

        A ProbingRectangle has one corner either at (0, ..., 0) or at
        (c_1, ..., c_N). This holds only for hop-level bounds (obtained
        without probing), but does not hold for balance bounds (obtained
        with probing). Hence, R_b is a Rectangle, whereas all others are
        ProbingRectangle's.
        """
        if not pss:
            self.R_h_l = ProbingRectangle(self, direction=dir0, bound=self.h_l)
            self.R_h_u = ProbingRectangle(self, direction=dir0, bound=self.h_u)
            self.R_g_l = ProbingRectangle(self, direction=dir1, bound=self.g_l)
            self.R_g_u = ProbingRectangle(self, direction=dir1, bound=self.g_u)
            self.R_b = Rectangle([b_l_i + 1 for b_l_i in self.b_l], self.b_u)
            self.S_F = self.S_F_generic(
                self.R_h_l, self.R_h_u, self.R_g_l, self.R_g_u, self.R_b
            )
        else:
            self.R_b = Rectangle([b_l_i + 1 for b_l_i in self.b_l], self.b_u)
            # optimal_h_l = max(self.h_l, sum(self.c) - self.g_u - 1)
            # optimal_h_u = min(self.h_u, sum(self.c) - self.g_l - 1)
            self.S_F = self.S_F_generic_pss(
                self.h_l, self.h_u, self.g_l, self.g_u, self.R_b
            )
        self.uncertainty = max(0, log2(self.S_F) - log2(self.granularity))
        assert all(
            -1 <= self.b_l[i] <= self.b_u[i] <= self.c[i]
            for i in range(len(self.c))
        ), self
        if not pss:
            assert -1 <= self.h_l < self.h <= self.h_u <= max(self.c), self
            assert -1 <= self.g_l < self.g <= self.g_u <= max(self.c), self
        else:
            assert -1 <= self.h_l < self.h <= self.h_u <= sum(self.c), self
            assert -1 <= self.g_l < self.g <= self.g_u <= sum(self.c), self
        # Assert that the true balances are inside F (as defined by the
        # current bounds)
        if not pss:
            b_inside_R_h_u = self.R_h_u.contains_point(self.b)
            b_inside_R_g_u = self.R_g_u.contains_point(self.b)
            b_inside_R_h_l = self.R_h_l.contains_point(self.b)
            b_inside_R_g_l = self.R_g_l.contains_point(self.b)
            # B must be within the upper bounds' rectangles
            assert b_inside_R_h_u, "\nB:\n" + "\n".join(
                [str(self.b), str(self.R_h_u)]
            )
            assert b_inside_R_g_u, "\nB:\n" + "\n".join(
                [str(self.b), str(self.R_g_u)]
            )
            # B must be outside the lower bounds' rectangles
            assert not b_inside_R_h_l, "\nB:\n" + "\n".join(
                [str(self.b), str(self.R_h_l)]
            )
            assert not b_inside_R_g_l, "\nB:\n" + "\n".join(
                [str(self.b), str(self.R_g_l)]
            )

        b_inside_R_b = self.R_b.contains_point(self.b)
        # B must be inside the current balance bounds rectangle
        assert b_inside_R_b, "\nB:\n" + "\n".join([str(self.b), str(self.R_b)])

    def reset_estimates(self, pss=False):
        """
        Set all variable hop parameters to their initial values. MUST be
        called on hop initialization and before running repeated probing
        on the same hops.
        """
        self.h_l = -1
        self.g_l = -1
        # NB: setting upper bound to max(self.c) (and not 0) if hop
        # can't forward is correct from the rectangle theory viewpoint
        if not pss:
            if self.can_forward(dir0):
                self.h_u = max(
                    [c for (i, c) in enumerate(self.c) if i in self.e[dir0]]
                )
            else:
                self.h_u = max(self.c)
        else:
            if self.can_forward(dir0):
                self.h_u = sum(
                    [c for (i, c) in enumerate(self.c) if i in self.e[dir0]]
                )
            else:
                self.h_u = sum(self.c)
        if not pss:
            if self.can_forward(dir1):
                self.g_u = max(
                    [c for (i, c) in enumerate(self.c) if i in self.e[dir1]]
                )
            else:
                self.g_u = max(self.c)
        else:
            if self.can_forward(dir1):
                self.g_u = sum(
                    [c for (i, c) in enumerate(self.c) if i in self.e[dir1]]
                )
            else:
                self.g_u = sum(self.c)

        self.b_l = [-1] * self.N
        self.b_u = [self.c[i] for i in range(len(self.c))]
        self.update_dependent_hop_properties(pss)

    def __str__(self):
        s = ""
        s += "Hop with properties:\n"
        s += "  channels: " + str(self.N) + "\n"
        s += "  capacities: " + str(self.c) + "\n"
        s += "  balances: " + str(self.b) + "\n"
        s += "  enabled in dir0: " + str(self.e[dir0]) + "\n"
        s += "  enabled in dir1: " + str(self.e[dir1]) + "\n"
        s += "  jammed in dir0: " + str(self.j[dir0]) + "\n"
        s += "  jammed in dir1: " + str(self.j[dir1]) + "\n"
        s += "  h if unjammed: " + str(self.h) + "\n"
        s += "  g if unjammed: " + str(self.g) + "\n"
        s += "  h_l: " + str(self.h_l) + "\n"
        s += "  h_u: " + str(self.h_u) + "\n"
        s += "  g_l: " + str(self.g_l) + "\n"
        s += "  g_u: " + str(self.g_u) + "\n"
        s += "  b_l: " + str(self.b_l) + "\n"
        s += "  b_u: " + str(self.b_u) + "\n"

        def effective_h(h):
            return h if self.can_forward(dir0) else 0

        def effective_g(g):
            return g if self.can_forward(dir1) else 0

        s += "Can forward in dir0 (effective h):\n"
        s += (
            "  "
            + str(effective_h(self.h_l + 1))
            + " -- "
            + str(effective_h(self.h_u))
            + "\n"
        )
        s += "Can forward in dir1 (effective g):\n"
        s += (
            "  "
            + str(effective_g(self.g_l + 1))
            + " -- "
            + str(effective_g(self.g_u))
            + "\n"
        )
        s += "Balance estimates:\n"
        s += (
            "  \n".join(
                [
                    "  " + str(self.b_l[i] + 1) + " -- " + str(self.b_u[i])
                    for i in range(len(self.b_l))
                ]
            )
            + "\n"
        )
        s += "Uncertainty: " + str(self.uncertainty) + "\n"
        return s

    def effective_vertex(self, direction, bound):
        """
        The coordinate of the _effective vertex_ corresponding to bound
        in direction is determined by:

        * the bound itself (bound = amount - 1);
        * if the i-th channel is enabled;
        * if the bound is lower than the i-th capacity.

        If bound <= i-th capacity and the i-th channel is enabled, the
        effective vertex' i-th coordinate equals the bound, otherwise it
        equals the i-th capacity.

        The corresponding ProbingRectangle is determined by the
        effective vertex and either [0, ... 0] (for dir0) or [c_1, ...
        c_N] (for dir1).

        Parameters:

        - direction: True if the bound corresponds to a probe in dir0,
          False otherwise
        - bound: equals to a - 1, where a is the probing amount
          (rationale: a probe of amount a cuts the points strictly less
          than a).

        Return:
        - eff_vertex: an N-element vector of coordinates of the
          effective vertex.
        """

        def effective_bound(bound, ch_i):
            # We're intentionally not accounting for jamming here. h and
            # g are "permanent" hop properties, assuming all channels
            # unjammed. For single-channel hops, h / g bounds are not
            # independent. Hence, it is sufficient for channel to be
            # enabled in one direction.
            if (
                ch_i in self.e[direction] or self.N == 1 or bound < 0
            ) and bound <= self.c[ch_i]:
                eff_bound = bound
            else:
                eff_bound = self.c[ch_i]
            # print("effective bound:", eff_bound)
            return eff_bound

        def effective_coordinate(bound, ch_i):
            return (
                effective_bound(bound, ch_i)
                if direction == dir0
                else self.c[ch_i] - effective_bound(bound, ch_i)
            )

        eff_vertex = [
            effective_coordinate(bound, ch_i) for ch_i in range(self.N)
        ]
        assert max(eff_vertex) <= max(self.c) + 1, (eff_vertex, max(self.c))
        # print("coordinates of effective vertex for bound = ", bound,
        # "in", ("dir0" if direction else "dir1"), ":", eff_vertex)
        return eff_vertex

    def S_F_generic(self, R_h_l, R_h_u, R_g_l, R_g_u, R_b):
        """
        Calculate S(F) determined by five rectangles (see paper for
        details). Note: the rectangles may not correspond to the current
        bounds.

        We use this function to calculate both: - the _actual_ S(F)
        (using the current rectangles from self. as arguments) - the
        _potential_ S(F) if we do a probe with amount a (when doing
        binary search for NBS amount selection)

        Parameters:
        - R_h_l: a rectangle defining the strict lower bound on h
        - R_h_u: a rectangle defining the non-strict upper bound on h
        - R_g_l: a rectangle defining the strict lower bound on g
        - R_g_u: a rectangle defining the non-strict upper bound on g
        - R_b : a rectangle defining the current knowledge about balance
          bounds (only if jamming)

        Return:
        - S_F: the number of points that:
          - belong to R_h_u and R_g_u and R_b
          - do NOT belong to R_h_l
          - do NOT belong to R_g_l

        The points in S_F are all possible positions of the true
        balances B.

        If jamming is used, we additionally intersect all rectangles
        with R_b. R_b reflects our current knowledge about individual
        balance bounds.

        """
        # Theoretically, we could intersect the final F with R_b, but we
        # can't do it easily because F may not be a rectangle. Instead,
        # we first intersect all four rectangles with R_b, and then
        # derive F as usual.
        R_u_u = R_h_u.intersect_with(R_g_u).intersect_with(R_b)
        R_u_l = R_h_u.intersect_with(R_g_l).intersect_with(R_b)
        R_l_u = R_h_l.intersect_with(R_g_u).intersect_with(R_b)
        R_l_l = R_h_l.intersect_with(R_g_l).intersect_with(R_b)
        """
        print("\nR_h_l:", R_h_l, R_h_l.S()) print("\nR_h_u:", R_h_u,
        R_h_u.S()) print("\nR_g_l:", R_g_l, R_g_l.S()) print("\nR_g_u:",
        R_g_u, R_g_u.S()) print("\nR_b:\n", R_b, R_b.S()) print("\nAfter
        intersecting with R_b:") print("\nR_u_u:", R_u_u, R_u_u.S())
        print("\nR_u_l:", R_u_l, R_u_l.S()) print("\nR_l_u:", R_l_u,
        R_l_u.S()) print("\nR_l_l:", R_l_l, R_l_l.S())
        """
        assert R_l_l.is_inside(R_u_u), self
        S_F = R_u_u.S() - R_u_l.S() - R_l_u.S() + R_l_l.S()
        # print(R_u_u.S(), "-", R_u_l.S(), "-", R_l_u.S(), "+",
        # R_l_l.S(), "=", S_F)
        assert S_F >= 0, self
        return S_F

    def S_F_generic_pss(self, h_l, h_u, g_l, g_u, R_b):
        """
        Return the count of the lattice points of the polytope defined
        by the the rectangle (or hyperbox in higher dimensions) R_b and
        the lower bound and higher bound of the forwarding abilities of
        the hop in both directions. h_l and g_l are strict, h_u and g_u
        non-strict.
        """

        if R_b.is_empty:
            return 0

        # LatteE is inclusive (non-strict) w.r.t. the polytope (the
        # lattice points in the facets count towards the sum of total
        # points contained), but h_l is strict. We make h_l non-strict:
        h_l += 1
        g_l += 1

        max_val_dir0 = sum(self.c[i] for i in self.e[dir0])
        max_val_dir1 = sum(self.c[i] for i in self.e[dir1])

        if h_l > max_val_dir0:
            raise ValueError(
                "Lower bound h_l cannot be larger than sum(capacities enabled in dir0) - 1"
            )

        if g_l > max_val_dir1:
            raise ValueError(
                "Lower bound g_l cannot be larger than sum(capacities enabled in dir1) - 1"
            )

        if h_u < 0 or g_u < 0:
            raise ValueError("Upper bound h_u or g_u cannot be smaller than 0")

        # If h_l and h_u are still at their initial values from hop
        # __init__, we don't need to do expensive lattice counts.
        if (
            h_l == 0
            and h_u == max_val_dir0
            and g_l == 0
            and g_u == max_val_dir1
        ):
            return R_b.S()

        # g_l and g_u are in dir1. We need to translate the values to dir0 values.
        g_l = max_val_dir1 - g_l
        g_u = max_val_dir1 - g_u

        # Count the number of dimensions.
        dimensions = len(R_b.l_vertex)

        # We wil describe the hyper-rectangle as a convex polytope,
        # using H-representation.

        # Create an identity array and its element-wise, numerical
        # negative. These are the vectors of the `m x n` matrix of the
        # H-representation, representing the facets of the
        # hyper-rectangle.
        id_arr = np.identity(dimensions, int)
        A = np.concatenate((id_arr, np.negative(id_arr)))

        # Add the vectors to the `m x n` matrix, representing the bounds

        # h_l:
        A = np.append(
            A,
            [
                list(
                    map(lambda i: 1 if i in self.e[dir0] else 0, range(self.N))
                )
            ],
            axis=0,
        )

        # h_u
        A = np.append(
            A,
            [
                list(
                    map(
                        lambda i: -1 if i in self.e[dir0] else 0, range(self.N)
                    )
                )
            ],
            axis=0,
        )

        # g_l:
        A = np.append(
            A,
            [
                list(
                    map(
                        lambda i: -1 if i in self.e[dir1] else 0, range(self.N)
                    )
                )
            ],
            axis=0,
        )
        # g_u
        A = np.append(
            A,
            [
                list(
                    map(lambda i: 1 if i in self.e[dir1] else 0, range(self.N))
                )
            ],
            axis=0,
        )

        # Create the `m x 1` matrix for the  H-representation.
        b = np.concatenate(
            (
                np.array([-x for x in R_b.l_vertex]),
                np.array(R_b.u_vertex),
            )
        )
        b = np.append(b, [-h_l])
        b = np.append(b, [h_u])
        b = np.append(b, [g_l])
        b = np.append(b, [-g_u])

        # Output A and b to file using LattE h-representation
        b = np.atleast_2d(b).T
        output = np.concatenate((b, A), axis=1)

        latte_file = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
        latte_file.writelines(
            " ".join([str(d) for d in np.shape(output)]) + "\n"
        )
        for row in output:
            latte_file.writelines(" ".join([str(x) for x in row]) + "\n")

        latte_file.flush()

        # Run LattE and capture output
        result = subprocess.run(
            ["../latte-distro/dest/bin/count", latte_file.name],
            capture_output=True,
            text=True,
        )

        # Close temporary file. This will automatically delete it.
        latte_file.close()

        if (
            result.returncode != 0
            and result.stderr.find("Empty polytope or unbounded polytope!")
            > -1
        ):
            return 0

        # Return the last line of stdout which contains the # of lattice points
        if len(result.stdout) > 0:
            return int(result.stdout.splitlines()[-1:][0])
        else:  # LattE sometimes gives back the result via stderr (don't ask why)
            last_line = result.stderr.splitlines()[-1:][0]
            last_word = last_line.split(" ")[-1:][0]
            return int(last_word.replace(".", ""))

    def S_F_a_expected(self, direction, a, pss=False):
        """
        Calculate the _potential_ S(F) if the probe of amount a fails
        ("area under the cut").

        Parameters:
        - direction: probe direction (dir0 / dir1)
        - a: the probe amount
        - pss: Payment Splitting and Switching assumed
        - success: mimic probe success (only for PSS)

        Return:
        - S_F_a: the number of points in S(F) "under the cut".
        """
        new_b_l = [0] * len(self.b_l)
        new_b_u = self.c.copy()
        # available channels are channels that are enabled and not
        # jammed
        available_channels = [
            i for i in self.e[direction] if i not in self.j[direction]
        ]
        jamming = len(self.j[dir0]) > 0 or len(self.j[dir1]) > 0
        if not pss:
            # mimic the scenario when probe fails
            if direction == dir0:
                new_R_h_u = (
                    self.R_h_u
                    if jamming
                    else ProbingRectangle(self, direction=dir0, bound=a - 1)
                )
                new_R_g_u = self.R_g_u
                for i in available_channels:
                    # probe failed => all available channels have
                    # insufficient balances
                    new_b_u[i] = min(new_b_u[i], a - 1)
            else:
                new_R_h_u = self.R_h_u
                new_R_g_u = (
                    self.R_g_u
                    if jamming
                    else ProbingRectangle(self, direction=dir1, bound=a - 1)
                )
                if len(available_channels) == 1:
                    # we can only update the lower bound if there is
                    # only one available channel and we know the probe
                    # went through this channel
                    new_b_l[available_channels[0]] = max(
                        new_b_l[available_channels[0]],
                        self.c[available_channels[0]] - a,
                    )
            new_R_b = Rectangle(new_b_l, new_b_u)
            S_F_a = self.S_F_generic(
                self.R_h_l, new_R_h_u, self.R_g_l, new_R_g_u, new_R_b
            )
        else:
            # PSS
            if direction == dir0:
                # new_h_u = a - 1
                # new_g_u = self.g_u
                for i in available_channels:
                    # probe failed => all available channels have
                    # insufficient balances
                    new_b_u[i] = min(new_b_u[i], a - 1)
                new_g_l = max(
                    sum(
                        self.c[j]
                        for j in available_channels
                        if j in self.e[not direction]
                    )
                    - a,
                    -1,
                )
                new_h_l = self.h_l
            else:
                # new_g_u = a - 1
                # new_h_u = self.h_u
                for i in available_channels:
                    # TODO: Seems to me you should use the real b_u
                    # here...
                    new_b_l[i] = max(
                        self.c[i]
                        - a
                        + sum(
                            self.c[j] - new_b_u[j]
                            for j in available_channels
                            if j != i
                        ),
                        -1,
                    )
                new_h_l = max(
                    sum(
                        self.c[j]
                        for j in available_channels
                        if j in self.e[not direction]
                    )
                    - a,
                    -1,
                )
                new_g_l = self.g_l
            new_R_b = Rectangle(new_b_l, new_b_u)
            S_F_a = self.S_F_generic_pss(
                new_h_l, self.h_u, new_g_l, self.g_u, new_R_b
            )
        # print("  expected area under the cut:", S_F_a, "(assuming
        # failed probe)")
        return S_F_a

    def worth_probing_h(self):
        # is there any uncertainty left about h that we resolve it
        # without jamming?
        # print(
        #     "Worth probing h (dir0)",
        #     str(self.can_forward(dir0) and self.h_u - self.h_l > 1),
        # )
        return self.can_forward(dir0) and self.h_u - self.h_l > 1

    def worth_probing_g(self):
        # is there any uncertainty left about g that we resolve it
        # without jamming?
        # print(
        #     "Worth probing g (dir1)",
        #     str(self.can_forward(dir1) and self.g_u - self.g_l > 1),
        # )
        return self.can_forward(dir1) and self.g_u - self.g_l > 1

    def worth_probing_h_or_g(self, direction):
        return (
            self.worth_probing_h()
            if direction == dir0
            else self.worth_probing_g()
        )

    def worth_probing_channel(self, i):
        # is it worth doing jamming-enhanced probing on this channel?
        return self.b_u[i] - self.b_l[i] > 1 and (
            self.can_forward(dir0) or self.can_forward(dir1)
        )

    def worth_probing(self):
        # is there any uncertainty left in the hop?
        return self.uncertainty > 0

    def next_a(self, direction, bs, jamming, pss=False):
        """
        Calculate the optimal (NBS) amount for probe in direction. The
        NBS amount shrinks S(F) by half. (In other words, the probe
        leaves S_F/2 under the cut.) We look for the NBS amount using a
        binary search: starting from the current bounds in the required
        direction, we choose a in the middle. We then check the area
        under the cut _if_ we probed with this amount. Depending on if
        S_F_a < S_F or S_F_a > S_F, we increase / decrease a.

        Parameters:
        - direction: dir0 or dir1

        Return:
        - a: the NBS amount, or None if the hop cannot forward
          in this direction
        """

        # def new_a(direction, start_a, a_l, a_u, pss, success=False):
        #     while True:
        #         S_F_a = self.S_F_a_expected(direction, start_a, pss, success)
        #         # print(a_l, a, a_u)
        #         if success:
        #             if S_F_a < S_F_half:
        #                 a_u = start_a
        #             else:
        #                 a_l = start_a
        #             new_a = (a_l + a_u + 1) // 2
        #         else:
        #             if S_F_a < S_F_half:
        #                 a_l = start_a
        #             else:
        #                 a_u = start_a
        #             new_a = (a_l + a_u + 1) // 2
        #         if new_a == start_a:
        #             break
        #         start_a = new_a
        #     return new_a

        S_F_half = max(1, self.S_F // 2)
        if not jamming:
            # only makes sense to send probes between current estimates
            # on h or g
            a_l, a_u = (
                (self.h_l + 1, self.h_u)
                if direction == dir0
                else (self.g_l + 1, self.g_u)
            )
        else:
            # individual balance bounds may be outside bounds for h / g
            # (those are bounds for maximums!)
            available_channels = [
                i for i in self.e[direction] if i not in self.j[direction]
            ]
            assert (
                len(available_channels) == 1
            ), "We only support probing one unjammed channel at a time"
            i = available_channels[0]
            a_l, a_u = (
                (self.b_l[i] + 1, self.b_u[i])
                if direction == dir0
                else (self.c[i] - self.b_u[i], self.c[i] - self.b_l[i] - 1)
            )
        a = (a_l + a_u + 1) // 2
        print("Binary Search amount: " + str(a))
        if not bs and not jamming:
            # we only do binary search over S(F) in pre-jamming probing
            # phase
            while True:
                S_F_a = self.S_F_a_expected(direction, a, pss)
                if S_F_a < S_F_half:
                    a_l = a
                else:
                    a_u = a
                new_a = (a_l + a_u + 1) // 2
                if new_a == a:
                    break
                a = new_a
            print("Optimal amount: " + str(a))
            # a = new_a(direction, a, a_l, a_u, pss, success)
            # if pss and (
            #     (a_u - a <= 1 and success) or (a - a_l <= 1 and not success)
            # ):
            #     # PSS amount can't give a lot of information, probably
            #     # because we assumed a fail, and the previous amount
            #     # succeeded or the other way around. Let's try assuming
            #     # a the opposite.
            #     a = new_a(direction, a, a_l, a_u, pss, not success)
            #     if (a_u - a <= 1 and not success) or (
            #         a - a_l <= 1 and success
            #     ):
            #         # PSS amount can't give a lot of information, either
            #         # way. Let's default back to binary search.
            #         a = (a_l + a_u + 1) // 2
            # # print("if a = ", a, ", then area under cut = ", S_F_a,
            # # ", need", S_F_half)
        assert a > 0
        return a

    def next_dir(
        self,
        bs,
        jamming,
        prefer_small_amounts=False,
        threshold_area_difference=0.1,
        pss=False,
    ):
        """
        Suggest the NBS direction for the next probe.

        Parameters:
        - bs: True if we do binary search only amounts;
          False if we use NBS amount choice
        - jamming: are we doing jamming-enhanced probing after regular
          probing
        - prefer_small_amounts: prioritize small amounts vs cutting S(F)
          in half more precisely
        - threshold_area_difference: the difference in S(F) that we
          neglect when choosing between two amounts
        - pss: assume Payment Split & Switch

        Return:
        - chosen_dir: the suggested direction
        """
        assert self.can_forward(dir0) or self.can_forward(dir1), self
        should_consider_dir0 = (
            jamming
            and self.can_forward(dir0)
            or not jamming
            and self.worth_probing_h()
        )
        should_consider_dir1 = (
            jamming
            and self.can_forward(dir1)
            or not jamming
            and self.worth_probing_g()
        )
        if not should_consider_dir0:
            chosen_dir = dir1
        elif not should_consider_dir1:
            chosen_dir = dir0
        else:
            # next_a for the case the probe fails
            a_dir0 = self.next_a(dir0, bs, jamming, pss)
            a_dir1 = self.next_a(dir1, bs, jamming, pss)

            if bs or prefer_small_amounts:
                # choose smaller amount: more likely to pass

                # if pss:
                #     a_smallest = min(
                #         a_dir0, a_dir1, a_dir0_success, a_dir1_success
                #     )
                #     switch = {
                #         a_dir0: dir0,
                #         a_dir1: dir1,
                #         a_dir0_success: dir0,
                #         a_dir1_success: dir1,
                #     }
                #     chosen_dir = switch.get(a_smallest)
                # else:
                #     chosen_dir = dir0 if a_dir0 < a_dir1 else dir1

                chosen_dir = dir0 if a_dir0 < a_dir1 else dir1

            else:
                # prefer amount that splits in half better
                S_F_half = max(1, self.S_F // 2)
                S_F_a_dir0 = self.S_F_a_expected(dir0, a_dir0, pss)
                S_F_a_dir1 = self.S_F_a_expected(dir1, a_dir1, pss)

                if (
                    abs(S_F_a_dir0 - S_F_a_dir1) / S_F_half
                    < threshold_area_difference
                ):
                    chosen_dir = dir0 if a_dir0 < a_dir1 else dir1
                else:
                    chosen_dir = (
                        dir0
                        if abs(S_F_a_dir0 - S_F_half)
                        < abs(S_F_a_dir1 - S_F_half)
                        else dir1
                    )
                # if pss:
                #     S_F_a_dir0_success = self.S_F_a_expected(
                #         dir0, a_dir0_success, pss, True
                #     )
                #     S_F_a_dir1_success = self.S_F_a_expected(
                #         dir1, a_dir1_success, pss, True
                #     )

                # if (
                #     abs(S_F_a_dir0 - S_F_a_dir1) / S_F_half
                #     < threshold_area_difference
                # ):
                #     chosen_dir = dir0 if a_dir0 < a_dir1 else dir1
                # elif (
                #     pss
                #  and abs(S_F_a_dir0_success - S_F_a_dir1_success) / S_F_half
                #     < threshold_area_difference
                # ):
                #     chosen_dir = (
                #         dir0 if a_dir0_success < a_dir1_success else dir1
                #     )
                # else:
                #     if pss:
                #         S_F_a_smallest = min(
                #             S_F_a_dir0 - S_F_half,
                #             S_F_a_dir1 - S_F_half,
                #             S_F_a_dir0_success - S_F_half,
                #             S_F_a_dir1_success - S_F_half,
                #         )
                #         switch = {
                #             S_F_a_dir0 - S_F_half: dir0,
                #             S_F_a_dir1 - S_F_half: dir1,
                #             S_F_a_dir0_success - S_F_half: dir0,
                #             S_F_a_dir1_success - S_F_half: dir1,
                #         }
                #         chosen_dir = switch.get(S_F_a_smallest)
                #     else:
                #         chosen_dir = (
                #             dir0
                #             if abs(S_F_a_dir0 - S_F_half)
                #             < abs(S_F_a_dir1 - S_F_half)
                #             else dir1
                #         )
        return chosen_dir

    def probe(self, direction, amount, pss=False):
        """
        Update the bounds as a result of a probe.

        Parameters:
        - direction: probe direction (dir0 or dir1)
        - amount: probe amount

        Return:
        - None (the current bounds are updated)
        """
        # print("doing probe", amount, "in", "dir0" if direction else
        # "dir1")
        jamming = len(self.j[dir0]) > 0 or len(self.j[dir1]) > 0
        # print("Are we jamming?", jamming)
        available_channels = [
            i for i in self.e[direction] if i not in self.j[direction]
        ]
        if jamming:
            # if we're jamming, we must jam all channels except one
            assert len(available_channels) <= 1

        def b_in_dir(i, direction):
            return self.b[i] if direction == dir0 else self.c[i] - self.b[i]

        probe_passed = (
            amount <= max(b_in_dir(i, direction) for i in available_channels)
            if not pss
            else amount
            <= sum(b_in_dir(i, direction) for i in available_channels)
        )
        if direction == dir0:
            # should only update if the amount is between current bounds
            # this is not always true for intermediary hops
            should_update_h = self.h_l < amount <= self.h_u
            if probe_passed:
                # print("probe passed in dir0") sic! lower bounds are
                # strict
                if should_update_h and not jamming:
                    # update hop-level lower bound
                    self.h_l = amount - 1
                    if not pss:
                        if len(self.e[dir0]) == 1:
                            # if only one channel is enabled, we can
                            # update this channel's lower bound

                            # FIXME: This assumption is incorrect. If
                            # only one channel with enough capacity
                            # (b_u) is enabled, you can update the
                            # bound.
                            self.b_l[self.e[dir0][0]] = max(
                                self.b_l[self.e[dir0][0]], self.h_l
                            )
                        if len(self.e[dir1]) > 0:
                            # if some channels are enabled in the
                            # opposite direction, update that upper
                            # bound
                            self.g_u = min(
                                self.g_u,
                                max(
                                    self.c[i] - self.b_l[i] - 1
                                    for i in self.e[dir1]
                                ),
                            )
                    else:
                        # PSS: Update a channel's lower bound if the
                        # amount is higher than the combined b_u of the
                        # other channels in this hop that are enabled in
                        # dir0
                        for i in self.e[dir0]:
                            self.b_l[i] = max(
                                self.b_l[i],
                                max(
                                    amount
                                    - sum(
                                        self.b_u[j]
                                        for j in self.e[dir0]
                                        if j != i
                                    )
                                    - 1,
                                    -1,
                                ),
                            )
                        if len(self.e[dir1]) > 0:
                            # if some channels are enabled in the
                            # opposite direction, update that upper
                            # bound
                            if all(
                                [chan in self.e[dir1] for chan in self.e[dir0]]
                            ):
                                # All channels enabled in dir0 are also
                                # enabled in dir1. We can use h_l to
                                # update g_u.
                                self.g_u = min(
                                    self.g_u,
                                    sum(self.c[i] for i in self.e[dir1])
                                    - self.h_l
                                    - 1,
                                )
                            else:
                                # e[dir0] is not a subset of e[dir1] so
                                # the best we can do is to update g_u
                                # with the b_l estimates
                                self.g_u = min(
                                    self.g_u,
                                    sum(
                                        self.c[i] - self.b_l[i] - 1
                                        for i in self.e[dir1]
                                    ),
                                )
                if jamming:
                    # if we're jamming, we can update the only unjammed
                    # channel's lower bound
                    self.b_l[available_channels[0]] = max(
                        self.b_l[available_channels[0]], amount - 1
                    )
            else:
                # print("probe failed in dir0")
                if should_update_h and not jamming:
                    # update hop-level upper bound
                    self.h_u = amount - 1
                    for i in self.e[dir0]:
                        if not pss:
                            # update all channels' upper bounds
                            self.b_u[i] = min(self.b_u[i], self.h_u)
                        else:
                            self.b_u[i] = min(
                                self.b_u[i],
                                min(
                                    (
                                        amount
                                        - sum(
                                            self.b_l[j] + 1
                                            for j in self.e[dir0]
                                            if j != i
                                        )
                                        - 1,
                                        self.c[i],
                                    )
                                ),
                            )
                    if len(self.e[dir1]) > 0:
                        if not pss:
                            # if some channels are enabled in the
                            # opposite direction, update their lower
                            # bound
                            self.g_l = max(
                                self.g_l,
                                min(
                                    self.c[i] - self.b_u[i] - 1
                                    for i in self.e[dir1]
                                ),
                            )
                        else:
                            # If we calculate [b_l_dir1 plus the
                            # b_u_dir1 of the other channels] for each
                            # channel involved in the payment in dir0
                            # that is also enabled in dir1, than g_l is
                            # the lowest result of that set.(We don't
                            # keep track of b_l_dir1, we only keep track
                            # of b_l and b_u wich are assumed to be in
                            # dir0). b_l_dir1 = c - b_u - 1 b_u_dir1 = c
                            # - b_l - 1

                            # if True in (
                            #     self.c[i] - self.b_u[i] > 0
                            #     for i in self.e[dir0]
                            #     if i in self.e[dir1]
                            # ):
                            #     self.g_l = max(
                            #         self.g_l,
                            #         min(
                            #             self.c[i]
                            #             - self.b_u[i]
                            #             - 1
                            #             + sum(
                            #                 self.c[j] - self.b_l[j] - 1
                            #                 for j in self.e[dir0]
                            #                 if j != i and j in self.e[dir1]
                            #             )
                            #             for i in self.e[dir0]
                            #             if i in self.e[dir1]
                            #         ),
                            #     )
                            self.g_l = max(
                                self.g_l,
                                sum(
                                    self.c[i]
                                    for i in self.e[dir0]
                                    if i in self.e[dir1]
                                )
                                - amount,
                            )
                if jamming:
                    # if we're jamming, we can update the only unjammed
                    # channel's upper bound
                    self.b_u[available_channels[0]] = min(
                        self.b_u[available_channels[0]], amount - 1
                    )
        else:
            should_update_g = self.g_l < amount <= self.g_u
            if probe_passed:
                # print("probe passed in dir1")
                if should_update_g and not jamming:
                    self.g_l = amount - 1
                    if not pss:
                        if len(self.e[dir1]) == 1:
                            self.b_u[self.e[dir1][0]] = min(
                                self.b_u[self.e[dir1][0]],
                                self.c[self.e[dir1][0]] - self.g_l - 1,
                            )
                        if len(self.e[dir0]) > 0:
                            self.h_u = min(
                                self.h_u,
                                max(self.b_u[i] for i in self.e[dir0]),
                            )
                    else:
                        # PSS: Update a channel's lower bound in dir1
                        # (b_l_dir1) if the amount is higher than the
                        # combined b_u_dir1 of the other channels in
                        # this hop that are enabled in dir1. Keep in
                        # mind we don't keep track of b_l_dir1 and
                        # b_u_dir1, but:
                        # b_l_dir1 = c - b_u - 1
                        # b_u_dir1 = c - b_l - 1
                        for i in self.e[dir1]:
                            self.b_u[i] = min(
                                self.b_u[i],
                                min(
                                    sum(
                                        self.c[j] - self.b_l[j] - 1
                                        for j in self.e[dir1]
                                        if j != i
                                    )
                                    + self.c[i]
                                    - amount,
                                    self.c[i],
                                ),
                            )
                        if len(self.e[dir0]) > 0:
                            if all(
                                [chan in self.e[dir0] for chan in self.e[dir1]]
                            ):
                                # All channels enabled in dir1 are also
                                # enabled in dir0. We can use g_l to
                                # update h_u.
                                self.h_u = min(
                                    self.h_u,
                                    sum(self.c[i] for i in self.e[dir0])
                                    - self.g_l
                                    - 1,
                                )
                            else:
                                # e[dir1] is not a subset of e[dir0] so
                                # the best we can do is to update h_u
                                # with the b_l estimates
                                self.h_u = min(
                                    self.h_u,
                                    sum(self.b_u[i] for i in self.e[dir0]),
                                )
                if jamming:
                    self.b_u[available_channels[0]] = min(
                        self.b_u[available_channels[0]],
                        self.c[available_channels[0]] - amount,
                    )
            else:
                # print("probe failed in dir1")
                if should_update_g and not jamming:
                    self.g_u = amount - 1
                    for i in self.e[dir1]:
                        if not pss:
                            self.b_l[i] = max(
                                self.b_l[i], self.c[i] - self.g_u - 1
                            )
                        else:
                            self.b_l[i] = max(
                                self.b_l[i],
                                max(
                                    self.c[i]
                                    + sum(
                                        self.c[j] - self.b_u[j]
                                        for j in self.e[dir1]
                                        if j != i
                                    )
                                    - amount,
                                    -1,
                                ),
                            )
                    if len(self.e[dir0]) > 0:
                        if not pss:
                            self.h_l = max(
                                self.h_l,
                                min(self.b_l[i] for i in self.e[dir0]),
                            )
                        else:
                            # If we calculate [b_l plus the b_u of the
                            # other channels] for each channel involved
                            # in the payment in dir1 that is also
                            # enabled in dir0, than h_l is the lowest
                            # result of that set.

                            # if True in (
                            #     self.b_l[i] > -1
                            #     for i in self.e[dir1]
                            #     if i in self.e[dir0]
                            # ):
                            #     self.h_l = max(
                            #         self.h_l,
                            #         min(
                            #             self.b_l[i]
                            #             + sum(
                            #                 self.b_u[j]
                            #                 for j in self.e[dir1]
                            #                 if j != i and j in self.e[dir0]
                            #             )
                            #             for i in self.e[dir1]
                            #             if i in self.e[dir0]
                            #         ),
                            #     )
                            self.h_l = max(
                                self.h_l,
                                sum(
                                    self.c[i]
                                    for i in self.e[dir1]
                                    if i in self.e[dir0]
                                )
                                - amount,
                            )

                if jamming:
                    self.b_l[available_channels[0]] = max(
                        self.b_l[available_channels[0]],
                        self.c[available_channels[0]] - amount,
                    )
        # print("after probe:", self.h_l, self.h_u, self.g_l, self.g_u)
        self.update_dependent_hop_properties(pss)
        if self.uncertainty == 0 and not pss:
            corner_points = self.get_corner_points()
            assert len(corner_points) <= 1
            if len(corner_points) == 1:
                p = corner_points[0]
                for i in range(self.N):
                    if self.worth_probing_channel(i):
                        self.b_l[i] = p[i] - 1
                        self.b_u[i] = p[i]
                if len(self.e[dir0]) > 0:
                    self.h_l = max(
                        self.h_l, min([p[i] for i in self.e[dir0]]) - 1
                    )
                    self.h_u = min(self.h_u, max([p[i] for i in self.e[dir0]]))
                if len(self.e[dir1]):
                    self.g_l = max(
                        self.g_l,
                        min([self.c[i] - p[i] for i in self.e[dir1]]) - 1,
                    )
                    self.g_u = min(
                        self.g_u, max([self.c[i] - p[i] for i in self.e[dir1]])
                    )
            else:
                print("Corners are not viable points, continue probing")
                pass
            self.update_dependent_hop_properties(pss)
        return probe_passed

    def analysis(self):
        """
        Returns what a perfect probe should return, based on the actual
        balances, h- and g-values. Where a probe doesn't have access to
        these values (we probe to *detect* these values!) the analysis
        does have access to these values. So we can calculate the exact
        data a perfect probe should return. The analysis assumes h,g and
        c are correct.
        """
        hop_info = dict()
        hop_info["h_l"] = self.h - 1
        hop_info["h_u"] = self.h
        hop_info["g_l"] = self.g - 1
        hop_info["g_u"] = self.g
        hop_info["b_l_dir0"] = [-1] * self.N
        hop_info["b_u_dir0"] = self.c.copy()
        hop_info["b_l_dir1"] = [-1] * self.N
        hop_info["b_u_dir1"] = self.c.copy()
        hop_info["b_l"] = [-1] * self.N
        hop_info["b_u"] = self.c.copy()
        hop_info["uncertainty"] = None

        def update_b_l_b_u(b_l, b_u, l, u, direction):
            if self.pss:
                largest_succesful_probe = l + 1
                smallest_failed_probe = u + 1
                # Use u to calculate b_u. h is the smallest failed probe
                # in the given direction.
                for i in self.e[direction]:
                    b_u[i] = min(
                        smallest_failed_probe
                        - sum(b_l[j] + 1 for j in self.e[direction] if j != i)
                        - 1,
                        self.c[i],
                        b_u[i],
                    )

                # Use l to calculate b_l. l + 1 is the largest
                # successful probe in the given direction.
                for i in self.e[direction]:
                    b_l[i] = max(
                        largest_succesful_probe
                        - sum(b_u[j] for j in self.e[direction] if j != i)
                        - 1,
                        -1,
                        b_l[i],
                    )
            else:
                if len(self.e[direction]) == 1:
                    # if only one channel is enabled, we can
                    # update this channel's lower bound

                    # FIXME: This assumption is incorrect. If
                    # only one channel with enough capacity
                    # (b_u) is enabled, you can update the
                    # bound.
                    b_l[0] = max(l, b_l[0])
                for i in self.e[direction]:
                    b_u[i] = min(b_u[i], u)

        def get_corner_points():
            R_u_u = R_h_u.intersect_with(R_g_u).intersect_with(R_b)
            R_u_l = R_h_u.intersect_with(R_g_l).intersect_with(R_b)
            R_l_u = R_h_l.intersect_with(R_g_u).intersect_with(R_b)
            ranges = [
                [R_u_u.l_vertex[i], R_u_u.u_vertex[i]]
                for i in range(len(R_u_u.l_vertex))
            ]
            points = []
            points_left = S_F
            for p in product(*ranges):
                if not R_u_l.contains_point(p) and not R_l_u.contains_point(p):
                    points.append(p)
                    points_left -= 1
                if points_left == 0:
                    break
            return points

        def worth_probing_channel(i):
            # is it worth doing jamming-enhanced probing on this channel?
            return hop_info["b_u"][i] - hop_info["b_l"][i] > 1 and (
                self.can_forward(dir0) or self.can_forward(dir1)
            )

        update_b_l_b_u(
            hop_info["b_l_dir0"],
            hop_info["b_u_dir0"],
            hop_info["h_l"],
            hop_info["h_u"],
            dir0,
        )

        for i, c in enumerate(self.c):
            hop_info["b_l_dir1"][i] = max(
                hop_info["b_l_dir1"][i], c - hop_info["b_u_dir0"][i] - 1
            )
            hop_info["b_u_dir1"][i] = min(
                hop_info["b_u_dir1"][i], c - hop_info["b_l_dir0"][i] - 1
            )

        update_b_l_b_u(
            hop_info["b_l_dir1"],
            hop_info["b_u_dir1"],
            hop_info["g_l"],
            hop_info["g_u"],
            dir1,
        )

        for i, c in enumerate(self.c):
            hop_info["b_l_dir0"][i] = max(
                hop_info["b_l_dir0"][i], c - hop_info["b_u_dir1"][i] - 1
            )
            hop_info["b_u_dir0"][i] = min(
                hop_info["b_u_dir0"][i], c - hop_info["b_l_dir1"][i] - 1
            )

        update_b_l_b_u(
            hop_info["b_l_dir0"],
            hop_info["b_u_dir0"],
            hop_info["h_l"],
            hop_info["h_u"],
            dir0,
        )

        for i, c in enumerate(self.c):
            hop_info["b_l_dir1"][i] = max(
                hop_info["b_l_dir1"][i], c - hop_info["b_u_dir0"][i] - 1
            )
            hop_info["b_u_dir1"][i] = min(
                hop_info["b_u_dir1"][i], c - hop_info["b_l_dir0"][i] - 1
            )

        for i, c in enumerate(self.c):
            hop_info["b_l"][i] = max(
                hop_info["b_l_dir0"][i], c - hop_info["b_u_dir1"][i] - 1
            )
            hop_info["b_u"][i] = min(
                hop_info["b_u_dir0"][i], c - hop_info["b_l_dir1"][i] - 1
            )

        R_b = Rectangle(
            [b_l_i + 1 for b_l_i in hop_info["b_l"]], hop_info["b_u"]
        )
        if not self.pss:
            R_h_l = ProbingRectangle(
                self, direction=dir0, bound=hop_info["h_l"]
            )
            R_h_u = ProbingRectangle(
                self, direction=dir0, bound=hop_info["h_u"]
            )
            R_g_l = ProbingRectangle(
                self, direction=dir1, bound=hop_info["g_l"]
            )
            R_g_u = ProbingRectangle(
                self, direction=dir1, bound=hop_info["g_u"]
            )
            S_F = self.S_F_generic(R_h_l, R_h_u, R_g_l, R_g_u, R_b)
        else:
            S_F = self.S_F_generic_pss(
                hop_info["h_l"],
                hop_info["h_u"],
                hop_info["g_l"],
                hop_info["g_u"],
                R_b,
            )

        hop_info["uncertainty"] = max(0, log2(S_F) - log2(self.granularity))

        if hop_info["uncertainty"] == 0 and not self.pss:
            corner_points = get_corner_points()
            assert len(corner_points) <= 1
            if len(corner_points) == 1:
                p = corner_points[0]
                for i in range(self.N):
                    if worth_probing_channel(i):
                        hop_info["b_l"][i] = p[i] - 1
                        hop_info["b_u"][i] = p[i]
                if len(self.e[dir0]) > 0:
                    hop_info["h_l"] = max(
                        hop_info["h_l"], min([p[i] for i in self.e[dir0]]) - 1
                    )
                    hop_info["h_u"] = min(
                        hop_info["h_u"], max([p[i] for i in self.e[dir0]])
                    )
                if len(self.e[dir1]):
                    hop_info["g_l"] = max(
                        hop_info["g_l"],
                        min([self.c[i] - p[i] for i in self.e[dir1]]) - 1,
                    )
                    hop_info["g_u"] = min(
                        hop_info["g_u"],
                        max([self.c[i] - p[i] for i in self.e[dir1]]),
                    )

        return hop_info
