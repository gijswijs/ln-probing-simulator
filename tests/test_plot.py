from plot import plot

NUM_CHANNELS_IN_TARGET_HOPS = [1, 2, 3]
Y_GAINS_LINES_VANILLA = [
    (
        [
            [1.0, 1.0],
            [0.7369380933924997, 0.8801330173326367],
            [0.6005750999438223, 0.5879422038536971],
        ],
        "Direct probing, non-pss",
        (0, ()),
        "blue",
    ),
    (
        [
            [1.0, 1.0],
            [0.9044341897216104, 0.6034647320602335],
            [0.47769200190624084, 0.6138007690327479],
        ],
        "Remote probing, non-pss",
        (0, (3, 1)),
        "red",
    ),
]
Y_GAINS_LINES_JAMMING = [
    (
        [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
        "Direct probing, non-pss",
        (0, ()),
        "blue",
    ),
    (
        [
            [0.8019811589232457, 0.7642346116526579],
            [0.8819068617882873, 1.0],
            [1.0, 1.0],
        ],
        "Remote probing, non-pss",
        (0, (3, 1)),
        "red",
    ),
]
Y_SPEED_LINES_VANILLA = [
    (
        [
            [1.0, 1.0],
            [0.7369380933924997, 0.8801330173326367],
            [0.6005750999438223, 0.5879422038536971],
        ],
        "Direct probing, non-pss",
        (0, ()),
        "blue",
    ),
    (
        [
            [1.0, 1.0],
            [0.9044341897216104, 0.6034647320602335],
            [0.47769200190624084, 0.6138007690327479],
        ],
        "Remote probing, non-pss",
        (0, (3, 1)),
        "red",
    ),
]
Y_SPEED_LINES_JAMMING = [
    (
        [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
        "Direct probing, non-pss",
        (0, ()),
        "blue",
    ),
    (
        [
            [0.8019811589232457, 0.7642346116526579],
            [0.8819068617882873, 1.0],
            [1.0, 1.0],
        ],
        "Remote probing, non-pss",
        (0, (3, 1)),
        "red",
    ),
]
X_LABEL = "\nNumber of channels in target hops\n"
FILENAME_GAINS = "gains_snapshot_test"
FILENAME_SPEED = "speed_snapshot_test"


def test_plot_gains():
    plot(
        x_data=NUM_CHANNELS_IN_TARGET_HOPS,
        y_data_lists=[Y_GAINS_LINES_VANILLA, Y_GAINS_LINES_JAMMING],
        x_label=X_LABEL,
        y_label="Information gain (share of initial uncertainty)\n",
        title="",
        filename=FILENAME_GAINS,
    )
    assert True


def test_plot_speed():
    plot(
        x_data=NUM_CHANNELS_IN_TARGET_HOPS,
        y_data_lists=[Y_SPEED_LINES_VANILLA, Y_SPEED_LINES_JAMMING],
        x_label=X_LABEL,
        y_label="Probing speed (bits / message)\n",
        title="",
        filename=FILENAME_SPEED,
    )
    assert True
