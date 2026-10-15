import unittest

import matplotlib.pyplot as plt

from plot_surface import build_dielectric_control_figure


class DielectricPlotLayoutTests(unittest.TestCase):
    def test_note_is_below_x_axis_label_without_overlap(self):
        figure, axis, note, _ = build_dielectric_control_figure()
        try:
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            label_box = axis.xaxis.label.get_window_extent(renderer)
            note_box = note.get_window_extent(renderer)
            self.assertGreater(label_box.y0 - note_box.y1, 4.0)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
