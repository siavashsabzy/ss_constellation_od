from astropy import units as u

from poliastro.earth import EarthSatellite
from poliastro.earth.plotting import GroundtrackPlotter
from poliastro.examples import iss
from poliastro.util import time_range

import plotly

from PySide6.QtWidgets import QApplication
from PySide6.QtWebEngineWidgets import QWebEngineView


def build_plot(fig):
    html = "".join(
        [
            "<html><body>",
            plotly.offline.plot(fig, output_type="div", include_plotlyjs="cdn"),
            "</body></html>",
        ]
    )
    return html


def create_plot():

    satellite_spacecraft = EarthSatellite(iss, None)
    t_span = time_range(iss.epoch - 1.5 * u.h, periods=150, end=iss.epoch + 1.5 * u.h)

    gp = GroundtrackPlotter()
    gp.update_layout(title="International Space Station groundtrack")

    # Plot previously defined EarthSatellite object
    gp.plot(
        satellite_spacecraft,
        t_span,
        label="Satellite",
        color="red",
        marker={"size": 10, "symbol": "triangle-right", "line": {"width": 1, "color": "black"}},
    )
    # For building geo traces
    import plotly.graph_objects as go

    # Faculty of Radiophysics and Computer Technologies coordinates
    STATION = [53.83821551524637, 27.476136409973797] * u.deg

    # Let us add a new trace in original figure
    gp.add_trace(
        go.Scattergeo(
            lat=STATION[0],
            lon=STATION[-1],
            name="Faculty of Radiophysics and Computer Technologies",
            marker={"color": "blue"},
        )
    )
    # Switch to three dimensional representation
    gp.update_geos(projection_type="orthographic")
    return gp.fig


def main():
    app = QApplication([])

    fig = create_plot()
    html = build_plot(fig)

    view = QWebEngineView()
    view.setHtml(html)
    view.resize(640, 480)
    view.show()

    app.exec()


if __name__ == "__main__":
    main()