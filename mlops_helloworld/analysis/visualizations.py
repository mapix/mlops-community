from typing import List


def draw_multiple_distribution(data: List[List[float]], systems: List[str]):
    import plotly.figure_factory as ff
    imgfile = "histogram.png"

    if len(data) > 0:
        fig = ff.create_distplot(data, systems)
        img = fig.write_image(imgfile)
    return fig, imgfile
