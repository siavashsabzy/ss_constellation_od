import numpy as np
from PIL import Image as PImage
import plotly.graph_objects as go
from plotly import tools
import plotly.offline
def create_rgb_surface(rgb_img, depth_img, depth_cutoff=20, **kwargs):
    rgb_img = rgb_img.swapaxes(0, 1)[:, ::-1]
    depth_img = depth_img.swapaxes(0, 1)[:, ::-1]
    eight_bit_img = PImage.fromarray(rgb_img).convert('P', palette='WEB', dither=None)
    idx_to_color = np.array(eight_bit_img.getpalette()).reshape((-1, 3))
    colorscale=[[i/255.0, "rgb({}, {}, {})".format(*rgb)] for i, rgb in enumerate(idx_to_color)]
    depth_map = depth_img.copy().astype('float')
    depth_map[depth_map<depth_cutoff] = np.nan
    return go.Surface(
        z=depth_map,
        surfacecolor=np.array(eight_bit_img),
        cmin=0, 
        cmax=255,
        colorscale=colorscale,
        showscale=False,
        **kwargs
    )

