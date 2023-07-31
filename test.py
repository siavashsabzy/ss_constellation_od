import pandas as pd
from bs4 import BeautifulSoup
import plotly.express as px
import numpy as np
from PIL import Image
import plotly.io as pio
from datetime import datetime
import plotly.graph_objs as go
import webbrowser
webbrowser.register('firefox',None,webbrowser.BackgroundBrowser("C:\\Program Files\\Mozilla Firefox\\firefox.exe"))
pio.renderers.default = "firefox"
def sphere(size, texture): 
    N_lat = int(texture.shape[0])
    N_lon = int(texture.shape[1])
    theta = np.linspace(0, 2*np.pi, N_lat)
    phi = np.linspace(0, np.pi, N_lon)
    
    # Set up coordinates for points on the sphere
    x0 = size * np.outer(np.cos(theta), np.sin(phi))
    y0 = size * np.outer(np.sin(theta), np.sin(phi))
    z0 = size * np.outer(np.ones(N_lat), np.cos(phi))
    
    # Set up trace
    return x0, y0, z0

texture = np.asarray(Image.open('earthmap1k.jpg'))
x, y, z = sphere(6378137.000, texture)
surf = go.Surface(x=x, y=y, z=z, surfacecolor=texture)

with open('S1A_OPER_AUX_RESORB_OPOD_20220103T002736_V20220102T202842_20220102T234612.EOF', 'r') as f:
    data = f.read()

Bs_data = BeautifulSoup(data, "xml")
b_unique = Bs_data.find_all('OSV')
position = pd.DataFrame(columns=['X', 'Y', 'Z'])

for i in range(len(b_unique)):
    s = b_unique[i].find_all("UTC")[0].contents[0]
    x = b_unique[i].find_all("X")[0].contents[0]
    y = b_unique[i].find_all("Y")[0].contents[0]
    z = b_unique[i].find_all("Z")[0].contents[0]
    a = datetime.strptime(s, "UTC=%Y-%m-%dT%H:%M:%S.%f")
    position.loc[a] = [float(x), float(y), float(z)]

fig = go.Figure(data=[surf])

fig.add_trace(
    go.Scatter3d(
        x=position['X'],
        y=position['Y'],
        z=position['Z'],
        mode='markers',
        marker=dict(
            size=2,
            color='red',
            opacity=0.8
        )
    )
)

fig.update_layout(scene=dict(aspectratio=dict(x=1, y=1, z=1)))

pio.show(fig)

