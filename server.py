def address(ticker, back='2024-01-01'):
    key = ''
    url = f'https://financialmodelingprep.com/api/v3/historical-price-full/{ticker}?from={back}&apikey={key}'
    return url

import requests
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ctypes

def Calibrate(close, dt):
    ror = close[1:]/close[:-1] - 1.0
    diff = []
    window = 30
    for i in range(window, len(ror)):
        hold = ror[i-window:i]
        diff.append(np.std(hold)**2)

    theta = np.mean(diff)
    sigma = np.std(diff)

    top = 0
    bot = 0

    for i in range(1, len(diff)):
        top += (diff[i] - theta)*(diff[i-1] - theta)
        bot += pow(diff[i] - theta, 2)

    kappa = -np.log(top/bot)/dt
        
    drift = np.mean(ror)
    S = close[-1]
    v0 = np.std(ror)**2

    return kappa, theta, sigma, S, drift, v0

cd = ctypes.c_double
ci = ctypes.c_int

models = ctypes.CDLL("./model.so")

models.Heston.argtypes = (
    cd,
    cd,
    cd,
    cd,
    cd,
    cd,
    cd,
    ci,
    ci
)

models.Heston.restype = cd

models.GBM.argtypes = (
    cd,
    cd,
    cd,
    cd,
    ci,
    ci
)

models.GBM.restype = cd

ticker = 'AMZN'

T = 1.0
N = 1500
P = 100
dt = T / N

data = requests.get(address(ticker)).json()
data = pd.DataFrame(data['historical'])
close = data[::-1]['adjClose'].values

kappa, theta, sigma, S, drift, v0 = Calibrate(close, dt)

heston, gbm = [], []
ux = []

SHeston = S
SGBM = S

days = 60
G = 100
adjClose = close[-G:]
adjX = list(range(G))

heston.append(SHeston)
gbm.append(SGBM)
ux.append(G-1)

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

for i in range(G, G+days+1):
    SHeston = models.Heston(kappa, theta, sigma, SHeston, drift, v0, dt, N, P)
    SGBM = models.GBM(SGBM, drift, v0, dt, N, P)

    ux.append(i)
    heston.append(SHeston)
    gbm.append(SGBM)

    ax.cla()
    ax.set_title("Geometric Browninan Motion vs. Heston Model")
    ax.set_xlabel("Time")
    ax.set_ylabel("Price")

    ax.plot(adjX, adjClose, color='red')
    ax.plot(ux, heston, color='blue')
    ax.plot(ux, gbm, color='black')
    plt.pause(0.1)

plt.show()


    
