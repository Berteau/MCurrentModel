from random import gauss

from pylab import *

import Neuron
from oldStuff import IonotropicConnection

# Time
tau = 0.1
tspan = arange(0,500, tau)
T1=len(tspan)/5

#Parameters
pyramidalParams = {}
pyramidalParams["C"]= 100
pyramidalParams["k"] = 3
pyramidalParams["v_r"] = -60
pyramidalParams["v_t"] = -50
pyramidalParams["v_peak"] = 50
# pyramidalParams["g_up"] = 3.0
# pyramidalParams["g_down"] = 5.0
pyramidalParams["a"] = 0.01
pyramidalParams["b"] = 5
pyramidalParams["c"] = -60
pyramidalParams["d"] = 400
pyramidalParams["tau"] = tau
pyramidalParams["type"] = "Pyramidal"
#

fsParams = {}
fsParams["C"]= 100
fsParams["k"] = 1
fsParams["v_r"] = -55
fsParams["v_t"] = -40
fsParams["v_peak"] = 25
fsParams["a"] = 0.15
fsParams["b"] = 8
fsParams["c"] = -55
fsParams["d"] = 200
fsParams["tau"] = tau
fsParams["type"] = "FS"

ltsParams = {}
ltsParams["C"]= 100
ltsParams["k"] = 1
ltsParams["v_r"] = -56
ltsParams["v_t"] = -42
ltsParams["v_peak"] = 40
ltsParams["a"] = 0.03
ltsParams["b"] = 8
ltsParams["c"] = -50
ltsParams["d"] = 200
ltsParams["tau"] = tau
ltsParams["type"] = "LTS"

#
# lowThresholdSpikingParams = {}
# params["a"]=0.02
# params["b"]=8
# params["c"]=-65
# params["d"]=6
# params["tau"] = tau
# params["V0"] = -56

popCount = 40

pyramidals = [Neuron.Neuron(pyramidalParams) for i in range(popCount)]
fastSpikings = [Neuron.Neuron(fsParams) for i in range(popCount)]
ltsNeurons = [Neuron.Neuron(ltsParams) for i in range(popCount)]
pyramidals2 = [Neuron.Neuron(pyramidalParams) for i in range(popCount)]
fastSpikings2 = [Neuron.Neuron(fsParams) for i in range(popCount)]
ltsNeurons2 = [Neuron.Neuron(ltsParams) for i in range(popCount)]
connections = []

# Set up connections in column one
for pyramidal in pyramidals:
    pyramidal.externalInput = gauss(550.0, 100.0)
    connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(2000 / popCount, 400 / popCount), pyramidal, pyramidal))
    for fastSpiking in fastSpikings:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(1200 / popCount, 2000 / popCount), pyramidal, fastSpiking))
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-200.0, 5.0), fastSpiking, pyramidal))

for ltsNeuron in ltsNeurons:
    for pyramidal in pyramidals:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(300.0, 50.0), pyramidal, ltsNeuron))
    for fastSpiking in fastSpikings:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-200.0, 50.0), fastSpiking, ltsNeuron))

# Set up connections in column two
for pyramidal2 in pyramidals2:
    pyramidal2.externalInput = gauss(950.0, 100.0)
    connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(50.0, 10.0), pyramidal2, pyramidal2))
    for fastSpiking2 in fastSpikings2:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(300.0, 50.0), pyramidal2, fastSpiking2))
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-200.0, 50.0), fastSpiking2, pyramidal2))

for ltsNeuron2 in ltsNeurons2:
    for pyramidal2 in pyramidals2:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(300.0, 50.0), pyramidal2, ltsNeuron2))
    for fastSpiking2 in fastSpikings2:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-200.0, 50.0), fastSpiking2, ltsNeuron2))

# # Set up long-range connections
for ltsNeuron in ltsNeurons:
    for fastSpiking2 in fastSpikings2:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-200.0, 50.0), ltsNeuron, fastSpiking2))
    for pyramidal2 in pyramidals2:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-50.0, 10.0), ltsNeuron, pyramidal2))

for ltsNeuron2 in ltsNeurons2:
    for fastSpiking in fastSpikings:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-200.0, 50.0), ltsNeuron2, fastSpiking))
    for pyramidal in pyramidals:
        connections.append(IonotropicConnection.IonotropicConnection(tau, 0, gauss(-50.0, 10.0), ltsNeuron2, pyramidal))


for t in range(len(tspan)):
    if t%10 == 0:
        print(t)
    for pyramidal in pyramidals:
        pyramidal.step()
    for pyramidal2 in pyramidals2:
        pyramidal2.step()
    for fastSpiking in fastSpikings:
        fastSpiking.step()
    for fastSpiking2 in fastSpikings2:
        fastSpiking2.step()
    for lts in ltsNeurons:
        lts.step()
    for lts2 in ltsNeurons2:
        lts2.step()
    # for connection in connections:
    #     connection.step()

pyramidalVoltages = [pyramidals[i].vv for i in range(len(pyramidals))]
fastSpikingVoltages = [fastSpikings[i].vv for i in range(len(fastSpikings))]
ltsVoltages = [ltsNeurons[i].vv for i in range(len(ltsNeurons))]
pyramidalVoltages2 = [pyramidals2[i].vv for i in range(len(pyramidals2))]
fastSpikingVoltages2 = [fastSpikings2[i].vv for i in range(len(fastSpikings2))]
ltsVoltages2 = [ltsNeurons2[i].vv for i in range(len(ltsNeurons2))]

# plot(tspan, pyramidals[0].vv)
# # plot(tspan, cell1.uu)
# title('Pyramidal Cell')

figure()
pcolor(pyramidalVoltages, vmin=-100, vmax=60)
colorbar()
title('Pyramidal Cells')

figure()
pcolor(pyramidalVoltages2, vmin=-100, vmax=60)
colorbar()
title('Pyramidal Cells 2')

# figure()
# plot(tspan, fastSpikings[0].vv)
# # plot(tspan, cell2.uu)
# title('FS Interneuron')

figure()
pcolor(fastSpikingVoltages, vmin=-100, vmax=60)
colorbar()
title('FS Interneurons')

figure()
pcolor(fastSpikingVoltages2, vmin=-100, vmax=60)
colorbar()
title('FS Interneurons 2')

figure()
pcolor(ltsVoltages, vmin=-100, vmax=60)
colorbar()
title('LTS Interneurons')

figure()
pcolor(ltsVoltages2, vmin=-100, vmax=60)
colorbar()
title('LTS Interneurons 2')

figure()
plot(tspan, fastSpikings[0].vv)
title('FSNeuron')

show()

