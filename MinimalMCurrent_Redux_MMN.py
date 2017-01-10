from __future__ import division
import dipoleWrapperNeuronActualDipole
import dipoleWrapperInterNeuron
import poissonSourceAxon
from pylab import *
from random import randint, gauss
import NeuronNeuronAxon
from scipy.stats import poisson
import pickle
from scipy.signal import butter, lfilter, freqz
from matplotlib.font_manager import FontProperties, findfont
import matplotlib.font_manager as fm
from mpl_toolkits.mplot3d import Axes3D


'''
Created on 1 juin 2014

@author: Berteau
'''

# Set up time
tau = 0.1
maxTime = 1.0
maxTimeInMS = maxTime * 1000
maxTimeStep = maxTime / tau
tspan = arange(0,maxTimeInMS/tau, 1)
T1=len(tspan)/5
T2 = 4*T1

# Neuron Parameters
neuronHeight = 0.0003
drivingOnset = 300
drivingOffset = 700
primingEOnset = 150
primingEOffset = 700
primingIOnset = 175
primingIOffset = 700
drivingWeight = 16  # 30.0
primingEWeight = 12  # 12
primingIWeight = 6   # 8
driveLambda = 0.01
primeELambda = 0.01
primeILambda = 0.01
inactivating = True
ACh = False
Retigabine = False
PCP = False
pharmas = 'Baseline'
# pharmas = 'PCP'
# pharmas = 'PCPRetigabine'
# pharmas = 'TargetPrimedBaseline'
# pharmas = 'ACh'
predictionInputCount = 200
sensoryInputCount = 200

targetPrimeELambda = 0.0000001
targetPrimeILambda = 0.0000001
comparatorTargetWeight = 1200

# Set up neurons
neuron = dipoleWrapperNeuronActualDipole.DipoleWrapperNeuron(inputs_a=[], inputs_b=[], outputs=[], externalInput_a=0.0, externalInput_b=0.0, distance = neuronHeight, position = [randint(1,100), randint(1,100), 0.0], debug=False, tau=tau, tspan=tspan, inactivating=inactivating)
targetNeuron = dipoleWrapperNeuronActualDipole.DipoleWrapperNeuron(inputs_a=[], inputs_b=[], outputs=[], externalInput_a=0.0, externalInput_b=0.0, distance=neuronHeight, position=[randint(1,100), randint(1,100), 0.0], debug=False, tau=tau, tspan=tspan, inactivating=inactivating)
# inhibitoryNeuron = dipoleWrapperInterNeuron.DipoleWrapperInterNeuron(inputs_a=[], inputs_b=[], outputs=[], externalInput_a=0.0, externalInput_b=0.0, distance=neuronHeight, position=[randint(1,100), randint(1,100), 0.0], debug=False, tau=tau, tspan=tspan, inactivating=inactivating)

if ACh:
    neuron.base.gdapt = 3.0
    neuron.apex.gdapt = 3.0
    targetNeuron.base.gdapt = 3.0
    targetNeuron.apex.gdapt = 3.0

if PCP:
    neuron.base.g_nmda = 1.5e-3 # Try setting this to -4 instead
    neuron.apex.g_nmda = 1.5e-3
    targetNeuron.base.g_nmda = 1.5e-3  # Try setting this to -4 instead
    targetNeuron.apex.g_nmda = 1.5e-3
    primingIWeight -= 4.0

if Retigabine:
    neuron.base.gdapt = 5.0
    neuron.apex.gdapt = 5.0
    targetNeuron.base.gdapt = 5.0
    targetNeuron.apex.gdapt = 5.0
    primingIWeight += 1.0

drivingConnections = []
primingEConnections = []
primingIConnections = []
targetConnections = []

# Set up driving sensory input
for n in range(sensoryInputCount):
    weight = drivingWeight*gauss(1.0, 0.0001) #Original weight 0.025
    drivingConnections.append(poissonSourceAxon.PoissonAxon(tau, driveLambda, drivingOnset, drivingOffset))
    neuron.addInput_b(drivingConnections[-1], weight)

# Set up priming prediction inputs
for n in range(predictionInputCount):
    weightE = primingEWeight*gauss(1.0, 0.0001)
    weightI = -1*primingIWeight*gauss(1.0, 0.0001)
    primingEConnections.append(poissonSourceAxon.PoissonAxon(tau, primeELambda, primingEOnset, primingEOffset))
    neuron.addInput_a(primingEConnections[-1], weightE)
    primingIConnections.append(poissonSourceAxon.PoissonAxon(tau, primeILambda, primingIOnset, primingIOffset))
    neuron.addInput_b(primingIConnections[-1], weightI)

# Set up Target connections
targetConnections.append(NeuronNeuronAxon.NeuronNeuronAxon(tau, [1.0], neuron, [targetNeuron], width=0.1, myelin=0.0))
targetNeuron.addInput_b(targetConnections[-1], comparatorTargetWeight*gauss(1.0, 0.0001))
for n in range(predictionInputCount):
    weightE = primingEWeight * gauss(1.0, 0.1)
    weightI = -1 * primingIWeight * gauss(1.0, 0.1);
    primingEConnections.append(poissonSourceAxon.PoissonAxon(tau, targetPrimeELambda, primingEOnset, primingEOffset))
    targetNeuron.addInput_a(primingEConnections[-1], weightE)
    primingIConnections.append(poissonSourceAxon.PoissonAxon(tau, targetPrimeILambda, primingIOnset, primingIOffset))
    targetNeuron.addInput_b(primingIConnections[-1], weightI)

# Set up records
KCNQOpen = zeros(maxTimeInMS/tau)
apexDV = zeros(maxTimeInMS/tau)
baseDV = zeros(maxTimeInMS/tau)
KCNQTau = zeros(maxTimeInMS/tau)
MCurrent = zeros(maxTimeInMS/tau)
Voltage = zeros(maxTimeInMS/tau)
hRecord = zeros(maxTimeInMS/tau)
electrodeTrace = zeros(maxTimeInMS/tau)
targetVoltage = zeros(maxTimeInMS/tau)
targetElectrode = zeros(maxTimeInMS/tau)

# Run the simluation
for time in range(len(tspan)):
    print "time:", time
    neuron.step(time)
    KCNQOpen[time] = neuron.base.z
    KCNQTau[time] = neuron.base.tau_z
    MCurrent[time] = neuron.base.I_M
    Voltage[time] = neuron.base.v
    apexDV[time] = neuron.apex.dv_temp
    baseDV[time] = neuron.base.dv_temp
    electrodeTrace[time] = (((neuron.base.IExternal - neuron.apex.IExternal) * neuronHeight * cos(0)) / (
        4 * pi * 8.854187817e-12 * (0.01) ** 2)) / 100000000000
    targetNeuron.step(time)
    targetVoltage[time] = targetNeuron.base.v
    targetElectrode[time] = (((targetNeuron.base.IExternal - targetNeuron.apex.IExternal) * neuronHeight * cos(0)) / (
        4 * pi * 8.854187817e-12 * (0.01) ** 2)) / 100000000000

    for connection in drivingConnections:
        connection.step()
    for connection in primingEConnections:
        connection.step()
    for connection in primingIConnections:
        connection.step()
    for connection in targetConnections:
        connection.step()

# Plot and save results
figDir = 'DissPlots'


# ## Electrode Trace
# figure()
# plot(electrodeTrace)
# title('Electric Field Trace from an electrode 1cm above the neuron')
# xlabel('Time in ms')
# ylabel('Voltage in mV')
# #        legend()
# # savefig(testTitle + '/' + runString + 'ElectrodeTrace.png')
# # savefig(testTitle + '/' + runString + 'ElectrodeTrace.svg')
#
# # Filter Electrode Trace

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

# Filter requirements.
order = 6
fs = 10000.0       # sample rate, Hz
cutoff = 14  # desired cutoff frequency of the filter, Hz

# Get the filter coefficients so we can check its frequency response.
b, a = butter_lowpass(cutoff, fs, order)

# # Plot the frequency response.
# w, h = freqz(b, a, worN=8000)
# plt.subplot(2, 1, 1)
# plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
# plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
# plt.axvline(cutoff, color='k')
# plt.xlim(0, 0.5*fs)
# plt.title("Lowpass Filter Frequency Response")
# plt.xlabel('Frequency [Hz]')
# plt.grid()

filteredElectrodeTracePrimed = butter_lowpass_filter(electrodeTrace, cutoff, fs, order)
filteredTargetTracePrimed = butter_lowpass_filter(targetElectrode, cutoff, fs, order)

# Set font properties
prop = fm.FontProperties(fname='C:\Users\stefa\OneDrive\Fonts\mreavessansreg\mreavessansreg\MrEavesSanR.ttf')

## Filtered Electrode Trace
figure()
plot(arange(0,maxTimeInMS,tau), filteredElectrodeTracePrimed, linewidth=1.5)
title('Electric Field: Primed Comparator', fontproperties=prop, fontsize=24)
xlabel('Time in ms', fontproperties=prop, fontsize=24)
ylabel('Voltage in mV', fontproperties=prop, fontsize=24)
a = gca()
ylim([-10,7])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_PrimedComparatorElectrodeTrace.png')
savefig(figDir + '/' + pharmas + '_PrimedComparatorElectrodeTrace.eps')

# fig = figure()
# ax = fig.gca(projection='3d')
# ax.plot(arange(0,maxTimeInMS,tau), Voltage, MCurrent)
# title("Membrane Voltage and M-Current over Time")
# xlabel("Time in ms")
# ylabel("Membrane Voltage in mV")
# savefig(figDir + '/' + pharmas + '_VoltageMCurrent.png')
# savefig(figDir + '/' + pharmas + '_VoltageMCurrent.eps')

figure()
plot(arange(0,maxTimeInMS,tau),MCurrent, linewidth=1.5)
title("Primed Comparator M-Current", fontproperties=prop, fontsize=24)
xlabel("Time in MS", fontproperties=prop, fontsize=24)
ylabel("Current in mA", fontproperties=prop, fontsize=24)
a = gca()
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_PrimedComparatorMCurrent.png')
savefig(figDir + '/' + pharmas + '_PrimedComparatorMCurrent.eps')

figure()
plot(arange(0,maxTimeInMS,tau),Voltage, linewidth=1.5)
title("Primed Comparator Voltage Trace", fontproperties=prop, fontsize=24)
xlabel("Time in MS", fontproperties=prop, fontsize=24)
ylabel("Voltage in mV", fontproperties=prop, fontsize=24)
a = gca()
ylim([-80,40])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)

savefig(figDir + '/' + pharmas + '_PrimedComparatorVoltage.png')
savefig(figDir + '/' + pharmas + '_PrimedComparatorVoltage.eps')

## Filtered Electrode Trace
figure()
plot(filteredTargetTracePrimed, linewidth=1.5)
title('Electric Field: Target (Primed Comparator)', fontproperties=prop, fontsize=24)
xlabel('Time in ms', fontproperties=prop, fontsize=24)
ylabel('Voltage in mV', fontproperties=prop, fontsize=24)
a = gca()
ylim([-10,7])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_PrimedTargetElectrodeTrace.png')
savefig(figDir + '/' + pharmas + '_PrimedTargetElectrodeTrace.eps')

figure()
plot(arange(0,maxTimeInMS,tau),targetVoltage, linewidth=1.5)
title("Target Voltage Trace (Primed Comparator)", fontproperties=prop, fontsize=24)
xlabel("Time in MS", fontproperties=prop, fontsize=24)
ylabel("Voltage in mV", fontproperties=prop, fontsize=24)
a = gca()
ylim([-80,40])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_PrimedTargetVoltage.png')
savefig(figDir + '/' + pharmas + '_PrimedTargetVoltage.eps')

# resolutionR = 40
# resolutionAngle = 360
#
# primeTime = int(185 / tau)
# driveTime = int(397 / tau)
# print primeTime
# print driveTime
#
# polarField = array([[(((apexDV[driveTime] - baseDV[driveTime])) * neuronHeight * cos(theta)) / (
# 4 * pi * 8.854187817e-12 * (r / 10) ** 2) for r in linspace(1, 2, resolutionR)] for theta in
#                     linspace(0, 2 * pi, resolutionAngle)])
# polarFieldPriming = array([[(((apexDV[primeTime] - baseDV[primeTime])) * neuronHeight * cos(theta)) / (
# 4 * pi * 8.854187817e-12 * (r / 10) ** 2) for r in linspace(1, 2, resolutionR)] for theta in
#                            linspace(0, 2 * pi, resolutionAngle)])
#
# rectField = zeros([resolutionR * 2, resolutionR * 2])
# rectFieldPriming = zeros([resolutionR * 2, resolutionR * 2])
#
# print size(rectField)
#
# for angle in range(resolutionAngle):
#     for r in range(resolutionR):
#         x = r * cos((angle / resolutionAngle) * 2 * pi)
#         y = r * sin((angle / resolutionAngle) * 2 * pi)
#         rectField[x + resolutionR, y + resolutionR] = polarField[angle, r]
#         rectFieldPriming[x + resolutionR, y + resolutionR] = polarFieldPriming[angle, r]
#
# X = arange(1, shape(rectFieldPriming)[0] + 1)
# Y = arange(1, shape(rectFieldPriming)[1] + 1)
#
# figure()
# pcolor(rectField)
# # pcolor(X, Y, rectField, vmin=-2.4e12, vmax=2.4e12)
# # contour(rectField, origin='lower')
# colorbar()
# title('Polar contour of the electric field during Driving')
# savefig(figDir + '/' + runString + 'DrivingField_CurrentClamp.png')
# savefig(figDir + '/' + runString + 'DrivingField_CurrentClamp.svg')
#
# figure()
# pcolor(rectFieldPriming)
# # pcolor(X, Y, rectFieldPriming, vmin=-2.4e12, vmax=2.4e12)
# # contour(rectFieldPriming, origin='lower')
# colorbar()
# title('Polar contour of the electric field during Priming')
# savefig(figDir + '/' + runString + 'PrimingField_CurrentClamp.png')
# savefig(figDir + '/' + runString + 'PrimingField_CurrentClamp.svg')
# #
#
# figure()
# plot(arange(0,maxTimeInMS,tau),Voltage)
# title("Neuron Voltage Trace")
# xlabel("Time in MS")
# ylabel("Voltage in mV")
# savefig(figDir+'/Voltage_'+runString+'.png')
# print "Saved:" + figDir+'/Voltage_'+runString+'.png'
#
# figure()
# plot(arange(0,maxTimeInMS,tau),KCNQOpen)
# title("Proportion of KCNQ Channels Open")
# xlabel("Time in MS")
# ylabel("Proportion (1.0=100%)")
# savefig(figDir+'/KCNQOpen_'+runString+'.png')
# print "Saved:" + figDir+'/KCNQOpen_'+runString+'.png'
#
# figure()
# plot(arange(0,maxTimeInMS,tau),KCNQTau)
# title("KCNQ Time Constant")
# xlabel("Time in MS")
# ylabel("Time Constant in MS")
# savefig(figDir+'/KCNQTau_'+runString+'.png')
# print "Saved:" + figDir+'/KCNQTau_'+runString+'.png'
#
# figure()
# plot(arange(0,maxTimeInMS,tau),MCurrent)
# title("MCurrent over Time")
# xlabel("Time in MS")
# ylabel("MCurrent in mA")
# savefig(figDir+'/MCurrent_'+runString+'.png')
# print "Saved:" + figDir+'/MCurrent_'+runString+'.png'
#
# figure()
# plot(arange(0,maxTimeInMS,tau),hRecord)
# title("Proportion of Non-Inactivated Na Channels")
# xlabel("Time in MS")
# ylabel("Proportion (1.0 = 100%)")
# savefig(figDir+'/hRecord_'+runString+'.png')
# print "Saved:" + figDir+'/hRecord_'+runString+'.png'

# show()


# Neuron Parameters
neuronHeight = 0.0003
drivingOnset = 300
drivingOffset = 700
primingEOnset = 150
primingEOffset = 700
primingIOnset = 175
primingIOffset = 700
driveLambda = 0.01
primeELambda = 0.000001
primeILambda = 0.000001
predictionInputCount = 200
sensoryInputCount = 200

# Set up neurons
neuron = dipoleWrapperNeuronActualDipole.DipoleWrapperNeuron(inputs_a=[], inputs_b=[], outputs=[], externalInput_a=0.0, externalInput_b=0.0, distance = neuronHeight, position = [randint(1,100), randint(1,100), 0.0], debug=False, tau=tau, tspan=tspan, inactivating=inactivating)
targetNeuron = dipoleWrapperNeuronActualDipole.DipoleWrapperNeuron(inputs_a=[], inputs_b=[], outputs=[], externalInput_a=0.0, externalInput_b=0.0, distance=neuronHeight, position=[randint(1,100), randint(1,100), 0.0], debug=False, tau=tau, tspan=tspan, inactivating=inactivating)


if ACh:
    neuron.base.gdapt = 3.0
    neuron.apex.gdapt = 3.0
    targetNeuron.base.gdapt = 3.0
    targetNeuron.apex.gdapt = 3.0

if Retigabine:
    neuron.base.gdapt = 5.0
    neuron.apex.gdapt = 5.0
    targetNeuron.base.gdapt = 5.0
    targetNeuron.apex.gdapt = 5.0

if PCP:
    neuron.base.g_nmda = 1.5e-3 # Try setting this to -4 instead
    neuron.apex.g_nmda = 1.5e-3
    targetNeuron.base.g_nmda = 1.5e-3  # Try setting this to -4 instead
    targetNeuron.apex.g_nmda = 1.5e-3


drivingConnections = []
primingEConnections = []
primingIConnections = []

# Set up driving sensory input
for n in range(sensoryInputCount):
    weight = drivingWeight*gauss(1.0, 0.0001) #Original weight 0.025
    drivingConnections.append(poissonSourceAxon.PoissonAxon(tau, driveLambda, drivingOnset, drivingOffset))
    neuron.addInput_b(drivingConnections[-1], weight)

# Set up priming prediction inputs
for n in range(predictionInputCount):
    weightE = primingEWeight*gauss(1.0, 0.0001)
    weightI = -1*primingIWeight*gauss(1.0, 0.0001);
    primingEConnections.append(poissonSourceAxon.PoissonAxon(tau, primeELambda, primingEOnset, primingEOffset))
    neuron.addInput_a(primingEConnections[-1], weightE)
    primingIConnections.append(poissonSourceAxon.PoissonAxon(tau, primeILambda, primingIOnset, primingIOffset))
    neuron.addInput_b(primingIConnections[-1], weightI)

# Set up Target connections
targetConnections.append(NeuronNeuronAxon.NeuronNeuronAxon(tau, [1.0], neuron, [targetNeuron], width=0.1, myelin=0.0))
targetNeuron.addInput_b(targetConnections[-1], comparatorTargetWeight*gauss(1.0, 0.0001))
for n in range(predictionInputCount):
    weightE = primingEWeight * gauss(1.0, 0.0001)
    weightI = -1 * primingIWeight * gauss(1.0, 0.0001);
    primingEConnections.append(poissonSourceAxon.PoissonAxon(tau, targetPrimeELambda, primingEOnset, primingEOffset))
    targetNeuron.addInput_a(primingEConnections[-1], weightE)
    primingIConnections.append(poissonSourceAxon.PoissonAxon(tau, targetPrimeILambda, primingIOnset, primingIOffset))
    targetNeuron.addInput_b(primingIConnections[-1], weightI)

# Set up records
KCNQOpen = zeros(maxTimeInMS/tau)
apexDV = zeros(maxTimeInMS/tau)
baseDV = zeros(maxTimeInMS/tau)
KCNQTau = zeros(maxTimeInMS/tau)
MCurrent = zeros(maxTimeInMS/tau)
Voltage = zeros(maxTimeInMS/tau)
hRecord = zeros(maxTimeInMS/tau)
electrodeTrace = zeros(maxTimeInMS/tau)
targetVoltage = zeros(maxTimeInMS/tau)
targetElectrode = zeros(maxTimeInMS/tau)

# Run the simluation
for time in range(len(tspan)):
    print "time:", time
    neuron.step(time)
    KCNQOpen[time] = neuron.base.z
    KCNQTau[time] = neuron.base.tau_z
    MCurrent[time] = -neuron.base.I_M
    Voltage[time] = neuron.base.v
    apexDV[time] = neuron.apex.dv_temp
    baseDV[time] = neuron.base.dv_temp
    electrodeTrace[time] = (((neuron.base.IExternal - neuron.apex.IExternal) * neuronHeight * cos(0)) / (
        4 * pi * 8.854187817e-12 * (0.01) ** 2)) / 100000000000
    targetNeuron.step(time)
    targetVoltage[time] = targetNeuron.base.v
    targetElectrode[time] = (((targetNeuron.base.IExternal - targetNeuron.apex.IExternal) * neuronHeight * cos(0)) / (
        4 * pi * 8.854187817e-12 * (0.01) ** 2)) / 100000000000

    for connection in drivingConnections:
        connection.step()
    for connection in primingEConnections:
        connection.step()
    for connection in primingIConnections:
        connection.step()
    for connection in targetConnections:
        connection.step()

filteredElectrodeTraceUnPrimed = butter_lowpass_filter(electrodeTrace, cutoff, fs, order)
filteredTargetTraceUnPrimed = butter_lowpass_filter(targetElectrode, cutoff, fs, order)

## Filtered Electrode Trace

prop = fm.FontProperties(fname='C:\Users\stefa\OneDrive\Fonts\mreavessansreg\mreavessansreg\MrEavesSanR.ttf')

fig = figure()
plot(arange(0,maxTimeInMS,tau), filteredElectrodeTraceUnPrimed, linewidth=1.5)
title('Electric Field: Unprimed Comparator', fontproperties=prop, fontsize=24)
xlabel('Time in ms', fontproperties=prop, fontsize=24)
ylabel('Voltage in mV', fontproperties=prop, fontsize=24)
a = gca()
ylim([-10,7])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_UnprimedComparatorElectrodeTrace.png')
savefig(figDir + '/' + pharmas + '_UnprimedComparatorElectrodeTrace.eps')


figure()
plot(arange(0,maxTimeInMS,tau), Voltage, linewidth=1.5)
title("Unprimed Comparator Voltage Trace", fontproperties=prop, fontsize=24)
xlabel("Time in MS", fontproperties=prop, fontsize=24)
ylabel("Voltage in mV", fontproperties=prop, fontsize=24)
a = gca()
a.set_ylim([-80, 40])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_UnprimedComparatorVoltage.png')
savefig(figDir + '/' + pharmas + '_UnprimedComparatorVoltage.eps')

figure()
plot(arange(0,maxTimeInMS,tau),MCurrent, linewidth=1.5)
title("Unprimed Comparator M-Current", fontproperties=prop, fontsize=24)
xlabel("Time in MS", fontproperties=prop, fontsize=24)
ylabel("Current in mA", fontproperties=prop, fontsize=24)
a = gca()
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_UnprimedComparatorMCurrent.png')
savefig(figDir + '/' + pharmas + '_UnprimedComparatorMCurrent.eps')

## Filtered Electrode Trace
figure()
plot(arange(0,maxTimeInMS,tau), filteredTargetTraceUnPrimed, linewidth=1.5)
title('Electric Field: Target (Unprimed Comparator)', fontproperties=prop, fontsize=24)
xlabel('Time in ms', fontproperties=prop, fontsize=24)
ylabel('Voltage in mV', fontproperties=prop, fontsize=24)
a = gca()
ylim([-10,7])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_UnprimedTargetElectrodeTrace.png')
savefig(figDir + '/' + pharmas + '_UnprimedTargetElectrodeTrace.eps')



figure()
plot(arange(0,maxTimeInMS,tau),targetVoltage, linewidth=1.5)
title("Target Voltage Trace (Unprimed Comparator)", fontproperties=prop, fontsize=24)
xlabel("Time in MS", fontproperties=prop, fontsize=24)
ylabel("Voltage in mV", fontproperties=prop, fontsize=24)
a = gca()
ylim([-80,40])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_UnprimedTargetElectrodeTrace.png')
savefig(figDir + '/' + pharmas + '_UnprimedTargetElectrodeTrace.eps')


differenceWave = filteredElectrodeTraceUnPrimed - filteredElectrodeTracePrimed
targetDifference = filteredTargetTraceUnPrimed - filteredTargetTracePrimed

figure()
plot(arange(0,maxTimeInMS,tau), differenceWave, linewidth=1.5)
title('Electric Field Difference Wave: Comparator', fontproperties=prop, fontsize=24)
xlabel('Time in ms', fontproperties=prop, fontsize=24)
ylabel('Voltage in mV', fontproperties=prop, fontsize=24)
a = gca()
ylim([-10,7])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_DifferenceWaveComparator.png')
savefig(figDir + '/' + pharmas + '_DifferenceWaveComparator.eps')

figure()
plot(arange(0,maxTimeInMS,tau), targetDifference, linewidth=1.5)
title('Electric Field Difference Wave: Target', fontproperties=prop, fontsize=24)
xlabel('Time in ms', fontproperties=prop, fontsize=24)
ylabel('Voltage in mV', fontproperties=prop, fontsize=24)
a = gca()
ylim([-10,7])
a.set_xticklabels(a.get_xticks(), fontproperties=prop, fontsize=20)
a.set_yticklabels(a.get_yticks(), fontproperties=prop, fontsize=20)
savefig(figDir + '/' + pharmas + '_DifferenceWaveTarget.png')
savefig(figDir + '/' + pharmas + '_DifferenceWaveTarget.eps')

show()