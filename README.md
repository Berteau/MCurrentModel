# MCurrentModel
A Simple Model of how input patterns can cause encoding switches in pyramidal cells
Based on Berteau, S., Mingolla, E., Bullock, D., "A Biophysical Model of the Mismatch Negativity".  
Paper in preperation.

Requirements: Scipy, MatPlotLib, Numpy

To Run: launch MinimalMCurrent_Redux_MMN.py

This file will caulculate the neural response for a L2/3 pyramidal cell exposed to a sensory input, 
both under expectational "priming" input and without it.  It will also calculate the quasi-static electric
fields near the neuron, and perform the subtraction necessecary to produce an MMN trace.

There are multiple pharmacological manipulations available as well.  By default, they are all commented out 
except the baseline, or no-pharma option.  To apply a manipulation, comment out line 49, and uncomment one 
of the following lines before running the program.
