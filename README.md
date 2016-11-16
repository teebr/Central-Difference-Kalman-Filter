# Central-Difference-Kalman-Filter
Implement a generic discrete-time central difference Kalman Filter (CDKF) in MATLAB/ Simulink.
This is an abstract class implemented as a matlab.System object. To use it for a specific application, you must inherit it and overwrite the stateFcn and outputFcn functions based on your specific model (these names can't change).
These take the respective forms:
function x_new = stateFcn(obj,x,u,w)
...

function y = outputFcn(obj,x,u,v)
...


Additional properties (e.g. model parameters) can be added using the adding a Properties definition block. Same applies to methods (Static methods might be useful for certain models) 

I have validated this and have had it running for a few months, including on a Speedgoat real-time target. However, I have only ever implemented SISO models. 

I have written a CDKF with adaptive output noise, and also a CDKF with box constraints on the state estimates, but there are a bit shaky - maybe I'll upload them one day.
