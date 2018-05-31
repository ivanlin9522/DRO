addpath C:\Users\Ivan\Desktop\robust_optimization\myExperiment; % adds the directory C:\XXX to the list of directories which Matlab "sees" (referred to as paths)
mlpath='C:\Users\Ivan\Desktop\robust_optimization\myExperiment' % The directory where mathlink.h is
mllib='C:\Users\Ivan\Desktop\robust_optimization\myExperiment\ml64i3m.lib' %The library ml32i3m.lib

%make command
command=sprintf('mex -D__STDC__ -I%s %s %s', mlpath, 'math.c', mllib);
%compile
eval(command)  