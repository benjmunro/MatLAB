function [time,force] = interpclass(dt)

% interpolating function for wind data

load('Wind.mat')         

time = [0:dt:Wind_record(1,end)];
force = interp1(Wind_record(1,:),Wind_record(2,:),time);