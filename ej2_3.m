% Script to perform excercise 2.3
clear;   close all;   clc;

%Define la función DIXMAANA
f=@(x) dixmaana(x);

% n en {240,960,3,9,15,51,87} 

x0 = 2*ones([n,1]);

