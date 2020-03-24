clc;
clear;
close all;

syms a1 a2 a3 a4 a5;
syms b1 b2 b3 b4;
syms x1 x2 x3 x4 x5;

x = [x1; x2; x3; x4; x5];

A = [a1 b1 0 0 0; 0 a2 b2 0 0; 0 0 a3 b3 0; 0 0 0 a4 b4; 0 0 0 0 a5];