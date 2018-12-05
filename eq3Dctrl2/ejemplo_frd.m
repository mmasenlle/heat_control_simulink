clear all; close all; clc;

G = tf(1,[1 2 1]);
w = logspace(-2,2,1000);
respG = squeeze(freqresp(G,w));
Gfrd = frd(respG,w)

bode(G,Gfrd)