clc; close all; clear all
nx = 41;
ny = 41;

a1 = dir('sol01_*.txt');
u = load(a1(1).name);
r0 = reshape(u(:,3),nx,ny);
x = reshape(u(:,1),nx,ny);
y = reshape(u(:,2),nx,ny);

csvwrite('csv/csv_0000.csv', u)

for i = 1:1:length(a1)
    r = load(a1(i).name);
%     r = reshape(r(:,3),nx,ny);
%     contour(x,y,r,linspace(-0.5,0.5,40))
    filenum = num2str(i, '%04.f');
    outstr = strcat('csv/csv_', filenum);
    outstr = strcat(outstr, '.csv');
    csvwrite(outstr, r)
end