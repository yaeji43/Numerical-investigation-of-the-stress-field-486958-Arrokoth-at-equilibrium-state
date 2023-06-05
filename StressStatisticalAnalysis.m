clear all;
close all;
clc;

%% ReadInput
%% part I
MV1 = load('Low_Resolution/1/node.txt')*1.e+3;
MF1 = load('Low_Resolution/1/elem.txt')+1;

lenE = length(MF1);
lenN = length(MV1);

N = [3,6];  % the number of columns for each section
H = [5,1];   % the number of header lines (note added in the trailer here
fid=fopen('THL_Output_.txt');
% node/ stress/ disp/ mass
% node/ stress/ den/ nodal/ mass


for i=1:length(N)
    fmt=repmat('%f',1,N(i));  % build the format string for the number columns
    a(i)=textscan(fid,fmt,'headerlines',H(i),'Delimiter',',','collectoutput',1); % read section
end
fclose('all');

MV = a{1};
Stress = a{2};

%% Stress distribution
angle = 35*pi/180;
DPS = zeros(lenN,1);
i1 = zeros(lenN,1);
j2 = zeros(lenN,1);
for i = 1:lenN
    s11 = Stress(i,1);
    s22 = Stress(i,2);
    s33 = Stress(i,3);
    s12 = Stress(i,4);
    s31 = Stress(i,5);
    s23 = Stress(i,6);

    i1(i) = (s11 + s22 + s33)/3;
    j2(i) = (1/6)*((s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2) + (s12^2 + s31^2 + s23^2);
    alpha = 2*sin(angle)/(sqrt(3)*(3-sin(angle)));
    DPS(i) = sqrt(3)*(3-sin(angle))*(alpha*3*i1(i) + sqrt(j2(i)))/(6*cos(angle));
end

for i = 1:lenN
    if DPS(i) < 0
        DPS(i) = 0;
    end
end





