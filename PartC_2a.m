clc
clear

R=100;
t=10;
P_1=40; %max pressure for cycling loading
P_2=80; %max pressure for short term increase
K_ic=90;

%cycling loading
sthet1=P_1*R/t; 
sax1=(P_1*R)/(2*t);

%short term increase
sthet2=P_2*R/t;
sax2=(P_2*R)/(2*t);

%determine failure mechanism
af1=(K_ic^2)/(pi*sax1^2)*1000; %leak before failure
af2=1000*(1/pi)*(K_ic/sax2)^2; %leak before failure
 
%cycle per mm of different crack lengths
dNda=[];
a0=0.1:0.1:t
for i=1:length(a0)
    dK=0.728*40*sqrt(pi*a0(i));
    dNda(i)=1/((10^-12)*(dK)^4);
end

%find prob of detection
pod=[];
a=0.1:0.1:10;

for j=1:length(a);
    pod(j)=fn_pod(1944637,a(j));
end

%probability of failure (1-pod)
pof=[];

for i=1:length(a)
    P_pod(i)=fn_pod(1944637,a(i));
    pof(i)=1-P_pod(i)
end

% plot(a0,dNda);
% ylabel('Cycle per mm');
% xlabel('Crack length (mm)')
% figure
% plot(a,pod);
% xlabel('Crack length (mm)');
% ylabel('Probability of detection');
% figure
% plot(a,pof);
% xlabel('Crack length (mm)');
% ylabel('Probability of failure');
% figure

%growth rate of crack

paris = [];
c_length = 0.0003;
C = 10^-12;
m = 4;
i = 1;
t = 0.01;
 
 
while c_length <= 0.01;
    Y = 0.728 + 0.373*(c_length/t)^2 - 0.029*(c_length/t)^4;
    K = Y*sax1*sqrt(pi*c_length);
    paris(i,1) = c_length; %crack radius
    paris(i,2) = K; 
    paris(i,3) = C*K^m; %dadN or crack growth rate
    c_length = c_length + paris(i,3);  %new length
    paris(i,4) = c_length; %new length 
    paris(i,5) = i; %initial length
    i = i+1;
    
end
% 
% dNda(:)=1/paris(:,3);
% crack(:)=(10^3).*paris(:,4);
% plot(crack,N);
% xlabel('Crack length (mm)');
% ylabel('dN/da (Cycles/mm)');
% xlim=([1 10]);


%% % Number of cycles from initial to failure

paris2 = [];
n2 = 1;
a = 1e-3; % initial crack length
a1 = 1;
af=10e-3;

while a <= af

    Y = 0.728 + (0.373*(a/0.01)^2) - (0.029*(a/0.01)^4);

    K = Y*sax1*sqrt(pi*a);

    paris2(n2,1) = a*1000; % crack radius

    paris2(n2,2) = 1 - fn_pod(1944637,a*1000); % prob of not detection cracks above this radius

    paris2(n2,3) = fn_pod(1944637,a*1000); % prob of detection
    
    paris2(n2,4) = n2;

    dadN(n2,1) = a*1000;

    dadN(n2,2) = C*(K^m); % cack growth rate at that crack radius

    a = a + dadN(n2,2);

    n2 = n2 + 1;


end

figure
plot(paris2(:,1),paris2(:,3))
xlabel('Crack Radius "a" [mm]')
ylabel('Probability of Detection')

figure
plot(paris2(:,1),paris2(:,2))
xlabel('Crack Length "a" [mm]')
ylabel('Probability of Failure')

%% Failure Variables

NC = length(dadN);
xaxis1(:,1) = 1:NC;

failure = [];

cycles = 5200:1:6200; % cycles per inspection
cycles = transpose(cycles);

insp_int = [];
a2 = 1;
failure = [];
IIT = {};

%% Finding Probability of undetection at inspection number

for j = 1:length(cycles)

    a2 = 1;

    cycle = cycles(j);

    for i = 1:length(paris2)
        
        insp_int(i,1) = xaxis1(i,1)/cycle;
    
        if rem(insp_int(i,1),1) == 0
    
            IIT{j,1}(a2,1) = insp_int(i,1); % inspection number
            IIT{j,1}(a2,2) = paris2(i,2); % prob. of undetection
            IIT{j,1}(a2,3) = paris2(i,1); % crack radius @ inspection

            a2 = a2 + 1;
    
        end

    end

end

%% Finding probability of failure for each cycle variable

for i = 1:length(IIT)

    for k = 2:length(IIT{i,1})

        IIT{i,1}(1,4) = IIT{i,1}(1,2);
    
        IIT{i,1}(k,4) = IIT{i,1}(k,2) * IIT{i,1}((k-1),4);

    end

    failure(i,1) = cycles(i); % cycles per iteration
    failure(i,2) = IIT{i,1}(end,4); % probability of failure

end

% Graph to show all variations of inspection
figure
scatter(failure(:,1),failure(:,2), 7)
hold on
yline(0.01)
xlabel('No. of Cycles per Inspection')
ylabel('Probability of Failure')
hold off
% 
% % Graph based on inspection for 8000:12000 cycles
% figure
% scatter(IIT{316,1}(:,1),IIT{316,1}(:,4), 7)
% hold on
% scatter(IIT{317,1}(:,1),IIT{317,1}(:,4), 7)
% yline(0.01)
% xlabel('Inspection Number')
% ylabel('Probability of Failure @ Inspection')
% hold off
% 
% % Graph based on inspection for 13000:16000 cycles
% figure
% scatter(IIT{361,1}(:,1),IIT{361,1}(:,4), 7)
% hold on
% scatter(IIT{362,1}(:,1),IIT{362,1}(:,4), 7)
% yline(0.01)
% xlabel('Inspection Number')
% ylabel('Probability of Failure @ Inspection')
% hold off
% 
% % Graph of Probability of non-detection vs number of cycles
% figure
% plot(paris2(:,4), paris2(:,2))
% xlabel('Cycle Number')
% ylabel('Probability of Non-detection')

















