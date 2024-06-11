import Usefulfunctions.*

% Define a radiostripe with a transmitter and nolinks links/repeaters.
nolinks=5;
RS=c_radiostripe(nolinks);


% generate a QAM signal, and use RRC pulse shaping
N=1e4; OS=5; % number of symbols, and oversampling factor
x=randconst(N,1);
x2=pulseshape(x,OS,0.1);

% Run signal over stripe. The input is a vector N*1. The output is a matrix
% N*(nolinks+1), where the first column is the output signal of the
% transmitter, and the rest of the columns are the outputs of each link.
Y=RS.run(x2);



%-----------------
% plots

figure(1)
spec(Y)
set(gcf,'position',[100 100 400 300])

RS.calibrate(x2,0); % calibrate all amplifiers to give an output power of approx 0 dBm
y2=RS.run(x2);
figure(2)
spec(y2)
set(gcf,'position',[100 100 400 300])

RS.calibrate(x2,-5); % calibrate all amplifiers to give an output power of approx -5 dBm
y3=RS.run(x2);
figure(3)
spec(y3)
set(gcf,'position',[100 100 400 300])

figure(4)
amam(x2,timealign(x2,Y(:,1)),x2,timealign(x2,Y(:,3)),x2,timealign(x2,Y(:,6)))
set(gcf,'position',[100 100 400 300])
figure(5)
amam(x2,timealign(x2,y3(:,1)),x2,timealign(x2,y3(:,3)),x2,timealign(x2,y3(:,6)))
set(gcf,'position',[100 100 400 300])
