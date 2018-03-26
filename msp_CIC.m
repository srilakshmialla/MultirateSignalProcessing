%% Multirate Signal Processing
%
% SRILAKSHMI ALLA
%
%
% San Diego State University
%
% Professor: fredric j.harris
%
% Design and Implementation of a 32-to-1 down sampling filter using CIC, 
% Polyphase and half band filters to process the output of Sigma-delta
% modulator and
% Up sampling a shaped input signal using a single 4-path filter, cascade of
% two 2-path filters and cascade of two true-half band 2-path FIR filters


%% Clearing Memory
clc;
clear all;
close all;
%%
% 
% We are going to compare the workload of a number of solutions to the same
% problem. The problem is to up-sample a shaped input signal by a fixed 
% factor of 4. The input signal is the output of a 41 tap square-root Nyquist
% filter with excess bandwidth of 0.25 operating at 2-samples per symbol. 
% The output is the 1-to-4 interpolated composite impulse response. 
% The spectral artifacts formed by the interpolators must be 60 dB below peak
% amplitude response. At the output of each interpolator chain we will pass
% the interpolated signal through its matched filter, a square-root Nyquist 
% filter with same excess bandwidth operating at 8-samples per symbol and 
% collect the sum of squares of every 8-th sample point to determine the 
% signal degradation relative to the levels obtained by correlating the 
% 8-sample per symbol shaping filter with itself. Remember to normalize 
% the filter response for peak time response of unity.   
%
%  
% The interpolator options we examine are:
% 
% a.A single 4-path FIR filter.
% b.A cascade of two 2-path FIR filters designed by the remez alghorithm.
% c.A cascade of two true-half band 2-path FIR filters (see true half band trick)
% 
% 
% 
%% 41 tap square root nyquist filter 
%
% Specifications
%
% Excess Bandwidth=0.25
% Fs= 2samples/symbol

h_nyq=0.5*sqrt_nyq_y2(2,0.25,10,1);  % 41 tap square root nyquist filter
x_nyq=[1 zeros(1,100)]; % Input to Nyquist filter
y_nyq=conv(h_nyq,x_nyq); % Output of Nyquist filter

figure(4);

% Impulse Response of output of nyquist filter(Input of Interpolator)

subplot(2,1,1)
plot(y_nyq/max(y_nyq),'linewidth',2)
grid on
title('Impulse Response of output of nyquist filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response of output of nyquist filter(Input of Interpolator)

subplot(2,1,2)
f_nyq=fft(y_nyq,4096);
fh_nyq=fftshift(20*log10(abs(f_nyq)));
plot((-0.5:1/4096:0.5-1/4096)*2,fh_nyq);
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'r','linewidth',2)
plot([-1 -0.625 -0.625],[-80 -60 0],'r','linewidth',2)
plot([+1 +0.625 +0.625],[-80 -60 0],'r','linewidth',2)
hold off;
grid on
axis([-1 1 -150 10]);
title('Frequency Response of output of nyquist filter')
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)')

%% Design of a FIR filter
 %
 % Specifications
 % Specifications of this filter depends on the specifications of the
 % nyquist filter
 % Passband=[0 0.625]
 % The passband of Nyquist filter is [0 0.5] but the question says it has
 % excess bandwidth of 0.25 so the pass band of this filter becomes [0
 % 0.625].
 % Stopband=[1.375 4]
 % Since we are sampling the signal,the second copy of the signal is
 % at 2.So the stopband starts from 2-0.635 which is 1.375 inorder to get a
 % neat one sample of signal.
 
fs=8;     % Sampling frequency
df=1.375-0.625; % difference between passband and stopband
NN = floor((fs/df)*(60/22));% NN gives order of filter


h_4p=remez(31,[0 0.625 1.375 4]/4,{'myfrf',[1 1 0 0]},[1 5]);


figure(5);

% Impulse Response of FIR filter

subplot(2,1,1)
plot(h_4p/max(h_4p),'linewidth',2)
grid on
title('Impulse Response of Prototype:FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response of output of nyquist filter(Input of Interpolator)

subplot(2,1,2)
f_4p=fft(h_4p,4096);
fh_4p=fftshift(20*log10(abs(f_4p)));
plot((-0.5:1/4096:0.5-1/4096)*8,fh_4p)
hold on;
plot([-0.625 -0.625 0.625 0.625],[-90 0 0 -90],'r','linewidth',2)
plot([-4 -1.375 -1.375],[-60 -60 0],'r','linewidth',2)
plot([+4 +1.375 +1.375],[-60 -60 0],'r','linewidth',2)
hold off;
grid on
axis([-4 4 -150 10]);
title('Frequency Response of Prototype: FIR filter')
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)')

%% Single 4-path FIR filter
% Interpolator

hh_4p=reshape(h_4p,4,8);
reg_hh_4p=zeros(1,8);


m_4p=0;

yy_4p=[y_nyq zeros(1,8) ];

for n=1:1:length(yy_4p)
    
    reg_hh_4p=[yy_4p(n) reg_hh_4p(1:7)];
    
    for k=1:4
        y_4p(m_4p+k)=reg_hh_4p*hh_4p(k,:)';
    end
    
       m_4p=m_4p+4;
end


figure(6);

% Impulse Response of FIR filter

subplot(2,1,1)
stem(y_4p,'linewidth',2)
grid on
title('Impulse Response of output of 4-path FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response of output of 4-path filter(Input of Interpolator)

subplot(2,1,2)
ff_4p=fft(y_4p,4096);
ffh_4p=fftshift(20*log10(abs(ff_4p)));
plot((-0.5:1/4096:0.5-1/4096)*8,ffh_4p);
hold on
plot((-0.5:1/4096:0.5-1/4096)*8,fh_4p,'r','linewidth',2)
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'k','linewidth',2)
plot([-4 -0.625 -0.625],[-60 -60 0],'k','linewidth',2)
plot([+4 +0.625 +0.625],[-60 -60 0],'k','linewidth',2)
hold off
axis([-4 4 -200 10]);
grid on
title('Frequency Response of output of 4-path FIR filter')
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)')

%% A cascade of two 2-path FIR filters designed by the remez alghorithm.

% Design of 2-path Prototype filter(1st in cascade)

h_21p=remez(14,[0 0.625 1.375 2]/2,[1 1 0 0],[1 15]);

figure(7);
% Impulse Response
subplot(2,1,1);
plot(h_21p);
grid on;
title('Impulse Response of 2-path Prototype FIR filter(1st in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response
subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*4,fftshift(20*log10(abs(fft(h_21p,2048)))))
hold on
plot([-0.625 -0.625 0.625 0.625],[-90 0 0 -90],'r','linewidth',2)
plot([-2 -1.375 -1.375],[-60 -60 0],'r','linewidth',2)
plot([+2 +1.375 +1.375],[-60 -60 0],'r','linewidth',2)
hold off
grid on
title('Frequency Response of 2-path Prototype FIR filter(1st in cascade)');
xlabel('Frequency (kHz)');
ylabel('Log Mag (dB)');

%% Upsampling by 2 fs=4

hh_21p=reshape(2*[h_21p 0],2,8);

reg_hh_21p=zeros(1,8);



c1=0;
for n=1:length(y_nyq)
    
    reg_hh_21p=[y_nyq(n) reg_hh_21p(1:7)];

    for r=1:2
        y_21p(c1+r)=reg_hh_21p*hh_21p(r,:)';
    end
    
    c1=c1+2;
end


figure(8)

% Impulse Response

subplot(2,1,1)
plot(y_21p)
grid on;
title('Impulse Response of 2-path FIR filter(1st in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*4,fftshift(20*log10(abs(fft(y_21p/sum(y_21p),2048)))))
hold on
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'r','linewidth',2)
plot([-2 -0.625 -0.625],[-60 -60 0],'r','linewidth',2)
plot([+2 +0.625 +0.625],[-60 -60 0],'r','linewidth',2)
hold off
grid on;
title('Frequency Response of 2-path FIR filter(1st in cascade)');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% Design of 2-path Prototype Filter (2nd in Cascade)

h_22p=remez(8,[0 0.3125 1.6875 2]/2,[1 1 0 0],[1 1]);

figure(9);

% Impulse Response

subplot(2,1,1);
plot(h_22p);
grid on;
title('Impulse Response of 2-path Prototype FIR filter(2nd in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*4,fftshift(20*log10(abs(fft(h_22p,2048)))))
hold on
plot([-0.3125 -0.3125 0.3125 0.3125],[-90 0 0 -90],'r','linewidth',2)
plot([-2 -1.6875 -1.6875],[-60 -60 0],'r','linewidth',2)
plot([+2 +1.6875 +1.6875],[-60 -60 0],'r','linewidth',2)
hold off
grid on;
title('Frequency Response of 2-path Prototype FIR filter(2nd in cascade)');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% Upsampling by 2 fs=8

hh_22p=reshape(2*[h_22p 0],2,5);

reg_hh_22p=zeros(1,5);


c2=0;
for n1=1:length(y_21p)
    
    reg_hh_22p=[y_21p(n1) reg_hh_22p(1:4)];

    for r=1:2
        y_22p(c2+r)=reg_hh_22p*hh_22p(r,:)';
    end
    
    c2=c2+2;
end


figure(10)

subplot(2,1,1)
plot(y_22p)
grid on;
title('Impulse Response of 2-path FIR filter(2nd in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*8,fftshift(20*log10(abs(fft(y_22p/sum(y_22p),2048)))))
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'k','linewidth',2)
plot([-4 -0.625 -0.625],[-60 -60 0],'k','linewidth',2)
plot([+4 +0.625 +0.625],[-60 -60 0],'k','linewidth',2)
hold off
grid on
axis([-4 4 -160 10]);
title('Frequency Response of 2-path FIR filter(2nd in cascade)');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');


%% 

NN1=8
h_21phb=remez(7,[0 0.3 0.5 0.5]/0.5,[1 1 0 0],[1 100]);

hh2(2:2:2*NN1)=h_21phb;
hh2=[hh2 0];
h_21p_hb=hh2;
h_21p_hb(NN1+1)=1;


figure(11);
% Impulse Response
subplot(2,1,1);
plot(h_21p_hb);
grid on;
title('Impulse Response of 2-path Prototype FIR filter(1st in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response
subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048),fftshift(20*log10(abs(fft(0.5*h_21p_hb,2048)))))
grid on
title('Frequency Response of 2-path Prototype FIR filter(1st in cascade)');
xlabel('Frequency (kHz)');
ylabel('Log Mag (dB)');


%% Upsampling by 2 fs=4

hh_21p_hb=reshape(2*[h_21p_hb 0],2,9);

reg_hh_21p_hb=zeros(1,9);



c1=0;
for n=1:length(y_nyq)
    
    reg_hh_21p_hb=[y_nyq(n) reg_hh_21p_hb(1:8)];

    for r=1:2
        y_21p_hb(c1+r)=reg_hh_21p_hb*hh_21p_hb(r,:)';
    end
    
    c1=c1+2;
end


figure(12)

%Impulse Response

subplot(2,1,1)
plot(y_21p_hb)
grid on;
title('Impulse Response of 2-path FIR filter(1st in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*4,fftshift(20*log10(abs(fft(y_21p_hb/sum(y_21p_hb),2048)))))
grid on;
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'k','linewidth',2)
plot([-2 -0.625 -0.625],[-60 -60 0],'k','linewidth',2)
plot([+2 +0.625 +0.625],[-60 -60 0],'k','linewidth',2)
hold off
title('Frequency Response of 2-path FIR filter(1st in cascade)');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% Design of 2-path Prototype Filter (2nd in Cascade)

NN2=8;

h_22phb=remez(7,[0 0.15625 0.5 0.5]/0.5,[1 1 0 0],[1 100]);

hh3(2:2:2*NN2)=h_22phb;
hh3=[hh3 0];
h_22p_hb=hh3;
h_22p_hb(NN2+1)=1;


figure(13);

% Impulse Response

subplot(2,1,1);
plot(h_22p_hb);
grid on;
title('Impulse Response of 2-path Prototype FIR filter(2nd in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048),fftshift(20*log10(abs(fft(0.5*h_22p_hb,2048)))))
grid on;
% hold on
% plot([-0.3125 -0.3125 0.3125 0.3125],[-90 0 0 -90],'r','linewidth',2)
% plot([-0.5 -1.6875 -1.6875],[-60 -60 0],'r','linewidth',2)
% plot([+2 +1.6875 +1.6875],[-60 -60 0],'r','linewidth',2)
% hold off

title('Frequency Response of 2-path Prototype FIR filter(2nd in cascade)');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% upsampling by 2 fs=8

hh_22p_hb=reshape(2*[h_22p_hb 0],2,9);

reg_hh_22p_hb=zeros(1,9);


c2=0;
for n1=1:length(y_21p_hb)
    
    reg_hh_22p_hb=[y_21p_hb(n1) reg_hh_22p_hb(1:8)];

    for r=1:2
        y_22p_hb(c2+r)=reg_hh_22p_hb*hh_22p_hb(r,:)';
    end
    
    c2=c2+2;
end


figure(14)

subplot(2,1,1)
plot(y_22p_hb);
grid on;
title('Impulse Response of 2-path FIR filter(2nd in cascade)');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*8,fftshift(20*log10(abs(fft(y_22p_hb/sum(y_22p_hb),2048)))))
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'k','linewidth',2)
plot([-4 -0.625 -0.625],[-60 -60 0],'k','linewidth',2)
plot([+4 +0.625 +0.625],[-60 -60 0],'k','linewidth',2)
hold off
grid on
axis([-4 4 -160 10]);
title('Frequency Response of 2-path FIR filter(2nd in cascade)');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%%
% We are to design and implement a 32-to-1 down sampling filter to process 
% the output of the sigma-delta modulator. The noise floor of the modulator is 100 dB
% below full scale and the noise suppression of the filters must match the 
% 100 db level.
% 
% We consider three filter options:
% a.	A 16-to-1 CIC filter followed by an inverse sinc correcting filter
% and then a half band filter 2-to-1 down sampling filter.
% b.	A cascade of an 8-path and 4-path polyphase filters.
% c.	A cascade of 5-half band filters.

%% Sigma delta modulator

xx=0.8*cos(2*pi*(0:4196)*0.005);    % Input sinusoid
 
reg1=0; reg2=0; reg3=0;
b1=1/128;

yy=zeros(1,length(xx));             % Initialize yy to 0
for nn=1:length(xx)
qq=sign(reg3);
   yy(nn)=qq;            % Output of Sigma Delta converter
   
   sm3=0.25*reg2+reg3         -qq;
   sm2=0.125*reg1+reg2-sm3*b1 -qq;
   sm1=xx(nn)+reg1            -qq;
   
   reg1=sm1;
   reg2=sm2;
   reg3=sm3;
end

% First 500 samples of input and output signals

figure(15);
subplot(2,1,1);
plot(xx(1:500));

hold on
plot(yy(1:500),'k')
hold off

xlabel('Time');
ylabel('Amplitude');
title('first 500 samples of input and output signals');

% Windowing

w=kaiser(4197,12)';
w=w/sum(w);
xx1=xx.*w;
yy1=yy.*w;

% Windowed Spectrum of the input and output time series

subplot(2,1,2);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(yy1,4197)))),'r')
hold on
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(xx1,4197)))),'k')
hold off
grid on
axis([-32 32 -110 10]);
xlabel('Frequency');
ylabel('Gain(db)');
title('Sigma-Delta modulator Output Spectrum');


%% Steps involved 
%
% The output of Sigma Delta Modulator is given to 16 to 1 CIC filter
% The output of CIC filter is given to Inverse Sinc correcting filter
% The output of Inverse Sinc correcting filter is given as input to 2 to 1
% downsampling half band filter.
%

%% Design of 16 to 1 CIC filter
%
% When building a CIC filter for downsampling, comb filter is placed at the
% output end of the filter.
%
% CIC filter generates sequence of 'M' ones,where M is sampling factor.
%
% Reason:When delta signal is given as input to integrator it generates
% multiple ones depending on the integrator (say infinite).When this is
% given as input to comb filter(which has one at 0 and -1 at M),we end up
% with sequence of ones.
%
% Usually, CIC filter is designed in multiple stages depending on
% operations.
% Boxcar filters when put together behave as a pretty good filter whereas
% boxcar operating as a individual filter is a poor filter.(Refer Pg:336)


%%
% CIC filter stages

cic1x=ones(1,16);
cic1=cic1x/sum(cic1x);
cic2=conv(cic1,cic1);
cic3=conv(cic2,cic1);
cic4=conv(cic3,cic1);

% Output of CIC filter at each stage

yx=filter(cic1,1,yy);  % Entire time series as input to first stage
yx1=filter(cic2,1,yy); % Entire time series as input to second stage
yx2=filter(cic3,1,yy); % Entire time series as input to third stage
y=filter(cic4,1,yy);   % Entire time series as input to fourth stage 


% Impulse Responses of Four stages of CIC filters,Length-16
% Refer to page 328

figure(16);

subplot(2,2,1)
stem(cic1x);
grid on;
title('Impulse Response of first CIC filter');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,2,2)
stem(cic2/max(cic2));
grid on;
title('Impulse Response of second CIC filter');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,2,3)
stem(cic3/max(cic3));
grid on;
title('Impulse Response of third CIC filter)');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,2,4)
stem(cic4/max(cic4));
grid on;
title('Impulse Response of fourth CIC filter)');
xlabel('Time Index');
ylabel('Amplitude');

%% Observations
%
% The successive outputs are rectangle,a triangle formed by convolving two
% rectangles,a piecewise quadratic formed by convolving three rectangles
% and piecewise cubic formed by convolving four rectangles.
%
%
%
% Frequency Responses of four stages,Length-16

figure(17);

subplot(2,2,1);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic1,4197)))),'r');
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(yx.*w,4197)))));
hold off;
axis([-32 32 -100 10]);
grid on;
xlabel('Frequency');
ylabel('Gain(db)');
title('Output of first CIC filter and CIC filter');


subplot(2,2,2);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic2,4197)))),'r');
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(yx1.*w,4197)))));
hold off;
axis([-32 32 -100 10]);
grid on;
xlabel('Frequency');
ylabel('Gain(db)');
title('Output of second CIC filter and CIC filter');


subplot(2,2,3);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic3,4197)))),'r');
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(yx2.*w,4197)))));
hold off;
axis([-32 32 -100 10]);
grid on;
xlabel('Frequency');
ylabel('Gain(db)');
title('Output of third CIC filter and CIC filter');


subplot(2,2,4);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic4,4197)))),'r');
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(y.*w,4197)))));
hold off;
axis([-32 32 -100 10]);
grid on;
xlabel('Frequency');
ylabel('Gain(db)');
title('Output of fourth CIC filter and CIC filter');

% Frequency Response of output of entire CIC(fourth stage)

figure(18);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic4,4197)))),'r');
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(y.*w,4197)))));
hold off;
axis([-32 32 -110 10]);
grid on;
xlabel('Frequency');
ylabel('Gain(db)');
title('Output of CIC filter and CIC filter');


%% Inverse sinc correcting filter
%
% Inverse sinc filter is used for compensating the roll off CIC filter in
% passband by lettig the CIC filter followed by symmetric FIR filter with a
% minimum order.
% (Refer to document you have)

cor=-0.346; % Value choosen by observation

% Coefficients of Compensation filter

com_fil=[cor 0 0 0 0 0 0 0 0  1 0 0 0 0 0 0 0 0  cor];
com_fil=com_fil/sum(com_fil);

cic_comp=conv(cic4,com_fil);


% Frequency Response of CIC filter(last stage) and Correcting filter

figure(19);
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic4,4197)))),'r');
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic_comp,4197)))));
hold off;
grid on;
axis([-0.5 0.5 -0.5 0.5]);
xlabel('Frequency');
ylabel('Gain(db)');
title('Zoomed Frequency Response of CIC filter and Correcting filter');


% Output of CIC Compensating filter and the filter

y2=filter(cic_comp,1,y);

figure(20);

plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(y2.*w,4197)))));
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(cic_comp,4197)))),'r');
hold off;
grid on;
axis([-32 32 -120 10]);
xlabel('Frequency');
ylabel('Gain(db)');
title(' Frequency Response of Correcting filter and its output');



%% 2 to 1 downsampling half band filter

% Filter design using Remez

h_2p=remez(18,[0 0.5 1.5 2]/2,{'myfrf',[1 1 0 0]},[1 50]);

figure(21);

% Impulse Response

subplot(2,2,1);
stem(h_2p);
grid on;
title('Impulse Response of 2-path Prototype FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(2,2,2)
plot((-0.5:1/4197:0.5-1/4197)*4,fftshift(20*log10(abs(fft(h_2p,4197)))))
grid on;
axis([-2 2 -150 10]);
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'r','linewidth',2)
plot([-2 -1.5 -1.5],[-100 -100 0],'r','linewidth',2)
plot([+2 +1.5 +1.5],[-100 -100 0],'r','linewidth',2)
hold off;
title('Frequency Response of 2-path Prototype FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

% Taking every 8th sample

yd=y2(1:16:length(y2));

y3=filter(h_2p,1,yd);



% Windowing

w1=kaiser(263,12)';
w1=w1/sum(w1);

% Impulse Response

subplot(2,2,3);
plot(y3);
grid on;
title('Impulse Response of output of 2:1 downsampler');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,2,4);
plot((-0.5:1/263:0.5-1/263)*4,fftshift(20*log10(abs(fft(y3.*w1,263)))),'r');
grid on;
hold on
plot((-0.5:1/4197:0.5-1/4197)*4,fftshift(20*log10(abs(fft(h_2p,4197)))))
hold off;
axis([-2 2 -120 10]);
title('Frequency Response of output of 2:1 downsampler');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% 2b Cascade of an 8-path and 4-path Polyphase filters

% Design of 8-path Polyphase filter

fs=64;

h_8p=remez(37,[0 0.5 7.5 8.5 15.5 16.5 23.5 24.5 31.5 32]/32,...
[1 1 0 0 0 0 0 0 0 0 ]);

figure(22);

% Impulse Response

subplot(3,2,1);
plot(h_8p);
grid on;
title('Impulse Response of 8-path Prototype FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(h_8p,4197)))))
grid on;
axis([-32 32 -150 10]);
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'r','linewidth',2)
plot([-8.5 -8.5 -7.5 -7.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([+8.5 +8.5 +7.5 +7.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([-16.5 -16.5 -15.5 -15.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([+16.5 +16.5 +15.5 +15.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([-24.5 -24.5 -23.5 -23.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([+24.5 +24.5 +23.5 +23.5],[-80 -100 -100 -80],'r','linewidth',2)
hold off;
axis([-32 32 -110 10]);
title('Frequency Response of 8-path Prototype FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

% Passing output of sigma delta modulator to 8-path Polyphase filter

yy_8p=filter(h_8p,1,yy);


subplot(3,2,2)
plot(yy_8p);
grid on;
title('Impulse Response of 8-path Polyphase FIR filter');
xlabel('Time Index');
ylabel('Amplitude');


subplot(3,2,[5 6])
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(yy_8p.*w,4197)))))
grid on;
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(h_8p,4197)))),'r');
hold off;
axis([-32 32 -110 10]);
title('Frequency Response of output of 8-path Polyphase FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');


%% 4-path Polyphase filter

% 4-path Prototype filter

h_4p=remez(40,[0 0.5 1.5 2.5 3.5 4]/4,[1 1 0 0 0 0],[1 50 50]);


figure(23);

% Impulse Response

subplot(3,2,1);
plot(h_4p);
grid on;
title('Impulse Response of 4-path Prototype FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*8,fftshift(20*log10(abs(fft(h_4p,4197)))))
grid on;
hold on;
plot([-0.5 -0.5 0.5 0.5],[-90 0 0 -90],'r','linewidth',2)
plot([-2.5 -2.5 -1.5 -1.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([+2.5 +2.5 +1.5 +1.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([-4 -4 -3.5 -3.5],[-80 -100 -100 -80],'r','linewidth',2)
plot([+4 +4 +3.5 +3.5],[-80 -100 -100 -80],'r','linewidth',2)
hold off;
title('Frequency Response of 4-path Prototype FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

% Taking every 8th sample


y_4p=yy_8p(1:8:length(yy_8p));
yy_4p=filter(h_4p,1,y_4p);


subplot(3,2,2)
plot(yy_4p);
grid on;
title('Impulse Response of 4-path Polyphase FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

wxx=kaiser(525,12)';
wxx=wxx/sum(wxx);


subplot(3,2,[5 6])
plot((-0.5:1/525:0.5-1/525)*8,fftshift(20*log10(abs(fft(yy_4p.*wxx,525)))));
hold on;
plot((-0.5:1/4197:0.5-1/4197)*8,fftshift(20*log10(abs(fft(h_4p,4197)))),'r');
hold off;
grid on;
axis([-4 4 -120 10]);
title('Frequency Response of 4-path Polyphase FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');



%% Cascade of 5 half band filters

% design of first half band filter

hh1_hb=remez(7,[0 0.5 31.5 32]/32,{'myfrf',[1 1 0 0]});

yy2=filter(hh1_hb,1,yy);

figure(24);

% Impulse Response

subplot(3,2,1);
plot(hh1_hb);
grid on;
title('Impulse Response of first half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');


% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(hh1_hb,4197)))))
grid on;
axis([-35 35 -120 5]);
title('Frequency Response of first half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

% Impulse Response

subplot(3,2,2);
plot(yy2);
grid on;
title('Impulse Response of output of first half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');


subplot(3,2,[5 6])
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(yy2.*w,4197)))))
hold on;
plot((-0.5:1/4197:0.5-1/4197)*64,fftshift(20*log10(abs(fft(hh1_hb,4197)))),'r')
hold off;
grid on;
axis([-32 32 -120 10]);
title('Frequency Response of first half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');


%% design of second half band filter

hh2_hb=remez(6,[0 0.5 15.5 16]/16,[1 1 0 0]);

yy2d=yy2(1:2:length(yy2));

yy3=filter(hh2_hb,1,yy2d);

figure(25);

% Impulse Response

subplot(3,2,1);
plot(hh2_hb);
grid on;
title('Impulse Response of second half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*32,fftshift(20*log10(abs(fft(hh2_hb,4197)))))
grid on;
axis([-16 16 -200 10])
title('Frequency Response of second half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

wp=kaiser(2099,12)';
wp=wp/sum(wp);

% Impulse Response

subplot(3,2,2);
plot(yy3);
grid on;
title('Impulse Response of output of second half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,2,[5 6])
plot((-0.5:1/2099:0.5-1/2099)*32,fftshift(20*log10(abs(fft(yy3.*wp,2099)))))
hold on;
plot((-0.5:1/4197:0.5-1/4197)*32,fftshift(20*log10(abs(fft(hh2_hb,4197)))),'r')
grid on;
axis([-16 16 -120 10])
title('Frequency Response of output of second half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% design of third half band filter


hh3_hb=remez(8,[0 0.5 7.5 8]/8,{'myfrf',[1 1 0 0]},[1 5]);

yy3d=yy3(1:2:length(yy3));

yy4=filter(hh3_hb,1,yy3d);

figure(26);

% Impulse Response

subplot(3,2,1);
plot(hh3_hb);
grid on;
title('Impulse Response of third half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*16,fftshift(20*log10(abs(fft(hh3_hb,4197)))))
axis([-8 8 -200 10])
grid on;
title('Frequency Response of third half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

% Impulse Response

subplot(3,2,2);
plot(yy4);
grid on;
title('Impulse Response of output of third half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');




wq=kaiser(1050,12)';
wq=wq/sum(wq);

subplot(3,2,[5 6])
plot((-0.5:1/1050:0.5-1/1050)*16,fftshift(20*log10(abs(fft(yy4.*wq,1050)))))

hold on;
plot((-0.5:1/4197:0.5-1/4197)*16,fftshift(20*log10(abs(fft(hh3_hb,4197)))),'r')
grid on;
axis([-8 8 -120 10]);
title('Frequency Response of output of third half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% design of fourth half band filter

hh4_hb=remez(12,[0 0.5 3.5 4]/4,{'myfrf',[1 1 0 0]},[1 2]);


yy4d=yy4(1:2:length(yy4));

yy5=filter(hh4_hb,1,yy4d);

figure(27);

% Impulse Response

subplot(3,2,1);
plot(hh4_hb);
grid on;
title('Impulse Response of fourth half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');


% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*8,fftshift(20*log10(abs(fft(hh4_hb,4197)))))
grid on;
axis([-4 4 -200 10])
title('Frequency Response of fourth half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');


% Impulse Response

subplot(3,2,2);
plot(yy5);
grid on;
title('Impulse Response of output of fourth half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');



wr=kaiser(525,12)';
wr=wr/sum(wr);

subplot(3,2,[5 6])
plot((-0.5:1/4197:0.5-1/4197)*8,fftshift(20*log10(abs(fft(hh4_hb,4197)))),'r')
hold on;
plot((-0.5:1/525:0.5-1/525)*8,fftshift(20*log10(abs(fft(yy5.*wr,525)))))
grid on;
axis([-4 4 -120 10])
title('Frequency Response of fourth half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

%% design of fifth half band filter

hh5_hb=remez(12+9,[0 0.5 1.5 2]/2,{'myfrf',[1 1 0 0]},[1 5]);


yy5d=yy5(1:2:length(yy5));

yy6=filter(hh5_hb,1,yy5d);

figure(28);

% Impulse Response

subplot(3,2,1)
plot(hh5_hb);
grid on;
title('Impulse Response of fifth half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');


% Frequency Response

subplot(3,2,[3 4])
plot((-0.5:1/4197:0.5-1/4197)*4,fftshift(20*log10(abs(fft(hh5_hb,4197)))))
grid on;
axis([-2 2 -120 10])
title('Frequency Response of fifth half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');

% Impulse Response

subplot(3,2,2)
plot(yy6);
grid on;
title('Impulse Response of output of fifth half band FIR filter');
xlabel('Time Index');
ylabel('Amplitude');

ws=kaiser(263,12)';
ws=ws/sum(ws);

subplot(3,2,[5 6])
plot((-0.5:1/4197:0.5-1/4197)*4,fftshift(20*log10(abs(fft(hh5_hb,4197)))),'r')
hold on;
plot((-0.5:1/263:0.5-1/263)*4,fftshift(20*log10(abs(fft(yy6.*ws,263)))))
grid on;
axis([-2 2 -120 10])
title('Frequency Response of output of fifth half band FIR filter');
xlabel('Frequency (kHz)')
ylabel('Log Mag (dB)');




