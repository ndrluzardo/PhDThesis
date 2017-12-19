% Trace Conditioned Inhibition
% Simulation of experiment 1 in: Williams, D. A., Todd, T. P., Chubala, C. M., & Ludvig, E. A. (2016). Intertrial unconditioned stimuli differentially impact trace conditioning. Learning & Behavior, 1?13. 

clear
close all

% Experiment consist of three groups:
% Embedded: US onset is 10 seconds before CS offset.
% Delay: US onset coincides with CS offset.
% Trace: US onset is 10 seconds after CS offset.
%
% Each group has 2 conditions:
% no-ITI US: There is no US in the ITI.
% ITI US: There are random presentations of USs (one per 30 seconds on average) in the ITI.
%
% CS lasts for 120 seconds, ITI is 340 seconds. Note that the paper used a
% non-uniform distribution of ITIs, but it is not clear what they mean by that.

%% 

%---cycles
cycle_num=400;
%---

%---CS and CX (context) duration in seconds
CS_dur=120;
CX_dur=460;
CG_dur=10; % trace gap
%---

%---parameters
h=0.01;
tau_x=50;
alpha_t=0.1;
alpha_E=0.07;
mu=1;
sigma=0.4;
m=0.15;
H=40;
Aini=1*10^(-3);
%---


%% No-ITI US condition, Delay conditioning group

%--Initialize structures
NoITIUS.Del.CS.A=zeros(cycle_num,1);
NoITIUS.Del.CS.A(1)=Aini;
NoITIUS.Del.CS.V=zeros(cycle_num,1);
NoITIUS.Del.CX.A_ITI=zeros(cycle_num,1); % context A for ITI
NoITIUS.Del.CX.A_CS=zeros(cycle_num,1); % context A for CS
NoITIUS.Del.CX.A_ITI(1)=Aini;
NoITIUS.Del.CX.A_CS(1)=Aini;
NoITIUS.Del.CX.V=zeros(cycle_num,1); % context V
NoITIUS.Del.CR=zeros(cycle_num,round((CS_dur+CG_dur)/h));
%---

%---Run script
TrInhNoITIUSDel
%---

%% No-ITI US condition, Embedded conditioning group

%---Initialize structures
NoITIUS.Emb.CS.A1=zeros(cycle_num,1);
NoITIUS.Emb.CS.A2=zeros(cycle_num,1);
NoITIUS.Emb.CS.A3=zeros(cycle_num,1);
NoITIUS.Emb.CS.A1(1)=Aini;
NoITIUS.Emb.CS.A2(1)=Aini;
NoITIUS.Emb.CS.A3(1)=Aini;

NoITIUS.Emb.CS.V=zeros(cycle_num,1);
NoITIUS.Emb.CS.V2=zeros(cycle_num,1);
NoITIUS.Emb.CX.A_ITI=zeros(cycle_num,1);
NoITIUS.Emb.CX.A_CS1=zeros(cycle_num,1);
NoITIUS.Emb.CX.A_CS2=zeros(cycle_num,1);
NoITIUS.Emb.CX.A_ITI(1)=Aini;
NoITIUS.Emb.CX.A_CS1(1)=Aini;
NoITIUS.Emb.CX.A_CS2(1)=Aini;
NoITIUS.Emb.CX.V=zeros(cycle_num,1);
NoITIUS.Emb.CR=zeros(cycle_num,round((CS_dur+CG_dur)/h));
%---

%---Run script
TrInhNoITIUSEmb
%---

% Problem here: there is a strong response during the gap between US and
% CSoffset. There should be a separate timer for this interval, and the CX
% should become inhibitory. The CX.V does not become inhibitory enough.

%% No-ITI US condition, Trace conditioning group

%---Initialize structures
NoITIUS.Tra.CS.A=zeros(cycle_num,1);
NoITIUS.Tra.CS.A(1)=Aini;
NoITIUS.Tra.CS.V=zeros(cycle_num,1);
NoITIUS.Tra.CX.A1=zeros(cycle_num,1);
NoITIUS.Tra.CX.A2=zeros(cycle_num,1);
NoITIUS.Tra.CX.A3=zeros(cycle_num,1);
NoITIUS.Tra.CX.A1(1)=Aini;
NoITIUS.Tra.CX.A2(1)=Aini;
NoITIUS.Tra.CX.A3(1)=Aini;
NoITIUS.Tra.CX.V=zeros(cycle_num,1);
NoITIUS.Tra.CG.A=zeros(cycle_num,1);
NoITIUS.Tra.CG.A(1)=Aini;
NoITIUS.Tra.CG.V=zeros(cycle_num,1);
NoITIUS.Tra.CR=zeros(cycle_num,round((CS_dur+CG_dur)/h));
%---

%---Run script
TrInhNoITIUSTra
%---

%% US-ITI condition, Delay conditioning

%---Initialize structures
ITIUS.Del.CS.A=zeros(cycle_num,1);
ITIUS.Del.CS.A(1)=Aini;
ITIUS.Del.CS.V=zeros(cycle_num,1);
ITIUS.Del.CX.A_ITI=zeros(cycle_num,1);
ITIUS.Del.CX.A_ITI(1)=Aini;
ITIUS.Del.CX.A_CS=zeros(cycle_num,1);
ITIUS.Del.CX.A_CS(1)=Aini;
ITIUS.Del.CX.V=zeros(cycle_num,1);
ITIUS.Del.CR=zeros(cycle_num,round((CS_dur+CG_dur)/h));
%---

%---Run script
TrInhITIUSDel
%---

%% US-ITI condition, Embedded conditioning

%---Initialize structures
ITIUS.Emb.CS.A_S1=zeros(cycle_num,1);
ITIUS.Emb.CS.A_S1(1)=Aini;
ITIUS.Emb.CS.A_S2=zeros(cycle_num,1);
ITIUS.Emb.CS.A_S2(1)=Aini;
ITIUS.Emb.CS.V=zeros(cycle_num,1);
ITIUS.Emb.CX.A_ITI=zeros(cycle_num,1);
ITIUS.Emb.CX.A_ITI(1)=Aini;
ITIUS.Emb.CX.A_S1=zeros(cycle_num,1);
ITIUS.Emb.CX.A_S1(1)=Aini;
ITIUS.Emb.CX.A_S2=zeros(cycle_num,1);
ITIUS.Emb.CX.A_S2(1)=Aini;
ITIUS.Emb.CX.V=zeros(cycle_num,1);
ITIUS.Emb.CR=zeros(cycle_num,round((CS_dur+CG_dur)/h));
%---

%---Run script
TrInhITIUSEmb
%---

%% US-ITI condition, Trace conditioning

%---Initialize structures
ITIUS.Tra.CS.A=zeros(cycle_num,1);
ITIUS.Tra.CS.A(1)=Aini;
ITIUS.Tra.CS.V=zeros(cycle_num,1);
ITIUS.Tra.CX.A_ITI=zeros(cycle_num,1);
ITIUS.Tra.CX.A_ITI(1)=Aini;
ITIUS.Tra.CX.A_CS=zeros(cycle_num,1);
ITIUS.Tra.CX.A_CS(1)=Aini;
ITIUS.Tra.CX.A_G=zeros(cycle_num,1);
ITIUS.Tra.CX.A_G(1)=Aini;
ITIUS.Tra.CX.V=zeros(cycle_num,1);
ITIUS.Tra.CG.A=zeros(cycle_num,1);
ITIUS.Tra.CG.A(1)=Aini;
ITIUS.Tra.CG.V=zeros(cycle_num,1);
ITIUS.Tra.CR=zeros(cycle_num,round((CS_dur+CG_dur)/h));
%---

%---Run script
TrInhITIUSTra
%---

%% close, save and clear
close all
save('TraceInhResults','NoITIUS','ITIUS','h')
clear
