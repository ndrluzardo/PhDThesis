% Trace Conditioned Inhibition
% US-ITI condition, Delay conditioning

close all

% Simulation is divided into the following timed intervals:
% 1 - The ITI;
% 2 - the CS;


%---DDM constants
N_CX_ITI=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX1
N_CX_CS=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX2
N_CS=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CS
%---

%---US distribution during ITI
probMatrix=rand(cycle_num,round((CX_dur+CG_dur-CS_dur)/h));
UStimes=(probMatrix<(1/(30/h)))*1;
%---

counterDDMCX_ITI=0;
counterDDMCX_CS=0;
counterDDMCS=0;

%---
ITIUS.Del.CX.A_ITI(:)=1/30;
ITIUS.Del.CX.V(:)=120/30;
%---


for trial=1:cycle_num
%     
%     %---ITI
%     
%     interval_length=round((CX_dur+CG_dur-CS_dur)/h); % 350 seconds
%     
%     %--initialize values for timer and stimuli representations
%     P_CX=zeros(1,interval_length);
%     x_CX=zeros(1,interval_length);
%     %--
%     
%     %--Obtain A and V values from memory
%     A_CX=ITIUS.Del.CX.A_ITI(trial); % A trial value for CX
    V_CX=ITIUS.Del.CX.V(trial); % V trial value for CX
%     %--
%     
    %--CR initialization
    CR=zeros(1,round((CS_dur+CG_dur)/h));
    %--

%     for t=1:interval_length
%         
%         counterDDMCX_ITI=counterDDMCX_ITI+1; % update counter for random process in DDM
%         
%         % min will take the minimum value: either DDM result or 3. This
%         % caps the value of integrator at 3.
%         P_CX(t+1)=min(DDM( P_CX(t), A_CX, h, m, N_CX_ITI(counterDDMCX_ITI) ), 3);
%         
%         % max ensures the minimum value the accumulator can reach is
%         % Aini. This avoids division by zero later.
%         P_CX(t+1)=max(P_CX(t+1), Aini);
%         
%         %---Stimuli representations (RBFs)
%         x_CX(t)=CStrace(P_CX(t+1),mu,sigma,tau_x,1,x_CX(t),h);
%         %---
%         
%         %---Slope Correction
%         A_CX=A_CX+UStimes(trial,t)*A_CX*alpha_t*(1-P_CX(t+1))/P_CX(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
%         %---
%         
%         %---V update
%         V_CX=(1-UStimes(trial,t))*V_CX+UStimes(trial,t)*RW(V_CX,alpha_E,x_CX(t),H,A_CX,P_CX(t+1));
%         %---
%         
%         %---P value
%         P_CX(t+1)=(1-UStimes(trial,t))*P_CX(t+1)+UStimes(trial,t)*Aini;
%         %---
%         
%         
%     end
%     
%     %---Slope Correction (for end of ITI)
%     A_CX=A_CX+A_CX*alpha_t*(1-P_CX(t+1))/P_CX(t+1);  
%     %---
%     
%     %---Move every A and V value 1 trial forward
%     ITIUS.Del.CX.A_ITI(trial+1)=A_CX; % A trial value for CX
%     %---
    
    %---End of ITI
    
    %---CS1
    
    interval_length=round((CS_dur)/h);
    
    %--initialize values for timer and stimuli representations
    P_CX=zeros(1,interval_length);
    P_CS=zeros(1,interval_length);
    
    x_CX=zeros(1,interval_length);
    x_CS=zeros(1,interval_length);
    %--
    
    %--Obtain A and V values from memory
    A_CX=ITIUS.Del.CX.A_CS(trial); % A trial value for CX
    A_CS=ITIUS.Del.CS.A(trial); % A trial value for CS
    V_CS=ITIUS.Del.CS.V(trial); % V trial value for CS
    %--
    
    for t=1:interval_length
        
        counterDDMCX_CS=counterDDMCX_CS+1;
        counterDDMCS=counterDDMCS+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_CX(t+1)=min(DDM( P_CX(t), A_CX, h, m, N_CX_CS(counterDDMCX_CS) ), 3);
        P_CS(t+1)=min(DDM( P_CS(t), A_CS, h, m, N_CS(counterDDMCS) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % Aini. This avoids division by zero later.
        P_CX(t+1)=max(P_CX(t+1), Aini);
        P_CS(t+1)=max(P_CS(t+1), Aini);
        
        %---Stimuli representations (RBFs)
        x_CX(t)=CStrace(P_CX(t+1),mu,sigma,tau_x,1,x_CX(t),h);
        x_CS(t)=CStrace(P_CS(t+1),mu,sigma,tau_x,1,x_CS(t),h);
        %---
        
        %---CR
        CR(t)=max(0,[x_CS(t) x_CX(t)]*[V_CS V_CX]');
        %---
        
    end
    
    %---Slope Correction
    A_CX=A_CX+A_CX*alpha_t*(1-P_CX(t+1))/P_CX(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_CS=A_CS+A_CS*alpha_t*(1-P_CS(t+1))/P_CS(t+1);
    %---
    
    %---V update
    V_CX=RW([V_CX V_CS],alpha_E,[x_CX(t) x_CS(t)],H,A_CX,P_CX(t+1));
    V_CS=RW([V_CS V_CX],alpha_E,[x_CS(t) x_CX(t)],H,A_CS,P_CS(t+1));
    %---
    
    %---Move every A and V value 1 trial forward
    ITIUS.Del.CX.A_CS(trial+1)=A_CX; % A trial value for CX
    ITIUS.Del.CS.A(trial+1)=A_CS; % A trial value for CS
    %ITIUS.Del.CX.V(trial+1)=V_CX; % V trial value for CX
    ITIUS.Del.CS.V(trial+1)=V_CS; % V trial value for CS
    ITIUS.Del.CR(trial,:)=CR;
    %---
   
end


%--calculate average response rate
MeanCR=mean(ITIUS.Del.CR((end-20):end,:),1);
%---

%--plots
% Associative strength CX, CS during CS

subplot(2,3,1)
plot(1:length(ITIUS.Del.CS.V),ITIUS.Del.CS.V,1:length(ITIUS.Del.CX.V),ITIUS.Del.CX.V)
legend('V CS','V CX')

subplot(2,3,2)
plot(1:length(ITIUS.Del.CX.A_ITI)-10,1./ITIUS.Del.CX.A_ITI(11:end))
legend('context time estimate: 30')

subplot(2,3,3)
plot(1:length(ITIUS.Del.CX.A_CS)-10,1./ITIUS.Del.CX.A_CS(11:end))
legend('context time estimate to reward: 120')

subplot(2,3,4)
plot(1:length(ITIUS.Del.CS.A)-10,1./ITIUS.Del.CS.A(11:end))
legend('CS time estimate to reward: 120')

subplot(2,3,5)
plot(MeanCR)

