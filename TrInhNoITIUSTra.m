% Trace Conditioned Inhibition
% No-ITI US condition, Trace conditioning group

close all

% Simulation is divided into three timed intervals:
% 1 - The ITI;
% 2 - the CS;
% 3 - the 10 second gap between CS offset and US.

%---DDM constants
N_CX1=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX1
N_CX2=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX2
N_CX3=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX3
N_CS=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CS
N_CG=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for gap
%---


%%

counterDDMCX1=0;
counterDDMCX2=0;
counterDDMCX3=0;
counterDDMCS=0;
counterDDMCG=0;

for trial=1:cycle_num
    
    %---ITI
    
    interval_length=round((CX_dur-CG_dur-CS_dur)/h); 
    
    
    %--initialize values for timer and stimuli representations
    P_CX=zeros(1,interval_length);
    x_CX=zeros(1,interval_length);
    %--
    
    %--Obtain A and V values from memory
    A_CX=NoITIUS.Tra.CX.A1(trial); % A trial value for CX
    V_CX=NoITIUS.Tra.CX.V(trial); % V trial value for CX
    %--
    
    for t=1:interval_length
        
        counterDDMCX1=counterDDMCX1+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_CX(t+1)=min(DDM( P_CX(t), A_CX, h, m, N_CX1(counterDDMCX1) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % Aini. This avoids division by zero later.
        P_CX(t+1)=max(P_CX(t+1), Aini);
        
        %---Stimuli representations (RBFs)
        x_CX(t)=CStrace(P_CX(t+1),mu,sigma,tau_x,1,x_CX(t),h);
        %---
        
    end
    
    %---Slope Correction
    A_CX=A_CX+A_CX*alpha_t*(1-P_CX(t+1))/P_CX(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_CX=RW(V_CX,alpha_E,x_CX(t),0,A_CX,P_CX(t+1));
    %---
    
    %---Move every A and V value 1 trial forward
    NoITIUS.Tra.CX.A1(trial+1)=A_CX; % A trial value for CX
    %---
    
    %---End of ITI
    
    %---CS
    
    interval_length=round(CS_dur/h);
    
    %--initialize values for timer and stimuli representations
    P_CX=zeros(1,interval_length);
    P_CS=zeros(1,interval_length);
    
    x_CX=zeros(1,interval_length);
    x_CS=zeros(1,interval_length);
    
    CR=zeros(1,round((CS_dur+CG_dur)/h));
    %--
    
    %--Obtain A and V values from memory
    A_CX=NoITIUS.Tra.CX.A2(trial); % A trial value for CX
    A_CS=NoITIUS.Tra.CS.A(trial); % A trial value for CS
    V_CS=NoITIUS.Tra.CS.V(trial); % V trial value for CS
    %--
    
    for t=1:interval_length
        
        counterDDMCX2=counterDDMCX2+1;
        counterDDMCS=counterDDMCS+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_CX(t+1)=min(DDM( P_CX(t), A_CX, h, m, N_CX2(counterDDMCX2) ), 3);
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
    V_CX=RW([V_CX V_CS],alpha_E,[x_CX(t) x_CS(t)],0,A_CX,P_CX(t+1));
    V_CS=RW([V_CS V_CX],alpha_E,[x_CS(t) x_CX(t)],0,A_CS,P_CS(t+1));
    %---
    
    %---Move every A and V value 1 trial forward
    NoITIUS.Tra.CX.A2(trial+1)=A_CX; % A trial value for CX
    NoITIUS.Tra.CS.A(trial+1)=A_CS; % A trial value for CS
    NoITIUS.Tra.CS.V(trial+1)=V_CS; % V trial value for CS
    %---
    %---End of CS
    
    %---Begin gap
    
    interval_length=round(CG_dur/h);
    
    %--initialize values for timer and stimuli representations
    P_CX=zeros(1,interval_length);
    P_CG=zeros(1,interval_length);

    x_CX=zeros(1,interval_length);
    x_CG=zeros(1,interval_length);
    %--
    
    %--Obtain A and V values from memory
    A_CX=NoITIUS.Tra.CX.A3(trial); % A trial value for CX
    A_CG=NoITIUS.Tra.CG.A(trial); % A trial value for CG
    V_CG=NoITIUS.Tra.CG.V(trial); % V trial value for CG
    %--
    
    for t2=1:interval_length
        
        counterDDMCX3=counterDDMCX3+1;
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_CX(t2+1)=min(DDM( P_CX(t2), A_CX, h, m, N_CX3(counterDDMCX3) ),3);
        P_CG(t2+1)=min(DDM( P_CG(t2), A_CG, h, m, N_CG(counterDDMCX3) ),3); 

        % max ensures the minimum value the accumulator can reach is
        % Aini. This avoids division by zero later.
        P_CX(t2+1)=max(P_CX(t2+1), Aini);
        P_CG(t2+1)=max(P_CG(t2+1), Aini);
        
        %---Stimuli representations (RBFs)
        x_CX(t2)=CStrace(P_CX(t2+1),mu,sigma,tau_x,1,x_CX(t2),h);
        x_CG(t2)=CStrace(P_CG(t2+1),mu,sigma,tau_x,1,x_CG(t2),h);
        %---
        
        %---CR
        CR(t+t2)=max(0,[x_CX(t2) x_CG(t2)]*[V_CX V_CG]');
        %---
        
    end
    
    %---Slope Correction
    A_CX=A_CX+A_CX*alpha_t*(1-P_CX(t2+1))/P_CX(t2+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_CG=A_CG+A_CG*alpha_t*(1-P_CG(t2+1))/P_CG(t2+1);
    %---
    
    %---V update
    V_CX=RW([V_CX V_CG],alpha_E,[x_CX(t2) x_CG(t2)],H,A_CX,P_CX(t2+1));
    V_CG=RW([V_CG V_CX],alpha_E,[x_CG(t2) x_CX(t2)],H,A_CG,P_CG(t2+1));
    %---
    
    %---Move every A and V value 1 trial forward
    NoITIUS.Tra.CX.A3(trial+1)=A_CX; % A trial value for CX
    NoITIUS.Tra.CX.V(trial+1)=V_CX; % V trial value for CX
    NoITIUS.Tra.CG.A(trial+1)=A_CG;
    NoITIUS.Tra.CG.V(trial+1)=V_CG;
    NoITIUS.Tra.CR(trial,:)=CR;
    %---
   
end


%--calculate average response rate
MeanCR=mean(NoITIUS.Tra.CR((end-40):end,:),1);
%---

%--plots
% Associative strength CX, CS during CS

subplot(2,3,1)
plot(1:length(NoITIUS.Tra.CS.V),NoITIUS.Tra.CS.V,1:length(NoITIUS.Tra.CX.V),NoITIUS.Tra.CX.V)
legend('V CS','V CX')

subplot(2,3,2)
plot(1:length(NoITIUS.Tra.CX.A1)-10,1./NoITIUS.Tra.CX.A1(11:end))
legend('CX time est. to CS: 340')

subplot(2,3,3)
plot(1:length(NoITIUS.Tra.CX.A2)-10,1./NoITIUS.Tra.CX.A2(11:end))
legend('CX time est. to CS off: 120')

subplot(2,3,4)
plot(1:length(NoITIUS.Tra.CS.A)-10,1./NoITIUS.Tra.CS.A(11:end))
legend('CS time est. to CS off: 120')

subplot(2,3,5)
plot(1:length(NoITIUS.Tra.CX.A3)-10,1./NoITIUS.Tra.CX.A3(11:end))
legend('CX time est. to US: 10')

subplot(2,3,6)
plot(MeanCR)
legend('Mean CR')