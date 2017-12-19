% Trace Conditioned Inhibition
% US-ITI condition, Embedded conditioning

close all

% Simulation is divided into three timed intervals:
% 1 - The ITI;
% 2 - the CS;
% 3 - the 10 second period between US and CS offset. 

%---DDM constants
N_CX1=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX1
N_CX2=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX2
N_CX3=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CX3
N_CS1=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CS1
N_CS2=normrnd(0,1,ceil((CX_dur/h)*cycle_num),1); % noise for CS2
%---

counterDDMCX1=0;
counterDDMCX2=0;
counterDDMCX3=0;
counterDDMCS1=0;
counterDDMCS2=0;

%---
ITIUS.Emb.CX.A_ITI(:)=1/30;
ITIUS.Emb.CX.V(:)=120/30;
%---

for trial=1:cycle_num
    
    %---ITI
    
%     interval_length=round((CX_dur+CG_dur-CS_dur)/h); % 350 seconds
%     
%     
%     %--initialize values for timer and stimuli representations
%     P_CX=zeros(1,interval_length);
%     x_CX=zeros(1,interval_length);
%     %--
%     
%     %--Obtain A and V values from memory
%     A_CX=NoITIUS.Emb.CX.A_ITI(trial); % A trial value for CX
    V_CX=ITIUS.Emb.CX.V(trial); % V trial value for CX
%     %--
%     
    %--CR initialization
    CR=zeros(1,round((CS_dur+CG_dur)/h));
    %--
%     
%     for t=1:interval_length
%         
%         counterDDMCX1=counterDDMCX1+1; % update counter for random process in DDM
%         
%         % min will take the minimum value: either DDM result or 3. This
%         % caps the value of integrator at 3.
%         P_CX(t+1)=min(DDM( P_CX(t), A_CX, h, m, N_CX1(counterDDMCX1) ), 3);
%         
%         % max ensures the minimum value the accumulator can reach is
%         % Aini. This avoids division by zero later.
%         P_CX(t+1)=max(P_CX(t+1), Aini);
%         
%         %---Stimuli representations (RBFs)
%         x_CX(t)=CStrace(P_CX(t+1),mu,sigma,tau_x,1,x_CX(t),h);
%         %---
%         
%         %---CR
%         if t<=round(CG_dur/h)
%             CR(end-round(CG_dur/h)+t)=max(0,x_CX(t)*V_CX);
%         end
%         %---
%         
%     end
%     
%     %---Slope Correction
%     A_CX=A_CX+A_CX*alpha_t*(1-P_CX(t+1))/P_CX(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
%     %---
%     
%     %---V update
%     V_CX=RW(V_CX,alpha_E,x_CX(t),0,A_CX,P_CX(t+1));
%     %---
%     
%     %---Move every A and V value 1 trial forward
%     NoITIUS.Emb.CX.A_ITI(trial+1)=A_CX; % A trial value for CX
%     %---
    
    %---End of ITI
    
    %---CS, 1st part
    
    interval_length=round((CS_dur-10)/h);
    
    %--initialize values for timer and stimuli representations
    P_CX=zeros(1,interval_length);
    P_CS=zeros(1,interval_length);
    
    x_CX=zeros(1,interval_length);
    x_CS=zeros(1,interval_length);
    %--
    
    %--Obtain A and V values from memory
    A_CX=ITIUS.Emb.CX.A_S1(trial); % A trial value for CX
    A_CS=ITIUS.Emb.CS.A_S1(trial); % A trial value for CS
    V_CS=ITIUS.Emb.CS.V(trial); % V trial value for CS
    %--
    
    for t=1:interval_length
        
        counterDDMCX2=counterDDMCX2+1;
        counterDDMCS1=counterDDMCS1+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_CX(t+1)=min(DDM( P_CX(t), A_CX, h, m, N_CX2(counterDDMCX2) ), 3);
        P_CS(t+1)=min(DDM( P_CS(t), A_CS, h, m, N_CS1(counterDDMCS1) ), 3);
        
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
    ITIUS.Emb.CX.A_S1(trial+1)=A_CX; % A trial value for CX
    ITIUS.Emb.CS.A_S1(trial+1)=A_CS; % A trial value for CS
    %---
    %---End of CS, 1st part
    
    %---Begin CS, 2nd part
    
%     interval_length=round(10/h);
%     
%     %--initialize values for timer and stimuli representations
%     P_CX=zeros(1,interval_length);
%     P_CS=zeros(1,interval_length);
%     
%     x_CX=zeros(1,interval_length);
%     x_CS=zeros(1,interval_length);
%     %--
%     
%     %--Obtain A and V values from memory
%     A_CX=ITIUS.Emb.CX.A_S2(trial); % A trial value for CX
%     A_CS=ITIUS.Emb.CS.A_S2(trial); % A trial value for CS
%     %--
%     
%     for t2=1:interval_length
%         
%         counterDDMCX3=counterDDMCX3+1;
%         counterDDMCS2=counterDDMCS2+1; % update counter for random process in DDM
%         
%         % min will take the minimum value: either DDM result or 3. This
%         % caps the value of integrator at 3.
%         P_CX(t2+1)=min(DDM( P_CX(t2), A_CX, h, m, N_CX3(counterDDMCX3) ), 3);
%         P_CS(t2+1)=min(DDM( P_CS(t2), A_CS, h, m, N_CS2(counterDDMCS2) ), 3);
%         
%         % max ensures the minimum value the accumulator can reach is
%         % Aini. This avoids division by zero later.
%         P_CX(t2+1)=max(P_CX(t2+1), Aini);
%         P_CS(t2+1)=max(P_CS(t2+1), Aini);
%         
%         %---Stimuli representations (RBFs)
%         x_CX(t2)=CStrace(P_CX(t2+1),mu,sigma,tau_x,1,x_CX(t2),h);
%         x_CS(t2)=CStrace(P_CS(t2+1),mu,sigma,tau_x,1,x_CS(t2),h);
%         %---
%         
%         %---CR
%         CR(t+t2)=max(0,[x_CS(t2) x_CX(t2)]*[V_CS V_CX]');
%         %---
%         
%     end
%     
%     %---Slope Correction
%     A_CX=A_CX+A_CX*alpha_t*(1-P_CX(t2+1))/P_CX(t2+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
%     A_CS=A_CS+A_CS*alpha_t*(1-P_CS(t2+1))/P_CS(t2+1);
%     %---
%     
%     %---V update
%     V_CX=RW([V_CX V_CS],alpha_E,[x_CX(t2) x_CS(t2)],0,A_CX,P_CX(t2+1));
%     V_CS=RW([V_CS V_CX],alpha_E,[x_CS(t2) x_CX(t2)],0,A_CS,P_CS(t2+1));
%     %---
    
    %---Move every A and V value 1 trial forward
    ITIUS.Emb.CX.A_S2(trial+1)=A_CX; % A trial value for CX
    ITIUS.Emb.CS.A_S2(trial+1)=A_CS; % A trial value for CS
%     NoITIUS.Emb.CX.V(trial+1)=V_CX; % V trial value for CX
    ITIUS.Emb.CS.V(trial+1)=V_CS; % V trial value for CS
    ITIUS.Emb.CR(trial,:)=CR;
    %---
   
end


%--calculate average response rate
MeanCR=mean(ITIUS.Emb.CR((end-40):end,:),1);
%---

%--plots
% Associative strength CX, CS during CS

subplot(2,4,1)
plot(1:length(ITIUS.Emb.CS.V),ITIUS.Emb.CS.V,1:length(ITIUS.Emb.CX.V),ITIUS.Emb.CX.V)
legend('V CS','V CX')

subplot(2,4,2)
plot(1:length(ITIUS.Emb.CX.A_ITI)-10,1./ITIUS.Emb.CX.A_ITI(11:end))
legend('CX time est. to CS: 350')

subplot(2,4,3)
plot(1:length(ITIUS.Emb.CX.A_S1)-10,1./ITIUS.Emb.CX.A_S1(11:end))
legend('CX time est. to reward: 110')

subplot(2,4,4)
plot(1:length(ITIUS.Emb.CS.A_S1)-10,1./ITIUS.Emb.CS.A_S1(11:end))
legend('CS time est. to reward: 110')

subplot(2,4,5)
plot(1:length(ITIUS.Emb.CX.A_S2)-10,1./ITIUS.Emb.CX.A_S2(11:end))
legend('CX time est. to CS off: 10')

subplot(2,4,6)
plot(1:length(ITIUS.Emb.CS.A_S2)-10,1./ITIUS.Emb.CS.A_S2(11:end))
legend('CS time est. to CS off: 10')

subplot(2,4,7)
plot(MeanCR)
xlim([0 130])