function [Xtr,ytr,wo,tr_acc,Xte,yte,te_acc,niter,tex] = uo_nn_solve_sequence(num_target, tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,icg,irc,nu,iheader)
%function uo_nn_solve_sequence(num_target, tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,icg,irc,nu,iheader)
%uo_nn_solve_sequence([1,1,2,2,2],1,123456,500,47904864,500,1,10^-6,5000,1,2,30,10^-3,0.01,0.45,1,2,2,1,0)
%0--------print the acutal parameters---------------------
%
fprintf('::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('           NN UNCONSTRAINED OPTIMIZATION SOLVE.\n')
fprintf('::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('    target sequence = %i \n',  num_target)

%
%1--------generate the training data set------------------------------
[Xtr,ytr] = uo_nn_sequence_dataset(tr_seed, tr_p, num_target, tr_freq);

%!!Ara es necessari analitzar 5 imatges de sopeton-> 5*35 pixels (175) o
%pas
%2--------Find the value of w minimizing the loss function------------
%definim funcions i varaibles
[n,p] = size(Xtr);
%w = zeros(35,1);
w = zeros(5*n,1);% sempre tindran 35 pixels(35 entrades)
Xseq = zeros(5*n,p);
for j=1:p
    if j==1
        Xseq(1:35,j)=Xtr(:,j);
    elseif j==2
        Xseq(36:70,j)=Xseq(1:35,j-1);
        Xseq(1:35,j)=Xtr(:,j);
    elseif j==3
        Xseq(71:105,j)=Xseq(36:70,j-1);
        Xseq(36:70,j)=Xseq(1:35,j-1);
        Xseq(1:35,j)=Xtr(:,j);
    elseif j==4
        Xseq(106:140,j)=Xseq(71:105,j-1);
        Xseq(71:105,j)=Xseq(36:70,j-1);
        Xseq(36:70,j)=Xseq(1:35,j-1);
        Xseq(1:35,j)=Xtr(:,j);
    else
        Xseq(141:175,j)=Xseq(106:140,j-1);
        Xseq(106:140,j)=Xseq(71:105,j-1);
        Xseq(71:105,j)=Xseq(36:70,j-1);
        Xseq(36:70,j)=Xseq(1:35,j-1);
        Xseq(1:35,j)=Xtr(:,j);
    end
end
yseq = zeros(1,p);
for j=1:p
    if j<4
        yseq(j)=0;
    elseif j>4 && ytr(j)==1 && ytr(j-1)==1 && ytr(j-2)==1 && ytr(j-3)==1 && ytr(j-4)==1
        yseq(j)=1;
    else
        yseq(j)=0;
    end  
end
%fprintf('%i,',Xseq(:,1));
sig = @(Xseq) 1./(1+exp(-Xseq));
y = @(Xseq,w) sig(w'*sig(Xseq));
L = @(w) norm(y(Xseq,w)-yseq)^2 + (la*norm(w)^2)/2;
gL = @(w) 2*sig(Xseq)*((y(Xseq,w)-yseq).*y(Xseq,w).*(1-y(Xseq,w)))'+la*w;
h = @(w) eye(n);
f= L(w);
g= gL(w);
%inicialitzem outputs
k = 0;
wk  = [];
alphak = [];
dk  = [];
gk  = [];
fk  = [];
Bbfgs = {};
tic;
while k < kmax && norm(gL(w)) >= epsG
    k = k + 1;
    wk = [wk,w];
    gk = [gk,g];
    fk = [fk,f];
    d=zeros(5*n,1);
    %% FIND D IN FUNCTION WITH THE DIFERENT METHODS
    if isd==1 
        %GRADIENT METHOD
        d = -g;
    elseif isd ==2
        break
    elseif isd == 3  
        % QNM/BFGS
        if k == 1
            Bbfgs{k} = eye(175);
            d = -g;      
        else
            sk = wk(:,k) - wk(:,k-1);      
            yk = gk(:,k) - gk(:,k-1);
            
            B = ((Bbfgs{k-1}*(sk*sk')*Bbfgs{k-1})/(sk'*Bbfgs{k-1}*sk));
            C = (yk*(yk)')/(yk'*sk);
            bbfgs =  Bbfgs{k-1} - B+ C;
            Bbfgs{k} = bbfgs;
            d = -inv(bbfgs)*gk(:,k);
        end
    end
    if g'*d >=0 %if the direction is not descendent break 
        dk=[dk,d];
        break
    end
    %%------------------------------------------
    %%FIND THE ALPHA
    %la funcio a minimitzar no es quadratica per tant es necessari fer
    %servir el BLS metode per trobar una alpha, no es pot trobar
    %l'exacta
    if k == 1
        almax = ialmax;
    else
        almax = 2*(fk(k)-fk(k-1))/(gk(:,k)'*d);  
        
    end
    %fprintf('-------------------almax = %3.1e\n', almax)
    [alpha,ioutBLS] = uo_BLSNW32(L,gL,w,d,almax,c1,c2,kmaxBLS,epsal);%uo_BLSNW32(f,g0,x0,d,alpham,c1,c2,maxiter,eps)
    % si iout_bls = 0 alpha succed
    %if (ioutBLS==1)
    %           fprintf('BLS reach max iterations(ioutBLS=1), may be dont fullfill the SWC\n')
    %else (ioutBLS==2)
    %            fprintf('BLS reach the optimality tolerance (ioutBLS=2), may be dont fullfill the SWC\n')
    %end
    %%--------------------------------------
    %% COMPUTE THE RESULTS ONCE WE'VE GOT ALPHA AND THE DIRECTION FOR NEXT ITERANTION
    alphak = [alphak,alpha];
    dk  = [dk,d];
    w   = w + alpha*d;
    f   = L(w);
    g   = gL(w);
end
%% OUTPUT THE FINAL RESULTS

niter = k;
wo = w;
wk = [wk,w];
fk = [fk,f];
gk = [gk,g];
%fprintf('k:%5i   gk(wk):%8.6e  %+4.1e \n', k, norm(gk(:,k)), la);    
%fprintf(' Optimal solution:\n');
%for i = 1:5*n
%   fprintf('final w(%i)  = [ %7.3f].\n ', i, w(i));
%end

if norm(g) <= eps
    iout = 0; 
elseif k == kmax 
    iout = 1;
elseif g'*d >= 0
    iout = 2;  
else
    iout = 3;  
end

tex=toc; %time to compute the optimal solution

   
%3--------------Calculate the accuracy training----------------------
tr_acc = uo_nn_accuracy(wo, Xseq, yseq);
fprintf('trainning accuracy: %i/n', tr_acc);


%4-----------------generate the test dataset------------------------
[Xte,yte] = uo_nn_sequence_dataset(te_seed, te_q, num_target, tr_freq);
Xteseq = zeros(5*n,1);
for j=1:p
    if j==1
        Xteseq(1:35,j)=Xte(:,j);
    elseif j==2
        Xteseq(36:70,j)=Xteseq(1:35,j-1);
        Xteseq(1:35,j)=Xte(:,j);
    elseif j==3
        Xteseq(71:105,j)=Xteseq(36:70,j-1);
        Xteseq(36:70,j)=Xteseq(1:35,j-1);
        Xteseq(1:35,j)=Xte(:,j);
    elseif j==4
        Xteseq(106:140,j)=Xteseq(71:105,j-1);
        Xteseq(71:105,j)=Xteseq(36:70,j-1);
        Xteseq(36:70,j)=Xteseq(1:35,j-1);
        Xteseq(1:35,j)=Xte(:,j);
    else
        Xteseq(141:175,j)=Xteseq(106:140,j-1);
        Xteseq(106:140,j)=Xteseq(71:105,j-1);
        Xteseq(71:105,j)=Xteseq(36:70,j-1);
        Xteseq(36:70,j)=Xteseq(1:35,j-1);
        Xteseq(1:35,j)=Xte(:,j);
    end
end
%fprintf('Xteseq(:,10)= %i \n',Xteseq(:,10));
%fprintf('Xte(:,10)= %i \n',Xte(:,10));
yteseq = zeros(1,p);
for j=1:p
    if j<4
        yteseq(j)=0;
    elseif j>4 && yte(j)==1 && yte(j-1)==1 && yte(j-2)==1 && yte(j-3)==1 && yte(j-4)==1
        yteseq(j)=1;
    else
        yteseq(j)=0;
    end  
end

%5-------------calculate the accuracity for this dataset------------
te_acc = uo_nn_accuracy(wo, Xteseq, yteseq);
fprintf('test accuracy: %i/n', te_acc);


%---------------Plot the sequence ---------------------------
%uo_nn_sequence_plot(Xteseq,yteseq,wo);