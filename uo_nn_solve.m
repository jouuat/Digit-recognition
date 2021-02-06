function [Xtr,ytr,wo,tr_acc,Xte,yte,te_acc,niter,tex] = uo_nn_solve(num_target, tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,icg,irc,nu,iheader)
%uo_nn_solve(8,0.5,123456,250,47904864,250,0,10^-6,5000,1,2,30,10^-3,0.01,0.45,3,2,2,1,0)
%0--------print the acutal parameters---------------------
%
fprintf('::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf('           NN UNCONSTRAINED OPTIMIZATION SOLVE.\n')
fprintf('::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
fprintf(' OPTIMIZATION PARAMETERS:\n')
fprintf('    lamba(add convexity)  : la = %3.1e\n', la)
fprintf('    Optimality tolerance  : epsG = %3.1e\n', epsG)
fprintf('    Max. iterations       : kmax = %i\n',  kmax)
fprintf('    Optimization method:\n')
if     (isd==1)
    fprintf('           Gradient Method (isd=1)\n')
elseif (isd==3)
    fprintf('           QNM/BFGS (isd=3)\n')
end
fprintf('    algo                       : ils = %5.3f\n', ils)
fprintf('    Max. steplength coef.      : ialmax = %5.3f\n', ialmax)
fprintf('    Max. iterations(BLS)       : kmaxBLS = %i\n', kmaxBLS)
fprintf('    SW1 condition, c1          : c1 = %3.1e\n', c1)
fprintf('    SW2 condition, c2          : c2 = %3.1e\n', c2)
fprintf('    Optimality tolerance(BLS)  : epsal = %3.1e\n', epsal)
fprintf('    search direction           : icg = %i\n', icg)
fprintf('    search direction           : irc = %i\n', irc)
fprintf('    search direction(coef. mult alpha when actualizes) : nu = %i\n',  nu)
fprintf('    Initial solution:\n');
%
%1--------generate the training data set------------------------------
[Xtr,ytr] = uo_nn_dataset(tr_seed, tr_p, num_target, tr_freq);


%2--------Find the value of w minimizing the loss function------------
%definim funcions i varaibles
[n,p] = size(Xtr);
%w = zeros(35,1);
w = zeros(n,1);% sempre tindran 35 pixels(35 entrades)
sig = @(Xtr) 1./(1+exp(-Xtr));
y = @(Xtr,w) sig(w'*sig(Xtr));
L = @(w) norm(y(Xtr,w)-ytr)^2 + (la*norm(w)^2)/2;
gL = @(w) 2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))'+la*w;
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
    d=zeros(n,1);
    %% FIND D IN FUNCTION WITH THE DIFERENT METHODS
    if isd==1 
        %GRADIENT METHOD
        d = -g;
    elseif isd ==2
        break
    elseif isd == 3  
        % QNM/BFGS
        if k == 1
            Bbfgs{k} = eye(n);
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
fprintf('k:%5i   gk(wk):%8.6e  %+4.1e \n', k, norm(gk(:,k)), la);    
fprintf(' Optimal solution:\n');
for i = 1:n
   fprintf('final w(%i)  = [ %7.3f].\n ', i, w(i));
end

if norm(g) <= eps
    iout = 0; 
elseif k == kmax 
    iout = 1;
    fprintf('reached the maximum iterations (%3d)\n', kmax)
elseif g'*d >= 0
    iout = 2;
    fprintf('!CONVERGENCE LOST! dk is not a descent direction: gk·dk = %+4.1e.\n', g'*d )   
else
    iout = 3;
    fprintf(' dk descent direction but !fk+1 > fk.!\n')   
end

tex=toc; %time to compute the optimal solution

   
%3--------------Calculate the accuracy training----------------------
tr_acc = uo_nn_accuracy(wo, Xtr, ytr);


%4-----------------generate the test dataset------------------------
[Xte,yte] = uo_nn_dataset(te_seed, te_q, num_target, tr_freq);


%5-------------calculate the accuracity for this dataset------------
te_acc = uo_nn_accuracy(wo, Xte, yte);


%---------------Plot the result (just in individuals runs)----------
%uo_nn_Xyplot(Xte,yte,wo);

