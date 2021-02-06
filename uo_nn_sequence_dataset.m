function [X,y] = uo_nn_sequence_dataset(seed, ncol, target, freq)
%uo_nn_sequence_dataset(123456, 2,[5,6,7], 0.5)
%fprintf('::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
%fprintf('                    DATA SET \n')
%fprintf('::::::::::::::::::::::::::::::::::::::::::::::::::::::\n')
%fprintf('    seed : %3.1e\n', seed)
%fprintf('    ncol  : %3.1e\n', ncol)
%fprintf('    target  : %3.1e\n', target)
%fprintf('    freq  :  %i\n',  freq)

%%
N = [
    0 0 1 0 0, 0 1 1 0 0, 0 0 1 0 0, 0 0 1 0 0, 0 0 1 0 0, 0 0 1 0 0, 0 1 1 1 0; % 1
    0 1 1 1 0, 1 0 0 0 1, 0 0 0 0 1, 0 0 0 1 0, 0 0 1 0 0, 0 1 0 0 0, 1 1 1 1 1; % 2
    0 1 1 1 0, 1 0 0 0 1, 0 0 0 0 1, 0 0 1 1 0, 0 0 0 0 1, 1 0 0 0 1, 0 1 1 1 0; % 3
    0 0 1 1 0, 0 1 0 1 0, 1 0 0 1 0, 1 0 0 1 0, 1 1 1 1 1, 0 0 0 1 0, 0 0 0 1 0; % 4
    1 1 1 1 1, 1 0 0 0 0, 1 1 1 1 0, 0 0 0 0 1, 0 0 0 0 1, 1 0 0 0 1, 0 1 1 1 0; % 5
    0 0 1 1 1, 0 1 0 0 0, 1 0 0 0 0, 1 1 1 1 0, 1 0 0 0 1, 1 0 0 0 1, 0 1 1 1 0; % 6
    1 1 1 1 1, 0 0 0 0 1, 0 0 0 0 1, 0 0 0 1 0, 0 0 1 0 0, 0 1 0 0 0, 1 0 0 0 0; % 7
    0 1 1 1 0, 1 0 0 0 1, 1 0 0 0 1, 0 1 1 1 0, 1 0 0 0 1, 1 0 0 0 1, 0 1 1 1 0; % 8
    0 1 1 1 0, 1 0 0 0 1, 1 0 0 0 1, 0 1 1 1 1, 0 0 0 0 1, 0 0 0 0 1, 0 1 1 1 0; % 9
    0 1 1 1 0, 1 0 0 0 1, 1 0 0 0 1, 1 0 0 0 1, 1 0 0 0 1, 1 0 0 0 1, 0 1 1 1 0; % 0
    ];
N = N';
%%
xval     = 10;
maxsigma = 1.; %dir com de borros les vols (0,0) sera maxim nitit
minsigma = 0.25;
if xval > 0
    N = xval*((N-1)+N);
end
if ~isempty(seed) rng(seed); end;
nT = size(target,2);%[n,p] = size(target);(numero de targets)
if nT < 5
    fprintf('Error missing targets (only %i targets)', nT);
    return
end
if nT > 5
    fprintf('Error too many targets (there are %i targets)', nT);
    return
end
%fprintf('target(1): %i /n',target(1));
nPixels = size(N,1);
k=0;
%% Generation
for j=1:ncol
    if rand() < (freq-0.1*nT)/(1-0.1*nT)%si freq es 0.5 "la mitat" de les vegades passara per aqui i com mes nT menys probabilitats tindra de passar per aqui
        i = target(randi([1,nT]));%si passa per aqui agafara un numero del target segur
    else
        i = randi(10);%si passa per aqui el numero escollit sera 100% random
    end
    X(:,j) = N(:,i);
    if isempty(target)%si no hi ha cap valor al target
        y(j) = i;
    else   
        if (k==0 && i == target(1))
            y(j) = 0;
            k=1;
        elseif (k==1 && i == target(1) && i~= target(2))
            y(j) = 0;
            k=1;
        elseif (k==1 && i == target(2))
            y(j) = 0;
            k=2;
        elseif (k==2 && i == target(3))
            y(j) = 0;
            k=3;
        elseif (k==2 && i ~= target(3) && i==target(2) && target(2)==target(1))
            y(j) = 0;
            k=2;
        elseif (k==3 && i == target(4))
            y(j) = 0;
            k=4;
        elseif (k==3 && i ~= target(4) && i==target(3) && target(3)==target(2) && target(2)==target(1))
            y(j) = 0;
            k=3;
        elseif (k==4 && i == target(5))
            y(j) = 1;
            y(j-1) = 1;
            y(j-2) = 1;
            y(j-3) = 1;
            y(j-4) = 1;
            k=0;  
        else
            y(j) = 0;
            k=0;
        end
    end
end
%% Blur
for j=1:ncol
    l = zeros(1,nPixels);
    sigma = minsigma + (maxsigma-minsigma)*rand();
    for k = 1:nPixels
        ii = randi([1 nPixels]);
        while l(ii)==1
            ii = randi([1 nPixels]);
        end
        l(ii) = 1;
        if X(ii,j) > 0
            X(ii,j) =  xval + sigma*xval*randn();
        else
            X(ii,j) = -(xval + sigma*xval*randn());
        end
    end
end
%disp('The matrix X is:') %Xtr->(35, 250)
%disp(X)
%disp('The matrix y is:')
%disp(y)
end
