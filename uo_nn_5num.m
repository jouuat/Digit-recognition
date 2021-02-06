function uo_nn_5_digits()

clear; clc
rep = 10;
total_acc = 0;
seed = 47904864;
freq = 1;
ncol = 1;
for k = 1:rep
    x = randi([1 9],[1 5]);
    te_acc = [];
    for i = 1:5
        target = x(i);%comprovarem cada digit ara
        % Agafarem els pesos que li pertany al target
        load("w5.mat");
        if target == 0
            w_i = w(:, 55) 
        else
            w_i = w(:, 1+(target-1)*6);
        end
        [Xte,yte] = uo_nn_dataset(seed, ncol, target, freq);
        te_acc_i = uo_nn_accuracy(w_i, Xte, yte);
        te_acc = [te_acc, te_acc_i];
    end
    % If all the accuracies are 100% (i.e., sum(te_acc) == 500)
    if sum(te_acc) == 500
        fprintf('set %i correctly predicted the digits\n', k)
        total_acc = total_acc + 1;
    else 
        fprintf('set %i UNcorrectly predicted\n', k) 
    end
end
total_percentage_accuracy = (total_acc/rep)*100;
fprintf('-----------------------------\n')
fprintf('a total accuracy of  %i \n', total_percentage_accuracy)
end
