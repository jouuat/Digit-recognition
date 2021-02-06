function accuracy = uo_nn_accuracy(w, xte, yte)
accuracy = 0;
n = length(yte);
sig = @(xte) 1./(1+exp(-xte));
y = @(xte,w) sig(w'*sig(xte));
for i=1:n
    if yte(i)==round(y(xte(:,i),w))
        accuracy = accuracy+1;
    end
end
accuracy = accuracy*(100/n);
fprintf('accuracy: %i\n', accuracy)
end