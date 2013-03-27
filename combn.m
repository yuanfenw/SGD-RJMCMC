function a = combn(s, inv_set)
% product of combinations (s(i) / inv_set(i))
a = 1;
for i=1:length(inv_set)
    a = a*prod((ones(1,inv_set(i))+(s(i)-1)) - [0:(inv_set(i)-1)]);
%     a = a*nchoosek(s(i), inv_set(i))*factorial(inv_set(i));
end
end