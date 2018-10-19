function H = lpfilter(type, M, N, DO)
% lpfilter
[U, V] = dftuv(M, N);

D = sqrt(U.^2 + V.^2);
switch type
    case 'gaussian'
        H = exp(-(D.^2)./(2*(DO^2)));
    case 'ideal'
        H = double(D<=DO);
    otherwise
        error('unknown filter type');
end
