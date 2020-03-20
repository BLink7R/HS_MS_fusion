function IMG = NormColor(IMG)
% NORMCOLOR
mini = min(IMG(:));
maxi = max(IMG(:));
IMG = (IMG-mini)/(maxi-mini);
end