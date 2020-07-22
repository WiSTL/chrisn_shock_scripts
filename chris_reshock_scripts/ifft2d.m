function [b] = ifft2d(a)
b = ifft(a, [], 2);
b = ifft(b, [], 1, 'symmetric');
end

