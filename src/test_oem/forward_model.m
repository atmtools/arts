function [ R, yi, Ki ] = forward_model( J, Hs, Q, R, x, iter )

    n = size(Hs,1);
    Ki = zeros(size(J));
    for i = 1:n
        Ki(i,:) = Ki(i,:) + 0.5 * (Hs(:,:,i) * x)';
    end
    Ki = Ki + J;
    yi = Ki * x;
    
end

