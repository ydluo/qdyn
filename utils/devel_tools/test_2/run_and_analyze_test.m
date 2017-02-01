% This is a template for QDYN tests

% 1. Run a qdyn example (a short one)
% ...

% 2. Compare the new output to a pre-computed reference output
%    For instance, compute the rms difference between two time series
% u = ... 
% load uref
% err = max(abs( u-uref )) ./ max(abs( uref ));
err = 1e-5;
test = err<0.01;

% 3. Write the result of the test
disp(['L1 misfit = ' num2str(err)])
disp(['Test = ' num2str(test) ])

% 4. Make pdf plots
%plot(t,u, t,uref, t,(uy-uyref)*10)
%print('test_1','-dpdf'))

