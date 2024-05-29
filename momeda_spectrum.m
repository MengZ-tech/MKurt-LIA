function [T MKurt f y T_best MKurt_best f_best y_best t] = momeda_spectrum(x,filterSize,window,range,plotMode)
    % MULTIPOINT OPTIMAL MINUMUM ENTROPY DECONVOLUTION ADJUSTED
    
    % Assign default values for inputs
    if( isempty(filterSize) )
        filterSize = 300;
    end
    if( isempty(plotMode) )
        plotMode = 0;
    end
    if( isempty(window) )
        window = ones(1,1);
    end
    if( isempty(range) )
        range = [5:0.05:300];
    end
    
    if( sum( size(x) > 1 ) > 1 )
        error('MOMEDA:InvalidInput', 'Input signal x must be 1d.')
    elseif(  sum(size(plotMode) > 1) ~= 0 )
        error('MOMEDA:InvalidInput', 'Input argument plotMode must be a scalar.')
    elseif( sum(size(filterSize) > 1) ~= 0 || filterSize <= 0 || mod(filterSize, 1) ~= 0 )
        error('MOMEDA:InvalidInput', 'Input argument filterSize must be a positive integer scalar.')
    elseif( sum(size(window) > 1) > 1 )
        error('MOMEDA:InvalidInput', 'Input argument window must be 1d.')
    elseif( min(range) <= length(window) )
        error('MOMEDA:InvalidInput', 'Range starting point must be larger than the length of the window.')
    elseif( filterSize >= length(x) )
        error('MOMEDA:InvalidInput', 'Input argument filterSize must be smaller than the length of input signal x.')
    end
    
    L = filterSize;
    x = x(:); % A column vector
    
    %%% Calculte X0 matrix
    N = length(x);
    X0 = zeros(L,N);
    
    for( l =1:L )
        if( l == 1 )
            X0(l,1:N) = x(1:N);
        else
            X0(l,2:end) = X0(l-1, 1:end-1);
        end
    end
    
                        % "valid" region only
    X0 = X0(:,L:N-1);   % y = f*x where only valid x is used
                        % y = Xm0'*x to get valid output signal
    
    autocorr = X0*X0';
    autocorr_inv = pinv(autocorr);
    
    % Built the array of targets impulse train vectors separated the by periods
    T = zeros(length(range),1);
    i = 1;
    t = zeros(N-L,length(range));
    for period = range
        points{i} = 1:period:(size(X0,2)-1);
        points{i} = round(points{i});
        t(points{i},i) = 1;
        T(i) = period;
        i = i + 1;
    end
    
    % Apply the windowing function to the target vectors
    t = filter(window, 1, t);
    
    % Calculate the spectrum of optimal filters
    f = autocorr_inv * X0 * t;

    % Calculate the spectrum of outputs
    y = X0'*f;
    
    % Calculate the spectrum of PKurt values for each output
    MKurt = mkurt(y,t);
    
    % Find the best match
    [MKurt_best index_max] = max(MKurt);
    T_best = T(index_max);
    f_best = f(:,index_max);
    y_best = y(:,index_max);
    
    % Plot the resulting spectrum
    if( plotMode > 0 )
        figure;
        plot(T,MKurt);
        ylabel('Multipoint Kurtosis')
        xlabel('Period (samples)');
        axis('tight')
        
        figure;
        subplot(3,1,1)
        plot(x)
        title('Input signal');
        xlabel('Sample number');
        
        subplot(3,1,2)
        plot(y_best)
        title(strcat(['Best output signal (period=', num2str(T_best), ')']));
        xlabel('Sample number');
        
        subplot(3,1,3)
        stem(f_best)
        title(strcat(['Best filter (period=', num2str(T_best), ')']));
        xlabel('Sample number');
    end
end

function [result] = mkurt(x,target)
    % This function simply calculates the summed kurtosis of the input
    % signal, x, according to the target vector positions.
    result = zeros(size(x,2),1);
    for i = 1:size(x,2)
        result(i) = ( (target(:,i).^4)'*(x(:,i).^4) )/(sum(x(:,i).^2)^2) * sum(abs(target(:,i)));
    end
end
