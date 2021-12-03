 function Answers = u1903643_lab3()
%% ES3C5 lab submission template
% Please DO NOT change this header
% Please DO NOT change any code in this top-level function
% Please DO NOT change function or subfunction arguments (Input OR Output)
% DO Change function name above to u<ID>_lab3(), where <ID> is your student #

% ES3C5 2021/2022 Lab 3
% Module Leader: Adam Noel

% Modify the SUBFUNCTIONS below with the code needed to determine or
% demonstrate the answers requested.

% See the Briefing Sheet for full instructions.

% Initialise answer structure. It will be the ONLY output argument and
% store the answers to all questions
Answers = [];

%% Template call
Q0();

%% Remaining Calls
Answers.Q1 = Q1Fun();
Answers.Q2 = Q2Fun();

end

%% Template Question hypotenuse length
% This is as a sample only to demonstrate what is expected
function c0 = Q0()
% Please DO NOT change function arguments (input OR output)
% Assign answer to c0 (double value)

% Define triangle lengths
a0 = 2; % 1st side
b0 = 1; % 2nd side

% Find length of hypothenuse
c0 = sqrt(a0^2 + b0^2); % Pythagorean theorem to find 3rd side

end
%%
%%%%%%%%%%%%%%%%%%% Start Modifying Below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%% Q1 Industrial Heating
function Q1 = Q1Fun()
% Please DO NOT change function arguments (input OR output)

    % DO NOT REMOVE THIS IF STATEMENT
    % WRITE YOUR CODE INSIDE THIS IF STATEMENT
    if exist('u1903643_lab3_signals.mat', 'file') == 2 % Update with your student ID
        load('u1903643_lab3_signals.mat', 'y1') % Update with your student ID
        
        Q1.Obs = []; % Q1a
        Q1.param = []; % Q1b
        Q1.yHat = []; % Q1c
        Q1.mse = []; % Q1d
        Q1.yFFT = []; % Q1f
        Q1.yHatFFT = []; % Q1f
        Q1.fRange = []; % Q1f
        
        % Q1a
        Ts1 = 0.025; %define interval 
        t = 0:Ts1:7.475; %vector of time corresponding to each sample 
        Q1.Obs = [(1-exp(-0.18*t))' ,(exp(-0.50*t).*sin(15.08*t))', (sin(41.05*t))']; %generation of observation matrix
        
        % Q1b
        Q1.param = Q1.Obs\y1; %column vector of estimates for unknown constants 
        A1 = Q1.param(1);
        B1 = Q1.param(2);
        C1 = Q1.param(3);
        %above positions extracted from Q1.param and assigned to
        %correspoding variables
        
        % Q1c
        Q1.yHat = (A1*(1-exp(-0.18*t)) + B1*exp(-0.50*t).*sin(15.08*t) + C1*sin(41.05*t))'; %predicted temperatures using model equation
        
        % Q1d
        Q1.mse = (1/length(y1))*sum((Q1.yHat - y1).^2); %calculation of mean square
       
        % Q1e
        figure; %new blank figure
        grid on; %applying grid in figure
        plot(t, y1, 'b'); %plot of y1 (actual temperatures) on y axis against t on x axis 
        hold on; %next plot will be in same figure as previous
        plot(t, Q1.yHat, 'r'); %plot of Q1.yHat (predicted temperatures) on y axis against t on x axis 
        legend('Actual', 'Predicted'); %labels of legend in figure
        xlabel('Time (s)'); %x-axis label
        ylabel('Temperature (°C)'); %y axis label
        title('Q1e ~ Comparison of Predictions with Actual Sensor Values ~ u1903643'); %figure title
        
        % Q1f
        % i
        Q1.yFFT = fft(y1); %fast fourier transform of y1 assigned to Q1.yFFT
        % ii
        Q1.yHatFFT = fft(Q1.yHat); %fast fourier transform of yHat assigned to Q1.yHatFFT
        % iii
        Q1.fRange = linspace(0, (1/Ts1)-1, 300); %vector of 300 equally spaced frequencies from 0 to (1/Ts1)-1 of  300 equally spaced points 
        % iv
        figure; 
        grid on;
        plot(Q1.fRange, abs(Q1.yFFT)); %fft plot of y1 on y axis against frequencies Q1.fRange on x axis 
        hold on;
        plot( Q1.fRange, abs(Q1.yHatFFT)); %fft plot of Q1.yHat on y axis against frequencies Q1.fRange on x axis 
        legend('Actual', 'Predicted');
        title('Q1f.iv ~ Comparison of FFT Outputs of Predictions with Actual Sensor Values ~ u1903643');
 
    else
        % LEAVE THESE DEFAULT ARGUMENTS EMPTY
        Q1.Obs = []; % Q1a
        Q1.param = []; % Q1b
        Q1.yHat = []; % Q1c
        Q1.mse = []; % Q1d
        Q1.yFFT = []; % Q1f
        Q1.yHatFFT = []; % Q1f
        Q1.fRange = []; % Q1f
    end

end

%% Q2 Communications Signal
function Q2 = Q2Fun()
% Please DO NOT change function arguments (input OR output)

    % DO NOT REMOVE THIS IF STATEMENT
    % WRITE YOUR CODE INSIDE THIS IF STATEMENT
    if exist('u1903643_lab3_signals.mat', 'file') == 2 % Update with your student ID
        load('u1903643_lab3_signals.mat', 'y2', 'x2') % Update with your student ID

        Q2.Obs = []; % Q2a
        Q2.W = []; % Q2c
        Q2.param = []; % Q2d
        Q2.yHat = []; % Q2e
        Q2.w = []; % Q2g
        Q2.var1 = []; % Q2g
        Q2.var2 = []; % Q2g
        
        % 2a
        A2_obs = x2; %assigning vector x2 to A2_obs
        B2_obs = [0; x2(1:end-1)]; %assigning 0 then vector x2 from first to the element before last to B2_obs
        C2_obs = [zeros(5, 1); x2(1:end-5)]; %assigning vector to C2_obs, first 5 elements being 0s and all elements of x2 excluding final 5 following
        D2_obs = [zeros(8, 1); x2(1:end-8)]; %assigning vector to D2_obs, first 8 elements being 0s and all elements of x2 excluding final 8 following
        Q2.Obs = [A2_obs B2_obs C2_obs D2_obs]; %vectors in the above lines vertically concatenated into Q2.Obs
        
        k2 = 3.9; %variance gain factor
        s2 = 54; %number of samples affected by k2
        unchanched_s = ones(176-s2, 1); %column vector of ones of 122 rows (176-s2) 
        increased_s = ones(s2, 1)/k2; %column vecto of ones of 54 rows then each element divided by k2
        Q2.W = [unchanched_s; increased_s]; %above column vectors vertically concatinated and assigned to Q2.W
        
        Q2.param = lscov(Q2.Obs,y2,Q2.W); %column vector of estimates for constants in eqution y2[n] in alphabetical order
        A2 = Q2.param(1);
        B2 = Q2.param(2);
        C2 = Q2.param(3);
        D2 = Q2.param(4);
        
        Q2.yHat = A2*Q2.Obs(:, 1) + B2*Q2.Obs(:, 2) + C2*Q2.Obs(:, 3) + D2*Q2.Obs(:, 4); %predicted values generated by multiplying constants by column in Q2.Obs 
        
        figure;
        grid on;
        n = 0:1:length(x2)-1; %sequence indeces  
        plot(n, y2, 'r'); %plot of y2 (sensor voltages) on y axis against n on x axis 
        hold on;
        plot(n, Q2.yHat, 'b'); %plot of Q2.yHat (predicted voltages) on y axis against n on x axis 
        legend('Actual', 'Predictions');
        xlabel('Sequence Index');
        ylabel('Voltage(V)');
        title('Q2f - Comparison of Predictions with Actual Voltages - u1903643');
        
        % 2g 
        Q2.w = y2 - Q2.yHat; %difference between predicted and sensor voltages 
        %vectors of mean as first element and standard deviation as second element using maximum likelihood estimates 
        mle_initial = mle(Q2.w(1:176-s2)); %mle_initial unaffected part of Q2.w by Q2.W
        mle_final = mle(Q2.w(end-s2+1:end)); %mle_initial unaffected part of Q2.w by Q2.W 
        %variance = standard deviation^2 calculated from avove vectors and
        %assigned to Q2.var1 and Q2.var2 respectively
        Q2.var1 = mle_initial(2)^2;  
        Q2.var2 = mle_final(2)^2;
        
    else
        % LEAVE THESE DEFAULT ARGUMENTS EMPTY
        Q2.Obs = []; % Q2a
        Q2.W = []; % Q2c
        Q2.param = []; % Q2d
        Q2.yHat = []; % Q2e
        Q2.w = []; % Q2g
        Q2.var1 = []; % Q2g
        Q2.var2 = []; % Q2g
    end

end