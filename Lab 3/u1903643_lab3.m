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
        Ts1 = 0.025;
        t = 0:Ts1:7.475;
        Q1.Obs = [(1-exp(-0.18*t))' ,(exp(-0.50*t).*sin(15.08*t))', (sin(41.05*t))'];
        
        % Q1b
        Q1.param = Q1.Obs\y1;
        A1 = Q1.param(1);
        B1 = Q1.param(2);
        C1 = Q1.param(3);
        
        % Q1c
        Q1.yHat = (A1*(1-exp(-0.18*t)) + B1*exp(-0.50*t).*sin(15.08*t) + C1*sin(41.05*t))';
        
        % Q1d
        Q1.mse = immse(Q1.yHat, y1);
        
        % Q1e
        figure;
        grid on;
        plot(t, y1, 'b');
        hold on;
        plot(t, Q1.yHat, 'r');
        legend('Actual', 'Predicted');
        xlabel('Time (s)');
        ylabel('Temperature (°C)');
        title('Q1e ~ Comparison of Predictions with Actual Sensor Values ~ u1903643');
        
        % Q1f
        % i
        Q1.yFFT = fft(y1);
        % ii
        Q1.yHatFFT = fft(Q1.yHat);
        % iii
        Q1.fRange = linspace(0, (1/Ts1)-1, 300);
        % iv
        figure;
        grid on;
        plot(Q1.fRange, abs(Q1.yFFT));
        hold on;
        plot( Q1.fRange, abs(Q1.yHatFFT));
        legend('Actual', 'Predicted');
        title('Q1f.iv ~ Comparison FFT Outputs of Predictions with Actual Sensor Values ~ u1903643');
 
        
        
        
        
        
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
        A2 = x2(1:end);
        B2 = [0; x2(1:end-1)];
        C2 = [zeros(5, 1); x2(1:end-5)];
        D2 = [zeros(8, 1); x2(1:end-8)];
        Q2.Obs = [A2 B2 C2 D2];
        
        k2 = 3.9;
        s2 = 54;
        unchanched_s = ones(176-s2, 1);
        increased_s = ones(54, 1)/k2; 
        Q2.W = [unchanched_s; increased_s];
        
    
        
        
        
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