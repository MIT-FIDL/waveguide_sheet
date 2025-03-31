clear all;clc

%% Initialize Workspace
Rf = 300e+3;
sensor_area = (2e-3)^2;     % m^2
gain = 2;
N_sensor_per_axis = 7;
ads1115_to_dB = @(signal,gain,R) 10*log10( 1e+6.*(4.096.*signal./(48848.*gain))./R );  % dB of lux, proportional to dB power

% initialize Ploting var
window_width_time = 20;     % sec
window_height = [-6 6];
plot_freq = 3; % plot every # data points
line_colors = [linspace(0.15,0.75,N_sensor_per_axis);linspace(0.7,0.8,N_sensor_per_axis);linspace(0.85,0.25,N_sensor_per_axis)];
line_colors = [line_colors [linspace(0.5,0.8,N_sensor_per_axis);linspace(0.25,0.15,N_sensor_per_axis);linspace(0.75,0.2,N_sensor_per_axis)]];

% data smoothing parameters
sampleSize = 4;  % Define the size of the moving window
dataBuffer = zeros(1, sampleSize);  % Pre-allocate buffer
smoothedDataX = [];  % Initialize an empty array to store smoothed data
smoothedDataY = [];  % Initialize an empty array to store smoothed data

% initialize serial port connection and figure
baud_rate = 9600;
total_sample_time = 3600;    % sec
dataLogger=serialport('COM5',baud_rate);
current_t = 0;
current_v = zeros(N_sensor_per_axis*2,1);
dBx = zeros(N_sensor_per_axis,0);
dBy = zeros(N_sensor_per_axis,0);
t = [];
it_count = 0;
figure;
ylim(window_height)
xlim([max(0,current_t-window_width_time), current_t+window_width_time*0.05])
grid on

% % filter design
% filter_order = 50;
% lp = designfilt('lowpassfir','FilterOrder',filter_order,'HalfPowerFrequency',3.5,'SampleRate',60);

%% Read data and plot in real time
while current_t < total_sample_time
    text = convertStringsToChars(readline(dataLogger));
    if size(text,2) > 0
        it_count = it_count + 1;
        ind = strfind(text,',');    % find index of comma
        current_t = str2double(text(1:ind(1)-1))*1e-6;
        for k = 1:N_sensor_per_axis*2
            current_v(k) = min(str2double(text(ind(k)+1:ind(k+1)-1)),48848);
        end
        dBx = [dBx ads1115_to_dB(current_v(1:2:end),gain,Rf)];
        dBy = [dBy ads1115_to_dB(current_v(2:2:end),gain,Rf)];
        t = [t current_t];
        
        dataBufferX = dBx(:,max([2,end-sampleSize]):end);
        dataBufferY = dBy(:,max([2,end-sampleSize]):end);
        smoothedDataX(:,end+1) = mean(dataBufferX,2);
        smoothedDataY(:,end+1) = mean(dataBufferY,2);        
        if it_count >= plot_freq %&& size(t,2)>filter_order*3
            hold off;
            for k = 1:N_sensor_per_axis
%                 pts = smoothdata(dBx(k,max(1,end-150):end),'movmean',1,'SamplePoints',t(max(1,end-150):end));
%                 plot(t(max(1,end-150):end),pts,'-','LineWidth',1.2,'Color',line_colors(:,k))
                plot(t,smoothedDataX(k,:),'-','LineWidth',1.2,'Color',line_colors(:,k))
%                 plot(t(max(1,end-150):end),dBx(k,max(1,end-150):end),'-','LineWidth',1.2,'Color',line_colors(:,k))
                hold on
            end
            for k = 1:N_sensor_per_axis
%                 plot(t(max(1,end-150):end),smoothdata(dBy(k,max(1,end-150):end),'movmean',0.6,'SamplePoints',t(max(1,end-150):end)),'-','LineWidth',1.2,'Color',line_colors(:,k+N_sensor_per_axis))
                plot(t,smoothedDataY(k,:),'-','LineWidth',1.2,'Color',line_colors(:,k+N_sensor_per_axis))
%                 plot(t(max(1,end-150):end),dBy(k,max(1,end-150):end),'-','LineWidth',1.2,'Color',line_colors(:,k+N_sensor_per_axis))
            end
%             window_height(1) = min(window_height(1),dBx(end)-0.25);
            ylim(window_height)
            xlim([max(0,current_t-window_width_time), current_t+window_width_time*0.05])
            xlabel('Time (sec)', 'FontSize', 14)
            ylabel('Illuminance (dB)', 'FontSize', 14)
            legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7',...
                   'y_1','y_2','y_3','y_4','y_5','y_6','y_7',...
                   'Location','northwest','NumColumns',2)
            pbaspect([1.6,1,1])
            grid on; drawnow
            it_count = 0;
        end
    end
end

%% Plot entire signal history
figure;
plot(t,dBx)
ylim(window_height)
xlim([min(t), max(t)])
xlabel('Time (sec)', 'FontSize', 16)
ylabel('Illuminance (dB)', 'FontSize', 16)
grid on