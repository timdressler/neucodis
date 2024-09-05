function [varargout] = pp_fastfft_TD(data,srate,channel,varargin)
%pp_fastfft_TD Performs fast fft on epoched data (total EROs)
%   my_fastfft2(data, srate, channel, 'wins', 2, 'overlap', 0.5, 'fig', 0)
try
    f = waitbar(0,'Please wait...');
    pause(.5)
    %set defaults
    wins = 2;
    overlap = 1-0.5;
    CHANI = channel;
    fig = 0;
    ylim_min = 0;
    ylim_max = 50;

    if nargin > 3 && all(isnumeric([srate, channel])) && isnumeric(data)
        for i = 1:2:length(varargin)
            if strcmp(varargin{i}, 'wins')
                wins = varargin{i+1};
            elseif strcmp(varargin{i}, 'overlap')
                overlap = 1-varargin{i+1};
            elseif strcmp(varargin{i}, 'fig')
                fig = varargin{i+1};
            elseif strcmp(varargin{i}, 'timevec')
                timevec = varargin{i+1};
            elseif strcmp(varargin{i}, 'ylim_min')
                ylim_min = varargin{i+1};
            elseif strcmp(varargin{i}, 'ylim_max')
                ylim_max = varargin{i+1};
            else
                warning('Unknown Input, Ignoring')
            end
        end
    else
        error('Not all mandatory Inputs specified')
    end
    waitbar(.33,f,'Starting fast fft');
    pause(1)

    %start fast fft
    win = ceil(wins*srate);
    fft_mat = [];
    curr_sam = 1;
    for epoch = 1:size(data,3)
        while size(data,2) >= curr_sam + win -1
            spex = abs(fft(data(CHANI,curr_sam:curr_sam+win-1,epoch)));
            spex = spex / length(spex);
            spex =  spex(1:length(spex)/2);
            spex(2:end) = spex(2:end) * 2;

            fft_mat(:,end+1,epoch) = spex';

            curr_sam = ceil(curr_sam +win*overlap);
        end
    end

    fft_mat = mean(fft_mat,3);

    spex3 =[];
    for s = 1:size(data,3)
        spex2 = abs(fft(data(channel,:,s)));
        spex2 = spex2 / length(spex2);
        spex2 =  spex2(1:length(spex2)/2);
        spex2(2:end) = spex2(2:end) * 2;
        spex3(:,:,s) = spex2;
    end
    freqvec2 = 0:1/(size(data,2)/srate):(srate/2)-1/(size(data,2)/srate);

    waitbar(.67,f,'Creating plot');
    pause(1)

    freqvec = 0:1/wins:srate/2-1/wins;
    if ~exist('timevec')
        timevec = linspace(0,length(data)/srate,size(fft_mat,2));
    end

    if fig == 1 || nargout == 0
        figure;
        subplot(1,5,2:5)
        imagesc(timevec,freqvec,fft_mat)
        xline(0,'--', 'LineWidth',2, 'Color','r')
        clim([0 max(fft_mat,[], 'all')])
        colormap jet
        colorbar
        axis xy
        ylim([ylim_min ylim_max])
        %%xlim([-500 700])
        xlabel('time [s]')

        subplot(1,5,1)
        plot(mean(spex3,3),freqvec2)
        ylim([ylim_min ylim_max])
        ylabel('frequency [Hz]')
        xlabel('Amplitude [uV]')
        set(gca, 'XDir', 'reverse')


        sgtitle(['Frequency Analysis for Channel ' num2str(channel)])
        set(gcf, "Position",[488   338   769   420])


    end

    switch nargout
        case 0
        case 1
            varargout{1} = fft_mat;
        case 2
            varargout{1} = fft_mat;
            varargout{2} = timevec;
        case 3
            varargout{1} = fft_mat;
            varargout{2} = timevec;
            varargout{3} = freqvec;
        otherwise
            error('Invalid number of outputs requested')
    end

    waitbar(1,f,'Finishing');
    pause(1)
    close(f)

catch
    varargout{1} = nan;
    varargout{2} = nan;
    varargout{3} = nan;
    warning('Unkown Error. Check Syntax and Data!')
    waitbar(1,f,'Exiting');
    pause(1)
    close(f)
end
end


