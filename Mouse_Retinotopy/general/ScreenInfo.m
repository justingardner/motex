classdef ScreenInfo
    %SCREENINFO the ScreenInfo object
    %   Detailed explanation goes here
    
    properties
        
        PixelSize = 0.0609;     % The size of the pixels in cm
        Xmax = 640;             % Total number of pixels, horizontal
        Ymax = 480;             % Total number of pixels, vertical
        FrameRate  = 124.8918;  % Refresh rate
        WhichScreen             % Screen number
        MonitorType = 'Dummy';  % Make and model of the monitor(s)
        MonitorSize             % Total size of monitor in cm, horizontal
        Dist = 64;              % Distance between observer and screen, cm
        PixelDepth = 8;         % The pixel depth
        windowPtr               % The pointer assigned to the window by Screen
        ScreenRect              % The rectangle assigned by Screen
        SyncSquare              % Object specifying the properties of the Sync Square
        WaveInfo                % Object specifying the properties of DAQ
        Calibration             % A struct with Calibration info (will be an object)      
    end
    
    methods
        
        function SI = ScreenInfo(RigInfo, LBflag)
            % ScreenInfo initializes the ScreenInfo object
            %
            % ScreenInfo(RigInfo)
            %
            % ScreenInfo(RigInfo,LBflag) lets you specify the
            % "Laura Busse" mode (default: 1, meaning that it is on)
            %
            % See also: RigInfoGet
            
            if nargin < 1
                return
            end
            
            if nargin<2
                LBflag = 1;
            end
            
            AssertOpenGL;
            
            SI.WhichScreen      = RigInfo.VsDisplayScreen;
            SI.MonitorType      = RigInfo.MonitorType;
            SI.MonitorSize      = RigInfo.MonitorSize;
            SI.Calibration      = Calibration; %#ok<CPROP>
            SI.Calibration.Directory   = RigInfo.VsHostCalibrationDir;
            SI.Dist             = RigInfo.DefaultMonitorDistance;
            
            % added by MC 2013-04-25:
            rect = RigInfo.VsDisplayRect; % defaults to [] i.e. whole screen
            
            SI.PixelDepth = 8;
            
            % suppress all the greetings and warnings
            oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel');
            oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings');

            Screen('CloseAll');
            
            % suppress all the greetings and warnings
            Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
            Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

            WaitSecs(0.5);
            % HACK TO TAKE CARE OF A FLAW IN 64 bit ver of Psychtoolbox (MC
            % 2013-01-23)
            if  strcmp(mexext,'mexw64')
                pixdepth = []; 
            else 
                pixdepth = SI.PixelDepth;
            end
              % the following line needs be enabled for proper linear superposition
            if LBflag
                [SI.windowPtr, SI.ScreenRect] = Screen('OpenWindow', SI.WhichScreen, [], rect, pixdepth, [], [], [], kPsychNeed16BPCFloat);
            else
                [SI.windowPtr, SI.ScreenRect] = Screen('OpenWindow', SI.WhichScreen, [], rect, pixdepth);
            end
            
            % the following will be overruled in various cases
            % Screen('BlendFunction', SI.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            % Screen('BlendFunction', ScreenInfo.windowPtr, GL_ONE, GL_ZERO);
            
            gray=GrayIndex(SI.WhichScreen);
            Screen('FillRect', SI.windowPtr, gray);
            Screen('Flip', SI.windowPtr);		% force gray screen
            
            
            % make a linear Clut (do this even though you will load the calibration later!!!)
            Screen('LoadNormalizedGammaTable', SI.WhichScreen, repmat((0:255)', 1, 3)/255);
%             % AS added next two lines 2012-08, cause calibration was never loaded before.
%             comment: I don't see this coming in anywhere in vs. 
%             SI = CalibrationLoad(SI);
%             Screen('LoadNormalizedGammaTable', SI.WhichScreen, SI.Calibration.monitorGamInv./255);
            
            SI.Xmax = RectWidth(SI.ScreenRect);
            SI.Ymax = RectHeight(SI.ScreenRect);
            
            SI.FrameRate =  1/Screen('GetFlipInterval',SI.windowPtr);
            % was FrameRate(SI.WhichScreen); but this occasionally got rid of
            % calibration, and needed to be flushed after changing framerate
            
            SI.PixelSize = SI.MonitorSize/SI.Xmax; % size of pixel
            
            SI.SyncSquare = RigInfo.SyncSquare;
            
            %% Deal with the WaveInfo
            
            SI.WaveInfo = RigInfo.WaveInfo;
            
            %% inform the user of how things are going
            
            fprintf('\n *** ScreenInfoInitialize *** \n');
            fprintf('VBLTimestampingMode is %d\n',...
                Screen('Preference', 'VBLTimestampingMode'));
            disp(['You are using a ' SI.MonitorType]);
            disp(['The refresh rate is ' num2str(SI.FrameRate,'%3.3f') ' Hz']);
            fprintf('The resolution is %dx%d pixels.\n',SI.Xmax,SI.Ymax);
            SI.WaveInfo.Describe;
            
        end % function ScreenInfo
        
        function [xPix, yPix] = Deg2PixCoord(SI,xDeg,yDeg)
            % Deg2PixCoord converts degrees coordinates into pixel coordinates
            %
            % [xPix, yPix] = ScreenInfo.Deg2PixCoord(xDeg,yDeg)
            % Rounds the results
            % Assumes you want to center things in the middle
            % Crops so you don't get out of the screen
            % 
            % See also Deg2Pix, Pix2Deg, and Deg2PixCirc
            
            xPixCtr = SI.Xmax/2;
            yPixCtr = SI.Ymax/2;
            
            xPix = 2 * SI.Dist/ SI.PixelSize * tan( pi/180 * xDeg/2 );
            yPix = 2 * SI.Dist/ SI.PixelSize * tan( pi/180 * yDeg/2 );
            
            xPix = round(xPixCtr + xPix);
            yPix = round(yPixCtr + yPix);
            
            xPix = max(1,xPix);
            yPix = max(1,yPix);
            
            xPix = min(xPix, SI.Xmax);
            yPix = min(yPix, SI.Ymax);
            
        end % function Deg2PixCoord
        
        function Pix = Deg2PixCirc(SI,Deg)
            % Deg2PixCirc converts degrees into pixels for a circular screen
            %
            % Pix = ScreenInfo.Deg2PixCirc(Deg) converts the degrees of
            % visual angle Deg into a round number of pixels Pix. This is
            % correct for a circular screen and obeys Deg2PixCirc(a+b) =
            % Deg2Pix(a)+Deg2Pix(b). If you want something appropriate for
            % a flat screen, see , use Deg2Pix.
            
            DegPerPix = 2* 180/pi * asin(SI.PixelSize/2 / SI.Dist);            
            Pix = round(Deg/DegPerPix);
            
        end % function Deg2Pix
        
        function Pix = Deg2Pix(SI,Deg)
            % Deg2Pix converts degrees into pixels for a flat screen
            %
            % Pix = ScreenInfo.Deg2Pix(Deg) converts the degrees of visual
            % angle Deg into a round number of pixels Pix. This is correct
            % for a flat screen and for objects centered in the middle of
            % teh screen. It uses tan (tangent) and has the
            % counterintuitive property that Deg2Pix(a+b) is not
            % Deg2Pix(a)+Deg2Pix(b). If you want that property, use
            % Deg2PixCirc.
            
            Pix = 2 * SI.Dist/ SI.PixelSize * tan( pi/180 * Deg/2 );
            
            Pix = round(Pix);
            
        end % function Deg2Pix
         function Deg = Pix2Deg(SI,Pix)
            % Pix2Deg converts pixels into degrees
            %
            % Deg = ScreenInfo.Pix2Deg(SI,Pix)
            %
            % See also Deg2Pix.
            
            Deg = 2 * 180/pi*atan(Pix / (2 * SI.Dist/ SI.PixelSize));
            
        end % function Pix2Deg
        
        function SI = CalibrationLoad(SI)
            % CalibrationLoad loads the calibration file and applies it
            
            SI.Calibration = Calibration.Load(SI); %#ok<PROP>
  
        end % function CalibrationLoad
        
        function CalibrationCheck(SI)
            % Check calibration automatically to ensure linearity
            
            Calibration.Check(SI); %#ok<PROP>  
            
        end
        
        function result = IsCalibrationStillOn(SI)
            % Checks whether the gamma table that we wanted is still on
            
            PresentGamma = Screen('ReadNormalizedGammaTable', SI.windowPtr);
            DesiredGamma = SI.Calibration.monitorGamInv / 255;	% corrected to have linear luminance
            MaxDifference = max(max(abs(PresentGamma-DesiredGamma)));
            result = MaxDifference<0.01;
            
        end % function IsCalibrationStillOn
    
        function FrameRate = RealFrameRate(SI)
            % RealFrameRate is for backward compatibility
            % MC 2011-06-15
            
            FrameRate = SI.FrameRate;
        end 
        
        function obj= loadCalibrationFromFile(obj,filepath)
            if exist(filepath,'file')
                % load the calibration data
                data = load(filepath);                                
                
                if isfield(data,'C') && isa(data.C,'Calibration')
                    obj.Calibration= data.C;
                elseif isfield(data,'Calibration') % backward compatibility
                    obj.Calibration.date = data.Calibration.date;
                    obj.Calibration.MonitorType = data.Calibration.ScreenInfo.MonitorType;
                    obj.Calibration.FrameRate = data.Calibration.ScreenInfo.FrameRate;
                    obj.Calibration.Directory = data.Calibration.ScreenInfo.CalibrationDir;
                    obj.Calibration.monitorGamInv = data.Calibration.monitorGamInv;                    
                else
                    error('invalid calibration file!');
                end
            
                if strcmpi(obj.MonitorType,obj.Calibration.MonitorType)
                    fprintf('Loading and applying calibration done on %s\n',obj.Calibration.date);
                else
                    warning('Inconsistent calibration file! loading a calibration made for %s on a %s!',...
                        obj.Calibration.MonitorType,obj.MonitorType);
                end
                if any(isnan(obj.Calibration.monitorGamInv))
                    error('unable to use the calibration! there are NaNs in inverse gamma function!')
                end                                        
                          
                Screen('LoadNormalizedGammaTable', obj.windowPtr, obj.Calibration.monitorGamInv / 255);
            else
                error('unable to find %s\n',filepath);
            end
        end
        
    end % Methods
    
end

