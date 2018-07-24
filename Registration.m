classdef Registration
    
    properties
        T % time vector in matlab's day format
        Tforms
        ref2d = imref2d([2048 2064]); 
    end
    
    methods
        function Reg = Registration(T,Tforms)
            if nargin==0
                return
            else
                Reg.T=T; 
                Reg.Tforms = Tforms; 
            end
        end
        function regimg = register(Reg,img,t,varargin)
            arg.resize=1; 
            arg.grouping = ones(size(t)); 
            arg = parseVarargin(varargin,arg); 
            sz=[size(img,1) size(img,2)]; 
            %Ryan: noticed an empty img can make size(img,3) == 1 not 0, creating
            %an error if t is 0 (which makes sense for an empty img)
            %however am not sure how we have an empty img and t
%             if isempty(img)
%                 siz = 0;
%             else
%                 siz = size(img,3);
%             end
%             assert(siz==numel(t),'Image Stack and time vector must of of the same size'); 
            assert(size(img,3)==numel(t),'Image Stack and time vector must of of the same size');
            assert(numel(Reg.T)==numel(Reg.Tforms),'Registration object not intialized propertly');
            regimg = zeros(size(img),'single'); 
            for i=1:numel(t)
                % to allows registration that within a specific acqusition
                % we introduce the arg,grouping. This is just a vector of
                % acq identity. The code that restricts the allowed
                % timepoints within the Reg.T to the ones that the time is
                % within the current group. 
                currentGroup = arg.grouping(i);
                allowedTimepoints = find(Reg.T>min(t(arg.grouping==currentGroup)) & Reg.T<=max(t(arg.grouping==currentGroup))); 
                [~,ix] = min((Reg.T(allowedTimepoints)-t(i)).^2);
                ix = allowedTimepoints(ix);
                if isempty(ix)
                    tform = Reg.Tforms(1);
                else
                    tform = Reg.Tforms(ix);
                end
                tform.T=ceil(tform.T*arg.resize); 
                regimg(:,:,i) = imwarp(img(:,:,i),tform,'OutputView',imref2d(sz));
            end
        end
        
        
        function plotAbsShift(Reg,varargin)
           arg.fig = []; 
           arg = parseVarargin(varargin,arg); 
           if isempty(arg.fig)
               arg.fig=figure; 
           end
           
           %%
           d=nan(numel(Reg.Tforms),1);
           xy=nan(numel(d),2); 
           for i=1:numel(Reg.Tforms)
               t=Reg.Tforms(i).T;
               xy(i,:)=t(end,1:2); 
               d(i)=sqrt(sum(t(end,1:2).^2)); 
           end
           
           %%
           figure(arg.fig)
           clf
           subplot(1,2,1)
           hold on
           clr=jet(size(xy,1)); 
           for i=2:size(xy,1)
               plot(xy(i-1:i,1),xy(i-1:i,2),'.-','color',clr(i,:))
           end
           xlabel('X correction (or maybe Y?) [pixels]')
           ylabel('Y correction (or maybe X?) [pixels]')
           subplot(1,2,2)
           t=Reg.T; 
           t=t-min(t); 
           t=t*24*60*60; 
           plot(t,d)
           xlabel('Time [min]')
           ylabel('Abs shift')
        end
    end
    
end