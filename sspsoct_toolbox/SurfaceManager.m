classdef SurfaceManager < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        surface = [];
        findex = [];
        
        dim
        
        weight = [];

    end
    
    methods
        function self = SurfaceManager(facedata)
            if nargin < 1
                %facedata = [];
            else
                self.surface = facedata.face;
                self.findex = facedata.index;
            end
        end
        
        function addline(self,line,frameindex)
            self.surface = cat(1,self.surface,gather(line));
            self.findex = cat(1,self.findex,gather(frameindex));
        end
        
        function line = get_line(self,frameindex,smooth)
            
            if nargin > 2
               tsurface = round(imgaussfilt(self.surface,smooth));
            else
               tsurface = self.surface;
            end
            
            if isempty(self.findex)
                line = tsurface(1,:);
            else
                i = find(self.findex == frameindex);
                line = tsurface(i,:);                  
            end
 
            
        end
        
        function set_image_size(self,dim)
            
            self.dim = dim;
            
            W = dim(2);
            H = dim(1);
            
            if isempty(self.surface)
                self.surface = ones(dim(3),dim(2));
            end
            
            trueW = length(self.surface(1,:));

            if W ~= trueW
                self.surface =  (round(interp1(1:trueW,self.surface',linspace(1,trueW,W),'linear','extrap')))';
            end
            
            trueD = length(self.surface(:,1));
            
            D = dim(3);
            if D ~= trueD
                self.surface =  (round(interp1(1:trueD,self.surface,linspace(1,trueD,D),'linear','extrap')));    
                self.findex = (round(interp1(1:trueD,self.findex,linspace(1,trueD,D),'linear','extrap')));
            end        
            
           self.surface(self.surface<1) = 1;
           self.surface(self.surface>H) = H;
           
           self.weight = zeros(length(self.surface(:,1)),1);
        end
        
        
        function roiP = line2roi(self,findex,key_alines)
                        
            Y = self.surface(findex,key_alines);
            roiP = [key_alines' Y'];
 
        end
        
%         function fit2D(self,key_alines,key_frames)
%             
%             [D,W] = size(self.surface);
%             
%             keysurface = self.surface(key_frames,key_alines);
% %             [Dkey,Wkey] = size(keysurface);
% 
%             [Xkey,Ykey] = meshgrid(key_alines,key_frames);  
%             sf = fit([Xkey(:),Ykey(:)],keysurface(:),'poly55');
%             [X,Y] = meshgrid(1:W,1:D);
%             z = sf(X,Y);
%             self.surface = round(z);
% 
%             self.surface(self.surface<1) = 1;
%             self.surface(self.surface>self.dim(1)) = self.dim(1);            
%             
%         end
        

        
        function weight_smooth(self,n,w)
            ss = [];
            flag = [];
            line1 = ones(1,length(self.surface(1,:)));
            line0 = zeros(1,length(self.surface(1,:)));
            
            
            [H,W] = size(self.surface);
            for index = 1:length(self.weight)
                if self.weight(index) == 1
                    aaa = repmat(self.surface(index,:),2*w+1,1);
                    flagaaa = repmat(line1,2*w+1,1);
                    
                    ss = cat(1,ss,aaa);
                    flag = cat(1,flag,flagaaa);
                    flag(end-w,:)=0;
                    
                else
                    ss = cat(1,ss,self.surface(index,:));
                    flag = cat(1,flag,line0);
                end
            end
            
            kenel = ones(n,5)/(n*5);
            ss = imfilter(ss,kenel,'replicate');
            ss(flag>0) = [];
            ss = reshape(ss,H,W);

            self.surface = round(ss);

            self.surface(self.surface<1) = 1;
            self.surface(self.surface>self.dim(1)) = self.dim(1);  

        end
        
        function med_filt(self,m)
            
            ss = [];
            flag = [];
            line1 = ones(1,length(self.surface(1,:)));
            line0 = zeros(1,length(self.surface(1,:)));

            [H,W] = size(self.surface);
            
            www = self.weight;
            www(1) = 1;
            www(end) = 1;
            
            for index = 1:length(www)
                if www(index) == 1
                    aaa = repmat(self.surface(index,:),2*m+1,1);
                    flagaaa = repmat(line1,2*m+1,1);
                    ss = cat(1,ss,aaa);
                    flag = cat(1,flag,flagaaa);
                    flag(end-m,:)=0;
                else
                    ss = cat(1,ss,self.surface(index,:));
                    flag = cat(1,flag,line0);
                end
            end

            ss = round(medfilt2(ss,[m 1],'symmetric'));
            ss(flag>0) = [];
            ss = reshape(ss,H,W);            

            self.surface = round(ss);
            self.surface(self.surface<1) = 1;
            self.surface(self.surface>self.dim(1)) = self.dim(1);  
        end
        
        function gauss_filt(self,sigma)
            self.surface = round(imgaussfilt(self.surface,sigma));
            self.surface(self.surface<1) = 1;
            self.surface(self.surface>self.dim(1)) = self.dim(1);              
        end

        function roi2line(self,roi,fi)
            
            self.weight(fi) = 1;
            
            [D,W] = size(self.surface);
            H = self.dim(1);
            line = (round(interp1(roi(:,1),roi(:,2),1:W,'spline','extrap')));
            self.surface(fi,:) = line';
            
            self.surface(self.surface<1) = 1;
            self.surface(self.surface>H) = H;

        end
        
        function roi2poly(self,roi,fi)
            
            self.weight(fi) = 1;
            
            [D,W] = size(self.surface);
            H = self.dim(1);
            
            pfit = fit(roi(:,1),roi(:,2),'poly5');

            line = round(pfit(1:W));
            self.surface(fi,:) = line';
            
            self.surface(self.surface<1) = 1;
            self.surface(self.surface>H) = H;

        end  
        
        function roi2spline(self,roi,fi)
            
            self.weight(fi) = 1;
            
            [D,W] = size(self.surface);
            H = self.dim(1);
            
            pfit = fit(roi(:,1),roi(:,2),'smoothingspline');

            line = round(pfit(1:W));
            self.surface(fi,:) = line';
            
            self.surface(self.surface<1) = 1;
            self.surface(self.surface>H) = H;

        end          
        
        
        
        
        function s = get_surface(self,fi)
            fi = round(fi);
            s = self.surface;
            s(fi,:) = max(s(:));
            return
        end
    end
end

