function [area,angle,frequency,countteeth,trilength,triratio,innerLength,toothList]=toothFinder(leafFile,showFlag)
%Locates and measures the teeth found at the margin of a leaf.
%Input: leafFile, a leaf object or a file containing a single leaf object. This
% object has fields 'imExt' = RGB image of single leaf; 'x' and 'y' = set of 
% Cartesian coordinates of boundary of leaf.
% showFlag, set to 1 to display various figures or 0 for no display (faster)
%
%Output:
%  area, total area of all teeth found
%  angle, mean angle of teeth tips
%  frequency, number of teeth per unit perimeter
%  countteeth, the total number of teeth found
%  trilength, total length of *outside* of all teeth (i.e. the outer two edges)
%  triratio, the mean ratio between the longer and shorter outer edges
%  innerLength, total length of 'inside' of all teeth (i.e. the single base edge)
%  toothList, a vector containing the area, tip angle and lengths of each tooth found
%
%Requires Matlab Image Processing Toolbox.
%
% Authors:
% Yin Hu, Jing Jin, David Corney
% University of Surrey, 2011-2012
% Contact d.corney@surrey.ac.uk

% v.1.0  March 2012

if ischar(leafFile)
    leafObj=load(leafFile);
else
    leafObj=leafFile;
end

imExt=leafObj.imExt;        %RGB image of leaf
x=leafObj.x;                %coordianates of leaf edge
y=leafObj.y;

if showFlag==1
    clf;
    imshow(imExt);hold on
    axis off;
end

%% Edge extraction
% color image change to grey imge and binary processing
graypic=(imExt(:,:,2)+imExt(:,:,3));              %pick out green & blue channels
level = graythresh(graypic);
BW = im2bw(graypic,level);      %to black & white
BW=~bwareaopen(~BW,6000);       %remove small black objects
BW=bwareaopen(BW,6000);         %remove small white objects
%for i=size(BW,1):-1:1
%    bknum=find(BW(i,:)==0);
%    if max(bknum)-min(bknum)<15
%        BW(i,bknum)=1;
%    end
%end
BW=~BW;

s=size(imExt);
bw=false(s(1),s(2));

%Interpolate large gaps in outline loaded from file
interThresh=150;
flagAdded=1;
while flagAdded==1
    flagAdded=0;
    for i=2:length(x)
        dsq=(x(i)-x(i-1))^2+(y(i)-y(i-1))^2;
        if dsq>interThresh
            interThresh=ceil(interThresh*.95);
            flagAdded=1;
            nx=0.5*(x(i)+x(i-1));       %add point half way between two distant (but adjacent) points
            ny=0.5*(y(i)+y(i-1));
            tx=[x(1:i-1);nx;x(i:end)];
            ty=[y(1:i-1);ny;y(i:end)];
            
            x=tx;
            y=ty;
        end
    end
end

%Find row and column min and max limits, with a boundary size of t
x=uint16(x);
y=uint16(y);
t=5;
for i=1:length(x)
    if x(i)-t<1
        rmin=1;
    else
        rmin=x(i)-t;
    end
    if x(i)+t>s(1)
        rmax=s(1);
    else
        rmax=x(i)+t;
    end
    if y(i)-t<1
        cmin=1;
    else
        cmin=y(i)-t;
    end
    if y(i)+t>s(2)
        cmax=s(2);
    else
        cmax=y(i)+t;
    end
    
    bw(rmin:rmax,cmin:cmax)=1;
end

bw = imfill(bw,'holes');

%Integrate both forms of binary image
f=find(bw==1 & BW==1);
BW(:,:)=1;
BW(f)=0;

BW=~BW;
areaList = regionprops(BW, 'area');
pixels = regionprops(BW, 'PixelIdxList');
maxSize=-1;index=0;

%Find biggest blob in image (just in case there are still more than one)
for i=1:length(areaList)
    if areaList(i).Area>maxSize
        maxSize=areaList(i).Area;
        index=i;
    end
end

BW(:,:)=1;
BW(pixels(index).PixelIdxList)=0;

BW=~BW;
%get central point
centerpos = regionprops(BW, 'centroid');
centerpt=centerpos.Centroid;


%invoke coutour_following from coutour_following.m to get edge image
[C1,re] = contour_following(BW);
if re==0
    BW=~BW;        %image may need inverting
    [C1,re] = contour_following(BW);
end

%Calculate centroid-contour distance vector
dist=sqrt((C1(:,2)-centerpt(2)).*(C1(:,2)-centerpt(2))+(C1(:,1)-centerpt(1)).*(C1(:,1)-centerpt(1)));

%Find maximally distant point from centre and re-order contour so it at start
[maxmid_v,maxmid_d]=max(dist);
disttemp = dist;
dist(maxmid_d-maxmid_d+1:length(dist)-maxmid_d+1)=disttemp(maxmid_d:length(dist));
dist(1+length(dist)-maxmid_d+1:length(dist))=disttemp(1:maxmid_d-1);
C1temp = C1;
C1(maxmid_d-maxmid_d+1:length(dist)-maxmid_d+1,:)=C1temp(maxmid_d:length(dist),:);
C1(1+length(dist)-maxmid_d+1:length(dist),:)=C1temp(1:maxmid_d-1,:);

%Extreme value calculation
%Finds non-trivial variations in centroid-contour distance vector
%Min=sinus, max=tooth tip
scaling=5;
IndMin=find(deviation(sign(deviation(dist,scaling)),scaling)>0)+scaling+1;
IndMax=find(deviation(sign(deviation(dist,scaling)),scaling)<0)+scaling+1;


%f=find(IndMin>length(dist));
%IndMin(f)=[];
%f=find(IndMax>length(dist));
%IndMax(f)=[];


%% Revise local optima
%Find local minimum near current sinus estimate
for i=1:length(IndMin)
    left = IndMin(i)-scaling;
    right= IndMin(i)+scaling;
    if left<=0
        left=1;
    end
    if right>length(dist)
        right=length(dist);
    end
    temp = dist(left:right);
    [~,mind] = min(temp);
    IndMin(i)=left+mind-1;
end

%Find local maximum near current tip estimate
for i=1:length(IndMax)
    left = IndMax(i)-scaling;
    right= IndMax(i)+scaling;
    if left<0
        left=1;
    end
    if right>length(dist)
        right=length(dist);
    end
    temp = dist(left:right);
    [~,maxd] = max(temp);
    IndMax(i)=left+maxd-1;
end


if showFlag==1
    figure(3);clf
    plot(1:length(dist),dist,'k');hold on
    axis square
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ht(1)=ylabel('distance from centroid');
    ht(2)=xlabel('distance around perimeter');
    %set(ht,'fontsize',14)
    %a=axis; a(2)=length(dist); axis(a);    
    plot(IndMax,dist(IndMax),'ko');
    plot(IndMin,dist(IndMin),'k+');
    drawnow    
end

%% Filter max and min value
%Find two local minima without a tooth between them, and delete "less minimal" one
tempmark=[];
for i=1:length(IndMin)-1
    left = IndMin(i);
    right= IndMin(i+1);
    f=find(IndMax>=left & IndMax<=right, 1); %no max value in the middle
    if isempty(f)
        %Delete "false sinus"
        if dist(left)<dist(right)
            tempmark=[tempmark i+1];
        else
            tempmark=[tempmark i];
        end
    end
end
IndMin(tempmark)=[];

%Find two local maxima without a sinus between them, and delete "less maximal" one
tempmark=[];
for i=1:length(IndMax)-1
    left = IndMax(i);
    right= IndMax(i+1);
    f=find(IndMin>=left & IndMin<=right, 1); %no max value in the middle
    if isempty(f)
        %Delete "false tooth tip"
        if dist(left)>dist(right)
            tempmark=[tempmark i+1];
        else
            tempmark=[tempmark i];
        end
    end
end
IndMax(tempmark)=[];

%% Find three teeth with two sinuses, and adjust location of middle tooth
tempIndMax = IndMax;
for i=2:length(tempIndMax)-1
    left = tempIndMax(i-1);
    mid = tempIndMax(i);
    right = tempIndMax(i+1);
    f1=find(IndMin>=left & IndMin<=mid);
    f2=find(IndMin>=mid & IndMin<=right);
    
    count = 1;
    temp=[];
    if ~isempty(f1) && ~isempty(f2)  %finding two min values
        for j=IndMin(f1(1)):IndMin(f2(1))
            x1 = C1(IndMin(f1(1)),1);
            y1 = C1(IndMin(f1(1)),2);
            x2 = C1(IndMin(f2(1)),1);
            y2 = C1(IndMin(f2(1)),2);
            x3 = C1(j,1);
            y3 = C1(j,2);
            temp(count)=pointline(x1,y1,x2,y2,x3,y3);
            count=count+1;
        end
        [~,maxd]=max(temp);
        IndMax(i)=IndMin(f1(1))+maxd-1;
    end
end



%% revise chain code
maxv=max(IndMax(2),IndMin(2));
len = length(C1);
C1=[C1;C1(1:maxv,:)];           %repeat first few contour points at end, to produce overlapping loop
IndMax=[IndMax; IndMax(1:2)+len];
IndMin=[IndMin; IndMin(1:2)+len];

%%  statistical leaf teeth and delete fake teeth

teethheight=0;
teethlength=0;
teethcount=0;
teethDims=[];

%First pass to calculate ranges
for i=1:length(tempIndMax)
    if i==1
        left = tempIndMax(length(tempIndMax));
        mid = tempIndMax(i);
        right = tempIndMax(i+1);
    elseif i==length(tempIndMax)
        left = tempIndMax(i-1);
        mid = tempIndMax(i);
        right = tempIndMax(1);
    else
        left = tempIndMax(i-1);
        mid = tempIndMax(i);
        right = tempIndMax(i+1);
    end
    f1=find(IndMin>=left & IndMin<=mid);
    f2=find(IndMin>=mid & IndMin<=right);
    if ~isempty(f1) && ~isempty(f2)
        x1 = C1(IndMin(f1(1)),1);       %(x1,y1)=first sinus
        y1 = C1(IndMin(f1(1)),2);
        x2 = C1(IndMin(f2(1)),1);       %(x2,y2)=second sinus
        y2 = C1(IndMin(f2(1)),2);
        x3 = C1(mid,1);                 %tooth tip
        y3 = C1(mid,2);
        
        %calculating area
        d=pointline(x1,y1,x2,y2,x3,y3);     %This is distance from (x3,y3) to line (x1,y1):(x2,y2)
        %Calculate area
        %a1=x1*y2+x2*y3+x3*y1;
        %a2=x2*y1+x3*y2+x1*y3;
        %polyarea=round(0.5*(a2-a1));       %assumes points are clockwise ordered
        
        ln=sqrt((x1-x2)^2+(y1-y2)^2);
        
        teethheight= teethheight+d;
        teethlength= teethlength+ln;
        teethcount = teethcount+1;
        teethDims  =[teethDims;d ln];
    end
end

teethheight=teethheight/teethcount;
if size(teethDims,1)==0
    teethHeightSD = 0;
    teethLengthSD = 0;
else
    hRange       = sort(teethDims(:,1));
    teethHeightSD= std(hRange);      %use mean+/- 2*s.d. later on as filter
    teethlength  = teethlength/teethcount;
    lRange       = sort(teethDims(:,2));
    teethLengthSD= std(lRange);
end

tempIndMax = IndMax;
area      =0;
countteeth=0;
leaflength=0;
angle     =0;
trilength =0;
triratio  =0;
toothList =[];

for i=1:length(tempIndMax)
    if i==1
        left = tempIndMax(length(tempIndMax));
        mid  = tempIndMax(i);
        right= tempIndMax(i+1);
    elseif i==length(tempIndMax)
        left = tempIndMax(i-1);
        mid  = tempIndMax(i);
        right= tempIndMax(1);
    else
        left = tempIndMax(i-1);
        mid  = tempIndMax(i);
        right= tempIndMax(i+1);
    end
    f1=find(IndMin>=left & IndMin<=mid);
    f2=find(IndMin>=mid & IndMin<=right);
    if ~isempty(f1) && ~isempty(f2)
        %Model each leaf as a triangle:
        x1 = C1(IndMin(f1(1)),1);        %(x1,y1)=first sinus (A)
        y1 = C1(IndMin(f1(1)),2);
        x2 = C1(IndMin(f2(1)),1);        %(x2,y2)=second sinus (B)
        y2 = C1(IndMin(f2(1)),2);
        x3 = C1(mid,1);                  %tip of tooth (C)
        y3 = C1(mid,2);
        
        v1=[x1-x3 y1-y3];
        v2=[x2-x3 y2-y3];
        v3=[x1-x2 y1-y2 ];
        
        AC=sqrt(sum(v1.^2));
        BC=sqrt(sum(v2.^2));
        AB=sqrt(sum(v3.^2));
        %re-order so that AB is longer than AC:
        if AC>AB
            ACt=AC; AC=AB; AB=ACt;      %for triRatio etc.
            v1t=v1; v1=v2; v2=v1t;      %used for angles
        end
        
        %Calculate tooth dimensions
        d =pointline(x1,y1,x2,y2,x3,y3);
        ln =sqrt((x1-x2)^2+(y1-y2)^2);
        
        %Calculate area, assuming tooth is a triangle
        a1=x1*y2+x2*y3+x3*y1;
        a2=x2*y1+x3*y2+x1*y3;
        polyarea=abs(0.5*(a2-a1));
        
        %size of incision as fraction of distance to blade centroid
        d1=sqrt( (x1-centerpt(1)).^2 + (y1-centerpt(2)).^2);
        d2=sqrt( (x2-centerpt(1)).^2 + (y2-centerpt(2)).^2);
        d3=sqrt( (x3-centerpt(1)).^2 + (y3-centerpt(2)).^2);
        incisionFrac=(d3-min(d1,d2))/min(d1,d2);
        if incisionFrac>.25         %It's a lobe and not a tooth (see "Manual of Leaf Architecture", Ellis et al., 2009)
            lobe=1;
        else
            lobe=0;
        end
        
        skipTooth=0;
        
        %teethLength = mean tooth length; ln=this tooth's length; teethLengthSD = (truncated) standard deviation
        
        %Skip over-sized teeth
        if ((ln-teethlength)>2*teethLengthSD || (d-teethheight)>2*teethHeightSD) %leaf apex is not tooth should be filtered
            skipTooth=1;
        end
        %Skip under-sized teeth
        if polyarea<=15
            skipTooth=1;        %ignore very small teeth
        end
        
        if lobe==1
            skipTooth=1;
        end
        if skipTooth==1
            if showFlag==1
                figure(2);
                line([y1,y2],[x1,x2],'color','m');
                for j=IndMin(f1(1)):IndMin(f2(1))
                    hl=line([C1(j,2),C1(j+1,2)],[C1(j,1),C1(j+1,1)]);
                    set(hl,'color','m');
                    if lobe==1
                        set(hl,'color','y')
                    end
                end
            end
            continue;
        end
        if showFlag==1
            figure(2);
            line([y1,y2],[x1,x2]);
            line([y1,y3],[x1,x3]);
            line([y3,y2],[x3,x2]);
            hold on;plot(y3,x3,'rx');
            drawnow
            %pause(0.01)
        end
        
        %trilength=trilength+AB+AC;
        trilength=trilength + AC+BC;
        if AC>0
            triratio=triratio + (AB/AC) ;
            %else degenerate case of collinearity, so no update to triratio
        end
        
        %ang = acos((v1*v2')/(sqrt(v1*v1')*sqrt(v2*v2')));
        %if (sqrt(v1*v1')*sqrt(v2*v2')) ==0
        %    ang=0;      %degenerate case of collinearity
        %end
        
        %Angle at tip:
        if (AC*BC)>0
            ang = acos( (AC^2 + BC^2 - AB^2) / (2 * AC * BC) );
        else
            ang = 0;       %degenerate case of collinearity
        end
        
        if skipTooth==0
            angle=angle+ang;            
            area=area+polyarea;
            %areaList=[areaList polyarea];   %list of individual tooth sizes
            thisTooth.area=polyarea;
            thisTooth.angle=ang;
            thisTooth.lenAB=AB;
            thisTooth.lenBC=BC;
            thisTooth.lenAC=AC;
            toothList=[toothList thisTooth];
            leaflength=leaflength+ln;       %add ln = length across base of teeth
            countteeth=countteeth+1;        %increment number of teeth counter
            if showFlag==1
                %fprintf('Tooth %03d\t Area %d\tLen,Height: %4.3f,%4.3f\t Cum area %d\n',countteeth,round(0.5*d*l),d,l,round(area));
                fprintf('Tooth %03d\t Area %d \n',countteeth,round(polyarea));
            end
        end
    end
end

%Calculate final scores, inc. calculating means from sums
innerLength =leaflength;
if countteeth>0
    angle=angle/countteeth;
else
    angle=0;
end

if leaflength>0
    frequency=countteeth/leaflength;
else
    frequency=0;
end

if countteeth>0
    triratio=triratio/countteeth;
else
    triratio=0;
end

end

%% subfunctions

function outdata=deviation(indata,scaling)
%Subsample a sequence at a given scale factor
outdata=zeros(size(indata));
for i=1:scaling:size(indata,1)-scaling
    outdata(i)=indata(i+scaling)-indata(i);
end
end

function [ d ] = pointline( x1,y1,x2,y2,x3,y3 )
%POINTLINE Distance from Point to edge
% Returns distance from point (x3,y3) to line through (x1,y1) and (x2,y2)

A=y2-y1;
B=-(x2-x1);
C=-(y2-y1)*x1+y1*(x2-x1);

if (abs(A)+abs(B))>0
    d=abs(A*x3+B*y3+C)/sqrt(A^2+B^2);
else
    d=0;        %degenerate case if points are collinear (-> area=0)
end

end



%% CONTOUR_FOLLOWING
function [C,re] = contour_following(BW)
% CONTOUR_FOLLOWING takes a binary array and returns the sorted row and
% column coordinates of contour pixels.
%
% C = CONTOUR_FOLLOWING(BW) takes BW as an input. BW is a binary array
% containing the image of an object ('1': foreground, '0': background). It
% returns a circular list (N x 2, C(1,:)=C(end,:)) of the
% (row,column)-coordinates of the object's contour, in the order of
% appearence (This function was inspired from the freeman contour coding
% algorithm).
%
% Note:
% - if the object is less than 3 pixels, CONTOUR_FOLLOWING sends back [0 0].
% - the algorithm is quite robust: the object can have holes, and can also
% be only one pixel thick in some parts (in this case, some coordinates
% pair will appear two times: they are counted "way and back").

% Based on code by Francois Mourougaya, via
% http://www.mathworks.com/matlabcentral/fileexchange/14947-contourfollowing

re=1;
[m,n]=size(BW);                                                            % getting the image height and width

Itemp=zeros(m+2,n+2);                                                      % we create a '0' frame around the image to avoid border problems
Itemp(2:(m+1),2:(n+1))=BW;
BW=Itemp;

BW = BW - imerode(BW,[0 1 0 ; 1 1 1 ; 0 1 0]);                             % gets the contour by substracting the erosion to the image
BW = bwmorph(BW,'thin',Inf);                                               % to be sure to have strictly 8-connected contour

if (sum(sum(BW))<3),                                                       % we consider that less than 3 pixels cannot make a contour
    C=[0 0];
    return;
end;

[row,col]=find(BW,1);                                                      % takes the first encountered '1' pixel as the starting point of the contour

MAJ=[6 6 0 0 2 2 4 4];                                                     % variable initialization
C=[0 0 ; 0 0];
k=0;
ended=0;
direction=4;
count1=0;
count2=0;
while(ended==0),
    k=k+1;
    found_next=0;
    count1=0;
    while(found_next==0)
        switch mod(direction,8),
            case 0,
                if (BW(row, col+1)==1),
                    row=row;
                    col=col+1;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 1;
                if (BW(row+1, col+1)==1),
                    row=row+1;
                    col=col+1;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 2;
                if (BW(row+1, col)==1),
                    row=row+1;
                    col=col;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 3;
                if (BW(row+1, col-1)==1),
                    row=row+1;
                    col=col-1;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 4;
                if (BW(row, col-1)==1),
                    row=row;
                    col=col-1;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 5;
                if (BW(row-1, col-1)==1),
                    row=row-1;
                    col=col-1;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 6;
                if (BW(row-1, col)==1),
                    row=row-1;
                    col=col;
                    C(k,:)=[row col];
                    found_next=1;
                end;
            case 7;
                if (BW(row-1, col+1)==1),
                    row=row-1;
                    col=col+1;
                    C(k,:)=[row col];
                    found_next=1;
                end;
                
        end
        
        if (found_next==0), direction=direction+1; end;
        count1=count1+1;
        if count1>100000;
            re=0;
            return;
        end
    end
    
    if(and((length(C)>3),(([C(1,:) C(2,:)]==[C((end-1),:) C(end,:)])))),
        ended=1;
    end;
    
    direction = MAJ((mod(direction,8)+1));
    count2=count2+1;
    if count2>10000;
        re=0;
        return;
    end
end

C=C(1:(end-1),:);                                                          % the first and last points in the list are the same (circular list)
C=C-1;                                                                     % to go back to the original coordinates (without the '0' frame)
end