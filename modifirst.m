function modifirst(X.tiff)
close all;  clc;
f = imread(X.tiff);
level = graythresh(uint8(f));
BW2 = im2bw(f,level);            %Binarization OTSU
figure,imshow(BW2),title('After Applying Dialation');
BW = imcomplement(BW2);          %inverting black and white
BW = bwareaopen(BW,40);          % Remove all object containing fewer than 40 pixels
[M,N]=size(BW)
%B=[0 1 0;1 1 1;0 1 0];
%BW1=imdilate(BW, B);
figure,imshow(BW1),title('After Applying Dialation');
%se = strel('disk',2);
%afterOpening = imopen(BW,se);
%figure, imshow(afterOpening,[]); title('After Applying Opening');
%figure,imshow(BW),title('After Applying Opening');

%.....Connected Component Labelling & region properties calculation........
[L,n] = bwlabeln(BW,4);
propied=regionprops(L,'BoundingBox', 'Area','PixelList','FilledImage');

%................Speckles and Ruling Line & Small CC removal...............
for i=1:length(propied)
    height(i)=propied(i).BoundingBox(:,4);         %CCs BB Height
    width(i)=propied(i).BoundingBox(:,3);          %CCs BB width
    % Area=propied(i).Area;                         %CC Area
    PL=propied(i).PixelList;                       %CC Pixel Lists
    [lenx,lengy]=size(PL);
    if (height(i)<=11 && height(i)/width(i) <= 0.2)||(width(i)<=5 && height(i)/width(i)<=20)
        for j=1:lenx
            xc=PL(j,1);
            yc=PL(j,2);
            BW(yc,xc)=0;         
        end
        
    end
end
%..........................................................................
%figure;imshow(BW);title('After removal');

%........ Mean Height and Width of all CC..................................
sortheight=sort(height, 'descend');
sortwidth=sort(width, 'descend');
fivepercent=ceil((n*5)/100);
sumheight1=0;
sumwidth1=0;
sumheight2=0;
sumwidth2=0;
heightcount=0;
widthcount=0;
for lp=1:fivepercent
    sumheight1=sumheight1+sortheight(lp);
end
for lp=1:fivepercent
    sumwidth1=sumwidth1+sortwidth(lp);
end

meantopheight=sumheight1/fivepercent;
meantopwidth=sumwidth1/fivepercent;
tenpercentheight=meantopheight*0.1;
tenpercentwidth=meantopwidth*0.1;
for i=1:n
    if sortheight(i)>=tenpercentheight
        sumheight2=sumheight2+sortheight(i);
        heightcount=heightcount+1;
    end
    if sortwidth(i)>=tenpercentwidth
        sumwidth2=sumwidth2+sortwidth(i);
        widthcount=widthcount+1;
    end
    
end
meanheightCC=sumheight2/heightcount;   %mean character height
meanWidthCC=sumwidth2/heightcount;     %mean character width

% half_height=max(height)-max(height)*0.05;
% half_width=max(width)-max(width)*0.05;
% kh=0;
% kw=0;
% hsum=0;
% wsum=0;
% for i=1:n
%     if (height(i)>=half_height)
%         hsum=hsum+height(i);
%         kh=kh+1;
%     end
%     if (width(i)>=half_width)
%         wsum=wsum+width(i);
%         kw=kw+1;
%     end
% end
% %kh;
% Mean_MaxheightCC=hsum/kh;
% Mean_MaxwidthCC=wsum/kw;
% CCTH=Mean_MaxheightCC-Mean_MaxheightCC*0.1;
% CCTW=Mean_MaxwidthCC-Mean_MaxwidthCC*0.1;
% kh=0;
% hsum=0;
% kw=0;
% wsum=0;
% for i=1:n
%     if (height(i)>=CCTH)
%         hsum=hsum+height(i);
%         kh=kh+1;
%     end
%     if (width(i)>=CCTW)
%         wsum=wsum+width(i);
%         kw=kw+1;
%     end
% end
% meanHeightCC=hsum/kh;          %Mean height using large component only
% meanwidthCC=wsum/kw;          % Mean width using long component only
h=mean(height);                %Mean height using all Connected component
D=ceil(meanheightCC/2);
% n
% sumheight1
% sumheight2
% meantopheight
% heightcount
% meanheightCC
% %meanHeightCC
% h
%
% sumwidth1
% sumwidth2
% meantopwidth
% heightcount
% meanWidthCC
%.............................Vertivcal Striping..........................%
No_Vertical_Strip=3;
chunk_width = floor(N/No_Vertical_Strip);   % vertical partiotioning
%HP=zeros(M);
% for k=1:M
%    for i=1:chunk_width
%       HP(k)=HP(k)+BW(k,i);
%  end
%end
chunk_no = 1;
chunk_pos(chunk_no) = 0;   %y-cordinate of first chunk
k = 1;
flag = 0;
for i = chunk_width:chunk_width:N
    chunk_no = chunk_no+1;
    chunk_pos(chunk_no) = i;       %y-cordinates of chunks
    if flag == 0
        chunk(k) = chunk_width;    %creating the vector for mat2cell
        k = k+1;
        flag =1;
    else
        chunk(k) = chunk_pos(chunk_no)-chunk_pos(chunk_no-1);
        k = k+1;
    end
end
% if chunk_width is not an integral multiple of y
if chunk_pos(chunk_no) < N
    chunk_pos(chunk_no+1) = N;
    chunk(k) = N-chunk_pos(chunk_no);
end
%................... Local minima Calculation..............................
y_pos = 1;
j = 1;
C = mat2cell(BW,M,chunk);
[~,n] = size(C) % n is the number of vertical chunks
%n=No_Vertical_Strip;
for l =1:n
    A = cell2mat(C(1,l));
    [u,v] = size(A);
    HP = sum(A,2);      %Hp is a column matrix
    
   
    %%..........................Mean Character Height .....................
    meanheightCC
    %figure;bar(HP);
    peak=max(HP)
    HPC=peak-HP;
    %figure;bar(HPC);
    [pks,locs]=findpeaks(HPC,'minpeakdistance',floor(1.5*meanheightCC));
%     for pk=1:size(locs)
%         locpk(l,pk)=locs(pk);
%     end
    locs
    %     [locpkx locpky]=size(locpk);
    %     for i=1:locpky
    %         count=0;
    %         sum1=0;
    %         for j=1:locpkx
    %             if locpk(j,i)~=0
    %                 sum1=sum1+locpk(j,i);
    %                 count=count+1;
    %             end
    %         end
    %         avgloc(i)=sum1/count;
    %     end
    %     for i=1:locpky-1
    %         diffloc(i)=avgloc(i+1)-avgloc(i);
    %     end
    %     meanloc=mean(diffloc);
    %......................................................................
    
    zerocount = 0;           %count No of consecutive rows with HP=0
    first = 0;
    i = 1;
    while i <= u    %can't use a for loop because it prevents change of index variable inside the loop
        while (HP(i) == 0) && (i <u)
            if first == 0
                first = i;
                pos_array(y_pos,j)=i;    %position of line boundary in each chunk (Column vector)
            end
            zerocount = zerocount+1;
            i = i+1;
        end
        if zerocount >= 1
            pos_array(y_pos,j)=first+floor(zerocount/2);
            j=j+1;
            first = 0;
            zerocount = 0;
        end
        i = i+1;
    end
    y_pos = y_pos+1;
    j = 1;
end
[p,q]=size(pos_array)     % size of array storing local minima
pos_array1=pos_array;

%..............Sorting Array..................................
for i=1:p
    for j=1:q-1
        if pos_array(i,j)>pos_array(i,j+1)
            temp=pos_array(i,j);
            pos_array(i,j)=pos_array(i,j+1);
            pos_array(i,j+1)=temp;
            
        end
        if pos_array(i,j)==0
            pos_array(i,j)=1;
        end
    end
end

%.............Line Segmentation Algo. Starts from Here .%..................
D=10;
rect= D*1.9;           %rectangle size almost double that of MeanMaxHeight
n=ceil(rect/2);
imgabove=0;            %Image Portion above line
imgbelow=0;            %Image Portion below line
z=1;A=zeros(q,chunk_width);
%imtool(BW);
figure,imshow(BW);title('After Local Search');
for localx=1:p
    for localy=1:q               %repeat until last local minima calculated
        y=pos_array(localx,localy);
        for x=z:(z+chunk_width)        %repeat till x become N alone x-axis
            for i=-n:n
                r=y-D+i;
                if r>=1 && r<=M
                    for j=-n:n
                        c=x+j;
                        if c>=1 && c<=N
                            imgabove=imgabove+BW(r,c);
                        end
                    end
                end
            end
            for i=-n:n
                r=y+D+i;
                if (r>=1 && r<=M)
                    for j=-n:n
                        c=x+j;
                        if(c>=1 && c<=N)
                            imgbelow=imgbelow+BW(r,c);
                        end
                    end
                end
            end
            %imgabove;
            %imgbelow;
            yk(localy,x)=y;
            if(imgabove < imgbelow)
                if y==0
                    yk(localy,x+1)=yk(localy,x)+1;
                    yk(localy,x)=1;
                elseif y==1
                    yk(localy,x+1)=yk(localy,x);
                else
                    yk(localy,x+1)=yk(localy,x)-1;
                end
                
            elseif(imgabove > imgbelow)
                if y>M-1
                    yk(localy,x+1)=yk(localy,x);
                else
                    yk(localy,x+1)=yk(localy,x)+1;
                end
            else
                yk(localy,x+1)=yk(localy,x);
            end
            imgabove=0;
            imgbelow=0;
            line([x,x+1],[yk(localy,x),yk(localy,x+1)],'LineWidth',2,'Color','red');
            y=yk(localy,x+1);
        end
    end
    z=floor(localx*chunk_width);
end
LINE=yk;
LINE1=LINE;
figure;imshow(BW);title('FROM LINE ARRAY');

%.......................Refinement.........................................
D=6;
rect= D*1.9;           %rectangle size almost double that of MeanMaxHeight
n=ceil(rect/2);
imgabove=0;            %Image Portion above line
imgbelow=0;            %Image Portion below line
for i=1:q
    for j=1:N-1
        if abs(LINE(i,j+1)-LINE(i,j))>=floor(meanheightCC*0.80)
            % if (LINE(i,j)~=LINE(i,j+1))&&(LINE(i,j)+1~=LINE(i,j+1))&&(LINE(i,j)-1~=LINE(i,j+1))
            y=LINE(i,j);
            x=j;
            for im=-n:n
                r=y-D+im;
                if r>=1 && r<=M
                    for jm=-n:n
                        c=x+jm;
                        if c>=1 && c<=N
                            imgabove=imgabove+BW(r,c);
                        end
                    end
                end
            end
            for im=-n:n
                r=y+D+im;
                if (r>=1 && r<=M)
                    for jm=-n:n
                        c=x+jm;
                        if(c>=1 && c<=N)
                            imgbelow=imgbelow+BW(r,c);
                        end
                    end
                end
            end
            if(imgabove < imgbelow)
                if y==0
                    LINE(i,j)=1;
                    LINE(i,j+1)=LINE(i,j);
                elseif y==1
                    LINE(i,j+1)=LINE(i,j);
                else
                    LINE(i,j+1)=LINE(i,j)-1;
                end
                
            elseif(imgabove > imgbelow)
                if y>M-1
                    LINE(i,j+1)=LINE(i,j);
                else
                    LINE(i,j+1)=LINE(i,j)+1;
                end
            else
                LINE(i,j+1)=LINE(i,j);
            end
            imgabove=0;
            imgbelow=0;
            
            %             pix=LINE(i,j);
            %             if pix>1 && pix<M-2
            %                 if BW(LINE(i,j)-1,j+1)==1 && BW(LINE(i,j)+1,j+1)==1 && BW(LINE(i,j),j+1)==1
            %                     LINE(i,j+1)=LINE(i,j);
            %                 elseif BW(LINE(i,j)-1,j+1)==1 && BW(LINE(i,j)+1,j+1)~=1 && BW(LINE(i,j),j+1)~=1
            %                     LINE(i,j+1)=LINE(i,j);
            %                 elseif BW(LINE(i,j)-1,j+1)~=1 && BW(LINE(i,j)+1,j+1)==1 && BW(LINE(i,j),j+1)~=1
            %                     LINE(i,j+1)=LINE(i,j);
            %                     %                 else
            %                     %                     LINE(i,j+1)=LINE(i,j);
            %                     % end
            %                     %elseif (pix>1 && pix<M-2)
            %                 elseif BW(LINE(i,j)+1,j+1)==1 && BW(LINE(i,j),j+1)==1 && BW(LINE(i,j)-1,j+1)~=1
            %                     LINE(i,j+1)=LINE(i,j)-1;
            %                 elseif  BW(LINE(i,j)+1,j+1)~=1 && BW(LINE(i,j),j+1)==1 && BW(LINE(i,j)-1,j+1)==1
            %                      LINE(i,j+1)=LINE(i,j)+1;
            %                 elseif BW(LINE(i,j)+1,j+1)==1 && BW(LINE(i,j),j+1)~=1
            %                     LINE(i,j+1)=LINE(i,j);
            %                 else
            %                     LINE(i,j+1)=LINE(i,j);
            %                 end
            %             else
            %                LINE(i,j+1)=LINE(i,j);
            %            end
        end
    end
end
LINE;
for i=1:q
    for j=1:N-1
        line([j,j+1],[LINE(i,j),LINE(i,j+1)],'LineWidth',2,'Color','red');
    end
end
%..........................................................................

%..............................Display Connected Component.................
% for n=1:size(propied,1)
%     rectangle('Position',propied(n).BoundingBox,'EdgeColor','g','LineWidth',2);
% end
...........................................................................
    
%.......................Display Line segment...............................
rowmin=[];%rowmax=[];
for r=1:q
    rmin(r)=LINE(r,1);
    rmax(r)=LINE(r,1);
    for s=1:N
        if(LINE(r,s)~=1)
            if rmin(r)>LINE(r,s)
                rmin(r)=LINE(r,s);
            end
            if rmax(r)<LINE(r,s)
                rmax(r)=LINE(r,s);
            end
        end
    end
    %rowmin(i)=rmin
    %rowmax(i)=rmax
end
rowmin=rmin;
rowmax=rmax;
%rowmin=min(LINE1,[],2)
%rowmax=max(LINE1,[],2)
text_line = cell(q-1,1);
for i=1:q-1
    % LINE(i)
    UPMIN=rowmin(i);
    LWMAX=rowmax(i+1);
    Line_Height=LWMAX-UPMIN;
    lineseg=zeros(Line_Height,N);
    for j=1:N
        yabove=LINE(i,j);
        ybelow=LINE(i+1,j);
        for k=1:Line_Height
            if (k-1)+UPMIN >=yabove && (k-1)+UPMIN<=ybelow
                lineseg(k,j)=BW(((k-1)+UPMIN),j);
            end
        end
    end
    text_line{i}=lineseg;
    VP=sum(lineseg,1);
    %figure;bar(VP);
    % figure;imshow(lineseg);title('Line Segment');
end
%..........................end display.....................................

%..........................Word segmentation...............................
x_pos=1;
for li=1:q-1
    VP=sum(text_line{li},1);
    zerocount = 0;         %count No of consecutive rows with VP=0
    first = 0;
    j=1;
    i = 1;
    while i < N    %can't use a for loop because it prevents change of index variable inside the loop
        while (VP(i) == 0) && (i <N)
            if first == 0
                first = i;
                %virtical_pos_array(x_pos,j)=i;    %position of line boundary in each chunk (Column vector)
            end
            zerocount = zerocount+1;
            i = i+1;
        end
        if zerocount >=20
            virtical_pos_array(x_pos,j)=first+floor(zerocount/2);
            j=j+1;
            first = 0;
            zerocount = 0;
        %end
        else
            zerocount = 0;
            first=0;
        end
        i = i+1;
    end
    x_pos = x_pos+1;
    j = 1;
    [~,vy]=size(virtical_pos_array);
    figure;imshow(text_line{li});title('Line Segment');
    [u,v]=size(text_line{li});
    for k=1:vy
        for j=1:u-1
           % virtical_pos_array(li,k)
           % j
            line([virtical_pos_array(li,k),virtical_pos_array(li,k)],[j,j+1],'LineWidth',2,'Color','yellow');
        end
    end
end
virtical_pos_array

%................................End Word segmentation......................................................

%...............................Store each Line into a CellArray............................................
% text_line = cell(q);
% for i=1:q-1
%     %for j=1:N
%             y=LINE(i,j);x=LINE(i+1,j);
%             pix=BW(y:x,:);
%             text_line{i}=pix;
%      %end
% end
% for i=1:q
%     %for j=1:N
%         fig/home/cseure;imshow(text_line{i});title('Line Segment');
%     %end
% end
%...........................................................................................................
end




