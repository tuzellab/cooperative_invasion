% Code to collate data and transform into radial and angular directions
% Written by Tüzel Lab
% Temple University, 2022

clear all
close all

% matrix index - 0: DMSO; 1:GM6001; 2:agarose 

[~, ~, raw] = xlsread('allfiles9_test.xlsx');
for idx = 1:numel(raw)  
    if isnumeric(raw{idx})
        raw{idx} = num2str(raw{idx});
    end
end

filemax=length(raw);

for k=1:filemax
    filename=raw{k,1};
    celltype=raw{k,2};
    location=raw{k,3};
    matrix=raw{k,4};
  
    [foldername,filenameshort,~] = fileparts(filename);
    idx=strfind(foldername,"Collagen");
    if(isempty(idx))
        idx=strfind(foldername,"Agarose");
    end
    rootfolder=extractBefore(foldername,idx);

    alldata=readtable(filename,'Sheet','Sheet1');

    alldata2=table2array(alldata);
    N=height(alldata);
    trackcount=alldata2(N,1);
    MM=max(alldata2(:,4));
    newdata=zeros(MM,2*trackcount+1);

    jold=1;
    r1=1;
    newdata(:,1)=1:MM;
    for i=1:N
        j=alldata2(i,1);
        if(j>jold) 
            r1=1;
            jold=j;
        end
        x=alldata2(i,5);
        y=alldata2(i,6);
        newdata(r1,2*j)=x;
        newdata(r1,2*j+1)=y;
        r1=r1+1;
    end

    maxtime=0;
    for i=2:2*trackcount+1
        if(nnz(newdata(:,i))>maxtime)
            maxtime=nnz(newdata(:,i));
        end
    end

    newdata2=zeros(maxtime,2*trackcount+1);
    newdata2=newdata(1:maxtime,1:2*trackcount+1);

    dummyidx = strfind(foldername,"Spheroid");
    spheroid=extractBetween(foldername,dummyidx+8,dummyidx+8);
    dummyidx = strfind(foldername,"Repeat");
    if(dummyidx)
        repeat=extractBetween(foldername,dummyidx+6,dummyidx+6);
    else
        repeat{1}='1';
    end
    
    if(matrix=='0') 
        cellfolder=sprintf('%sDMSO/cell%s',rootfolder,celltype);
    else if(matrix=='1') 
        cellfolder=sprintf('%sGM6001/cell%s',rootfolder,celltype);
    else if(matrix=='2') 
        cellfolder=sprintf('%sAgarose/cell%s',rootfolder,celltype);
        end
        end
    end

    if(location=='0') 
        locationcellfolder=sprintf('%s/internal',cellfolder);
    else
        locationcellfolder=sprintf('%s/interface',cellfolder);
    end

    
    interfacefolder=sprintf('%s/data_rotated_r%s_s%s_c%s_l%s_m%s',locationcellfolder,repeat{1},spheroid{1},celltype,location,matrix);

    if ~exist(interfacefolder, 'dir')
        mkdir(interfacefolder);
    else
        rmdir(interfacefolder, 's');
        mkdir(interfacefolder);
    end

    for i=1:trackcount
        filename1=sprintf('%s/file_r%03d.dat',interfacefolder,i);
        t=newdata2(1:length(nonzeros(newdata2(:,2*i))),1);
        x=nonzeros(newdata2(:,2*i));
        y=nonzeros(newdata2(:,2*i+1));
        x0=x(1);
        y0=y(1);
        x0=0; y0=0;
        xx=x-x0;
        yy=y-y0;
        plot(xx,yy);
        hold on;
        clear B;
        rr=sqrt(xx.*xx+yy.*yy);
        theta=atan2(yy,xx);
        B(:,1)=t;
        B(:,2)=rr;
        B(:,3)=theta;
      %  dlmwrite(filename1,B,'delimiter','\t','precision',3);

    end
  %  pause;
    hold off;  
    disp(interfacefolder)
    
end

   
   
