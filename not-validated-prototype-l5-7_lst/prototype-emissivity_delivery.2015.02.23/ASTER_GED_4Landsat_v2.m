function [eL7,nL7] = ASTER_GED_4Landsat_v2(dataloc,resfolder,resplot, vers, Lat,Lon)


% Find all ASTER GED 1x1 degree tiles necessary to cover Landsat scene
maxlat = ceil(max(Lat(:)));
minlat = floor(min(Lat(:)));

lats = maxlat:-1:minlat+1;

maxlon = ceil(max(Lon(:)));
minlon = floor(min(Lon(:)));

lons = minlon:maxlon-1;

latlim = minlat:maxlat;
lonlim = minlon:maxlon;

latlen = length(latlim);
lonlen = length(lonlim);
latins = zeros(1,(latlen-1)*(lonlen-1));
lonins = zeros(1,(latlen-1)*(lonlen-1));
c=1;
for i = 1:latlen-1
    for j = 1:lonlen-1
        latins(c) = latlim(i)+1;
        lonins(c) = lonlim(j);
        c=c+1;
    end
end

e6t = cell(1,length(lats));
nvt = cell(1,length(lats));

for r = 1:length(lats)      % Loop over all ASTER GED tiles
    
    latin = lats(r);
    
    e6 = cell(1,length(lons));
    nv = cell(1,length(lons));
   
    for c = 1:length(lons);
        
        lonin = lons(c);
        
        
        if abs(lonin)<100 && abs(lonin)>=10
            
            if lonin<0
                fileopen = [dataloc,'\AST_L3_',resfolder,'\asterged.v',vers,'.',num2str(latin),'.-0',num2str(abs(lonin)),'.',resplot,'.h5'];
            else
                fileopen = [dataloc,'\AST_L3_',resfolder,'\asterged.v',vers,'.',num2str(latin),'.0',num2str(lonin),'.',resplot,'.h5'];
            end
            
        elseif abs(lonin)<10
            
            if lonin<0
                fileopen = [dataloc,'\AST_L3_',resfolder,'\asterged.v',vers,'.',num2str(latin),'.-00',num2str(abs(lonin)),'.',resplot,'.h5'];
            else
                fileopen = [dataloc,'\AST_L3_',resfolder,'\asterged.v',vers,'.',num2str(latin),'.00',num2str(lonin),'.',resplot,'.h5'];
            end
            
        else
            fileopen = [dataloc,'\AST_L3_',resfolder,'\asterged.v',vers,'.',num2str(latin),'.',num2str(lonin),'.',resplot,'.h5'];
        end
        
        % ASTER Emissivities
        try
            hinfo = hdf5info(fileopen);
            
            % Read in Latitude
            %Lat{c}  = double(hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1)))';
            %Lon{c}  = double(hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(2)))';
            
            % Read in emissivity info
            Emis_mean = double(hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1)));
            
            
            % band 13
            e13_mean = Emis_mean(:,:,4)./1000;
            e13_mean = e13_mean';
            
            % Set ASTER fill values (-9999 and 0) to nan
            ast_fill = (e13_mean<=0);
            e13_mean(ast_fill) = nan;
            
            % band 14
            e14_mean = Emis_mean(:,:,5)./1000;
            e14_mean = e14_mean';
            
            % Set ASTER fill values (-9999 and 0) to nan
            ast_fill = (e14_mean<=0);
            e14_mean(ast_fill) = nan;
            
            % Estimate Landsat band6 emissivity using regression analysis
            e6{c} = 0.44.*e13_mean + 0.4.*e14_mean + 0.156;
            
            NDVI = double(hdf5read(hinfo.GroupHierarchy.Groups(5).Datasets(1)));
            NDVI = NDVI'./100;
            nv{c} = NDVI;
            
            
        catch ME  % Most likely ocean
            
            
            e13_mean = repmat(0.99,[1000 1000]);
            e14_mean = repmat(0.99,[1000 1000]);
            
            e6{c} = 0.44.*e13_mean + 0.4.*e14_mean + 0.156;
            
            NDVI = zeros(1000,1000);
            nv{c} = NDVI;
            
        end
        
    end
    
   e6t{r} = cat(2,e6{1:length(e6)});
    nvt{r} = cat(2,nv{1:length(nv)});
    
end

eL7_v3 = cat(1,e6t{1:length(e6t)});
nN_v3 = cat(1,nvt{1:length(nvt)});



%% Gridding onto Landsat scene

% Define ASTER GED lat/lons
latind = maxlat:-0.01:minlat;
lonind = minlon:0.01:maxlon;


latind = latind(1:length(latind)-1);
lonind = lonind(1:length(lonind)-1);

[C,rs1] = min(abs(latind - max(Lat(:))));
[C,rs2] = min(abs(latind - min(Lat(:))));
[C,cs1] = min(abs(lonind - min(Lon(:))));
[C,cs2] = min(abs(lonind - max(Lon(:))));

e_sub = eL7_v3(rs1:rs2,cs1:cs2);
n_sub = nN_v3(rs1:rs2,cs1:cs2);

% Rescale by factor of 2
e_sub = imresize(e_sub,1/2,'bicubic');
n_sub = imresize(n_sub,1/2,'bicubic');

e_sub = imresize(e_sub,1/2,'bicubic');
n_sub = imresize(n_sub,1/2,'bicubic');

sze = size(e_sub);

% ASTER lat/lon coords to cover Landsat scene
[Late_sub,Lone_sub] = meshgrat([latind(rs2) latind(rs1)],[lonind(cs1) lonind(cs2)],[sze(1),sze(2)]);
Late_sub = flipud(Late_sub);

LatL = imresize(Lat,1/2,'bicubic'); LonL = imresize(Lon,1/2,'bicubic');

szL = size(Late_sub);
szn = szL(1)*szL(2);
x = reshape(Late_sub,1,szn);
y = reshape(Lone_sub,1,szn);
z1 = reshape(e_sub,1,szn);
z2 = reshape(n_sub,1,szn);
 
% Interpolate ASTER GED grid onto Landsat pixels
F1 = TriScatteredInterp(x',y',z1');
F2 = TriScatteredInterp(x',y',z2');

sfull = size(Lat);

eres = F1(LatL,LonL);
eL7 = imresize(eres,sfull,'bicubic');

nres = F2(LatL,LonL);
nL7 = imresize(nres,sfull,'bicubic');




%% Snow and veg lab spectra for Landsat

% sensor = 'LANDSAT7';
%
% % Spectral response functions
% filename = 'F:\Documents and Settings\ghulley\My Documents\MATLAB\Lab_spectra\Landsat\misc_tir_wavnum.txt';
% fid = fopen(filename);
% blah = fgets(fid);
% blah = fgets(fid);
% blah = fgets(fid);
%
% a = fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g %g', [13 inf]); a=a';
% waveb = a(:,1);
%
% if strcmp(sensor,'LANDSAT5')==1
%     rsrb5 = a(:,7);
%     rsrb = rsrb5(270:501);   % 10-13 micron
% elseif strcmp(sensor,'LANDSAT7')==1
%     rsrb7 = a(:,8);
%     rsrb = rsrb7(270:501);   % 10-13 micron
% end
%
% waveb = 10000./waveb(270:501);   % 10-13 micron
%
% fclose(fid);
%
% snowmode = ['fine  ';'medium';'coarse';'ice   '];
% modes = cellstr(snowmode);
%
% count=1;
% for i = 2:2 %length(modes)
%
%     labin = ['F:\Documents and Settings\ghulley\My Documents\MATLAB\Lab_spectra\ASTER\Speclib\Speclib_total\JHU\becknic\water\txt\',modes{i},'.txt'];
%     fid = fopen(labin,'r');
%
%     % Skip header
%     for j = 1:26
%         headg = fgets(fid);
%     end
%
%     a = fscanf(fid, '%g %g', [2 inf]); a=a';
%
%     ref = a(:,2)./100;
%     wave = a(:,1);
%     fclose(fid);
%
%     rneg = ref<0;
%     ref(rneg) = 0;
%
%     % Sometimes ASTlib data is unsorted and there are repeated wavelengths
%     ma = [wave ref];
%     mas = sortrows(ma);
%     wsns = mas(:,1);
%     Rs = mas(:,2);
%
%     c=1;
%     while c==1
%
%         M = mode(wsns);
%         mf = find(wsns==M);
%         if length(mf)>1
%             wsns(mf(2:length(mf)))=[];
%             Rs(mf(2:length(mf)))=[];
%         else
%             c=2;
%         end
%     end
%
%     if i==length(modes)
%         wave_ice = wsns;
%         emis_ice = 1-Rs;
%     else
%
%         waveout(:,count) = wsns;
%         emisout(:,count) = 1-Rs;
%
%     end
%     count=count+1;
%
%
%     e_int = interp1(waveout(:,1),emisout(:,1),waveb);
%     esnow = sum(rsrb.*e_int)/sum(rsrb)
%
%
% end
%
%
%
% filenames1 = 'F:\Documents and Settings\ghulley\My Documents\MATLAB\Lab_spectra\ASTER\Speclib\Conifer.txt';
% filenames2 = 'F:\Documents and Settings\ghulley\My Documents\MATLAB\Lab_spectra\MODIS\Pine_old.txt';
%
% fid1 = fopen(filenames1);fid2 = fopen(filenames2);
%
% a1 = fscanf(fid1, '%g %g', [2 inf]);a2 = fscanf(fid2, '%g %g', [3 inf]);
%
% a1=a1';a2=a2';
% emiscon = 1 - (a1(:,2)/100);emisconm = a2(:,3);
% wavecon = a1(:,1);waveconm = a2(:,1);
%
% fclose(fid1);fclose(fid2);
%
% e_int1 = interp1(wavecon,emiscon,waveb);
% e_int2 = interp1(waveconm,emisconm,waveb);
% eintc = [e_int1 e_int2];
% e_int = mean(eintc,2);
% eveg = sum(rsrb.*e_int)/sum(rsrb)
%
% % figure()
% % plot(wavecon,emiscon,'r-');hold on
% % plot(waveconm,emisconm,'b-')
% % xlim([8 13])
% % ylim([0.95 1])










