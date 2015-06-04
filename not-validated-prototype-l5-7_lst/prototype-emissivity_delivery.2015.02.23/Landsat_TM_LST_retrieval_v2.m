%% Calculate Landsat LST given RIT atmospheric correction using ASTER GED for emissivity correction

clc
clear
fclose('all');

%% Check to see how many files in each Landsat path/row/year combo

% paths = 38:45;
%
% data = struct('Name', 'Landsat5');
% years = 1982:2011;
% for y = 1:length(years)
%     n = ['y',num2str(years(y))];
%     data.(n) = 0;
% end
%
% for j = 1:length(paths)
%
%     rows = dir(['X:\',num2str(paths(j)),'\']);
%     if isempty(rows)
%         continue
%     end
%
%     for i = 1:length(rows);
%         if strcmp(rows(i).name,'.') || strcmp(rows(i).name,'..')
%             continue
%         end
%
%         yrs = dir(['X:\',num2str(paths(j)),'\',rows(i).name,'\']);
%         if isempty(yrs)
%             continue
%         end
%
%         for k = 1:length(yrs);
%             if strcmp(yrs(k).name,'.') || strcmp(yrs(k).name,'..')
%                 continue
%             end
%
%             files = dir(['X:\',num2str(paths(j)),'\',rows(i).name,'\',yrs(k).name,'\']);
%
%             if length(files)>2
%
%                 data.(['y',yrs(k).name]) = data.(['y',yrs(k).name]) + length(files)-2;
%
%             end
%
%         end
%     end
% end


%% Inputs

% Landsat sensor
%sensor = 'LANDSAT7'; 
% Landsatdir = '\\clearlake\landsat-etm-collection\north_america';
% sres = 60;

sensor = 'LANDSAT5';
Landsatdir = '\\monolake\collections\primary\landsat\tm\north_america';
sres = 120;

% Output L2 data to geotiff (no/yes)
outputdata = 'no'; 
 
% ASTER GED info
vers = '003';  
resfolder = '100m';     
resplot = '0001';        
ASTERdir = '\\monolake\collections\primary\aster_ged\hdf5';  % ASTER GED location

% Temporary folder to unzip Landsat data
tfolder = 'C:\Users\ghulley\Documents\MATLAB\temp_landsat\';

% MODTRAN module
mmod = 'C:\Users\ghulley\Documents\MATLAB\MODTRAN5_2';

% Open Brightness temperature LUT
if strcmp(sensor,'LANDSAT7')
    
    % Compute using BT LUT-----------------------------------------------------
    fid = fopen('C:\Users\ghulley\Documents\MATLAB\BT_LUT\L7\L7.txt');
    BT = fscanf(fid, '%g %g', [2 inf]); BT=BT';
    fclose(fid);
    
elseif strcmp(sensor,'LANDSAT5')
    
    % Compute using BT LUT-----------------------------------------------------
    fid = fopen('C:\Users\ghulley\Documents\MATLAB\BT_LUT\L5\L5.txt');
    BT = fscanf(fid, '%g %g', [2 inf]); BT=BT';
    fclose(fid);
    
end


% Temp Landsat output folder
LSTfolder = 'C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\';
L1Bfolder = 'C:\Users\ghulley\Documents\MATLAB\temp_landsat\L1B\';

% Constants
c1 = 3.7418e-22;
c2 = 0.014388;
wave = 10.8686*10^-6;

K1 = 666.09;
K2 = 1282.71;
c_sky = [0.0194 0.5469 0.0254];


%% Loop files to process

% Death Valley: pathL 40, row 35, yr 2011, file 8
% Sierras: pathL 42, row 34, yr 2010, file 14

for j = 1:1 %length(paths)
    
    pathL = '39';
    %pathL = num2str(paths(j));
    
    rows = dir([Landsatdir,'\',pathL,'\']);
    
    for i = 1:1 %length(rows);
        
        rname = '37';
        %rname = rows(i).name;
        
        if strcmp(rname,'.') || strcmp(rname,'..')
            continue
        end
        
        yrs = dir([Landsatdir,'\',pathL,'\',rname,'\']);
        
        for k = 1:1 %length(yrs);
            
            yr = '2010';
            %year = yrs(k).name;
            
            if strcmp(yr,'.') || strcmp(yr,'..')
                continue
            end
            
            files = dir([Landsatdir,'\',pathL,'\',rname,'\',yr,'\*.gz']);
            
            if isempty(files)
                continue
            end
            
            for p = 8:8 %length(files);
                
                fileopen = [Landsatdir,'\',pathL,'\',rname,'\',yr,'\',files(p).name]
                
                % Untar to temporary local folder
                untar(fileopen,tfolder)
                
                %% Process Thermal band 6 data
                disp('process band 6')
                
                filest = dir([tfolder,'*_B60.TIF']);
                filem = dir([tfolder,'*MTL.txt']);
                
                f = filest(1).name;
                year = str2double(f(13:16));
                mth = str2double(f(17:18));
                day = str2double(f(19:20));
                
                
                fileband6 = [tfolder,filest(1).name];
                file_met = [tfolder,filem(1).name];
                
                % Get Tiff info
                info = geotiffinfo(fileband6);
                
                % Get TOA radiance and scale
                Qcal6 = geotiffread(fileband6);
                Qcal6 = double(Qcal6);
                
                szero = Qcal6==0;
                Qcal6(szero) = nan;
                
                % Get radiance scaling factors
                fid = fopen(file_met,'r');
                fs = fscanf(fid,'%c',inf);
                fclose(fid);
                
                %                 for t = 1:21
                %                     fget = fgets(fid);
                %                 end
                %
                %                 obst = fgets(fid);
                %                 for g = 1:length(obst)
                %                     if obst(g)=='='
                %                         obs_hr = str2double(obst(g+2:g+3));
                %                         break
                %                     end
                %                 end
                %                 fclose(fid);
                
                gf = strfind(fs,'TIME');
                obs_hr = str2double(fs(gf+18:gf+19));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                gf = strfind(fs,'LMAX_BAND6');
                Lmax = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND6');
                Lmin = str2double(fs(gf(1)+13:gf(1)+18));
                
                L = ((Lmax-Lmin)/(Qcalmax-Qcalmin))*(Qcal6-Qcalmin) + Lmin;
                Ls = size(L);
                
                L = single(L);
                
                clear Qcal6
                
                %% Geocoordinates

                disp('Read in Geocoords')
                lat1 = [info.CornerCoords.Lat(1) info.CornerCoords.Lat(2); info.CornerCoords.Lat(4) info.CornerCoords.Lat(3)];
                lon1 = [info.CornerCoords.Lon(1) info.CornerCoords.Lon(2); info.CornerCoords.Lon(4) info.CornerCoords.Lon(3)];
                
                % Expand lat/lon coords to full size
                [Xi,Yi] = meshgrid(1:2,1:2);
                [XI,YI] = meshgrid(1:1/Ls(2):2,1:1/Ls(1):2);
                latint = interp2(Xi,Yi,lat1,XI,YI);
                lonint = interp2(Xi,Yi,lon1,XI,YI);
                
                LatLsatf = latint(1:Ls(1),1:Ls(2));
                LonLsatf = lonint(1:Ls(1),1:Ls(2));
                
                % Resize Lat/Lon to 50 km for atmos correction
                sfac = 50000;
                gs1 = round((Ls(2)*sres)/sfac);
                gs2 = round((Ls(1)*sres)/sfac);
                
                % Expand lat/lon coords to 5 km grid size
                [Xi,Yi] = meshgrid(1:2,1:2);
                [XI,YI] = meshgrid(1:1/gs1:2,1:1/gs2:2);
                latint = interp2(Xi,Yi,lat1,XI,YI);
                lonint = interp2(Xi,Yi,lon1,XI,YI);
                
                LatLsat = latint(1:gs2,1:gs1);
                LonLsat = lonint(1:gs2,1:gs1);
                szl = size(LatLsat);
                
                %     LatLsatf = single(LatLsatf);
                %     LonLsatf = single(LonLsatf);
                
                clear latint lonint Xi Yi XI YI
                
                %% Atmospheric correction with NCEP data
                disp('Atmospheric correction')
                
                % Global Coordinates for NCEP data
                [latnc,lonnc] = meshgrat([-90 90],[0 359],[181 360]);
                latnc = flipud(latnc);
                
                prod = 'gdas1';
                
                % Find closest NCEP time to sensor observation time
                NCEPobs = ['00'; '06'; '12'; '18'];
                NCEP_utc = [0 6 12 18];
                NCEPobs = cellstr(NCEPobs);
                tdiff = obs_hr - NCEP_utc;
                [Cx,Ix] = min(abs(tdiff));
                
                time = NCEPobs{Ix};
                
                monthin = num2str(mth);
                din = num2str(day);
                if mth<10
                    monthin = ['0',monthin];
                end
                
                if day<10
                    din = ['0',din];
                end
                
                yr = num2str(year);
                gribloc = ['X:\',yr,'\',monthin,'\gdas1.PGrbF00.',yr(3:4),monthin,din,'.',time,'z'];
              
                
                [Tni,RHni,O3ni,gphni,P,QA,PWV] = NCEP_4Landsat(gribloc,latnc,lonnc,LatLsat,LonLsat);

                if QA==0
                    error('Oops, NCEP data did not read in correctly')
                end
                
                DEM_atm = GTOPO_get_Lsat(LatLsat,LonLsat);
                
                [pathrad,t1] = modtran_4landsat(sensor,Tni,RHni,O3ni,gphni,P,DEM_atm,mmod);
                
                
                sky = c_sky(1) + c_sky(2).*pathrad + c_sky(3).*pathrad.^2;
                
                % MODTRAN Sky Irradiance
                %sky = modtran_4landsat_sky(sensor,Tni,RHni,O3ni,gphni,P,DEM,mmod);
                
                
                % Resize to Landsat coords
                pathr = imresize(pathrad,[Ls(1) Ls(2)],'bicubic');
                t1r = imresize(t1,[Ls(1) Ls(2)],'bicubic');
                skyr = imresize(sky,[Ls(1) Ls(2)],'bicubic');
                PWVr = imresize(PWV,[Ls(1) Ls(2)],'bicubic');
                
                %% Get reflectances
                
          %----------------------------------------------------------
                % Get band 1 Tiff info
                %info = geotiffinfo(file_Lsat_3);
                
                disp('process visible bands')
                
                % Get TOA radiance and scale
                filest = dir([tfolder,'*_B10.TIF']);
                fileband1 = [tfolder,filest(1).name];
                Qcal = geotiffread(fileband1);
                Qcal = double(Qcal);
                
                gf = strfind(fs,'LMAX_BAND1');
                Lmax1 = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND1');
                Lmin1 = str2double(fs(gf(1)+13:gf(1)+18));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                L1 = ((Lmax1-Lmin1)/(Qcalmax-Qcalmin))*(Qcal-Qcalmin) + Lmin1;
                L1s = size(L1);
                L1 = single(L1);
                
                s1zero = Qcal==0;
                L1(s1zero) = nan;
                
                clear Qcal s1zero
                
                %----------------------------------------------------------
                % Get band 2 Tiff info
                %info = geotiffinfo(file_Lsat_3);
                
                % Get TOA radiance and scale
                filest = dir([tfolder,'*_B20.TIF']);
                fileband2 = [tfolder,filest(1).name];
                Qcal = geotiffread(fileband2);
                Qcal = double(Qcal);
                
                gf = strfind(fs,'LMAX_BAND2');
                Lmax2 = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND2');
                Lmin2 = str2double(fs(gf(1)+13:gf(1)+18));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                L2 = ((Lmax2-Lmin2)/(Qcalmax-Qcalmin))*(Qcal-Qcalmin) + Lmin2;
                L2s = size(L2);
                L2 = single(L2);
                
                s2zero = Qcal==0;
                L2(s2zero) = nan;
                
                clear Qcal s2zero
                
                %----------------------------------------------------------
                % Get band 3 Tiff info
                %info = geotiffinfo(file_Lsat_3);
                
                % Get TOA radiance and scale
                filest = dir([tfolder,'*_B30.TIF']);
                fileband3 = [tfolder,filest(1).name];
                Qcal = geotiffread(fileband3);
                Qcal = double(Qcal);
                
                gf = strfind(fs,'LMAX_BAND3');
                Lmax3 = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND3');
                Lmin3 = str2double(fs(gf(1)+13:gf(1)+18));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                L3 = ((Lmax3-Lmin3)/(Qcalmax-Qcalmin))*(Qcal-Qcalmin) + Lmin3;
                L3s = size(L3);
                L3 = single(L3);
                
                s3zero = Qcal==0;
                L3(s3zero) = nan;
                
                clear Qcal s3zero
                
                %----------------------------------------------------------
                % Get band 4 Tiff info
                filest = dir([tfolder,'*_B40.TIF']);
                fileband3 = [tfolder,filest(1).name];
                %info = geotiffinfo(fileband3);
                
                % Get TOA radiance and scale
                Qcal = geotiffread(fileband3);
                Qcal = double(Qcal);
                
                gf = strfind(fs,'LMAX_BAND4');
                Lmax4 = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND4');
                Lmin4 = str2double(fs(gf(1)+13:gf(1)+18));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                L4 = ((Lmax4-Lmin4)/(Qcalmax-Qcalmin))*(Qcal-Qcalmin) + Lmin4;
                L4s = size(L4);
                L4 = single(L4);
                
                s4zero = Qcal==0;
                L4(s4zero) = nan;
                
                clear Qcal s4zero
                
                %----------------------------------------------------------
                % Get band 5 Tiff info
                
                % Get TOA radiance and scale
                filest = dir([tfolder,'*_B50.TIF']);
                fileband5 = [tfolder,filest(1).name];
                Qcal = geotiffread(fileband5);
                Qcal = double(Qcal);
                
                gf = strfind(fs,'LMAX_BAND5');
                Lmax5 = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND5');
                Lmin5 = str2double(fs(gf(1)+13:gf(1)+18));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                L5 = ((Lmax5-Lmin5)/(Qcalmax-Qcalmin))*(Qcal-Qcalmin) + Lmin5;
                L5s = size(L5);
                L5 = single(L5);
                
                s5zero = Qcal==0;
                L5(s5zero) = nan;
                
                clear Qcal s5zero
                
                %----------------------------------------------------------
                % Get band 7 Tiff info
                
                % Get TOA radiance and scale
                filest = dir([tfolder,'*_B70.TIF']);
                fileband7 = [tfolder,filest(1).name];
                Qcal = geotiffread(fileband7);
                Qcal = double(Qcal);
                
                gf = strfind(fs,'LMAX_BAND7');
                Lmax7 = str2double(fs(gf(1)+13:gf(1)+18));
                
                gf = strfind(fs,'LMIN_BAND7');
                Lmin7 = str2double(fs(gf(1)+13:gf(1)+18));
                
                Qcalmin = 1;
                Qcalmax = 255;
                
                L7 = ((Lmax7-Lmin7)/(Qcalmax-Qcalmin))*(Qcal-Qcalmin) + Lmin7;
                L7s = size(L7);
                L7 = single(L7);
                
                s7zero = Qcal==0;
                L7(s7zero) = nan;
                
                clear Qcal s7zero
                
                
                %----------------------------------------------------------
                
                disp('compute reflectances')
                
                % Compute day of year, given acquisition date
                yd = dayofyear(year,mth,day);
                
                % Earth eccentricity
                e = 0.016710219;
                
                % Earth-sun distance in AU
                d = 1/((1 + e*cos((yd-4)*2*pi/365.25))/(1-e^2));
                
                % Exoatmospheric radiances (from Markham paper)
                e1 = 1957;
                e2 = 1826;
                e3 = 1554;
                e4 = 1036;
                e5 = 215;
                e7 = 80.67;
                
                
                gf = strfind(fs,'SUN_ELEVATION');
                solzen = str2double(fs(gf(1)+16:gf(1)+26));
                
                % Calculate TOA planetary albedo/reflectance
                r1 = (L1.*d^2/(e1*cos(solzen*pi/180)))*pi;
                clear L1
                r2 = (L2.*d^2/(e2*cos(solzen*pi/180)))*pi;
                clear L2
                r3 = (L3.*d^2/(e3*cos(solzen*pi/180)))*pi;
                clear L3
                r4 = (L4.*d^2/(e4*cos(solzen*pi/180)))*pi;
                clear L4
                r5 = (L5.*d^2/(e5*cos(solzen*pi/180)))*pi;
                clear L5
                r7 = (L7.*d^2/(e7*cos(solzen*pi/180)))*pi;
                clear L7
                
                % NDVI
                %NDVIS = 1.5.*(r4-r3)./(r4+r3+0.5);
                NDVI = (r4-r3)./(r4+r3);
                
                % Compute snow percentage on scene
                NDSI = (r2 - r5)./(r2 + r5);
                
                
                % Brightness temperature band 6
                Tb6 = K2./log((K1./L)+1);
                
                
                
                %% Cloud detection
                
                %     disp('Cloud processing')
                %
                %     r1 = imresize(r1,1/2,'nearest');
                %     r2 = imresize(r2,1/2,'nearest');
                %     r3 = imresize(r3,1/2,'nearest');
                %     r4 = imresize(r4,1/2,'nearest');
                %     r5 = imresize(r5,1/2,'nearest');
                %
                %     rs = size(r1);
                %     %DEM = GTOPO_get_Lsat(LatLsat,LonLsat);
                %     DEM = imresize(DEM,rs,'bicubic');
                %     Tb6 = imresize(Tb6,rs,'nearest');
                %
                %     s1 = rs(1); s2 = rs(2);
                %
                %     [CM,CM1,Tavg,Tmax,Tmin,Tstd,Skew,desratio,snowper,Cloud_percent_tot,Cloud_percent_12] = Landsat_CMASK_v2(r1,r2,r3,r4,r5,Tb6,s1,s2,DEM);
                %
                %
                %     CM = logical(CM);
                %
                %         figure()
                %         imagesc(CM)
                %
                %         break
                %
                %% Get ASTER GED Emissivity
                disp('Get Emissivity...')
                
                try
                    
                    %[eLn,tN,NDVI_ASTER] = NAALSED_4landsat_v2(LatLsatf,LonLsatf,NAALSEDloc,seas,season,naalvers,resfolder,resplot);
                    
                    [eLandsat,NDVI_ASTER] = ASTER_GED_4Landsat_v2(ASTERdir,resfolder,resplot, vers, LatLsatf,LonLsatf);
                    
                catch ME
                    disp('Problem getting emissivity')
                    continue
                end
                
                %     figure()
                %     imagesc(eLn)
                %     caxis([0.9 1])
                
                % Replace negative values with zero
                nneg = NDVI<0;
                NDVI(nneg) = 0;
                
                nneg2 = NDVI_ASTER<0;
                NDVI_ASTER(nneg) = 0;
                
                % Normalize NDVI by max value
                NDVI = NDVI./max(NDVI(:));
                NDVI_ASTER = NDVI_ASTER./max(NDVI_ASTER(:));
                
                % Calculate fractional vegetation cover
                fv_L = 1-((max(NDVI(:))-NDVI)./(max(NDVI(:))-min(NDVI(:))));
                fv_A = 1-((max(NDVI_ASTER(:))-NDVI_ASTER)./(max(NDVI_ASTER(:))-min(NDVI_ASTER(:))));
                
                % Snow index
                snow = NDSI>0.4;
                
                % Fill in holes between snow pixels
                se = strel('disk',10);
                snow = imclose(snow,se);
                
                % Adjust ASTER emissivity for vegetation and snow
                if strcmp(sensor,'LANDSAT7')
                    
                    easter_soil = (eLandsat - 0.975.*fv_N3)./(1-fv_N3);
                    easter_mod = 0.9848.*fv_L + easter_soil.*(1-fv_L);
                    easter_mod(snow) = 0.9869;  % Medium snow
                    
                    % ** NOTE: 0.975 emis for ASTER veg calculated from mean pixels above NDVI of 0.9
                    % Landsat veg and snow emis values calculated from lab
                    % spectra (ASTER Spectral library)
                    
                elseif strcmp(sensor,'LANDSAT5')
                    
                    easter_soil = (eLandsat - 0.975.*fv_N3)./(1-fv_N3);
                    easter_mod = 0.9851.*fv_L + easter_soil.*(1-fv_L);
                    easter_mod(snow) = 0.9851;   % Medium snow
                    
                end
                
                clear snow easter_soil fv_N3 fv_L NDSI NDVI NDVI_ASTER tN nneg nneg2
                clear r1 r2 r3 r4 r5
                
               
                %% LST calculation
                
                % Compute Surface radiance
                surfrad = (L - pathr)./t1r;
                
                %     if surfrad==0
                %         disp('Zero Radiance value')
                %         continue
                %     end
                
                
                R = surfrad - (1-easter_mod).*skyr;
                
                Rz = R<=0;
                R(Rz) = nan;
                
                % Account for emissivity to get Planck emitted radiance
                Re = R./easter_mod;
                szR = size(R);
                
                clear surfrad R Rz
                
                % Use Brightness temperature LUT to get skin temperature
                RadBT = BT(:,2);
                TBT = BT(:,1);
                
                Rint = reshape(Re,1,szR(1).*szR(2));
                LST = interp1(RadBT,TBT,Rint);
                LST = reshape(LST,szR);
               
                
                %wave = 10.8686*10^-6;
                %LST2 = (c2/wave).*(log((c1*easter_mod)./(pi.*R.*wave^5) + 1)).^-1;
                
                %atmostest = skyr(I,J)/(pi*L(I,J));
                %atmostest = skyi(I,J);
                
                %     figure()
                %     imagesc(LST2-LST)
                %     caxis([240 300])
                %
                %     figure()
                %     imagesc(LST)
                %     caxis([240 300])
                %
                %     break
                %     figure()
                %     imagesc(L)
                
                %geotiffwrite('C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\test.tif',TL,info);
                
                %% Output data in GEOtiff if necessary


                if strcmp(dataoutput,'yes')

                disp('Outputting Data')
                
                TL = zeros(szR(1),szR(2),4);
                TL(:,:,1) = L;
                TL(:,:,2) = LST;
                TL(:,:,3) = easter_mod;
                TL(:,:,4) = eLn;
                
                
                % Delete files in current temp Landsat folder
                cd('C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\')
                system(['del ' '*.tif']);
                system(['del ' 'terminal.txt']);
                system(['del ' 'profile.txt']);
                
                geotiffwrite('C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\LST_Emis.tif', TL, Rinfo,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
                
                clear TL R eLn Rint LST skyr L pathr tr L LatLsatf LonLsatf
                
                cd('C:\Users\ghulley\Documents\MATLAB\temp_landsat\');
                
                gzipname = [regexprep(fname,'_LST.tar.gz',''),'_LST_Emis'];
                tarname = [gzipname,'.tar'];
                tar(tarname,'C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\*');
                gzip('*.tar',loc);
                
                try
                    system(['del ' 'LE*']);
                catch ME
                    system(['del ' 'LT*']);
                end
                
                cd('C:\Users\ghulley\Documents\MATLAB\temp_landsat\L1B\')
                system(['del ' '*.tif']);
                system(['del ' '*.txt']);
                system(['del ' '*.jpg']);
                
                cd('C:\Users\ghulley\Documents\MATLAB\work')
                
                
                end
                
            end
            
        end
    end
    
end

%info = geotiffinfo('C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\test.tif');
%   [G,R] = geotiffread('C:\Users\ghulley\Documents\MATLAB\temp_landsat\LST\test.tif');

break
%% Read in data

clc
clear

% Desert Rock
sitename = 'Desert_Rock_NV';
siteoutname = 'DesertRock';
pr = '40_35';

% Bondville
sitename = 'Bondville_IL';
siteoutname = 'Bondville';
pr = '22_32';

% Table Mountain
sitename = 'Boulder_CO';
siteoutname = 'Boulder';
pr = '33_32';
%
% % Goodwin Creek
% sitename = 'Goodwin_Creek_MS';
% siteoutname = 'GoodwinCreek';
% pr = '22_36';
%
% Penn State
sitename = 'Penn_State_PA';
siteoutname = 'PennState';
pr = '16_32';
%
% % Sioux Falls
% sitename = 'Sioux_Falls_SD';
% siteoutname = 'SiouxFalls';
% pr = '29_30';

% % Fort Peck
% sitename = 'Fort_Peck_MT';
% siteoutname = 'FortPeck';
% pr = '35_26';


% File location
filename = [floc,':\Sites\',sitename,'\',siteoutname,'_',pr,'_L5.csv'];

[import, result]= readtext(filename,',', '', '', '');

fieldnames = import(1,:);
nFields = length(fieldnames);

fieldnames = regexprep(fieldnames,' ','');

% Create rudimentary struct with only one field
data= struct('NameOfImportedCSV', filename);

% Read in all fields
for i=1:nFields
    data.(fieldnames{i}) = import(2:end, i);
end

% Tsr = data.real_lst;
% pwv = data.pwv;
% BTdiff = data.lever_val;
w = data.L5_LST;

c=1;
for i = 1:length(w)
    
    Ts_L5(c) = data.L5_LST{i};
    Ts_SR(c) = data.SRAD_LST{i};
    
    Tb6m(c) = data.TB6{i};
    att(c) = data.atmostest{i};
    cloud(c) = data.cloud{i};
    
    c=c+1;
    
end


Ts_diff = Ts_L5-Ts_SR;
Tc = Ts_diff<-5;
Ts_L5(Tc) = [];
Ts_SR(Tc) = [];
Ts_diff(Tc) = [];
att(Tc) = [];
Tb6m(Tc) = [];
cloud(Tc) = [];


% Stats
L5_meandiffn = round2(mean(Ts_diff),0.01)
L5_rmsen = round2(norm(Ts_diff)./sqrt(length(Ts_diff)),0.01)

% Get R^2 statistics
Xd = [ones(size(Ts_diff))' Ts_SR'];
[b,bint,r,rint,stats] = regress(Ts_L5',Xd);
L5_R2 = round2(stats(1),0.001)


% % Plots
% sp=5;
% dlims_d = [floor(min(Ts_L5(:)))-sp ceil(max(Ts_L5(:)))+sp];
%
% figure()
% plot(Ts_SR,Ts_L5,'ro','Markerfacecolor','r','Markersize',5)
% line(dlims_d,dlims_d,'Color','k')
% xlim(dlims_d); ylim(dlims_d)
% xlabel('SURFRAD LST [K]','Fontsize',12)
% ylabel('Landsat5 LST [K]','Fontsize',12)
% set(gca,'Fontsize',12)
% title(siteoutname,'Fontsize',14)
% text(0.05,0.95,['Bias = ',num2str(L5_meandiffn),' K'],'Units','normalized','Fontsize',12)
% text(0.05,0.87,['RMSE = ',num2str(L5_rmsen),' K'],'Units','normalized','Fontsize',12)
% text(0.05,0.79,['r^{2} = ',num2str(L5_R2)],'Units','normalized','Fontsize',12)
% text(0.05,0.70,['N = ',num2str(length(Ts_diff))],'Units','normalized','Fontsize',12)
%
% set(gcf,'Color','White')
% export_fig(['Q:\Sites\',sitename,'\',siteoutname,'_LST_scatter'],'-jpeg','-r250')
%
%
% figure()
% plot(abs(Ts_diff),att,'bs','Markerfacecolor','b','Markersize',5)
% xlabel('LST error (L5 - SURFRAD) [K]','Fontsize',12)
% ylabel('Ld/L','Fontsize',12)
% set(gca,'Fontsize',12)
% title(siteoutname,'Fontsize',14)
% ylim([0 0.3])
%
% set(gcf,'Color','White')
% export_fig(['Q:\Sites\',sitename,'\',siteoutname,'_LSTerror_Atmostest'],'-jpeg','-r250')










