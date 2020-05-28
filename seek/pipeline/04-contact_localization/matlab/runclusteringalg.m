clear
close all
clc

%% Setup Paths and Directories for Registration Results and Raw EEG data_examples
addpath("/home/adam2392/Documents/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Documents/MATLAB/spm12");
% addpath("../../.lib/fieldtrip-20181108");

% run setup of global variables to make tool GUI work
ft_defaults

% neuroimaging output data dir
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/outputfiles/';
RESULTS_DIR = '/Users/adam2392/Dropbox/phd_research/data/neuroimaging_results/freesurfer_output/outputfiles';
% RESULTS_DIR = '/Users/adam2392/Downloads/neuroimaging_results/';
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/";

% subj id to analyze
subjID = 'la03';

% results directory
SUBJDIR = fullfile(RESULTS_DIR, subjID);

% output filepath for your electrode localization
COREGISTRATION_DIR = fullfile(RESULTS_DIR, subjID, 'coregistration')
elec_coords_filepath = fullfile(COREGISTRATION_DIR, [subjID '_elec_initialize.mat']);

%% SUBJECT
% load brainmask data
% filename  = fullfile(datapath,sprintf('EFRI_%02d',subjID),'mri','orig.mgz') ;
% orig = MRIread(filename) ;
% orig.name = sprintf('EFRI_%02d',subjID) ;
% 
% % load CT data
% filename = fullfile(datapath,sprintf('EFRI_%02d',subjID),'ct','postopCT.mgz') ;
% % filename = fullfile(datapath,sprintf('EFRI_%02d',subjID),'elec_recon','flirt','ctINt1.mgz') ;
% CT = MRIread(filename) ;
% CT.name  = sprintf('EFRI_%02d',subjID) ;
% 
% % load brainmask data
% filename  = fullfile(datapath,sprintf('EFRI_%02d',subjID),'mri','brainmask.mgz') ;
% brainmask = MRIread(filename) ;
% brainmask.name = sprintf('EFRI_%02d',subjID) ;
% 
% % load brainmask in CT data
% filename = fullfile(datapath,sprintf('EFRI_%02d',subjID),'elec_recon','flirt','brainmask_in_CT.mgz') ;
% brainmask_in_CT = MRIread(filename) ;
% brainmask_in_CT.name = sprintf('EFRI_%02d',subjID) ;
% registration matrix
% filename = fullfile(datapath,sprintf('EFRI_%02d',subjID),'elec_recon','flirt','ct2t1.dat') ;
% fid = fopen(filename,'r') ;
% A = textscan(fid,'%f %f %f %f',4,'HeaderLines',4) ;
% R = cell2mat(A) ;
% filename = fullfile(datapath,sprintf('EFRI_%02d',subjID),'elec_recon','flirt','ct2t1.mat') ;
% ss = load(filename,'-ascii') ;

% ctimgfile = fullfile(strcat('CT_IN_T1_mgz.nii.gz'))
ctimgfile = fullfile(strcat('CT.nii.gz'))
t1imgfile = fullfile(strcat('T1.nii.gz'))
brainmaskfile = fullfile('brainmask.mgz');

orig = ft_read_mri(fullfile(SUBJDIR, t1imgfile));
% read in CT in original format
CT = ft_read_mri(fullfile(SUBJDIR, ctimgfile));
brainmask = ft_read_mri(fullfile(SUBJDIR, brainmaskfile));

% params.letters  = {'I','T','O','J','Q','B','E','U','C','F','X','P'} ;
params.num      = [10 ,10 ,10 ,10 ,10, 10, 10, 10, 10, 10, 10, 10 ] ;
params.th = 0.75
% params.D  = 3.5

subjID    = 3 ;
% % % % params.letters  =  ;
% % % params.th = 0.75
% % % params.depth.NumElec = 12 ;
% % % params.depth.letter  = {'I','T','O','J','Q','B','E','U','C','F','X','P'} ;
% % params.depth.num     = [10 ,10 ,10 ,10 ,10, 10, 10, 10, 10, 10, 10, 10 ] ;
% params.D  = 1.5+2.5+2.5 %5.5%3.5

%% Get MRI, CT and Brainmask Anatomy Matrices
X = orig.anatomy       ;
X = X(1:end,1:end-1,1:end-2) ;

B = brainmask.anatomy ;
B = B(1:end,1:end-1,1:end-2) ;

% B_in_CT = brainmask_in_CT.vol ;
% B_in_CT = B_in_CT(1:end,1:end-1,1:end-2) ;

Y = CT.anatomy       ;
Y = Y(1:end,1:end-1,1:end-2) ;

r_MRI = 0:size(X,1)-1; %
c_MRI = 0:size(X,2)-1; %
s_MRI = 0:size(X,3)-1; %

r_CT = 0:size(Y,1)-1; %
c_CT = 0:size(Y,2)-1; %
s_CT = 0:size(Y,3)-1; %

% [R_MRI,C_MRI,S_MRI] = ndgrid(r_MRI,c_MRI,s_MRI) ;

% [f,v_CRS_MRI] = isosurface(c_MRI,r_MRI,s_MRI,double(B>0)) ;
% 
% % v_CRS_CT = inv( inv(orig.tkrvox2ras) * inv(R) * CT.tkrvox2ras ) * [ v_CRS_MRI.' ; ones(1,size(v_CRS_MRI,1))] ;
% v_CRS_CT = inv(CT.tkrvox2ras) * R * orig.tkrvox2ras  * [ v_CRS_MRI.' ; ones(1,size(v_CRS_MRI,1))] ;
% v_CRS_CT = v_CRS_CT(1:3,:).' ;
% 
% v_RAS_CT = R * orig.tkrvox2ras  * [ v_CRS_MRI.' ; ones(1,size(v_CRS_MRI,1))] ;
% v_RAS_CT = v_RAS_CT(1:3,:).' ;


% figure
% s = slice(c_MRI,r_MRI,s_MRI,X,150,90,100) ;
% for i=1:length(s)
%     s(i).EdgeColor='none';
% end
% xlabel('c')
% ylabel('r')
% zlabel('s')
% daspect(1./orig.volres)
% hold on
% x0   = inv(orig.vox2ras) * [0;0;0;1] ;
% bCRS = inv(orig.vox2ras) * [eye(3);ones(1,3)] - repmat(x0,[1,3]) ;
% vR = bCRS(1:3,1) ;
% vA = bCRS(1:3,2) ;
% vS = bCRS(1:3,3) ;
% quiver3(x0(1),x0(2),x0(3),vR(1),vR(2),vR(3),100,'-r')
% quiver3(x0(1),x0(2),x0(3),vA(1),vA(2),vA(3),100,'-g')
% quiver3(x0(1),x0(2),x0(3),vS(1),vS(2),vS(3),100,'-b')
% 
% % p = patch('Vertices',v_CRS_MRI,'Faces',f) ;
% % p.EdgeColor = 'none';
% % p.FaceColor = 'yellow'
% 
% colormap(gray)
% 
% 
% figure
% s = slice(c_MRI,r_MRI,s_MRI,B,150,90,100) ;
% for i=1:length(s)
%     s(i).EdgeColor='none';
% end
% xlabel('c')
% ylabel('r')
% zlabel('s')
% daspect(1./brainmask.volres)
% hold on
% x0   = inv(brainmask.vox2ras) * [0;0;0;1] ;
% bCRS = inv(brainmask.vox2ras) * [eye(3);ones(1,3)] - repmat(x0,[1,3]) ;
% vR = bCRS(1:3,1) ;
% vA = bCRS(1:3,2) ;
% vS = bCRS(1:3,3) ;
% quiver3(x0(1),x0(2),x0(3),vR(1),vR(2),vR(3),100,'-r')
% quiver3(x0(1),x0(2),x0(3),vA(1),vA(2),vA(3),100,'-g')
% quiver3(x0(1),x0(2),x0(3),vS(1),vS(2),vS(3),100,'-b')
% 
% colormap(gray)
% % p = patch('Vertices',v_CRS_MRI,'Faces',f) ;
% % p.EdgeColor = 'none';
% % p.FaceColor = 'yellow'
% 
% 
% figure
% s = slice(c_CT,r_CT,s_CT,Y,150,90,100) ;
% for i=1:length(s)
%     s(i).EdgeColor='none';
% end
% xlabel('c')
% ylabel('r')
% zlabel('s')
% daspect(1./CT.volres)
% hold on
% x0   = inv(CT.vox2ras) * [0;0;0;1] ;
% bCRS = inv(CT.vox2ras) * [eye(3);ones(1,3)] - repmat(x0,[1,3]) ;
% vR = bCRS(1:3,1) ;
% vA = bCRS(1:3,2) ;
% vS = bCRS(1:3,3) ;
% quiver3(x0(1),x0(2),x0(3),vR(1),vR(2),vR(3),100,'-r')
% quiver3(x0(1),x0(2),x0(3),vA(1),vA(2),vA(3),100,'-g')
% quiver3(x0(1),x0(2),x0(3),vS(1),vS(2),vS(3),100,'-b')
% colormap(gray)
% 
% p = patch('Vertices',v_CRS_CT,'Faces',f) ;
% p.EdgeColor = 'none';
% p.FaceColor = 'yellow';
% 
% 
% figure
% p = patch('Vertices',v_RAS_CT,'Faces',f) ;
% p.EdgeColor = 'none';
% p.FaceColor = 'yellow';
% daspect([1 1 1])
% 
% 
% figure
% s = slice(c_CT,r_CT,s_CT,Y,90,90,100) ;
% for i=1:length(s)
%     s(i).EdgeColor='none';
% end
% xlabel('c')
% ylabel('r')
% zlabel('s')
% daspect(1./CT.volres)
% hold on
% x0   = inv(CT.vox2ras) * [0;0;0;1] ;
% bCRS = inv(CT.vox2ras)*[eye(3);ones(1,3)] - repmat(x0,[1,3]) ;
% vR = bCRS(1:3,1) ;
% vA = bCRS(1:3,2) ;
% vS = bCRS(1:3,3) ;
% quiver3(x0(1),x0(2),x0(3),vR(1),vR(2),vR(3),100,'-r')
% quiver3(x0(1),x0(2),x0(3),vA(1),vA(2),vA(3),100,'-g')
% quiver3(x0(1),x0(2),x0(3),vS(1),vS(2),vS(3),100,'-b')
% 
% colormap(gray)


%
% % filename = fullfile(datapath,sprintf('EFRI_%02d',subjID),'elec_recon','flirt','ct2t1.mat') ;
% % load(filename,'-ASCII')
% % 
% % ct2t1

%% Create XYZ Coordinates of the CT Scan
[C_CT,R_CT,S_CT] = ndgrid(c_CT,r_CT,s_CT) ;
CRS_CT = [C_CT(:),R_CT(:),S_CT(:)];

% ras_CT = CT.vox2ras * [CRS_CT.';ones(1,size(CRS_CT,1))] ;
ras_CT = CT.transform * [CRS_CT.'; ones(1, size(CRS_CT, 1))];
ras_CT = ras_CT(1:3,:).' ;

X_CT = reshape(ras_CT(:,1),size(C_CT)) ;
Y_CT = reshape(ras_CT(:,2),size(C_CT)) ;
Z_CT = reshape(ras_CT(:,3),size(C_CT)) ;


%% Find electrodes
figure
subplot(2,1,1)
histogram(Y(:))
set(gca,'XLim',[0 255])

subplot(2,1,2)
histogram(Y(:)/max(Y(:)))
set(gca,'XLim',[0 1])

% tmp = Y/255 ;

% % % tmp( tmp < 0.5 )    = NaN ;
% % % tmp( B_in_CT == 0 ) = NaN ;
% % 
% % yyy = tmp(~isnan(tmp)) ;
% % yyy = sort(yyy(:))     ;
% % 
% % [N,edges] = histcounts(yyy(:),'Normalization','cdf') ;
% % 
% % % [f,x] = ecdf(yyy(:)) ;
% % 
% % figure
% % plot(1/2*(edges(1:end-1) + edges(2:end)),N)
% % 

type = 'remove_if_partially_outside' ;

thvec  = 0.5:0.005:1      ;
NumObj = nan(size(thvec)) ;
CC = cell(size(thvec))    ;
for i = 1:length(thvec)
    
    i/length(thvec)
    CC{i} = bwconncomp(Y/255 > thvec(i)) ;
    
    for k = 1:CC{i}.NumObjects
        
        k/CC{i}.NumObjects ;
        
        tmp  = CC{i}.PixelIdxList{k} ;

%         IN = (B_in_CT(tmp) > 0) ;
%         
%         % keep inside
%         if all(~IN) % all outside brain
%             tmp = [] ;
%         elseif any(~IN)
%             switch type
%                 case 'keep_if_partially_outside'
%                     % nothing
%                 case 'crop_if_partially_outside'
%                     tmp = tmp(IN) ;
%                 case 'remove_if_partially_outside'
%                     tmp = [];
%             end
%         end
        
        % store new cluster
        CC{i}.PixelIdxList{k} = tmp ;
        
    end
    iE = cellfun(@isempty,CC{i}.PixelIdxList)  ;
    CC{i}.PixelIdxList(iE) = [];
    CC{i}.NumObjects = length(CC{i}.PixelIdxList) ;
        
    NumObj(i) = CC{i}.NumObjects ;
end

edges = (0:max(NumObj)+1)-0.5 ;
NN    = histcounts(NumObj,edges) ;

[NNmax,indmax] = max(NN) ;

NNthres = mean(edges(indmax+[0,1]))

figure
histogram(NumObj,50:250)

% ind = find( diff(NumObj) > 0 , 1 , 'last' ) 

ind1 = find( NumObj > sum(params.num) , 1 , 'last' ) 
ind2 = find( NumObj == NNthres , 1 , 'first' ) 

ttt = [sum(params.num),NNthres]
iii = [ind1,ind2]
cmap_th = distinguishable_colors(length(iii)) ;


figure
subplot(3,1,1)
plot(thvec,NumObj)
hold on
for j = 1:length(iii)
    plot(thvec,ttt(j)*ones(size(thvec)),'-','Color',cmap_th(j,:))
    plot(thvec(iii(j)),NumObj(iii(j)),'o','Color',cmap_th(j,:))
end
subplot(3,1,2)
plot(thvec,zeros(size(thvec)),':k')
hold on
plot(thvec(2:end),diff(NumObj))
subplot(3,1,3)
plot(thvec(3:end),diff(NumObj,2))



figure
for j = 1:length(iii) % 1:length(thvec)
    
    i = iii(j) ;
    
    % subplot(3,3,i)
            
    hold on
    zzz = zeros(size(Y)) ;
    for k = 1:CC{i}.NumObjects% min(5,CCtmp.NumObjects)
        zzz(CC{i}.PixelIdxList{k}) = 1 ;
    end
    p = patch(isosurface(c_CT,r_CT,s_CT,zzz,0.5)) ;
    p.EdgeColor = 'none' ;
    p.FaceColor = cmap_th(j,:) ;
    xlabel('c')
    ylabel('r')
    zlabel('s')
    daspect(1./CT.volres)
    set(gca,'XLim',c_CT([1 end]))
    set(gca,'YLim',r_CT([1 end]))
    set(gca,'ZLim',s_CT([1 end]))
    
end

% hax   = findobj(gcf,'Type','Axes') ;
% hlink = linkprop(hax,{'CameraUpVector','CameraPosition','CameraTarget'}) ;
% setappdata(gcf,'StoreTheLink',link)

% % % 
% % % 
% % % keyboard
% % % 
% % % 
% % % 
% % % 
% % % 
% % % thvec = 0.6:0.05:0.9 %thvec([ind-1 ind ind+1]) % [0.875 0.88 0.885]%0.6:0.1:0.9;
% % % cmap_th = distinguishable_colors(length(thvec)) ;
% % % 
% % % figure
% % % 
% % % % % p = patch(isosurface(c_CT,r_CT,s_CT,double(B_in_CT > 0),0.5)) ;
% % % % % p.EdgeColor = 'none' ;
% % % % % p.FaceColor = 'black' ;
% % % % % p.FaceAlpha = 0.1 ;
% % % hold on
% % % for i = 1:length(thvec)
% % %     
% % %     i/length(thvec)
% % %         
% % %     zzz = (tmp > thvec(i) & B_in_CT > 0) ;
% % %     p = patch(isosurface(c_CT,r_CT,s_CT,zzz,0.5)) ;
% % %     p.EdgeColor = 'none' ;
% % %     p.FaceColor = cmap_th(i,:) ;
% % %     
% % % end
% % % legend(string(thvec))
% % % s = slice(c_CT,r_CT,s_CT,Y,median(c_CT),median(r_CT),median(s_CT)) ;
% % % for i=1:length(s)
% % %     s(i).EdgeColor='none';
% % % end
% % % 
% % % keyboard
% % % 
% % % 
% % % figure
% % % for i = 1:length(thvec)
% % %     
% % %     i/length(thvec)
% % %     
% % %     subplot(3,3,i)
% % %     
% % %     CCtmp = bwconncomp(tmp > thvec(i) & B_in_CT > 0 ) ;
% % %     
% % %     cmap_G = distinguishable_colors(CCtmp.NumObjects) ;
% % %     
% % %     hold on
% % %     for k = 1:CCtmp.NumObjects% min(5,CCtmp.NumObjects)
% % %         zzz = zeros(size(tmp)) ;
% % %         zzz(CCtmp.PixelIdxList{k}) = 1 ;
% % %         
% % %         p = patch(isosurface(c_CT,r_CT,s_CT,zzz,0.5)) ;
% % %         p.EdgeColor = 'none' ;
% % %         p.FaceColor = cmap_G(k,:) ;
% % %     end
% % %     xlabel('c')
% % %     ylabel('r')
% % %     zlabel('s')
% % %     daspect(1./CT.volres)
% % %     set(gca,'XLim',c_CT([1 end]))
% % %     set(gca,'YLim',r_CT([1 end]))
% % %     set(gca,'ZLim',s_CT([1 end]))
% % %         
% % % end
% % % 
% % % hax   = findobj(gcf,'Type','Axes') ;
% % % hlink = linkprop(hax,{'CameraUpVector','CameraPosition','CameraTarget'}) ;
% % % % setappdata(gcf,'StoreTheLink',link)
% % % 
% % % keyboard
% % % 
% % % 
% % % figure
% % % subplot(2,1,1)
% % % histogram(tmp(:))
% % % set(gca,'XLim',[0 255])
% % % 
% % % subplot(2,1,2)
% % % histogram(tmp(:)/max(tmp(:)))
% % % set(gca,'XLim',[0 1])
% % % 
% % % 
% % % figure
% % % histogram(xxx(~isnan(xxx)),'Normalization','cdf','DisplayStyle','stairs')
% % % 
% % % 
% % % 
% % % CC = bwconncomp((Y/255) > params.th) ;
% % % 
% % % % remove clustered_voxels that are outside brain shell
% % % type = 'crop_if_partially_outside' ;
% % % 
% % % for k = 1:CC.NumObjects
% % %     
% % %     k/CC.NumObjects
% % %     
% % %     tmp  = CC.PixelIdxList{k} ;
% % %     QPTS = CRS_CT(tmp,:) ;
% % %     IN   = inpolyhedron(f,v_CRS_CT,QPTS) ;
% % %     
% % %     % keep inside
% % %     if all(~IN)
% % %         tmp = [] ;
% % %     elseif any(~IN)
% % %         switch type
% % %             case 'keep_if_partially_outside'
% % %                 % nothing
% % %             case 'crop_if_partially_outside'
% % %                 tmp = tmp(IN) ;
% % %             case 'remove_if_partially_outside'
% % %                 tmp = [];
% % %         end
% % %     end
% % % 
% % %     % store new cluster
% % %     CC.PixelIdxList{k} = tmp ;
% % %     
% % % end
% % % iE = cellfun(@isempty,CC.PixelIdxList)  ;
% % % CC.PixelIdxList(iE) = [];
% % % CC.NumObjects = length(CC.PixelIdxList) ;


CC = CC{ind2} ;

cmap_G = distinguishable_colors(CC.NumObjects) ;

% % figure
% % p = patch('Vertices',v_RAS_CT,'Faces',f) ;
% % p.EdgeColor = 'none';
% % p.FaceColor = 'yellow';
% % p.FaceAlpha = 0.1
% % 
% % hold on
% % zzz = zeros(size(Y)) ;
% % for k = 1:CC.NumObjects% min(5,CCtmp.NumObjects)
% %     zzz(CC.PixelIdxList{k}) = 1 ;
% % end
% % 
% % [ff,vv] = isosurface(c_CT,r_CT,s_CT,zzz,0.5) ;
% % 
% % vv = CT.tkrvox2ras * [ vv.' ; ones(1,size(vv,1))] ;
% % vv = vv(1:3,:).' ;
% % 
% % 
% % p = patch('Vertices',vv,'Faces',ff) ;
% % p.EdgeColor = 'none' ;
% % p.FaceColor = 'blue'%cmap_G(k,:) ;
% % xlabel('r')
% % ylabel('a')
% % zlabel('s')
% % daspect([1 1 1])
% % % set(gca,'XLim',c_CT([1 end]))
% % % set(gca,'YLim',r_CT([1 end]))
% % % set(gca,'ZLim',s_CT([1 end]))



% compute distance between cluster
D = NaN(CC.NumObjects) ;
for j = 1:CC.NumObjects
    for k = 1:CC.NumObjects
        if j < k
            indj = CC.PixelIdxList{j};
            indk = CC.PixelIdxList{k};
                        
            d_tmp = pdist2(ras_CT(indj,:),ras_CT(indk,:),'euclidean') ;
            d_min = min(d_tmp(:)) ;
            
            D(j,k) = d_min;
            D(k,j) = d_min;
        end
    end
end
dMin = min(D) ;



figure
histogram(dMin)

D(logical(eye(size(D)))) = 0 ;

figure
pcolor(D)


Dvec = 0:0.1:20 ;
for i = 1:length(Dvec)
    D_sat = D ;
    D_sat(D>Dvec(i)) = 0;
    G    = graph(D_sat) ;
    bins = conncomp(G)  ;
    
    Bmax(i) = max(bins(:)) ;
end


figure
plot(Dvec,Bmax)



return


% % 
% % % remove point that are far from the rest
% % iE = dMin > 10 ;
% % CC.PixelIdxList(iE) = [];
% % CC.NumObjects = length(CC.PixelIdxList) ;
% % dMin(iE) = [];

cmap_G = distinguishable_colors(CC.NumObjects) ;


keyboard


figure
hold on
s = slice(c_CT,r_CT,s_CT,Y,median(c_CT),median(r_CT),median(s_CT)) ;
for i=1:length(s)
    s(i).EdgeColor='none';
end
xlabel('c')
ylabel('r')
zlabel('s')
daspect(1./CT.volres)
meanCRS = zeros(CC.NumObjects,3) ;
for i = 1:CC.NumObjects
    ind = CC.PixelIdxList{i} ;
    
    siz = size(Y)   ;
    [rr,cc,ss] = ind2sub(siz,ind) ;
    cc = cc-1;
    rr = rr-1;
    ss = ss-1;
    tmpCRS = [cc,rr,ss] ;
        
    meanCRS(i,:) = mean( tmpCRS , 1 ) ;
    
    cmap = jet(256)  ;
    cmin = 0;
    cmax = max(dMin) ;
    m = size(cmap,1) ;
    CData = dMin(i) ;
    
    idx = min(m,round((m-1)*(CData-cmin)/(cmax-cmin))+1) ;
    
    scatter3(meanCRS(i,1),meanCRS(i,2),meanCRS(i,3),[],cmap(idx,:),'.')
        
    % text(meanCRS(i,1),meanCRS(i,2),meanCRS(i,3),sprintf('%03d %03d %03d',meanCRS(i,1),meanCRS(i,2),meanCRS(i,3)),'Fontsize',10,'Color','red')
    
end

colormap(gray)

TT = zeros(size(Y));
for ic = 1:CC.NumObjects
    TT(CC.PixelIdxList{ic}) = ic ;
end


sz = size(Y);
for l = 1:size(Y,3)
    
    100*l / size(Y,3)
    
    if any(any(TT(:,:,l)>0))
        
        
        filename = sprintf('EFRI_%02d_z_%03d',subjID,l) ;
        savename = fullfile( '_fig' , filename ) ;
        
        h_fig = figure ;
        
        imshow(Y(:,:,l) / 255)
        hold on
        
        for ic = 1:CC.NumObjects
            
            % ic/CC.NumObjects
            
            I = false(size(Y)) ;
            I(CC.PixelIdxList{ic}) = true ;
            
            if any(any(I(:,:,l),1),2)
                r_CC = round(mean(R_CT(I))) ;
                c_CC = round(mean(C_CT(I))) ;
                s_CC = round(mean(S_CT(I))) ;
                
                col   = cat(3,repmat(cmap_G(ic,1),sz(1),sz(2)),repmat(cmap_G(ic,2),sz(1),sz(2)),repmat(cmap_G(ic,3),sz(1),sz(2))) ;
                h_col = imshow( col ) ;
                set(h_col,'AlphaData',I(:,:,l))
                str = sprintf('%d,%d,%d',c_CC,r_CC,s_CC) ;
                plot(c_CC,r_CC,'.w')
                h_txt = text(c_CC,r_CC-5,str,'Color','w','FontSize',5) ;
                set(h_txt,'Rotation',90)
            end
            
        end
        
        daspect([1 1 1])
        title(sprintf('th=%0.3f, z=%03d',params.th,l))
        
        axis on
        print( h_fig , savename ,'-djpeg','-r300')
        
        close
    end
end







return


% % figure
% % s = slice(c_MRI,r_MRI,s_MRI,brainmask.vol,150,90,100) ;
% % for i=1:length(s)
% %     s(i).EdgeColor='none';
% % end
% % xlabel('c')
% % ylabel('r')
% % zlabel('s')
% % daspect(1./brainmask.volres)
% % hold on
% % for i = 1:CC.NumObjects
% %     ind = CC.PixelIdxList{i} ;
% %     
% %     % tmpCRS = CRS_CT(ind,:) 
% %     
% %     siz = size(Y)   ;
% %     [rr,cc,ss] = ind2sub(siz,ind) ;
% %     cc = cc-1;
% %     rr = rr-1;
% %     ss = ss-1;
% %     tmpCRS = [cc,rr,ss] ;
% %         
% %     tmpCRS = mean( tmpCRS , 1 ) ;
% %         
% %     tmp1 = inv(orig.tkrvox2ras) * inv(R) * CT.tkrvox2ras * [ tmpCRS(:) ; 1] ;
% % 
% %     scatter3(tmp1(1),tmp1(2),tmp1(3),'.m')
% %     %scatter3(tmp2(1),tmp2(2),tmp2(3),'.g')
% %     
% % end
% % hold on
% % x0   = inv(brainmask.vox2ras) * [0;0;0;1] ;
% % bCRS = inv(brainmask.vox2ras) * [eye(3);ones(1,3)] - repmat(x0,[1,3]) ;
% % vR = bCRS(1:3,1) ;
% % vA = bCRS(1:3,2) ;
% % vS = bCRS(1:3,3) ;
% % quiver3(x0(1),x0(2),x0(3),vR(1),vR(2),vR(3),100,'-r')
% % quiver3(x0(1),x0(2),x0(3),vA(1),vA(2),vA(3),100,'-g')
% % quiver3(x0(1),x0(2),x0(3),vS(1),vS(2),vS(3),100,'-b')
% % 
% % 
% % colormap(gray)





return

%%

% % % Rescale image to start at 0.
% % minPixel = min(CT.vol(:));
% % maxPixel = max(CT.vol(:));
% % m = 1/(maxPixels - minPixel) ;
% % W_CT = imlincomb(m, double(CT.vol), -(m * minPixel));


for l = 100%1:CT.volsize(3)
    
    figure
    imshow( CT.vol(:,:,l) / 255 )
    colorbar
    
end

%%

CRS_CT         = [149;320;106;1] ;
tkrRAS      = CT.tkrvox2ras * CRS_CT ;
tkrRAS_int1 = inv(R) * CT.tkrvox2ras * CRS_CT ;

% % 
% % test    = ct2t1 * CRS_CT


%%

close all;

params.savepath = fullfile( '_fig' , sprintf('EFRI%02d',subjID) ) ;
params.R = R 


% elec_X = find_elec_X_CT(CT,params)

save( sprintf('EFRI%02d-electrodes.mat',subjID) , 'elec_X' )

%%

function elec_X = find_elec_X_CT(CT,params)

% params
th   = params.th %0.75
% maxD = 5
%

X = CT.vol / 255 ;

% X = X(1:256,:,:) ;


% Find convex hull of brain
h   = 0.01 ;
x_0 = find_grey_peak(X,h) ;
x_1 = x_0 - h/2;
x_2 = x_0 + h/2;

figure
histogram(X(:),0:h:1)
hold on
plot(x_1*[1 1],get(gca,'YLim'))
plot(x_2*[1 1],get(gca,'YLim'))

B = double(x_1 <= X & X <= x_2) ;
C = bwconncomp(B) ;
P = cellfun(@(x) length(x)/prod(C.ImageSize) , C.PixelIdxList) ;

[~,ind] = max( P ) ;

brainmask = false(size(X)) ;
brainmask(C.PixelIdxList{ind}) = true ;
tmp = brainmask ;

% SE = strel('sphere',5) ;
% tmp = imerode(tmp,SE) ;
% tmp = imdilate(tmp,SE) ;

%


for l = 100%1:CT.volsize(3)
    
    figure
    subplot(2,2,1)
    imshow( ind2rgb(CT.vol(:,:,l) + 1 , gray(256)) )
    hold on
    contour( tmp(:,:,l) , 0.5*[1 1] , '-r' )
    keyboard
    
end


r_CT  = 0:size(X,1)-1; %
c_CT  = 0:size(X,2)-1; %
s_CT  = 0:size(X,3)-1; %

[R_CT,C_CT,S_CT] = ndgrid(r_CT,c_CT,s_CT);

RCS = [R_CT(:),C_CT(:),S_CT(:)]       ;
rcs = RCS .* CT.volres                ;
V   = RCS(brainmask(:),:)             ;
F   = convhull(V(:,1),V(:,2),V(:,3))  ;

% % ind = find(X(:) > 0.5);
% % IN  = inpolyhedron(F,V, RCS(ind,:)) ;
% % ind = ind(IN);
% %
% % figure
% % histogram(X(ind),0:h:1)



% Find contact candidates (bright point)
CC = bwconncomp(X > th) ;


% remove clustered_voxels that are outside brain shell
for k = 1:CC.NumObjects
    
    k/CC.NumObjects
    
    tmp  = CC.PixelIdxList{k} ;
    QPTS = RCS(tmp,:) ;
    IN   = inpolyhedron(F,V,QPTS) ;
    
    % keep inside
    tmp  = tmp(IN) ;
    
    % remove is partially outside
    %     if any(~IN)
    %         tmp = [];
    %     end
    
    % store new cluster
    CC.PixelIdxList{k} = tmp ;
    
end

iE = cellfun(@(x) isempty(x) , CC.PixelIdxList) ;
CC.PixelIdxList(iE) = [];
CC.NumObjects = length(CC.PixelIdxList) ;


save
close all

% plot
dirname = fullfile( params.savepath , sprintf('CT-col-th-%0.3f',th) ) ;
if not(exist( dirname , 'dir' ))
    mkdir( dirname )
end

cmap_G = distinguishable_colors(CC.NumObjects) ;

sz = size(X);

for l = 1:size(X,3)

    100*l / size(X,3)

    filename = sprintf('z-%03d',l) ;
    savename = fullfile( dirname , filename ) ;

    h_fig = figure ;

    imshow(X(:,:,l))
    hold on

    for ic = 1:CC.NumObjects

        % ic/CC.NumObjects

        I = false(sz) ;
        I(CC.PixelIdxList{ic}) = true ;

        if any(any(I(:,:,l),1),2)
            r_CC = round(mean(R_CT(I))) ;
            c_CC = round(mean(C_CT(I))) ;
            s_CC = round(mean(S_CT(I))) ;

            col   = cat(3,repmat(cmap_G(ic,1),sz(1),sz(2)),repmat(cmap_G(ic,2),sz(1),sz(2)),repmat(cmap_G(ic,3),sz(1),sz(2))) ;
            h_col = imshow( col ) ;
            set(h_col,'AlphaData',I(:,:,l))
            str = sprintf('%d,%d,%d',c_CC,r_CC,s_CC) ;
            plot(c_CC,r_CC,'.w')
            h_txt = text(c_CC,r_CC-5,str,'Color','w','FontSize',5) ;
            set(h_txt,'Rotation',90)
        end

    end

    daspect([1 1 1])
    title(sprintf('th=%0.3f, z=%03d',th,l))

    axis on
    print( h_fig , savename ,'-djpeg','-r300')

    close

end

keyboard




% compute distance between cluster
D = zeros(CC.NumObjects) ;
for j = 1:CC.NumObjects
    for k = 1:CC.NumObjects
        
        indj = CC.PixelIdxList{j};
        indk = CC.PixelIdxList{k};
        
        d_tmp = pdist2(rcs(indj,:),rcs(indk,:),'euclidean') ;
        d_tmp = min(d_tmp(:)) ;
        
        D(j,k) = d_tmp;
        D(k,j) = d_tmp;
        
    end
end


% RCS_mean = cellfun( @(x) mean(RCS(x,:),1) , CC.PixelIdxList , 'UniformOutput',false) ;
RCS_mean = cellfun( @(x) median(RCS(x,:),1) , CC.PixelIdxList , 'UniformOutput',false) ;
RCS_mean = RCS_mean(:) ;
RCS_mean = cell2mat(RCS_mean) ;
rcs_mean = RCS_mean .* CT.volres ;

filename = sprintf('%s-scatter3-th-%0.3f',CT.name,th) ;
savename = fullfile( dirname , filename ) ;

close all

h_fig = figure;
h = slice(r_CT,c_CT,s_CT,permute(X,[2,1,3]),100,100,100) ;
set(h,'EdgeColor','none')
set(h,'FaceAlpha',0.25)
colormap(gray)
hold on
% p = patch('Faces',F,'Vertices',V) ;
% p.FaceColor = 'red' ;
% p.EdgeColor = 'none';
% set(p,'FaceAlpha',0.5)
scatter3(RCS_mean(:,1),RCS_mean(:,2),RCS_mean(:,3),'filled')
daspect( 1./CT.volres )
saveas( h_fig , [savename,'.fig'] ,'fig' )



figure
imagesc(D)


% create graph (connected n)
% % D_sat = D ;
% % D_sat(D>maxD) = 0;
% % G    = graph(D_sat) ;

close all

%Dmin = 1.5 ;
%Dmax = 1.5+2*2.5 ;

Dvec = params.D ;% 3.5%Dmin:0.5:Dmax ;
D_sat = D ;
D_sat(D>Dvec) = 0;
G    = graph(D_sat) ;

bins = conncomp(G)  ;

%%
cmap_bins = distinguishable_colors(max(bins)) ;

figure;
h = slice(r_CT,c_CT,s_CT,permute(X,[2,1,3]),round(r_CT(end)/2),round(c_CT(end)/2),round(s_CT(end)/2)) ;
set(h,'EdgeColor','none')
set(h,'FaceAlpha',0.25)
colormap(gray)
hold on
scatter3(RCS_mean(:,1),RCS_mean(:,2),RCS_mean(:,3),[],cmap_bins(bins,:),'filled')

%%

new_bins = NaN(size(bins)) ;

for i = 1:max(bins)
    
    ind_bins = bins==i ;
    if sum(ind_bins) > 5
        
        % find line going through points
        I_bins = cat(1,CC.PixelIdxList{ind_bins}) ;
        X_bins = rcs(I_bins,:) ;
        [p0,s] = lineReg(X_bins)    ; 
        
        D_el = cellfun(@(x) min( dist2Line(rcs(x,:),p0,s) ) , CC.PixelIdxList ) ;
        
        ind_bins_new = D_el < 2 ;
        
        figure(1)
        hold on
        histogram(D_el,0:0.25:max(D_el))
        
        if any(not(isnan(new_bins(ind_bins_new))))
            error('trying to assign already assigned...')
        end
        
        new_bins(ind_bins_new) = i ;
                
        % find distance between line and all points
        % M0M1 = 
        
    end
    
end

[C,ia,ic] = unique( new_bins ) ;
ic(isnan(new_bins)) = min( ic(isnan(new_bins)) ) ;
bins = ic ;


CC.ElecName = zeros(size(CC.PixelIdxList)) ;

for i = 1:max(bins)-1
    
    ind_bins = find( bins==i ) ;
    
    % find line going through points
    % I_bins   = cat(1,CC.PixelIdxList{ind_bins}) ;
    % X_bins   = RCS(I_bins,:) ;
    % avg_bins = mean(RCS,1) ;
    
    RCS_bins_mean = RCS_mean(ind_bins,:) ;
    switch sign( mean(RCS_bins_mean(:,2),1) - (c_CT(end)-c_CT(1))/2 )
        case -1
            [~,isrt] = sort(RCS_bins_mean(:,2),'descend') ;
        case  1
            [~,isrt] = sort(RCS_bins_mean(:,2),'ascend') ;
    end
    
    
    for j = 1:length(ind_bins)
        j
        isrt(j)
        ind_bins(isrt(j))
        
        CC.ElecName( ind_bins(isrt(j)) ) = j ;
    end
   
end


%%

cmap_bins = distinguishable_colors(max(bins)) ;

%% 

for i = 1:max(bins)
    
    ind_bins = find( bins==i ) ;
    
    RCS_bins_mean = RCS_mean(ind_bins,:) ;
    side = sign( mean(RCS_bins_mean(:,2),1) - (c_CT(end)-c_CT(1))/2 ) ;
    switch side
        case -1
            [~,isrt] = sort(RCS_bins_mean(:,2),'descend') ;
        case  1
            [~,isrt] = sort(RCS_bins_mean(:,2),'ascend') ;
    end    
    % RCS_bins_mean(isrt,:) ;
    
    elec(i).PixelIdxList = CC.PixelIdxList(ind_bins(isrt))  ;
    elec(i).CRS_CT   = RCS_mean(ind_bins(isrt),[2,1,3]) ;
    elec(i).num   = (1:length(ind_bins)).'           ;
    elec(i).side  = side                     ;
        
    tmp1 = CT.tkrvox2ras * [ elec(i).CRS_CT.' ; ones(1,size(elec(i).CRS_CT,1))] ;
    elec(i).tkrRAS = tmp1(1:3,:).' ;

    tmp2 = inv(params.R) * CT.tkrvox2ras * [ elec(i).CRS_CT.' ; ones(1,size(elec(i).CRS_CT,1))] ;
    elec(i).tkrRAS_inT1 = tmp2(1:3,:).' ;

end

%%

dist  = arrayfun(@(x) norm( x.CRS_CT(1,:) .* CT.volres([2,1,3]) ) , elec ) ;
[~,I] = sort(dist) ;
elec  = elec(I) ;

for i = 1:length(elec)
    if isfield(params,'letters')
        elec(i).let = params.letters{i} ;
        elec(i).max = params.num(i)     ;
    else
        elec(i).let = num2roman(i)      ;
        elec(i).max = elec(i).num(end) ;
    end
end

%%

% % % filename = sprintf('scatter3-th-%0.3f-D-%0.3f',th,Dvec) ;
% % % savename = fullfile( dirname , filename ) ;
% % % 
% % % %     h_fig = figure;
% % % %     scatter3(RCS_mean(:,1),RCS_mean(:,2),RCS_mean(:,3),[],cmap_bins(bins,:))
% % % 
% % % h_fig = figure;
% % % h = slice(r_CT,c_CT,s_CT,permute(X,[2,1,3]),round(r_CT(end)/2),round(c_CT(end)/2),round(s_CT(end)/2)) ;
% % % set(h,'EdgeColor','none')
% % % set(h,'FaceAlpha',0.25)
% % % colormap(gray)
% % % hold on
% % % % p = patch('Faces',F,'Vertices',V) ;
% % % % p.FaceColor = 'red' ;
% % % % p.EdgeColor = 'none';
% % % % set(p,'FaceAlpha',0.5)
% % % scatter3(RCS_mean(:,1),RCS_mean(:,2),RCS_mean(:,3),[],cmap_bins(bins,:),'filled')
% % % for i = 1:length(elec)
% % %     
% % %     text(elec(i).CRS_CT(end,2),elec(i).CRS_CT(end,1),elec(i).CRS_CT(end,3), num2roman(i) )
% % %     
% % % end
% % % 
% % % daspect( 1./CT.volres )
% % % 
% % % set(gca,'XDir','reverse')
% % % 
% % % saveas( h_fig , [savename,'.fig'] ,'fig' )


%%

filename = sprintf('%s-scatter3-th-%0.3f-D-%0.3f',CT.name,th,Dvec) ;
if isfield(params,'letters')
    filename = [filename,'-letters'] ;
else
    filename = [filename,'-romans']  ;
end
savename = fullfile( dirname , filename ) ;

%     h_fig = figure;
%     scatter3(RCS_mean(:,1),RCS_mean(:,2),RCS_mean(:,3),[],cmap_bins(bins,:))

h_fig = figure;
h = slice(r_CT,c_CT,s_CT,permute(X,[2,1,3]),round(r_CT(end)/2),round(c_CT(end)/2),round(s_CT(end)/2)) ;
set(h,'EdgeColor','none')
set(h,'FaceAlpha',0.25)
colormap(gray)
hold on
% p = patch('Faces',F,'Vertices',V) ;
% p.FaceColor = 'red' ;
% p.EdgeColor = 'none';
% set(p,'FaceAlpha',0.5)
scatter3(RCS_mean(:,1),RCS_mean(:,2),RCS_mean(:,3),[],cmap_bins(bins,:),'filled')
for i = 1:length(elec)
    
    text(elec(i).CRS_CT(end,2),elec(i).CRS_CT(end,1),elec(i).CRS_CT(end,3), elec(i).let )
    % elec(i)
    
    
    tCRS = elec(i).CRS_CT( 1 ,:) ;
    eCRS = elec(i).CRS_CT(end,:) ;
    
    tRAS = CT.tkrvox2ras * [tCRS(:);1] ;
    tRAS = tRAS(1:3).' ;
    eRAS = CT.tkrvox2ras * [eCRS(:);1] ;
    eRAS = eRAS(1:3).' ;
    
    oRAS = eRAS - tRAS ;
    oRAS = oRAS / norm(oRAS) ;
    
    DD = 3.5 ;
    NN = elec(i).max  ;
    
    RAS = tRAS  + (0:NN-1).' * DD * oRAS ;
    RAS = [RAS.';ones(1,size(RAS,1))]    ;
    CRS_CT = inv( CT.tkrvox2ras ) * RAS     ;
    CRS_CT = CRS_CT(1:3,:).' ;
    
    scatter3(CRS_CT(:,2),CRS_CT(:,1),CRS_CT(:,3),[],'r','+')
           
%     for j = 1:length(ind_bins)
%         % text(RCS_mean(j,1),RCS_mean(j,2),RCS_mean(j,3),sprintf('%d',CC.ElecName(j)))
%         if CC.ElecName(ind_bins(j))==max( CC.ElecName(ind_bins) )
%             
%             text(RCS_mean(ind_bins(j),1),RCS_mean(ind_bins(j),2),RCS_mean(ind_bins(j),3), num2roman( i ) )
%             
%         end
%     end
    
end

daspect( 1./CT.volres )

set(gca,'XDir','reverse')

saveas( h_fig , [savename,'.fig'] ,'fig' )

%%

for i = 1:length(elec)
    
    dd = diff( elec(i).tkrRAS , 1 , 1) ;
    
    nn = nan(size(dd,1),1) ;
    for j = 1:size(dd,1)
        nn(j) = norm(dd(j,:)) ;
    end
    
    nn
end


keyboard



%%


% % ENTRY_CRS = arrayfun(@(x) x.CRS_CT(end,:) , elec ,'UniformOutput',false) ;
% % ENTRY_CRS = cell2mat( ENTRY_CRS(:) ) ;
% % 
% % [~,isrt_C] = sort(ENTRY_CRS(:,1)) ;
% % rk_C = 1:length(elec);
% % rk_C(isrt_C) = rk_C  ;
% % 
% % [~,isrt_R] = sort(ENTRY_CRS(:,2)) ;
% % rk_R = 1:length(elec);
% % rk_R(isrt_R) = rk_R  ;
% % 
% % [~,isrt_S] = sort(ENTRY_CRS(:,3)) ;
% % rk_S = 1:length(elec);
% % rk_S(isrt_S) = rk_S  ;
% % 
% % XX  = cell(length(elec))
% % for i = 1:length(elec)
% %     XX{rk_S(i),rk_R(i)} = elec(i).let ;
% % end

% % XX = cell(length(elec)*ones(1,2)) ;
% % for i = 1:length(elec)
% %     XX{isrt_R(i),isrt_S(i)} = elec(i).let ;
% % end
    
% XX = squeeze( any(XX,2) )


keyboard

%%

filename = sprintf('%s-scatter3-th-%0.3f-D-%0.3f-CRS_CT',CT.name,th,Dvec) ;
if isfield(params,'letters')
    filename = [filename,'-letters'] ;
else
    filename = [filename,'-romans']  ;
end
savename = fullfile( dirname , filename ) ;
fileID = fopen([savename,'.txt'],'w') ;
fprintf(fileID,'# Name\tTarget_C\tTarget_R\tTarget_S\tEntry_C\tEntry_R\tEntry_S\tNum_contacts\n');
for i = 1:length(elec)
    fprintf(fileID,'%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t% 2d\n',...
        elec(i).let,...
        elec(i).CRS_CT( 1 ,1),elec(i).CRS_CT( 1 ,2),elec(i).CRS_CT( 1 ,3),...
        elec(i).CRS_CT(end,1),elec(i).CRS_CT(end,2),elec(i).CRS_CT(end,3),...
        elec(i).max);
end
fclose(fileID);

filename = sprintf('%s-scatter3-th-%0.3f-D-%0.3f-tkrRAS',CT.name,th,Dvec) ;
if isfield(params,'letters')
    filename = [filename,'-letters'] ;
else
    filename = [filename,'-romans']  ;
end
savename = fullfile( dirname , filename ) ;
fileID = fopen([savename,'.txt'],'w') ;
fprintf(fileID,'# Name\tTarget_x\tTarget_y\tTarget_z\tEntry_x\tEntry_y\tEntry_z\tNum_contacts\n');
for i = 1:length(elec)
    fprintf(fileID,'%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t% 2d\n',...
        elec(i).let,...
        elec(i).tkrRAS( 1 ,1),elec(i).tkrRAS( 1 ,2),elec(i).tkrRAS( 1 ,3),...
        elec(i).tkrRAS(end,1),elec(i).tkrRAS(end,2),elec(i).tkrRAS(end,3),...
        elec(i).max);
end
fclose(fileID);

filename = sprintf('%s-scatter3-th-%0.3f-D-%0.3f-tkrRAS_inT1',CT.name,th,Dvec) ;
if isfield(params,'letters')
    filename = [filename,'-letters'] ;
else
    filename = [filename,'-romans']  ;
end
savename = fullfile( dirname , filename ) ;
fileID = fopen([savename,'.txt'],'w') ;
fprintf(fileID,'# Name\tTarget_x\tTarget_y\tTarget_z\tEntry_x\tEntry_y\tEntry_z\tNum_contacts\n');
for i = 1:length(elec)
    fprintf(fileID,'%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t% 2d\n',...
        elec(i).let,...
        elec(i).tkrRAS_inT1( 1 ,1),elec(i).tkrRAS_inT1( 1 ,2),elec(i).tkrRAS_inT1( 1 ,3),...
        elec(i).tkrRAS_inT1(end,1),elec(i).tkrRAS_inT1(end,2),elec(i).tkrRAS_inT1(end,3),...
        elec(i).max);
end
fclose(fileID);


%%

elec_X = elec ;

return

close all

cmap_G = cmap_bins(bins,:) ;

for l = 1:size(X,3)
    
    100*l / size(X,3)
    
    filename = sprintf('z-%03d',l) ;
    savename = fullfile( dirname , filename ) ;
    
    h_fig = figure ;
    
    imshow(X(:,:,l))
    hold on
    
    for ic = 1:CC.NumObjects
        
        % ic/CC.NumObjects
        
        I = false(sz) ;
        I(CC.PixelIdxList{ic}) = true ;
        
        if any(any(I(:,:,l),1),2)
            r_CC = round(mean(R_CT(I))) ;
            c_CC = round(mean(C_CT(I))) ;
            s_CC = round(mean(S_CT(I))) ;
            
            col   = cat(3,repmat(cmap_G(ic,1),sz(1),sz(2)),repmat(cmap_G(ic,2),sz(1),sz(2)),repmat(cmap_G(ic,3),sz(1),sz(2))) ;
            h_col = imshow( col ) ;
            set(h_col,'AlphaData',I(:,:,l))
            str = sprintf('%d,%d,%d',c_CC,r_CC,s_CC) ;
            plot(c_CC,r_CC,'.w')
            h_txt = text(c_CC,r_CC-5,str,'Color','w','FontSize',5) ;
            set(h_txt,'Rotation',90)
        end
        
    end
    
    axis on
    
    daspect([1 1 1])
    title(sprintf('th=%0.3f, D=%0.3f, z=%03d',th,Dvec,l))
    drawnow
    
    
    print( h_fig , savename ,'-djpeg','-r300')
    
    
    close
    
end

keyboard









% %
% % MMMM = 10
% %         DD = D ;
% %         DD(D > MMMM) = NaN ;
% %         DD( D==0 ) = NaN ;
% %
% %         figure
% %         histogram(DD(:),0:0.5:MMMM)
% %
% %         Dmin = 1.5
% %         Dmax = 1.5+2*2.5 ;
% %         ind = Dmin < DD(:) & DD(:) < Dmax ;
% %
% %         figure
% %         histogram(DD(ind),Dmin:0.5:Dmax)


keyboard


elec_X = [];

end


function x_0 = find_grey_peak(X,h)

edges     = 0:h:1;
[N,edges] = histcounts(X(:),edges) ;
x         = edges(1:end-1) + h/2 ;

idx_crop = x > 0.1 ;
x_crop   = x(idx_crop) ;
N_crop   = N(idx_crop) ;

[~,ind] = max(N_crop) ;
x_0     = x_crop(ind) ;

end


function [p0,d] = lineReg(points)

avg = mean(points,1) ;
substracted=bsxfun(@minus,points,avg);
[~,~,V] = svd(substracted);
direction=V(:,1);

p0 = avg      ;
d  = direction;

end

function D = dist2Line(M,p0,s)

D = arrayfun(@(x) norm( cross(M(x,:)-p0,s) ) / norm(s) , 1:size(M,1) ) ;

end

