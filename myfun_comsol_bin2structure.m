% Convert 2D binary image structure to microstructure parameters for FEA input. 

% Assuming fillers are ellipses and domain window is square. 

% Input: 
% Imagefile: MAT file with matrix 'img_out' that contains binary image. Black (0) is matrix
% ActualLength: side length [pixel]
% dimension to pixel: ratio of physical size to pix
% vf expt: labeld VF

% Output: 
% Image file saved in the same folder containing 
% 1. img_para [um]. 
% 3. dimensionX [um], side length of the simulation domain 

function [] = myfun_comsol_bin2structure(imagefile, ActualLength, dimension_to_pixel,  vf_expt)

% =============================
isReScale = 1;

% Add buffer around sides of the square, so that fillers with added interphase do not overlap with boundary 
CutSide = 0.1; % fraction of side length that is cut from the ends 
RemainSide = 1-2*CutSide;

% =============================

% Obtain dimension to pixel convertion and simulation box side length.
disp(['Ratio of physical length to pixel: ',num2str(dimension_to_pixel),'nm-per-pixel'])
dimensionX = ActualLength*dimension_to_pixel*1e-3; % unit in um
dimensionY = dimensionX;

load(imagefile)
image = double(img_out);
% Label clusters in image
img = bwlabel(image);
ClusterNo = max(max(img));
Center_list = [];   % List all the cluster centers coordinates in a N*2 matrix
MajorAxis =[];
MinorAxis = [];
Angle = zeros(ClusterNo, 1);
cc = regionprops(img,'Centroid'); % cc is the cluster center
cra= regionprops(img,'MajorAxisLength');
crb= regionprops(img,'MinorAxisLength');

for ii = 1:1:ClusterNo
    Center_list = [Center_list; cc(ii).Centroid];
    MajorAxis = [MajorAxis; cra(ii).MajorAxisLength];
    MinorAxis = [MinorAxis; crb(ii).MinorAxisLength];
    [cy,cx] = find(img==ii); % pixel coordinates in each cluster
    PixelNo = length(cy); % number of pixels in current cluster
    % Find orientation angle of ellipse that encloses each cluster
    dist = 0;
    for j1=1:PixelNo
        for j2=1:PixelNo
            distnew = (cy(j1)-cy(j2))^2+(cx(j1)-cx(j2))^2;
            if distnew > dist
                dist = distnew;
                c1y = cy(j1); c2y = cy(j2);
                c1x = cx(j1); c2x = cx(j2);
            end
        end
    end
    Angle(ii) = atan((c1y-c2y)/(c1x-c2x));
end
MajorAxis = MajorAxis/2;
MinorAxis = MinorAxis/2;

% Convert to white - matrix, black - filler
I_bright = zeros(length(image), length(image));
for i = 1:length(image)
    for j = 1:length(image)
        if image(i,j) > 0
            I_bright(i,j) = 0;
        else
            I_bright(i,j) = 1;
        end
    end
end


% Return original X and Y coordinates of filler centers
x0 = Center_list(:,1);
y0 = Center_list(:,2);

disp(['Original number of clusters: ',num2str(ClusterNo)])

% Seed ellipses at central 90% by 90% area
scaling = RemainSide*dimensionX/ActualLength;
LongAxis = scaling*MajorAxis;
ShortAxis = scaling*MinorAxis;
PosX = CutSide*dimensionX + scaling*x0;
PosY = CutSide*dimensionY + scaling*y0;

NewClusterNo    = ClusterNo;
NewLongAxis     = LongAxis(1:ClusterNo);
NewShortAxis    = ShortAxis(1:ClusterNo);
NewPosX         = PosX(1:ClusterNo);
NewPosY         = PosY(1:ClusterNo);
NewAngle        = Angle(1:ClusterNo);

EllipseMatrix=zeros(NewClusterNo,5);
EllipseMatrix(:,1)=NewAngle*180/3.1416; % Angles in degree in COMSOL build
EllipseMatrix(:,2)=NewLongAxis;
EllipseMatrix(:,3)=NewShortAxis;
EllipseMatrix(:,4)=NewPosX;
EllipseMatrix(:,5)=NewPosY;

disp(['Number of clusters in FEA geometry: ',num2str(NewClusterNo)])
ActualVF = 3.1416*EllipseMatrix(:,2)'*EllipseMatrix(:,3)/(dimensionX*dimensionY);
disp(['Actual VF in simulation window: ',num2str(ActualVF)])

% Correct long/short axes to match apparent VF with labeled
if isReScale == 1
    ReScale = sqrt(vf_expt/ActualVF);
    EllipseMatrix(:,2) = EllipseMatrix(:,2)*ReScale;
    EllipseMatrix(:,3) = EllipseMatrix(:,3)*ReScale;
end

CorrectedVF = 3.1416*EllipseMatrix(:,2)'*EllipseMatrix(:,3)/(dimensionX*dimensionY);
disp(['Corrected VF in simulation window: ',num2str(CorrectedVF)])

img_para = EllipseMatrix;
save([imagefile,'_2D_structure_output'],'img_para', 'dimensionX')
end