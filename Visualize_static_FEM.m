function Visualize_static_FEM

clear all;
close all;
clc;

%% ReadInput
%% part I
MV1 = load('node.txt')*1.e+3;
MF1 = load('elem.txt')+1;

lenE = length(MF1);
lenN = length(MV1);

N = [3,6];  % the number of columns for each section
H = [5,1];   % the number of header lines (note added in the trailer here
fid=fopen('THL_Output_Low_3.txt');
% node/ stress/ disp/ mass
% node/ stress/ den/ nodal/ mass


for i=1:length(N)
    fmt=repmat('%f',1,N(i));  % build the format string for the number columns
    a(i)=textscan(fid,fmt,'headerlines',H(i),'Delimiter',',','collectoutput',1); % read section
end
fclose('all');

MV = a{1};
Stress = a{2};
% NDPY = a{4};

%% Stress distribution
angle = 35*pi/180;
DPS = zeros(lenN,1);
i1 = zeros(lenN,1);
j2 = zeros(lenN,1);
for i = 1:lenN
    s11 = Stress(i,1);
    s22 = Stress(i,2);
    s33 = Stress(i,3);
    s12 = Stress(i,4);
    s31 = Stress(i,5);
    s23 = Stress(i,6);

    i1(i) = (s11 + s22 + s33)/3;
    j2(i) = (1/6)*((s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2) + (s12^2 + s31^2 + s23^2);
    alpha = 2*sin(angle)/(sqrt(3)*(3-sin(angle)));
    DPS(i) = sqrt(3)*(3-sin(angle))*(alpha*3*i1(i) + sqrt(j2(i)))/(6*cos(angle));
%     if DPS(i) < 0
%         DPS(i) = 0;
%     end
end

%% Density distribution
[F,J,T] = boundary_faces(MF1);
% % MV = MV./1.0e+3;
V = MV;
T = MF1;
[U,G,J,~] = slice_tets(V,T,[0 1 0 1]);
str_slice = sort_data(MF1,J,U,i1);
ystar_slice = sort_data(MF1,J,U,DPS);

figure()
h = colorbar;
patch('Vertices',MV,'Faces',F,'FaceVertexCData',i1/1.e+3,'FaceColor','interp','EdgeColor','none');
axis tight
axis equal
axis on
grid on
lighting none
zlim([floor(min(MV(:,3))) max(MV(:,3))])
xlim([min(MV(:,1)) max(MV(:,1))])
ylabel(h, 'Normal Stress (KPa)');
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
set(gca,'fontsize', 20);
view(0,0)

figure()
h = colorbar;
patch('Vertices',U,'Faces',G,'FaceVertexCData',str_slice/1.e+3,'FaceColor','flat','EdgeColor','none');
axis tight
axis equal
grid on
lighting none
ylabel(h, 'Normal Stress (KPa)');
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
set(gca,'fontsize', 20);
view(0,0)

figure()
h = colorbar;
patch('Vertices',MV,'Faces',F,'FaceVertexCData',DPS/1.e+3,'FaceColor','interp','EdgeColor','none');
axis tight
axis equal
axis on
grid on
lighting none
zlim([floor(min(MV(:,3))) max(MV(:,3))])
xlim([min(MV(:,1)) max(MV(:,1))])
ylabel(h, 'Y^* (KPa)');
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
set(gca,'fontsize', 20);
view(0,0)

figure()
h = colorbar;
patch('Vertices',U,'Faces',G,'FaceVertexCData', ystar_slice/1.e+3,'FaceColor','flat','EdgeColor','none');
axis tight
axis equal
grid on
lighting none
ylabel(h, 'Y^* (KPa)');
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
set(gca,'fontsize', 20);
view(0,0)
end

function data = sort_data(MF,J,U,NDPY)

data = zeros(length(J),1);
for i = 1:length(J)
    order = MF(J(i),:);
    s1 = NDPY(order(1));
    s2 = NDPY(order(2));
    s3 = NDPY(order(3));
    s4 = NDPY(order(4));
    avgV = (s1+s2+s3+s4)/4.;
    data(i) = avgV;
end
end

function [F,J,K] = boundary_faces(T)
  % BOUNDARY_FACES Determine boundary faces of tetrahedra stored in T
  %
  % F = boundary_faces(T)
  %
  % Input:
  %   T  #F by simplex-size list of elements, if simplex-size is 4 (tetrahedra)
  %     then output are boundary faces (triangles), if simplex-size is 3
  %     (trianlges) then this is currently just wrapping `outline`
  % Output:
  %   F  list of boundary faces, n by 3, where n is the number of boundary faces
  %   J  list of indices into T, n by 1
  %   K  list of indices revealing across from which vertex is this face
  %

  if isempty(T)
    F = zeros(0,size(T,2)-1);
    J = [];
    K = [];
    return;
  end

  ss = size(T,2);
  switch ss
  case 2
    % Elements are edges, boundary "facets" are vertices
    F = find(sparse(T,1,1,max(T(:)),1)==1);
    assert(nargout<=1);
  case 3
    % Elements are triangles, boundary "facets" are edges
    F = outline(T);
    [~,I] = ismember(F,[T(:,[2 3]);T(:,[3 1]);T(:,[1 2])],'rows');
    J = mod(I-1,size(T,1))+1;
    K = floor((I-1)/size(T,1))+1;
  case 4
    % get all faces
    allF = [ ...
      T(:,[4 2 3]); ...
      T(:,[3 1 4]); ...
      T(:,[2 4 1]); ...
      T(:,[1 3 2])];
    % sort rows so that faces are reorder in ascending order of indices
    sortedF = sort(allF,2);

    % This is a wild hack but shaves nearly a 2x speed up. Convert the rows to
    % integers before calling unique
    if all(sortedF(:))>=0
      sortedF = uint64(sortedF);
      max_f = max(sortedF(:))+1;
      max_value = intmax('uint64');
      if max_f*max_f*max_f < max_value;
        sortedF = [sortedF(:,1)+sortedF(:,2)*max_f+sortedF(:,3)*max_f^2];
      elseif max_f*max_f < max_value
        sortedF = [sortedF(:,1)+sortedF(:,2)*max_f sortedF(:,3)];
      end
    end

    % determine uniqueness of faces
    [u,m,n] = unique(sortedF,'rows');
    % determine counts for each unique face
    counts = accumarray(n(:), 1);
    % extract faces that only occurred once
    I = m(counts == 1);
    L = I-1;
    J = mod(L,size(T,1))+1;
    K = floor(L/size(T,1))+1;
    F = allF(I,:);
  end
end

function [U,G,J,BC] = slice_tets(V,T,plane,varargin)
  % SLICE_TETS Slice through a tet mesh (V,T) along a given plane (via its
  % implicit equation).
  %
  % [U,G] = slice_tets(V,T,plane)
  % [U,G,J,BC] = slice_tets(V,T,plane,'ParameterName',parameter_value, ...)
  %
  % Inputs:
  %   V  #V by 3 list of tet mesh vertices
  %   T  #T by 4 list of tet indices into V
  %   plane  list of 4 coefficients in the plane equation: [x y z 1]'*plane = 0
  %     or
  %   S  #V by 1 list of values per vertex
  % Outputs:
  %   U  #U by 3 list of triangle mesh vertices along slice
  %   G  #G by 3 list of triangles indices into U
  %   J  #G list of indices into T revealing which tet this face came from
  %   BC  #U by #V list of barycentric coordinates (or more generally: linear
  %     interpolation coordinates) so that U = BC*V
  %
%   Example:
%     % Tet mesh in (V,T)
%     F = boundary_faces(T);
%     % Solve poisson equation for interesting function inside
%     L = cotmatrix(V,T);
%     M = massmatrix(V,T);
%     b = unique(F);
%     int = setdiff(1:size(V,1),b);
%     H = zeros(size(V,1),1);
%     H(int) = (-L(int,int))\(M(int,int)*ones(numel(int),1));
%     clf;
%     t = tsurf(F,V,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.2,'EdgeAlpha',0.2);
%     hold on;
%       s = tsurf([1 1 1],V,'EdgeColor','none',fphong);
%       BB = bounding_box(V(:,1:2));
%       BB = BB([1 2 4 3 1],:);
%       p = plot3(BB(:,1),BB(:,2),min(V(:,3))*[1;1;1;1;1],'-','LineWidth',3);
%     hold off;
%     caxis([min(H) max(H)]);
%     axis equal;
%     for z = linspace(min(V(:,3)),max(V(:,3)))
%       [U,G,J,BC] = slice_tets(V,T,[0 0 1 -z]);
%       set(s,'Vertices',U,'Faces',G,'CData',BC*H);
%       p.ZData = z*[1;1;1;1;1];
%       drawnow;
%     end
  %
  % Example:
  %   % slice a triangle mesh to produce an oriented curve
  %   plane = [0 0 1 -40];
  %   [U,E,J] = slice_tets(V,F(:,[1 2 3 3]),plane);
  %   [U,~,I] = remove_duplicate_vertices(U,eps);
  %   E(E(:,2) == E(:,1),2) = E(E(:,2) == E(:,1),3);
  %   E = E(:,1:2);
  %   N = normals(V,F(J,:));
  %   N = N-sum(N.*plane(1:3),2).*plane(1:3);
  %   M = U(E(:,2),:)-U(E(:,1),:);
  %   Q = [1 plane(1:3)].*[cos(pi/4) sin(pi/4)*[1 1 1]];
  %   W = quatmultiply(quatmultiply(Q,[zeros(size(M,1),1) M]),Q.*[1 -1 -1 -1]);
  %   W = W(:,2:4);
  %   R = sign(sum(W.*N,2))>0;
  %   E(R,:) = fliplr(E(R,:));
  %

  flipped_order = flipped_tet_orders();

  function [U,G] = one_below(V,T,IT)
    [sIT,sJ] = sort(IT,2);
    sT = T(sub2ind(size(T),repmat(1:size(T,1),size(T,2),1)',sJ));
    U = [repmat(sT(:,1),3,1) reshape(sT(:,2:4),size(sT,1)*3,1)];
    G = bsxfun(@plus,1:size(sT,1),[0;1;2]*size(sT,1))';
    flip = ismember(sJ,flipped_order,'rows');
    G(flip,:) = fliplr(G(flip,:));
  end

  function [U,G] = two_below(V,T,IT)
    [sIT,sJ] = sort(IT,2);
    sT = T(sub2ind(size(T),repmat(1:size(T,1),size(T,2),1)',sJ));
    U =  ...
        [repmat(sT(:,1),2,1) reshape(sT(:,3:4),size(sT,1)*2,1); ...
         repmat(sT(:,2),2,1) reshape(sT(:,3:4),size(sT,1)*2,1)];
    G = [ ...
      bsxfun(@plus,1:size(sT,1),[0;1;3]*size(sT,1))'; ...
      bsxfun(@plus,1:size(sT,1),[0;3;2]*size(sT,1))'];
    flip = ismember([sJ;sJ],flipped_order,'rows');
    G(flip,:) = fliplr(G(flip,:));
  end

  % default values
  manifold = true;
  construct_BC = nargout >= 4;

  % Map of parameter names to variable names
  % 'Manifold' is kept for legacy reason, but is no longer needed (always
  % "manifold")
  params_to_variables = containers.Map( ...
    {'Manifold'}, ...
    {'manifold'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  % Homogeneous coordinates
  if numel(plane) == 4
    IV = sum(bsxfun(@times,[V ones(size(V,1),1)],plane),2);
  else
    IV = plane;
  end
  IT = IV(T);
  IT = reshape(IT,size(T));

  I13 = sum(IT<0,2) == 1;
  [U13,G13] = one_below(V,T(I13,:),IT(I13,:));
  I31 = sum(IT>0,2) == 1;
  [U31,G31] = one_below(V,T(I31,:),-IT(I31,:));
  I22 = sum(IT<0,2) == 2;
  [U22,G22] = two_below(V,T(I22,:),IT(I22,:));

  U = [U13;U31;U22];

  % sort edges in U
  sU = sort(U,2);
  % find unique edges in U
  [E,uI,uJ] = unique(sU,'rows');
  % Values at each corner of each unique edge
  IE = IV(E);
  % Sort endpoints of each edge by value
  [sIE,sJ] = sort(IE,2);
  sE = E(sub2ind(size(E),repmat(1:size(E,1),size(E,2),1)',sJ));
  % Lambda from smallest to largest endpoint
  lambda = sIE(:,2)./(sIE(:,2)-sIE(:,1));
  % Vertex position on each unique edge
  U = V(sE(:,1),:).*lambda+ V(sE(:,2),:).*(1-lambda);
  G = [G13;size(U13,1)+[fliplr(G31);size(U31,1)+[G22;]]];
  G = uJ(G);

  if construct_BC
    BC = sparse( ...
      repmat(1:size(sE,1),2,1)', ...
      sE, ...
      [lambda 1-lambda], ...
      size(sE,1),size(V,1));
  end
  J = [find(I13);find(I31);repmat(find(I22),2,1)];

end


function flipped_order = flipped_tet_orders()
  % FLIPPED_TET_ORDERS
  %
  % flipped_order = flipped_tet_orders()
  %
  % Outputs:
  %   flipped_order  20 by 4 list of tet index orders that are flipped inside
  %     out (negative volume)
  %
  flipped_order = [ ...
      4 3 1 2
      4 2 3 1
      4 1 2 3
      3 4 2 1
      3 2 1 4
      3 1 4 2
      2 4 1 3
      2 3 4 1
      2 1 3 4
      1 4 3 2
      1 3 2 4
      1 2 4 3];
end
