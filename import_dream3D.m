
function ebsd = import_dream3D(fname, varargin)

    data = h5load(fname);
    d = data.DataContainers.ImageDataContainer.CellData;      

    % Import rotations
    EulerAngles = d.EulerAngles;
    phi1 = permute(EulerAngles(:,:,:,1),[3 2 1]);
    Phi  = permute(EulerAngles(:,:,:,2),[3 2 1]);
    phi2 = permute(EulerAngles(:,:,:,3),[3 2 1]);
    rot = rotation.byEuler(phi1(:),Phi(:),phi2(:));
    
    % Import phases
    Phases = permute(d.Phases,[3 2 1]);
    phaseId = Phases(:);
    
    % Import Image Quality
    ImageQuality = permute(d.ImageQuality,[3 2 1]);
    iq = ImageQuality(:);

    % % Import Mask Property
    % maskProp = permute(d.Mask,[3 2 1]);
    % mask = maskProp(:);


    % Create EBSD object in MTEX
    
    %Define crystal symmetry
    CS = crystalSymmetry('m-3m','mineral','\gamma phase');

    CSList = {'notIndexed', CS};
    % CSList = {CS};

    %Create manual xyz data
    [x_size, y_size, z_size] = size(squeeze(Phases));
    [x, y, z] = meshgrid(1:y_size,1:x_size,1:z_size);

    pos = vector3d(x(:),y(:),z(:));
    
    prop = struct;
    prop.iq = iq;
    % prop.mask = mask;
    prop.x = x;
    prop.y = y;
    prop.z = z;
    unitCell = 1/2 * vector3d([1 1 -1 -1 1 1 -1 -1],...
                                       [1 -1 -1 1 1 -1 -1 1],...
                                       [1 1 1 1 -1 -1 -1 -1]);
    ebsd = EBSD3(pos,rot,phaseId,CSList,prop,unitCell);
    ebsd.unitCell = unitCell;

%     ebsd = EBSD(pos,rot,phaseId,CSList,prop);
%     ebsd.unitCell = unitCell;
end

function data=h5load(filename, path)
    %
    % data = H5LOAD(filename)
    % data = H5LOAD(filename, path_in_file)
    %
    % Load data in a HDF5 file to a Matlab structure.
    %
    % Parameters
    % ----------
    %
    % filename
    %     Name of the file to load data from
    % path_in_file : optional
    %     Path to the part of the HDF5 file to load
    %
    
    % Author: Pauli Virtanen <pav@iki.fi>
    % This script is in the Public Domain. No warranty.
    
    if nargin > 1
      path_parts = regexp(path, '/', 'split');
    else
      path = '';
      path_parts = [];
    end
    
    loc = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    try
      data = load_one(loc, path_parts, path);
      H5F.close(loc);
    catch exc
      H5F.close(loc);
      rethrow(exc);
    end
end

function data=load_one(loc, path_parts, full_path)
    % Load a record recursively.
    
    while ~isempty(path_parts) & strcmp(path_parts{1}, '')
      path_parts = path_parts(2:end);
    end
    
    data = struct();
    
    num_objs = H5G.get_num_objs(loc);
    
    % 
    % Load groups and datasets
    %
    for j_item=0:num_objs-1,
      objtype = H5G.get_objtype_by_idx(loc, j_item);
      objname = H5G.get_objname_by_idx(loc, j_item);
      
      if objtype == 0
        % Group
        name = regexprep(objname, '.*/', '');
      
        if isempty(path_parts) | strcmp(path_parts{1}, name)
          if ~isempty(regexp(name,'^[a-zA-Z].*'))
	    group_loc = H5G.open(loc, name);
	    try
	      sub_data = load_one(group_loc, path_parts(2:end), full_path);
	      H5G.close(group_loc);
	    catch exc
	      H5G.close(group_loc);
	      rethrow(exc);
	    end
	    if isempty(path_parts)
	      data = setfield(data, erase(name," "), sub_data);
	    else
	      data = sub_data;
	      return
	    end
          end
        end
       
      elseif objtype == 1
        % Dataset
        name = regexprep(objname, '.*/', '');
      
        if isempty(path_parts) | strcmp(path_parts{1}, name)
          if ~isempty(regexp(name,'^[a-zA-Z].*'))
	    dataset_loc = H5D.open(loc, name);
	    try
	      sub_data = H5D.read(dataset_loc, ...
	          'H5ML_DEFAULT', 'H5S_ALL','H5S_ALL','H5P_DEFAULT');
	      H5D.close(dataset_loc);
	    catch exc
	      H5D.close(dataset_loc);
	      rethrow(exc);
	    end
	    
	    sub_data = fix_data(sub_data);
	    
	    if isempty(path_parts)
	      data = setfield(data, erase(name," "), sub_data);
	    else
	      data = sub_data;
	      return
	    end
          end
        end
      end
    end
    
    % Check that we managed to load something if path walking is in progress
    if ~isempty(path_parts)
      error(sprintf('Path "%s" not found in the HDF5 file', full_path));
    end
end

function data=fix_data(data)
    % Fix some common types of data to more friendly form.
    
    if isstruct(data)
      fields = fieldnames(data);
      if length(fields) == 2 & strcmp(fields{1}, 'r') & strcmp(fields{2}, 'i')
        if isnumeric(data.r) & isnumeric(data.i)
          data = data.r + 1j*data.i;
        end
      end
    end
    
    if isnumeric(data) & ndims(data) > 1
      % permute dimensions
      data = permute(data, fliplr(1:ndims(data)));
    end
end