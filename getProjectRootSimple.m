function root = getProjectRootSimple()
% Returns <root> assuming this file lives in <root>/code/
fp = mfilename('fullpath');      % .../code/getProjectRootSimple.m
codeDir = fileparts(fp);         % .../code
[root, folder] = fileparts(codeDir);
assert(strcmpi(folder, 'code'), 'getProjectRootSimple.m must be in <root>/code/');
end
