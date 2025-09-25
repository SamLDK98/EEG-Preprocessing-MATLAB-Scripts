function root = getProjectRootSimple()
% Returns <root> assuming this file lives in <root>/code/
fp = mfilename('fullpath');     
codeDir = fileparts(fp);         
[root, folder] = fileparts(codeDir);
assert(strcmpi(folder, 'code'), 'getProjectRootSimple.m must be in <root>/code/');
end
