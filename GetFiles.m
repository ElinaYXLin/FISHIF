%% baseline (get name from file)
function [names,class_num] = GetFiles(files)

size0 = size(files);
length = size0(1);
names = files(3:length);
class_num = size(names);
end