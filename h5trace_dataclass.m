classdef h5trace_dataclass < handle
    % Class to act as a pointer to the data so that we don't always have to
    % read the file.
    % Additionally the data can be filtered, to define the filter use the
    % method define_filter(function_handle)
    properties
        filename = [];
        dataset_path = [];
        filter_function = @(x)x; %equivalent to no filtering
    end
    methods
        function obj = h5trace_dataclass(filename,dataset_path)
            obj.filename = filename;
            obj.dataset_path = dataset_path;
        end
        function info = get_h5info(obj)
            info = h5info(obj.filename, obj.dataset_path);
        end
        function out = len(obj)
            info = obj.get_h5info;
            out = info.Dataspace.Size;
        end
        function define_filter(obj,filter_function)
            obj.filter_function = filter_function;
        end
        function argout = subsref(obj, subs)
            if strcmp(subs(1).type, '()') %bypass builtin method
                if isempty(subs.subs) || strcmp(subs(1).subs,':')
                    argout = obj.filter_function(h5read(obj.filename, obj.dataset_path));
                else
                    n = ceil(min(subs.subs{1}));
                    nn = ceil(max(subs.subs{1}));
                    argout = obj.filter_function(h5read(obj.filename, obj.dataset_path,n,nn-n+1));
                end
            else
                builtin('subsref',obj,subs); % use the builtin method...
            end
            
        end
        function ind = end(obj,k,n)
            ind = obj.len();
        end
    end
end