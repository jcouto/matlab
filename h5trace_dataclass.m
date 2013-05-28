classdef h5trace_dataclass < handle
    properties
        filename = [];
        dataset_path = [];
    methods
        function varargout = data_class(obj, filename,dataset_path)
            obj.filename = filename;
            obj.dataset_path = dataset_path;
        end
    end
end