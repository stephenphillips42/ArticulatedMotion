function value = definedOrDefault_long(name,default,args_in)
% DEFINEDORDEFAULT_LONG takes in arguments and returns default if name is.
% If default is empty, returns error if value with name not provided
% not in the args
    ind = (find(strcmp(args_in,name),1));
    if isempty(ind) && ~isempty(default)
        value = default;
    elseif isempty(ind) && isempty(default)
        error(['Argument ' name ' required'])
    else
        value = args_in{ind+1};
    end
end
