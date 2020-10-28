function mystruct = obj2struct(objArray,fieldNames)
    M=numel(objArray);
    N=numel(fieldNames);
    
    for m=1:M
        for n=1:N
            field=fieldNames{n};
            mystruct(m).(field)=objArray{m}.(field);
        end
    end
    mystruct=reshape(mystruct,size(objArray));
end