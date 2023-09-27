function index = ArrayValueToIndex(value, array)

size_l = numel(array);

index = 1;
 
while( index < size_l && (array(index) > value || value >= array(index+1)))
    index = index + 1;
end

if array(index) == value
    bool = 1;
end

if (index > size_l && ~bool)
    error('no matching value');
end