function mid_points_range = mid_points(range)

mid_points_range = zeros(1,length(range)-1);
for jj =1:length(range)-1
    mid_points_range(jj) = (range(jj) + range(jj+1)) / 2;
end
end