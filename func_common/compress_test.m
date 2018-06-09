fields = fieldnames(data);
data_com = struct;
for i = 3:length(fields)
    field = fields(i);
    data_com.(field{1}) = compress_matrix(data.(field{1}),10);
    
    figure(i)
    plot(data.time,data.(field{1}),'.','markersize',2); hold on
    plot(data_com.time, data_com.(field{1}),'.','markersize',5); 
    title(field{1});
end
