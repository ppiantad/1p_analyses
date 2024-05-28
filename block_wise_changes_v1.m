% Pre_RDT_RM & RDT D1 matched
% RDT D1
% collectionTime
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 1
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3

collect_all_possible = sum([sum(respClass_all_array{1, 1} ==1), sum(respClass_all_array{1, 2} ==1)])
collect_conserved = respClass_all_array{1, 1}==1 & respClass_all_array{1, 2}==1;
collect_conserved_sum = sum(collect_conserved)
collect_original_sum = sum(respClass_all_array{1, 1}==1)
collect_remapped = respClass_all_array{1, 1} ~=1 & respClass_all_array{1, 2}==1;
collect_remapped_sum = sum(collect_remapped)
not_collect_sum = sum(respClass_all_array{1, 1} ~=1)

collect_lost = respClass_all_array{1, 1} ==1 & respClass_all_array{1, 2} ~=1;
collect_lost_sum = sum(collect_lost)
collect_all_possible = sum([sum(respClass_all_array{1, 1} ==1), sum(respClass_all_array{1, 2} ==1)])
collect_conserved_percent = collect_conserved_sum/neuron_num
collect_remapped_percent = collect_remapped_sum/neuron_num
collect_lost_percent = collect_lost_sum/neuron_num

inner_pie = [collect_conserved_sum/neuron_num,...
            
            collect_remapped_sum/neuron_num,...
            
            collect_lost_sum/neuron_num,...
           
            (neuron_num- [collect_conserved_sum+collect_remapped_sum+collect_lost_sum]) / neuron_num];


labels = ["conserved", "re-mapped", "lost", "not sig responsive"];

figure; donutchart(inner_pie, labels, 'InnerRadius', 0.5)