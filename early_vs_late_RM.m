neuron_num = size(RM_D1_vs_Pre_RDT_RM_ensemble_logical_indices, 2);

RM_early_indices = RM_D1_vs_Pre_RDT_RM_ensemble_logical_indices;
RM_late_indices = RM_D1_vs_Pre_RDT_RM_Pre_RDT_RM_ensemble_logical_indices;

early_RM_true_neutral = all(RM_early_indices == 0, 1);
late_RM_true_neutral = all(RM_late_indices == 0, 1);

early_RM_responsive = any(RM_early_indices == 1, 1);
late_RM_responsive = any(RM_late_indices == 1, 1);




early_RM_active = any(RM_early_indices(1:3, :) == 1, 1);
sum(early_RM_active)


pre_choice_active_stable = RM_early_indices(1, :) == 1 & RM_late_indices(1, :) == 1;
sum(pre_choice_active_stable)
post_choice_active_stable = RM_early_indices(2, :) == 1 & RM_late_indices(2, :) == 1;
sum(post_choice_active_stable)
consumption_active_stable = RM_early_indices(3, :) == 1 & RM_late_indices(3, :) == 1;
sum(consumption_active_stable)

pre_choice_inhibited_stable = RM_early_indices(4, :) == 1 & RM_late_indices(4, :) == 1;
sum(pre_choice_inhibited_stable)
post_choice_inhibited_stable = RM_early_indices(5, :) == 1 & RM_late_indices(5, :) == 1;
sum(post_choice_inhibited_stable)
consumption_inhibited_stable = RM_early_indices(6, :) == 1 & RM_late_indices(6, :) == 1;
sum(consumption_inhibited_stable)

neutral_early_to_pre_choice_late = early_RM_true_neutral == 1 & (RM_late_indices(1, :) == 1 | RM_late_indices(4, :) == 1);
sum(neutral_early_to_pre_choice_late)

neutral_early_to_post_choice_late = early_RM_true_neutral == 1 & (RM_late_indices(2, :) == 1 | RM_late_indices(5, :) == 1);
sum(neutral_early_to_post_choice_late)

neutral_early_to_consumption_late = early_RM_true_neutral == 1 & (RM_late_indices(3, :) == 1 | RM_late_indices(6, :) == 1);
sum(neutral_early_to_consumption_late)

pre_choice_early_to_post_choice_late = RM_early_indices(1, :) == 1 & (RM_late_indices(2, :) == 1 | RM_late_indices(5, :) == 1);
sum(pre_choice_early_to_post_choice_late)

pre_choice_early_to_consumption_late = RM_early_indices(1, :) == 1 & (RM_late_indices(3, :) == 1 | RM_late_indices(6, :) == 1);
sum(pre_choice_early_to_consumption_late)

post_choice_early_to_pre_choice_late = RM_early_indices(2, :) == 1 & (RM_late_indices(1, :) == 1 | RM_late_indices(4, :) == 1);
sum(post_choice_early_to_pre_choice_late)

post_choice_early_to_consumption_late = RM_early_indices(2, :) == 1 & (RM_late_indices(3, :) == 1 | RM_late_indices(6, :) == 1);
sum(post_choice_early_to_consumption_late)

consumption_early_to_pre_choice_late = RM_early_indices(3, :) == 1 & (RM_late_indices(1, :) == 1 | RM_late_indices(4, :) == 1);
sum(consumption_early_to_pre_choice_late)

consumption_early_to_post_choice_late = RM_early_indices(3, :) == 1 & (RM_late_indices(2, :) == 1 | RM_late_indices(5, :) == 1);
sum(consumption_early_to_post_choice_late)


%%
neutral_early_to_responsive_late = early_RM_true_neutral == 1 & late_RM_responsive == 1;
sum(neutral_early_to_responsive_late)

stable_responsive = early_RM_responsive == 1 & late_RM_responsive == 1;
sum(stable_responsive)

early_active = any(RM_early_indices(1:3, :) == 1, 1);
sum(early_active )

late_active = any(RM_late_indices(1:3, :) == 1, 1);
sum(late_active)

early_inhibited = any(RM_early_indices(4:6, :) == 1, 1);
sum(early_inhibited)

late_inhibited = any(RM_late_indices(4:6, :) == 1, 1);
sum(late_inhibited)


early_rm_responsive_to_late_rm_neutral = early_RM_responsive == 1 & late_RM_true_neutral == 1;
sum(early_rm_responsive_to_late_rm_neutral)

stable_active = early_active == 1 & late_active == 1; 
sum(stable_active)

stable_inhibited = early_inhibited == 1 & late_inhibited == 1; 
sum(stable_inhibited)

neutral_early_to_neutral_late = early_RM_true_neutral == 1 & late_RM_true_neutral == 1;
sum(neutral_early_to_neutral_late)

inner_pie = [sum(neutral_early_to_responsive_late)/neuron_num,...
            
            sum(early_rm_responsive_to_late_rm_neutral)/neuron_num,...
            
            sum(stable_responsive)/neuron_num,...
           
            sum(neutral_early_to_neutral_late)/neuron_num];


labels = ["neutral to responsive", "responsive to neutral", "stable", "neutral to neutral"];

figure; donutchart(inner_pie, labels, 'InnerRadius', 0.5)

%%

early_active = any(RM_early_indices(1:3, :) == 1, 1);
sum(early_active )

late_active = any(RM_late_indices(1:3, :) == 1, 1);
sum(late_active)

stable_active = early_active == 1 & late_active == 1; 
sum(stable_active)

stable_active_pre_choice = RM_early_indices(1, :) == 1 & RM_late_indices(1, :) == 1;
sum(stable_active_pre_choice)

stable_active_post_choice = RM_early_indices(2, :) == 1 & RM_late_indices(2, :) == 1;
sum(stable_active_post_choice)

stable_active_consumption = RM_early_indices(3, :) == 1 & RM_late_indices(3, :) == 1;
sum(stable_active_consumption)

stable_active_within_category = sum(stable_active_pre_choice) + sum(stable_active_post_choice) + sum(stable_active_consumption);
sum(stable_active_within_category)

neutral_early_to_neutral_late = early_RM_true_neutral == 1 & late_RM_true_neutral == 1;
sum(neutral_early_to_neutral_late)

neutral_early_to_active_late = early_RM_true_neutral == 1 & late_active == 1;
sum(neutral_early_to_active_late)

early_active_to_neutral_late = early_active == 1 & late_RM_true_neutral == 1;
sum(early_active_to_neutral_late)

inhibited_early_to_active_late = early_inhibited == 1 & late_active == 1; 
sum(inhibited_early_to_active_late)


early_inhibited_pre_choice_to_late_inhibited_pre_choice = RM_early_indices(4, :) == 1 & RM_late_indices(4, :) == 1;
sum(early_inhibited_pre_choice_to_late_inhibited_pre_choice)

early_active_pre_choice_to_late_inhibited_pre_choice = RM_early_indices(1, :) == 1 & RM_late_indices(4, :) == 1;
sum(early_active_pre_choice_to_late_inhibited_pre_choice)

early_active_pre_choice_to_late_inhibited_post_choice = RM_early_indices(1, :) == 1 & RM_late_indices(5, :) == 1;
sum(early_active_pre_choice_to_late_inhibited_post_choice)

early_active_pre_choice_to_late_inhibited_consumption = RM_early_indices(1, :) == 1 & RM_late_indices(6, :) == 1;
sum(early_active_pre_choice_to_late_inhibited_consumption)




inner_pie = [sum(neutral_early_to_active_late)/sum(late_RM_active),...
                                    
            sum(stable_active_within_category)/sum(late_RM_active),...

            sum(inhibited_early_to_active_late)/sum(late_RM_active)];


labels = ["neutral to active", "stably active",  "inhibited early to active late"];

figure; donutchart(inner_pie, labels, 'InnerRadius', 0.5)


%%
late_RM_active = any(RM_late_indices(1:3, :) == 1, 1);
sum(late_RM_active)



