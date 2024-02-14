
%TrialFilter: respClass_all_array{1,1} = 'REW', 1.2, respClass_all_array{1,2} = 'REW', 0.3
large_reward_active_exclusive = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
large_reward_active_exclusive_sum = sum(large_reward_active_exclusive);
small_reward_active_exclusive = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
small_reward_active_exclusive_sum = sum(small_reward_active_exclusive);


large_reward_inhibited_exclusive = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 3;
large_reward_inhibited_exclusive_sum = sum(large_reward_inhibited_exclusive);
small_reward_inhibited_exclusive = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 2;
small_reward_inhibited_exclusive_sum = sum(small_reward_inhibited_exclusive);

large_reward_neutral_and_small_reward_neutral = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3;
large_reward_neutral_and_small_reward_neutral_sum = sum(large_reward_neutral_and_small_reward_neutral);


%TrialFilter: respClass_all_array{1,3} = 'REW', 1.2, 'BLOCK', 1
large_reward_exclusive_block_1 = large_reward_active_exclusive == 1 & respClass_all_array{1,3} == 1;
large_reward_exclusive_block_1_sum = sum(large_reward_exclusive_block_1);
%TrialFilter: respClass_all_array{1,4} = 'REW', 1.2, 'BLOCK', 2
large_reward_exclusive_block_2 = large_reward_active_exclusive == 1 & respClass_all_array{1,4} == 1;
large_reward_exclusive_block_2_sum = sum(large_reward_exclusive_block_2);
%TrialFilter: respClass_all_array{1,5} = 'REW', 1.2, 'BLOCK', 3
large_reward_exclusive_block_3 = large_reward_active_exclusive == 1 & respClass_all_array{1,5} == 1;
large_reward_exclusive_block_3_sum = sum(large_reward_exclusive_block_3);

large_activated_exclusive_all_blocks = respClass_all_array{1,3} == 1 & respClass_all_array{1,4} == 1 & respClass_all_array{1,5} == 1;
large_activated_exclusive_all_blocks_sum = sum(large_activated_exclusive_all_blocks);