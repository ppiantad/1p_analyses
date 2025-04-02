


non_shk_only_BLA = 100 - percent_shk_only_BLA;

non_shk_only_bla_nacs = 100 - percent_shk_only_bla_nacs;


outer_pie = [percent_shk_only_BLA non_shk_only_BLA];
inner_pie = [percent_shk_only_bla_nacs non_shk_only_bla_nacs];

C = {...
    inner_pie,... % Inner to outer layer
    outer_pie};

% Spider plot
nested_pie(C,'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');

figure; donutchart(outer_pie, 'InnerRadius', 0.5, 'ExplodedWedges', [1 2])

%%


non_consum_only_BLA = 100 - percent_consum_only_BLA;

non_consum_only_bla_nacs = 100 - percent_consum_only_bla_nacs;


outer_pie = [percent_consum_only_BLA non_consum_only_BLA];
inner_pie = [percent_consum_only_bla_nacs non_consum_only_bla_nacs];

C = {...
    inner_pie,... % Inner to outer layer
    outer_pie};

% Spider plot
nested_pie(C,...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');

%%

non_shk_only_BLA = 100 - percent_shk_only_BLA;

non_shk_only_bla_nacs = 100 - percent_shk_only_bla_nacs;


outer_pie = [ non_shk_only_BLA percent_shk_only_BLA];
inner_pie = [ non_shk_only_bla_nacs percent_shk_only_bla_nacs];

C = {...
    outer_pie,... % Inner to outer layer
    inner_pie };

% Spider plot
nested_pie(C,'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');

figure; donutchart(outer_pie, 'InnerRadius', 0.5, 'ExplodedWedges', [1 2])

%%


non_consum_only_BLA = 100 - percent_consum_only_BLA;

non_consum_only_bla_nacs = 100 - percent_consum_only_bla_nacs;


outer_pie = [ non_consum_only_BLA percent_consum_only_BLA];
inner_pie = [ non_consum_only_bla_nacs percent_consum_only_bla_nacs];

C = {...
    outer_pie,... % Inner to outer layer
     inner_pie};

% Spider plot
nested_pie(C,...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');

%%

prechoice_b1_percent_bla = sum(prechoice_block_1_bla)/size(prechoice_block_1_bla, 2)*100;
prechoice_b1_percent_bla_nacs = sum(prechoice_block_1_bla_nacs)/size(prechoice_block_1_bla_nacs, 2)*100;

non_prechoice_b1_percent_bla = 100 - prechoice_b1_percent_bla;

non_prechoice_b1_percent_bla_nacs = 100 - prechoice_b1_percent_bla_nacs;


outer_pie = [ non_prechoice_b1_percent_bla prechoice_b1_percent_bla];
inner_pie = [ non_prechoice_b1_percent_bla_nacs prechoice_b1_percent_bla_nacs];

C = {...
    outer_pie,... % Inner to outer layer
     inner_pie};

% Spider plot
nested_pie(C,...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');



%%

postchoice_b1_percent_bla = sum(postchoice_reward_block_1_bla)/size(postchoice_reward_block_1_bla, 2)*100;
postchoice_b1_percent_bla_nacs = sum(postchoice_reward_block_1_bla_nacs)/size(postchoice_reward_block_1_bla_nacs, 2)*100;

non_postchoice_reward_block_1_bla = 100 - postchoice_b1_percent_bla;

non_postchoice_reward_block_1_bla_nacs = 100 - postchoice_b1_percent_bla_nacs;


outer_pie = [ non_postchoice_reward_block_1_bla postchoice_b1_percent_bla];
inner_pie = [ non_postchoice_reward_block_1_bla_nacs postchoice_b1_percent_bla_nacs];

C = {...
    outer_pie,... % Inner to outer layer
     inner_pie};

% Spider plot
nested_pie(C,...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');
%%

collect_b1_percent_bla = sum(collect_block_1_bla)/size(collect_block_1_bla, 2)*100;
collect_b1_percent_bla_nacs = sum(collect_block_1_bla_nacs)/size(collect_block_1_bla_nacs, 2)*100;

non_collect_b1_percent_bla = 100 - collect_b1_percent_bla;

non_collect_b1_percent_bla_nacs = 100 - collect_b1_percent_bla_nacs;


outer_pie = [ non_collect_b1_percent_bla collect_b1_percent_bla];
inner_pie = [ non_collect_b1_percent_bla_nacs collect_b1_percent_bla_nacs];

C = {...
    outer_pie,... % Inner to outer layer
     inner_pie};

% Spider plot
nested_pie(C,...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

% Title
title('Nested Pie Chart');

