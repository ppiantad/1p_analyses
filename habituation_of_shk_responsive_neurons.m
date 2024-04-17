% first vs last SHK
% first load a dataset with the following data (created by
% eventRelatedActivity):

% assuming input data are:
%pre-choice 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1 choiceTime -4 to 0
%post-choice 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1 choiceTime 0 to 2
%consumption 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1 collectionTime 1 to 3
%shk 'SHK', 1 choiceTime 0 to 2


array_to_examine = 1;

%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        first_shk_data(zz,:) = zall_array{array_to_examine, qq}(1,:);
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end

%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        last_shk_data(zz,:) = zall_array{array_to_examine, qq}(end,:);
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end

%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        middle_shk_data(zz,:) = zall_array{array_to_examine, qq}(round(end/2),:);
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end


first_shk_data_mean = mean(first_shk_data);
first_shk_data_sem = nansem(first_shk_data);

middle_shk_data_mean = mean(middle_shk_data);
middle_shk_data_sem = nansem(middle_shk_data);

final_shk_data_mean = mean(last_shk_data);
final_shk_data_sem = nansem(last_shk_data);
%%

figure;
shadedErrorBar(ts1, first_shk_data_mean, first_shk_data_sem, 'lineProps', {'color', batlowW(iter+50,:), 'LineWidth', 2});
hold on; shadedErrorBar(ts1, middle_shk_data_mean, middle_shk_data_sem, 'lineProps', {'color', batlowW(iter+75,:), 'LineWidth', 2});
hold on; shadedErrorBar(ts1, final_shk_data_mean, final_shk_data_sem, 'lineProps', {'color', batlowW(iter+100,:), 'LineWidth', 2});
legend({'first shock', 'middle shock', 'final shock'}, 'Location','northwest');

%% for analyzing SHOCK TEST DATA
% first run eventRelatedActivity with SHK, 1 as the filter! 
array_to_examine = 1;

% if using shock test, uncomment the below line:
shk_activated = respClass_all_array{1, 1} == 1;

% use this array if you want to use all neurons
% shk_activated = ones(1, size(respClass_all_array{1, 1}, 2));


%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        first_five_shk_data(zz,:) = mean(zall_array{array_to_examine, qq}(2:6,:));
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end

%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        second_five_shk_data(zz,:) = mean(zall_array{array_to_examine, qq}(7:11,:));
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end

%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        third_five_shk_data(zz,:) = mean(zall_array{array_to_examine, qq}(12:16,:));
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end


%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        fourth_five_shk_data(zz,:) = mean(zall_array{array_to_examine, qq}(17:21,:));
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end

%create a dataset where only shk responsive cells are stored for each mouse
zz=0;
for qq = 1:size(zall_array, 2)
    if shk_activated(qq) == 1
        zz = zz+1;
        fifth_five_shk_data(zz,:) = mean(zall_array{array_to_examine, qq}(22:26,:));
        corresponding_mouse(zz) = mouse_cells(array_to_examine, qq);
    end
end

first_five_shk_data_mean = mean(first_five_shk_data);
first_five_shk_data_sem = nansem(first_five_shk_data);

second_five_shk_data_mean = mean(second_five_shk_data);
second_five_shk_data_sem = nansem(second_five_shk_data);

third_five_shk_data_mean = mean(third_five_shk_data);
third_five_shk_data_sem = nansem(third_five_shk_data);

fourth_five_shk_data_mean = mean(fourth_five_shk_data);
fourth_five_shk_data_sem = nansem(fourth_five_shk_data);

fifth_five_shk_data_mean = mean(fifth_five_shk_data);
fifth_five_shk_data_sem = nansem(fifth_five_shk_data);

%%

figure;
shadedErrorBar(ts1, first_five_shk_data_mean, first_five_shk_data_sem, 'lineProps', {'color', batlowW(iter+50,:), 'LineWidth', 2});
hold on; shadedErrorBar(ts1, second_five_shk_data_mean, second_five_shk_data_sem, 'lineProps', {'color', batlowW(iter+75,:), 'LineWidth', 2});
hold on; shadedErrorBar(ts1, third_five_shk_data_mean, third_five_shk_data_sem, 'lineProps', {'color', batlowW(iter+100,:), 'LineWidth', 2});
hold on; shadedErrorBar(ts1, fourth_five_shk_data_mean, fourth_five_shk_data_sem, 'lineProps', {'color', batlowW(iter+125,:), 'LineWidth', 2});
hold on; shadedErrorBar(ts1, fifth_five_shk_data_mean, fifth_five_shk_data_sem, 'lineProps', {'color', batlowW(iter+150,:), 'LineWidth', 2});
legend({'0.02 to 0.10 mA', '0.12 to 0.20 mA', '0.21 to 0.30 mA', '0.31 to 0.40 mA', '0.41 to 0.50 mA'}, 'Location','northwest');


