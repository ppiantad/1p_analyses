data = [excited_to_excited_sum excited_to_inhibited_sum excited_to_neutral_sum]

pie(data)
figure; pie(data); labels = {'large excited & small excited: ' + string(excited_to_excited_sum(1)), 'large excited + small inhibited: ' + string(excited_to_inhibited_sum(1)), 'large excited + small neutral: ' + string(excited_to_neutral_sum(1))};
legend(labels)

%%
data = [inhibited_to_excited_sum inhibited_to_inhibited_sum inhibited_to_neutral_sum]
figure; pie(data); labels = {'large inhibited & small excited: ' + string(inhibited_to_excited_sum(1)), 'large inhibited + small inhibited: ' + string(inhibited_to_inhibited_sum(1)), 'large inhibited + small neutral: ' + string(inhibited_to_neutral_sum(1))};
legend(labels)

%%
data = [neutral_to_excited_sum neutral_to_inhibited_sum neutral_to_neutral_sum]
figure; pie(data); labels = {'large neutral & small excited: ' + string(neutral_to_excited_sum(1)), 'large neutral + small inhibited: ' + string(neutral_to_inhibited_sum(1)), 'large neutral + small neutral: ' + string(neutral_to_neutral_sum(1))};
legend(labels)








% Use below if SMALL was the first event filtered on, and LARGE was the
% second
%%
data = [excited_to_excited_sum excited_to_inhibited_sum excited_to_neutral_sum]

pie(data)
figure; pie(data); labels = {'small excited & large excited: ' + string(excited_to_excited_sum(1)), 'small excited + large inhibited: ' + string(excited_to_inhibited_sum(1)), 'small excited + large neutral: ' + string(excited_to_neutral_sum(1))};
legend(labels)

%%
data = [inhibited_to_excited_sum inhibited_to_inhibited_sum inhibited_to_neutral_sum]
figure; pie(data); labels = {'small inhibited & large excited: ' + string(inhibited_to_excited_sum(1)), 'small inhibited + large inhibited: ' + string(inhibited_to_inhibited_sum(1)), 'small inhibited + large neutral: ' + string(inhibited_to_neutral_sum(1))};
legend(labels)

%%
data = [neutral_to_excited_sum neutral_to_inhibited_sum neutral_to_neutral_sum]
figure; pie(data); labels = {'small neutral & large excited: ' + string(neutral_to_excited_sum(1)), 'small neutral + large inhibited: ' + string(neutral_to_inhibited_sum(1)), 'small neutral + large neutral: ' + string(neutral_to_neutral_sum(1))};
legend(labels)