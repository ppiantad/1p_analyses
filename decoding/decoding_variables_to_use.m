% ALL COMPARISONS

neuron_subgroup = {...
    'all';
    'all';
    'all';
    'all';
    'pre-choice';
    'pre-choice';
    'pre-choice';
    'pre-choice';
    'post-choice reward';
    'post-choice reward';
    'post-choice reward';
    'post-choice reward';
    'consumption';
    'consumption';
    'consumption';
    'consumption';
    'true neutral';
    'true neutral';
    'true neutral';
    'true neutral';
    };

shuffle_confirm = [...
    0;
    0;
    0;
    1;
    0;
    0;
    0;
    1;
    0;
    0;
    0;
    1;
    0;
    0;
    0;
    1;
    0;
    0;
    0;
    1;
    ];
ca_period_to_use = [...
    -4 0;
    0 2;
    1 3;
    1 3;
    -4 0;
    0 2;
    1 3;
    1 3;
    -4 0;
    0 2;
    1 3;
    1 3;
    -4 0;
    0 2;
    1 3;
    1 3;
    -4 0;
    0 2;
    1 3;
    1 3;
    ];
current_epoc = {...
    'choiceTime';
    'choiceTime';
    'collectionTime';
    'collectionTime';
    'choiceTime';
    'choiceTime';
    'collectionTime';
    'collectionTime';
    'choiceTime';
    'choiceTime';
    'collectionTime';
    'collectionTime';
    'choiceTime';
    'choiceTime';
    'collectionTime';
    'collectionTime';
    'choiceTime';
    'choiceTime';
    'collectionTime';
    'collectionTime';
    };


%% SUBSET OF COMPARISONS

neuron_subgroup = {...

    'pre-choice';

    'pre-choice';

    'post-choice reward';

    'post-choice reward';

    'consumption';
    'consumption';

    'true neutral';
    'true neutral';
    };

shuffle_confirm = [...

    0;

    1;

    0;

    1;

    0;
    1;
    0;

    1;
    ];
ca_period_to_use = [...

    -4 0;

    1 3;

    0 2;

    1 3;

    1 3;
    1 3;

    1 3;
    1 3;
    ];
current_epoc = {...

    'choiceTime';

    'collectionTime';

    'choiceTime';

    'collectionTime';

    'collectionTime';
    'collectionTime';

    'collectionTime';
    'collectionTime';
    };
