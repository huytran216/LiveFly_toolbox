feature_description = {
    'ON/OFF', ... 1
    'Relative First Activation time', ... 2
    'Relative Last Activation time', ... 3
    'Relative Total active duration', ... 4
    'Relative Total spot duration of existence ', ... 5
    'Relative Time to max intensity', ... 6
    'Brightest spot intensity', ... 7
    'Mean Intensity produced (absolute time) (including non-expressing nuclei and non-expressing frame)', ... 8
    'Total Intensity produced (including non-expressing nuclei and non-expressing frame)', ... 9
    'Total interphase duration', ... 10
    'Absolute First Activation time', ... 11
    'Absolute Last Activation time', ... 12
    'Absolute Total active duration', ... 13
    'Absolute Total spot duration of existence', ... 14
    'Absolute Time to max intensity', ... 15
    'Mean Spot Intensity (including non-expressing nuclei and non-expressing frame)', ... 16
    'Mean Spot Intensity (only expressing nuclei (full traces) including non-expressing frames)', ... 17
    'Mean Spot Intensity (only expressing nuclei (trimmed traces) including non-expressing frames)', ... 18
    'Mean Spot Intensity (only expressing nuclei and expressing frames)', ... 19
    'Mean Maximum Spot Intensity (only expressing nuclei and expressing frames)', ... 20
    'Mean PSpot (including non-expressing nuclei)', ... 21
    'Mean PSpot (only expressing nuclei (full traces))', ... 22
    'Mean PSpot (only expressing nuclei (trimmed traces))', ... 23
    'Mean Spot Intensity (after 1st spot appearance, 1st burst only, with a time limit)', ... 24
    'Mean Pspot (after 1st spot appearance, 1st burst only, with a time limit)', ... 25
    'Mean Spot Intensity (after 1st spot appearance, within a time window)', ... 26
    'Mean Pspot (after 1st spot appearance, within a time window)', ... 27
    'Number of burst (whole trace)', ... 28
    'Number of burst (within a window)' ... 29
};

feature_label={'ON',... 1
    't_{init}%',...2
    't_{end}%',... 3
    't_{active}%',...4
    't_{exist}%',...5
    't_{max}%',...6
    'max Ispot',...7
    '\mu{I}_{active}',...8
    '\Sigma{I}',...9
    'tphase',...10
    't_{init}',...11
    't_{end}',...12
    't_{active}',...13
    't_{exist}',...14
    't_{max}',...15
    '\mu{I}',...16
    '\mu{I}_{ONfull}',...17
    '\mu{I}_{ONtrim}',...18
    '\mu{I}_{spot}',...19
    '\mu{I}_{spot,max}',...20
    'P_{Spot}',...21
    'P_{Spot,ONfull}',...22
    'P_{Spot,ONtrim}',...23
    '\mu_{I}_{burst}',...24
    'P_{Spot}_{burst}',...25
    '\mu_{I}_{ONwindow}',...26
    'P_{Spot}_{ONwindow}',...27
    'N_{burst}',...28
    'N_{burst}_{window}',...29
    };

feature_distribution={'binomial',... 1
    'gamma',...2
    'gamma',... 3
    'gamma',...4
    'gamma',...5
    'gamma',...6
    'gamma',...7
    'gamma',...8
    'gamma',...9
    'gamma',...10
    'gamma',...11
    'gamma',...12
    'gamma',...13
    'gamma',...14
    'gamma',...15
    'gamma',...16
    'gamma',...17
    'gamma',...18
    'gamma',...19
    'gamma',...20
    'beta',...21
    'beta',...22
    'beta',...23
    'gamma',...24
    'beta',...25
    'gamma',...26
    'beta',...27
    'poisson',...28
    'poisson',...29
    };

feature_Tsensitive=[0,... 1
    0,...2
    0,... 3
    0,...4
    0,...5
    0,...6
    0,...7
    0,...8
    1,...9
    1,...10
    1,...11
    1,...12
    1,...13
    1,...14
    0,...15
    0,...16
    0,...17
    0,...18
    0,...19
    1,...20
    0,...21
    0,...22
    0,...23
    0,...24
    0,...25
    0,...26
    0,...27
    1,...28
    1,...29
    ];

Nfea=numel(feature_label);

feature_info = struct;
for i=1:Nfea
    feature_info(i).label = feature_label{i};
    feature_info(i).description = feature_description{i};
    feature_info(i).distribution = feature_distribution{i};
    feature_info(i).Tsensitive = feature_Tsensitive(i);
end
save('feature_label_moreinfo.mat','feature_Tsensitive','feature_label','feature_description','feature_distribution','Nfea','feature_info');