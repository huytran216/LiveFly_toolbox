    %       % 1. ON/OFF
    %       % 2. Relative First Activation time
    %       % 3. Relative Last Activation time
    %       % 4. Relative Total active duration
    %       % 5. Relative Total spot duration of existence 
    %       % 6. Relative Time to max intensity
    %       % 7. Brightest spot intensity
    %       % 8. Mean Intensity produced (absolute time) (including non-expressing nuclei and non-expressing frame)
    %       % 9. Total Intensity produced (including non-expressing nuclei and non-expressing frame)
    %       % 10. Total interphase duration
    %       % 11. Absolute First Activation time
    %       % 12. Absolute Last Activation time
    %       % 13. Absolute Total active duration
    %       % 14. Absolute Total spot duration of existence
    %       % 15. Absolute Time to max intensity
    %       % 16. Mean Spot Intensity (including non-expressing nuclei and non-expressing frame)
    %       % 17. Mean Spot Intensity (only expressing nuclei (full traces) including non-expressing frames)
    %       % 18. Mean Spot Intensity (only expressing nuclei (trimmed traces) including non-expressing frames)
    %       % 19. Mean Spot Intensity (only expressing nuclei and expressing frames)
    %       % 20. Mean PSpot (including non-expressing nuclei)
    %       % 21. Mean PSpot (only non-expressing nuclei (full traces))
    %       % 22. Mean PSpot (only non-expressing nuclei (trimmed traces))
    %       % 23. Mean Spot Intensity (after 1st spot appearance, 1st burst only, with a time limit)
    %       % 24. Mean Pspot (after 1st spot appearance, 1st burst only, with a time limit)
    %       % 25. Mean Spot Intensity (after 1st spot appearance, within a time window)
    %       % 26. Mean Pspot (after 1st spot appearance, within a time window)
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
    'P_{Spot}',...20
    'P_{Spot,ONfull}',...21
    'P_{Spot,ONtrim}',...22
    '\mu_{I}_{burst}',...23
    'P_{Spot}_{burst}',...24
    '\mu_{I}_{ONwindow}',...25
    'P_{Spot}_{ONwindow}',...26
    };
Nfea=numel(feature_label);
save('feature_label.mat','feature_label','Nfea');