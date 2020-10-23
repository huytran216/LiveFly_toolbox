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
    %       % 20. Mean Maximum Spot Intensity (only expressing nuclei and expressing frames)
    %       % 21. Mean PSpot (including non-expressing nuclei)
    %       % 22. Mean PSpot (only expressing nuclei (full traces))
    %       % 23. Mean PSpot (only expressing nuclei (trimmed traces))
    %       % 24. Mean Spot Intensity (after 1st spot appearance, 1st burst only, with a time limit)
    %       % 25. Mean Pspot (after 1st spot appearance, 1st burst only, with a time limit)
    %       % 26. Mean Spot Intensity (after 1st spot appearance, within a time window)
    %       % 27. Mean Pspot (after 1st spot appearance, within a time window)
    %       % 28. Number of burst (whole trace)
    %       % 29. Number of burst (within a window)
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
Nfea=numel(feature_label);
save('feature_label.mat','feature_label','Nfea');