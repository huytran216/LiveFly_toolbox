mov_list = {'\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200721_GF1\table_summary\',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200721_GF2\table_summary\',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200723_GF1\table_summary\',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200723_GF2\table_summary\',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200723_GF3\table_summary',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200723_GF4\table_summary\',...
    };

for i=1:numel(mov_list)
    LiveFly_MITOSIS(mov_list{i},1,true);    
end