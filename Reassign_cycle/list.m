mov_list = {'\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200205_GF1\table_summary\',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200205_GF2\table_summary\',...
    '\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\200205_GF4\table_summary\'
    };

for i=1:numel(mov_list)
    LiveFly_MITOSIS(mov_list{i},1,true);    
end