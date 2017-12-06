function txt = myupdatefcn(~,event_obj)
    pos = get(event_obj,'Position');
    txt={['Frame: ' num2str(round(pos(2)))],['Cell: ' num2str(round(pos(3)))]};
end