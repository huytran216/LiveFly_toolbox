function newfilename=conv_tiff(fileName,framenum)
    info=imfinfo([fileName '.tif']);
    fp = fopen([fileName '.tif'], 'rb');
    for cnt = 1:framenum
        if cnt==1
            writemode='overwrite';
        else
            writemode='append';
        end
        fseek(fp, info(cnt).StripOffsets, 'bof');
        imData = uint8(fread(fp, [info(cnt).Width info(cnt).Height], 'uint8', 0, 'ieee-be'))';
        imwrite(imData,[fileName '_interleaved.tif'],'tif','compression','none','writemode',writemode);
    end
    fclose(fp);
    newfilename=[fileName '_interleaved'];
end