ss = 1;

ABET_removeheader = ABETfile(2:end,:);

tbl_ABET = cell2table(ABET_removeheader);
tbl_ABET.Properties.VariableNames = ABETfile(1,:);

ABET_TTL_Times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'TTL #1'));

for qq = 1:size(ABET_TTL_Times,1)
    if ABET_TTL_Times(qq) == ABET_TTL_Times(qq+1)
        ABET_TTL_Times(qq+1) = [];
    end
end

trial_ref = 1

for gg = 1:size(tbl_ABET, 1)
    if strcmp(tbl_ABET.Item_Name(gg),'TTL #1')
        tbl_ABET.ABET_Trial_Ref(gg) = trial_ref
        trial_ref = trial_ref + 1
    end
end
    

    
tbl_ABET.Evnt_Time < ABET_TTL_Times(2) & tbl_ABET.Evnt_Time > ABET_TTL_Times(1)


for qq = 1:size(tbl_ABET,1)
    if strcmp(tbl_ABET.Item_Name(qq),'TTL #1')
        tbl_ABET.ABET_Trial_Ref(qq) = 0;
        for zz = 1:size(tbl_ABET,1)
            
    end
end

end

for zz = 1:size(ABET_TTL_Times,1)
    for qq = 1:size(tbl_ABET,1)
        if tbl_ABET.Evnt_Time(qq) == ABET_TTL_Times(zz)
            tbl_ABET.ABET_Trial_Ref(qq) = ABET_TTL_Times(zz);
        end
    end
end


for zz = 2:size(ABETfile,1)
    if regexp(ABETfile{zz, 4},'TTL #1')
        ABETfile{zz,19} = 0;
    end
end


        for qq = zz+1:size(ABETfile,1)
            ABETfile{qq,19} = ABETfile{zz,1}-ABETfile{qq,1}
            qq=qq+1;
        end
    zz=zz+1;
end
end
end
% for zz = 2:size(ABETfile,1)
%     if regexp(ABETfile{zz, 4},'TTL #1')
%         ABETfile{zz,19} = ttl_filtered(ss)
%         ABETfile{zz+1,19} = ttl_filtered(ss
%         
%     end
% end