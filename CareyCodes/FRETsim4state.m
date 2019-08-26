function traj = FRETsim4state(trajtime,val0,val1,val2,val3,tau01,tau10,tau02,tau20,tau12,tau21,tau23,tau32)

state = 0; % Initial state

traj = [];

% Build time array
while length(traj) < trajtime
    if state == 0
        temp = -tau01*log(1 - rand(1,1));
        temp2 = -tau02*log(1 - rand(1,1));
        if temp <= temp2
            trajconc = ones(ceil(temp),1)*val0;
            state = 1;
        else
            trajconc = ones(ceil(temp2),1)*val0;
            state = 2;
        end
        traj = [traj;trajconc];
    elseif state == 1
        temp = -tau10*log(1 - rand(1,1));
        temp2 = -tau12*log(1 - rand(1,1));
        if temp <= temp2
            trajconc = ones(ceil(temp),1)*val1;
            state = 0;
        else
            trajconc = ones(ceil(temp2),1)*val1;
            state = 2;
        end
        traj = [traj;trajconc];
    elseif state == 2
        temp = -tau21*log(1 - rand(1,1));
        temp2 = -tau20*log(1 - rand(1,1));
        temp3 = -tau23*log(1 - rand(1,1));
        if temp <= temp2
            if temp <= temp3
                trajconc = ones(ceil(temp),1)*val2;
                state = 1;
            else
                trajconc = ones(ceil(temp3),1)*val2;
                state = 3;
            end
        else
            if temp2 <= temp3
                trajconc = ones(ceil(temp2),1)*val2;
                state = 0;
            else
                trajconc = ones(ceil(temp3),1)*val2;
                state = 3;
            end
        end
        traj = [traj;trajconc];
    else
        temp = -tau32*log(1 - rand(1,1));
        trajconc = ones(ceil(temp),1)*val3;
        state = 2;
        traj = [traj;trajconc];
    end
end

traj = traj(1:trajtime);

end