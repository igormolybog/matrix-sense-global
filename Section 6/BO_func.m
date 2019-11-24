function [min_ball, min_sphere] = BO_func(caseN, R, omega)
%BO_FUNC Summary of this function goes here
%   Detailed explanation goes here
    addpath('../util/');
    
    A = casetoA(caseN);
    H = A.'*A;
    [S,W] = get_structure(H);
    
    [ballvars, spherevars] = getvars(caseN, R);
    
    if caseN == 9
        obj_ball = @(t) fun9_ball(t, S, W, omega);
        obj_sphere = @(t) fun9_sphere(t, R, S, W, omega);
    elseif caseN == 14
        obj_ball = @(t) fun14_ball(t, S, W, omega);
        obj_sphere = @(t) fun14_sphere(t, R, S, W, omega);
    elseif caseN == 30
        obj_ball = @(t) fun30_ball(t, S, W, omega);
        obj_sphere = @(t) fun30_sphere(t, R, S, W, omega);
    else
        'ERROR: caseN must take values 9, 14 or 30. BO_func got: '
        caseN
    end
    
    sprintf('Working on case%d, radius %f...', caseN, R)    
    
    result_ball = bayesopt(obj_ball, ballvars, 'PlotFcn', {},'Verbose',0);
    result_sphere = bayesopt(obj_sphere, spherevars, 'PlotFcn', {},'Verbose',0); 
    
    min_ball = result_ball.MinObjective;
    min_sphere = result_sphere.MinObjective;
    
    sprintf('solved. ball: %f, sphere: %f', min_ball, min_sphere)
end

function [bvars, svars] = getvars(caseN, R)
% generating variables for both LMI_P (bvars)
% and LMI_P_dB (svars)
    [rxball, rzball] = ballRvars(R);
    [rxsphere, rzsphere] = ballRvars(R); % rxsphere is unused
    [phixball, phizball] = phivars(caseN);
    [phixsphere, phizsphere] = phivars(caseN);
    bvars = [rxball phixball rzball phizball];
    svars = [phixsphere rzsphere phizsphere];
end

function [rx, rz] = ballRvars(R)
    rx = optimizableVariable('x1',[0,R]);
    rz = optimizableVariable('z1',[0,R]);
end


function result = fun30_ball(t, S, W, omega)
    
    rx = t.x1;
    rz = t.z1;
    phix = [t.x2;t.x3;t.x4;t.x5;t.x6;t.x7;t.x8;t.x9;t.x10;t.x11;t.x12;t.x13;t.x14;...
        t.x15;t.x16;t.x17;t.x18;t.x19;t.x20;t.x21;t.x22;t.x23;t.x24;t.x25;t.x26;t.x27;t.x28;t.x29;t.x30].';
    phiz = [t.z2;t.z3;t.z4;t.z5;t.z6;t.z7;t.z8;t.z9;t.z10;t.z11;t.z12;t.z13;t.z14;...
        t.z15;t.z16;t.z17;t.z18;t.z19;t.z20;t.z21;t.z22;t.z23;t.z24;t.z25;t.z26;t.z27;t.z28;t.z29;t.z30].';
    
    xcart = npol2ncart(rx, phix).';
    zcart = npol2ncart(rz, phiz).';
    
    result = LMI_P(xcart+omega, zcart+omega, S, W);
end

function result = fun30_sphere(t, R, S, W, omega)
    
    rz = t.z1;
    phix = [t.x2;t.x3;t.x4;t.x5;t.x6;t.x7;t.x8;t.x9;t.x10;t.x11;t.x12;t.x13;t.x14;...
        t.x15;t.x16;t.x17;t.x18;t.x19;t.x20;t.x21;t.x22;t.x23;t.x24;t.x25;t.x26;t.x27;t.x28;t.x29;t.x30].';
    phiz = [t.z2;t.z3;t.z4;t.z5;t.z6;t.z7;t.z8;t.z9;t.z10;t.z11;t.z12;t.z13;t.z14;...
        t.z15;t.z16;t.z17;t.z18;t.z19;t.z20;t.z21;t.z22;t.z23;t.z24;t.z25;t.z26;t.z27;t.z28;t.z29;t.z30].';

    xcart = npol2ncart(R, phix).';
    zcart = npol2ncart(rz, phiz).';
    
    result = LMI_P_dB(xcart+omega, zcart+omega, S, W, omega);
end

function result = fun14_ball(t, S, W, omega)

    rx = t.x1;
    rz = t.z1;
    phix = [t.x2;t.x3;t.x4;t.x5;t.x6;t.x7;t.x8;t.x9;t.x10;t.x11;t.x12;t.x13;t.x14].';
    phiz = [t.z2;t.z3;t.z4;t.z5;t.z6;t.z7;t.z8;t.z9;t.z10;t.z11;t.z12;t.z13;t.z14].';

    xcart = npol2ncart(rx, phix).';
    zcart = npol2ncart(rz, phiz).';
    
    result = LMI_P(xcart+omega, zcart+omega, S, W);
end

function result = fun14_sphere(t, R, S, W, omega)
    
    rz = t.z1;
    phix = [t.x2;t.x3;t.x4;t.x5;t.x6;t.x7;t.x8;t.x9;t.x10;t.x11;t.x12;t.x13;t.x14].';
    phiz = [t.z2;t.z3;t.z4;t.z5;t.z6;t.z7;t.z8;t.z9;t.z10;t.z11;t.z12;t.z13;t.z14].';

    xcart = npol2ncart(R, phix).';
    zcart = npol2ncart(rz, phiz).';
    
    result = LMI_P_dB(xcart+omega, zcart+omega, S, W, omega);
end


function result = fun9_ball(t, S, W, omega)
    
    rx = t.x1;
    rz = t.z1;
    phix = [t.x2;t.x3;t.x4;t.x5;t.x6;t.x7;t.x8;t.x9].';
    phiz = [t.z2;t.z3;t.z4;t.z5;t.z6;t.z7;t.z8;t.z9].';

    xcart = npol2ncart(rx, phix).';
    zcart = npol2ncart(rz, phiz).';
    
    result = LMI_P(xcart+omega, zcart+omega, S, W);
    
end

function result = fun9_sphere(t, R, S, W, omega)
    
    rz = t.z1;
    phix = [t.x2;t.x3;t.x4;t.x5;t.x6;t.x7;t.x8;t.x9].';
    phiz = [t.z2;t.z3;t.z4;t.z5;t.z6;t.z7;t.z8;t.z9].';

    xcart = npol2ncart(R, phix).';
    zcart = npol2ncart(rz, phiz).';
    
    result = LMI_P_dB(xcart+omega, zcart+omega, S, W, omega);
    
end

% getting variables for optimization with BiasOpt tool
function [phix, phiz] = phivars(caseN)
    if caseN == 9
        x2 = optimizableVariable('x2',[0, pi]);
        x3 = optimizableVariable('x3',[0, pi]);
        x4 = optimizableVariable('x4',[0, pi]);
        x5 = optimizableVariable('x5',[0, pi]);
        x6 = optimizableVariable('x6',[0, pi]);
        x7 = optimizableVariable('x7',[0, pi]);
        x8 = optimizableVariable('x8',[0, pi]);
        x9 = optimizableVariable('x9',[0, 2*pi]);

        z2 = optimizableVariable('z2',[0, pi]);
        z3 = optimizableVariable('z3',[0, pi]);
        z4 = optimizableVariable('z4',[0, pi]);
        z5 = optimizableVariable('z5',[0, pi]);
        z6 = optimizableVariable('z6',[0, pi]);
        z7 = optimizableVariable('z7',[0, pi]);
        z8 = optimizableVariable('z8',[0, pi]);
        z9 = optimizableVariable('z9',[0, 2*pi]);
        phix = [x2,x3,x4,x5,x6,x7,x8,x9];
        phiz = [z2,z3,z4,z5,z6,z7,z8,z9];
    elseif caseN == 14
        x2 = optimizableVariable('x2',[0, pi]);
        x3 = optimizableVariable('x3',[0, pi]);
        x4 = optimizableVariable('x4',[0, pi]);
        x5 = optimizableVariable('x5',[0, pi]);
        x6 = optimizableVariable('x6',[0, pi]);
        x7 = optimizableVariable('x7',[0, pi]);
        x8 = optimizableVariable('x8',[0, pi]);
        x9 = optimizableVariable('x9',[0, pi]);
        x10 = optimizableVariable('x10',[0, pi]);   
        x11 = optimizableVariable('x11',[0, pi]);
        x12 = optimizableVariable('x12',[0, pi]);
        x13 = optimizableVariable('x13',[0, pi]);
        x14 = optimizableVariable('x14',[0, 2*pi]);

        z2 = optimizableVariable('z2',[0, pi]);
        z3 = optimizableVariable('z3',[0, pi]);
        z4 = optimizableVariable('z4',[0, pi]);
        z5 = optimizableVariable('z5',[0, pi]);
        z6 = optimizableVariable('z6',[0, pi]);
        z7 = optimizableVariable('z7',[0, pi]);
        z8 = optimizableVariable('z8',[0, pi]);
        z9 = optimizableVariable('z9',[0, pi]);
        z10 = optimizableVariable('z10',[0, pi]);
        z11 = optimizableVariable('z11',[0, pi]);
        z12 = optimizableVariable('z12',[0, pi]);
        z13 = optimizableVariable('z13',[0, pi]);
        z14 = optimizableVariable('z14',[0, 2*pi]);
        
        phix = [x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14];
        phiz = [z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14];
    elseif caseN == 30
        x2 = optimizableVariable('x2',[0, pi]);
        x3 = optimizableVariable('x3',[0, pi]);
        x4 = optimizableVariable('x4',[0, pi]);
        x5 = optimizableVariable('x5',[0, pi]);
        x6 = optimizableVariable('x6',[0, pi]);
        x7 = optimizableVariable('x7',[0, pi]);
        x8 = optimizableVariable('x8',[0, pi]);
        x9 = optimizableVariable('x9',[0, pi]);
        x10 = optimizableVariable('x10',[0, pi]);
        x11 = optimizableVariable('x11',[0, pi]);
        x12 = optimizableVariable('x12',[0, pi]);
        x13 = optimizableVariable('x13',[0, pi]);
        x14 = optimizableVariable('x14',[0, pi]);
        x15 = optimizableVariable('x15',[0, pi]);
        x16 = optimizableVariable('x16',[0, pi]);
        x17 = optimizableVariable('x17',[0, pi]);
        x18 = optimizableVariable('x18',[0, pi]);
        x19 = optimizableVariable('x19',[0, pi]);
        x20 = optimizableVariable('x20',[0, pi]);
        x21 = optimizableVariable('x21',[0, pi]);
        x22 = optimizableVariable('x22',[0, pi]);
        x23 = optimizableVariable('x23',[0, pi]);
        x24 = optimizableVariable('x24',[0, pi]);
        x25 = optimizableVariable('x25',[0, pi]);
        x26 = optimizableVariable('x26',[0, pi]);
        x27 = optimizableVariable('x27',[0, pi]);
        x28 = optimizableVariable('x28',[0, pi]);
        x29 = optimizableVariable('x29',[0, pi]);
        x30 = optimizableVariable('x30',[0, 2*pi]);

        z2 = optimizableVariable('z2',[0, pi]);
        z3 = optimizableVariable('z3',[0, pi]);
        z4 = optimizableVariable('z4',[0, pi]);
        z5 = optimizableVariable('z5',[0, pi]);
        z6 = optimizableVariable('z6',[0, pi]);
        z7 = optimizableVariable('z7',[0, pi]);
        z8 = optimizableVariable('z8',[0, pi]);
        z9 = optimizableVariable('z9',[0, pi]);
        z10 = optimizableVariable('z10',[0, pi]);
        z11 = optimizableVariable('z11',[0, pi]);
        z12 = optimizableVariable('z12',[0, pi]);
        z13 = optimizableVariable('z13',[0, pi]);
        z14 = optimizableVariable('z14',[0, pi]);
        z15 = optimizableVariable('z15',[0, pi]);
        z16 = optimizableVariable('z16',[0, pi]);
        z17 = optimizableVariable('z17',[0, pi]);
        z18 = optimizableVariable('z18',[0, pi]);
        z19 = optimizableVariable('z19',[0, pi]);
        z20 = optimizableVariable('z20',[0, pi]);
        z21 = optimizableVariable('z21',[0, pi]);
        z22 = optimizableVariable('z22',[0, pi]);
        z23 = optimizableVariable('z23',[0, pi]);
        z24 = optimizableVariable('z24',[0, pi]);
        z25 = optimizableVariable('z25',[0, pi]);
        z26 = optimizableVariable('z26',[0, pi]);
        z27 = optimizableVariable('z27',[0, pi]);
        z28 = optimizableVariable('z28',[0, pi]);
        z29 = optimizableVariable('z29',[0, pi]);
        z30 = optimizableVariable('z30',[0, 2*pi]);
        
        phix = [x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30];
        phiz = [z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,z21,z22,z23,z24,z25,z26,z27,z28,z29,z30];
    else
        'ERROR: caseN must take values 9, 14 or 30. we got: '
        caseN
    end
end
