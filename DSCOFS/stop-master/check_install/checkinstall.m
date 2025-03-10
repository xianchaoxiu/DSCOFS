try
    FOForth_checkinstall
    fprintf('FoForth is installed correctly.\n')
catch errmsg
    warning('FOForth may not be installed currectly. \n')
end

try
    PCAL_checkinstall
    fprintf('PCAL is installed correctly.\n')
catch errmsg
    warning('PCAL may not be installed currectly. \n ')
end


try
    PenCF_checkinstall
    fprintf('PenCF is installed correctly.\n')
catch errmsg
    warning('PenCF may not be installed currectly. \n')
end


try
    PenCS_checkinstall
    fprintf('PenCS is installed correctly.\n')
catch errmsg
    warning('PenCS may not be installed currectly. \n')
end


try
    PenCPG_checkinstall
    fprintf('PenCPG is installed correctly.\n')
catch errmsg
    warning('PenCPG may not be installed currectly. \n')
end


try
    ProxOrth_checkinstall
    fprintf('ProxOrth is installed correctly.\n')
catch errmsg
    warning('ProxOrth may not be installed currectly. \n')
end



try
    SLPG_checkinstall
    fprintf('SLPG is installed correctly.\n')
catch errmsg
    warning('SLPG may not be installed currectly. \n')
end


try
    SLPG_l21_checkinstall
    fprintf('SLPG_l21 is installed correctly.\n')
catch errmsg
    warning('SLPG_l21 may not be installed currectly. \n')
end


try
    SLPG_smooth_checkinstall
    fprintf('SLPG_smooth is installed correctly.\n')
catch errmsg
    warning('SLPG_smooth may not be installed currectly. \n')
end

