function [cost] = Costfunction_assist(beta, tau_data, fun_data, config)

    fit_para = config.fit_para;
    value_0 = fit_para(:,3);
    %value_0 = config.fit_para(:,3)';
    NVars = size(config.fit_para,1);
    kz = config.kz;
    kr = config.kr;
    G = config.G;
    vhc = config.vhc;
    d = config.d;
    for index = 1:1:NVars
        switch config.fit_para(index,2)
            case 1
                kz(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 2
                kr(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 3
                vhc(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 4
                d(config.fit_para(index,1)) = value_0(index)*beta(index);
            case 5
                G(config.fit_para(index,1)) = value_0(index)*beta(index);
        end
    end
    % equivalent parameter
    if isfield(config, 'equal_para')
        for index = 1:1:size(config.equal_para,1)
            switch config.equal_para(index,3)
                case 1
                    kz(config.equal_para(index,2)) = value_0(config.equal_para(index,1))*beta(config.equal_para(index,1));
                case 2
                    kr(config.equal_para(index,2)) = value_0(config.equal_para(index,1))*beta(config.equal_para(index,1));
                case 3
                    vhc(config.equal_para(index,2)) = value_0(config.equal_para(index,1))*beta(config.equal_para(index,1));
                case 4
                    d(config.equal_para(index,2)) = value_0(config.equal_para(index,1))*beta(config.equal_para(index,1));
                case 5
                    G(config.equal_para(index,2)) = value_0(config.equal_para(index,1))*beta(config.equal_para(index,1));
            end
        end
    end
    NLayer = size(config.kz, 2);
    for index = 1:1:NLayer
        if config.iso(index) == 1
            kr(index) = kz(index);
        end
    end
    [cost] = Costfunction(kz,kr,G,d,vhc,config.w,tau_data,fun_data, config);
end