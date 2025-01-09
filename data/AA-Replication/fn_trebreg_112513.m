function [beta,beta_eth,r2,resid,r2_within] = fn_trebreg_112513(y,x,x2,FEi,FEj)

    % getting the sample down to finite numbers
        sample = isfinite(y) & isfinite(x);
        y = y(sample);
        x = x(sample);
        FEi = FEi(sample);
        FEj = FEj(sample);
        for j =1:size(x2,2)
            tempa = x2(:,j);
            temp(:,j) = tempa(sample);
        end      
        x2 = temp;
        
    % getting the fixed effect variables
        dummy1 = ((ones(size(x))*(unique(FEi)'))==(FEi*ones(1,size(unique(FEi),1))));
        dummy2 = ((ones(size(x))*(unique(FEj)'))==(FEj*ones(1,size(unique(FEj),1))));
        
    % doing the regression (getting rid of multicollinearities by using the
    % generalized inverse)
        X = [x,ones(size(x,1),1),x2,dummy1(:,2:size(dummy1,2)),dummy2(:,2:size(dummy2,2))];
        warning('off','all')
        temp_beta = (X'*X)\(X'*y);
        [temp1,temp2,resid,temp4,stats] = regress(y,X);
        warning('on','all')
        beta = temp_beta(1);
        beta_eth = temp_beta(3:2+size(x2,2)); % the coefficient on the other covariates
        r2 = stats(1);
        
    % getting the within r-2
    
        % de-fixed effecting all the variables
            temp_X = double([dummy1(:,2:size(dummy1,2)),dummy2(:,2:size(dummy2,2))]);
            temp_beta = ((temp_X')*temp_X)\(temp_X'*x);
            x_dm = x - temp_X*temp_beta;
            x2_dm = NaN(size(x2));
            for j =1:size(x2,2)
                temp_beta = (temp_X'*temp_X)\(temp_X'*x2(:,j));
                x2_dm(:,j) = x2(:,j) - temp_X*temp_beta;;
            end
            temp_beta = (temp_X'*temp_X)\(temp_X'*y);
            y_dm = y - temp_X*temp_beta;
            
            % doing the de-meaned regression
                [temp1,temp2,temp3,temp4,tempstats] = regress(y_dm,[x_dm,x2_dm,ones(size(x_dm))]);
                r2_within = tempstats(1);
