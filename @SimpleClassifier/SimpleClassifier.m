classdef SimpleClassifier
    
    %Classifies a dataset with two labels [0 1]
    %NB: for this classification the 0 label is considered the true label
    
    properties(SetAccess=private)
    
        %input for the whole class (at constructor)
        dataset
        label
        title1 % A unique identifies for the classification (will also be legend name)
        pat %patient name
        dat %dataset name
        
        %changing with each classification
        label_clas
        Conf
        Acc
        Se
        Sp
        pos
        
        SeSp %TPR and FPR
        AUC
        CI
        Confvec %vector of the Confustion matrix
        
    end

    methods
        
        %Constructor
        function S=SimpleClassifier(dataset, label, title1)
            S.dataset=dataset;
            S.label=label;
            S.title1=title1;
            
            disp(title1)
        end

        
        %TEST1: split into training and test sets
        function S=clas_split(S, split, clas_type, prior)
            
            if nargin <4
                prior=[0.5 0.5];
            end
            
            label_train=S.label(1:split-1);
            label_test=S.label(split:end);
            train=S.dataset(1:split-1,:);
            test=S.dataset(split:end, :);
            
            [S.label_clas,~,S.pos]=classify(test, train, label_train, clas_type, prior);
            
            [S.Conf S.Acc]=S.conf_mat(label_test, S.label_clas);
            S.label=label_test; 
            disp('Warning: labels permanently divided. Create new object for rerun')
        end
        
        
        % TEST2: leave-one-out cross validation
        function S=clas_LOOCV(S, clas_type, prior)            
            if nargin <3
                prior.prob=[0.5 0.5];
                prior.group=[0 1];
            end
            
            
            S.label_clas=zeros(size(S.label));
            for ii=1:length(S.label)
                aa=false(size(S.label));
                aa(ii)=1;

                label_train=S.label(~aa);
                label_test=S.label(aa);
                train=S.dataset(~aa, :);
                test=S.dataset(aa, :);
                
                % Output
                [S.label_clas(ii),~,S.pos(ii,:)]=classify(test, train, label_train, clas_type, prior);
                
            end
            
            [S.Conf S.Acc]=S.conf_mat(S.label, S.label_clas);

        end
        
        %using matlabs built in roc curve generator
        function S=roc_curve_perf_pos(S)
                
                range1=0:0.1:1;
                
                [xx, yy, ~, S.AUC]=perfcurve(S.label, S.pos(:,1), logical(0),...
                    'NBOOT',1000, 'BootType', 'bca', 'TVals', range1);
                S.SeSp(:, 1)=yy(:,1);
                S.SeSp(:,2)=xx(:,1);
                S.CI=yy(:,1)-yy(:,2);  
        end
        
        
        %Display Confusion matrix
        function S=disp_conf(S)
            disp('Confusion matrix')
            formatspec='TN | %d %d | FP\nFN | %d %d | TP\n';
            fprintf(formatspec, S.Conf(1,1), S.Conf(1,2), S.Conf(2,1), S.Conf(2,2));            
        end
        
        %Multiple rocs
        function plot(varargin)
%             varargin=[{S} varargin];
            figure('Units', 'pixels', 'Position', [100 100 600 600]),
            clf;
%             hTitle=title1('ROC with 95% confidence intervals for detection of TB cases (Bootstrap 1000 using BCA)');
            hXLabel=xlabel('1-Specificity');
            hYLabel=ylabel('Sensitivity');
            cols1=([[0 206 209]/255; [20 50 170]/255; [112 138 144]/255; [138 43 226]/255]);
            marker1='sdvsdvsdvsdvsdvsdv';
            
            hold on
            for ii=1:size(varargin, 2)

                x=varargin{ii}.SeSp(:,2);
                y1=varargin{ii}.SeSp(:,1);
                c1=y1-varargin{ii}.CI;
                c2=y1+varargin{ii}.CI;
                e1=varargin{ii}.CI;
                
                h=errorbar(x, y1, e1, [marker1(ii) '--']);
                set(h, 'MarkerFaceColor', cols1(mod(ii,size(cols1,1))+1,:), 'MarkerEdgeColor', [.2 .2 .2], 'LineWidth', 0.9);
                set(h, 'MarkerSize', 7, 'Color', [.2 .2 .2]);
                
                %title
                varargin{ii}.AUC=round(varargin{ii}.AUC*100)/100;
                
                if length(varargin{ii}.AUC)==1
                    leg1{ii}=[varargin{ii}.title1 ' AUC: ' sprintf('%g', varargin{ii}.AUC(1))];
                elseif length(varargin{ii}.AUC)==3
                    leg1{ii}=[varargin{ii}.title1 ' AUC: ' sprintf('%g', varargin{ii}.AUC(1))...
                        ' (' sprintf('%g', varargin{ii}.AUC(2)) '-' sprintf('%g', varargin{ii}.AUC(3)) ')'];

                end
                
                
            end
            %plot([0 1], [0 1])
            hold off
            hLegend=legend(leg1);
            axis([0 1 0 1])
            axis square;
            
            set(gca, 'FontName', 'Helvetica');
            set(hLegend, 'FontSize', 12, 'Location', 'South', 'Box', 'off'); 
            set([hXLabel, hYLabel], 'FontSize', 12);
%             set(hTitle, 'FontSize', 10)
            
            set(gcf, 'PaperPositionMode', 'auto');
            print -depsc2 finalPlot1.eps
        end
        
    end
    
    methods(Static)
        
        %confidence
        function [C A]=conf_mat(truelab, claslab)
            %TP FN
            %FP TN
            C=[sum( (truelab==0) & (claslab==0) ) sum( (truelab==0) & (claslab==1));...
               sum( (truelab==1) & (claslab==0)) sum( (truelab==1) & (claslab==1))];
           
            A=sum(truelab==claslab)/length(truelab);
            
        end
        
         %sensitivity and specificity
        function [Sens Spec]=sens_spec(Conf)
            %Sensitivity=TP/(TP+FN)
            Sens=Conf(1,1)/(Conf(1,1)+Conf(1,2));
            %Specificity=TN/(TN+FP)
            Spec=Conf(2,2)/(Conf(2,2)+Conf(2,1));
            
        end
        
        %display sensitivty and specificity
        function disp_sens_spec(Sensitivity, Specificity)
            fprintf('Sensitivity %.2f and Specificity %.2f \n', Sensitivity, Specificity)
        end
        
        function [AUC]=area_under_curve(SeSp)
            %trapezoidal rule
            %1/2*f1*f1*deltaX
            delX=SeSp(1,2);
            AUC=0.5*(SeSp(1,1)+0)*delX;
            for jj=1:length(SeSp)-1
                delX=SeSp(jj+1,2)-SeSp(jj,2);
                AUC=AUC + 0.5*(SeSp(jj, 1)+SeSp(jj+1,1))*delX;
            end
            delX=1-SeSp(end,2);
            AUC=AUC + 0.5*(SeSp(end,1) + 1)*delX;
            
        end 

    end
end
