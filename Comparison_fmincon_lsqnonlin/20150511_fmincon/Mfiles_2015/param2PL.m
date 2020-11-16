function PL=param2PL(param,Data)
for TSi=1:length(Data)
    A_list{TSi}=unique(Data(TSi).Stimuli);
    A_list{TSi}=A_list{TSi}(:);
    TS(TSi,:)=[Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    
    len_ind{TSi}=zeros(size(A_list{TSi}));
    R1{TSi}=zeros(size(A_list{TSi}));
    for Aj=1:length(A_list{TSi})
        ind_list=find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001); % to prevent possible small difference in the collect amplitudes
        len_ind{TSi}(Aj)=length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj)=sum(Data(TSi).Responses(ind_list)); % the number of trils with 'detected' response
%         [phat{TSi}(Aj),pci{TSi}(Aj,:)] = binofit(R1{TSi}(Aj),length(ind_list));
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;
[F]=nln_likelihood_component_J(param,TS,A_list,R1,len_ind);
PL=2*sum(F.^2);