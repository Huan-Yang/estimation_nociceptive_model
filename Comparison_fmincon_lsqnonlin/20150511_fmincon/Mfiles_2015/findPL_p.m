function F=findPL_p(p,parami,psub,current_PL,qPLthres,Data)

pall(parami)=p;
pall(setdiff([1:6],parami))=psub;
funhPL=@param2PL;
F=feval(funhPL,pall,Data)-(current_PL+qPLthres);