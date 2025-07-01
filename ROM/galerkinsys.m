function gsys = galerkinsys(~,X,L,QQ,F,Re);
%nmodes = length(X);
XX = X*X.'; XX = XX(:);

Xp = 1/Re*L*X + QQ*XX + F; Xp = full(Xp);

%for i=1:nmodes
%    for j=1:nmodes
%        for k=1:nmodes
%            Xp(i) = Xp(i)+Q(i,j,k)*X(j)*X(k);
%            Qaux = squeeze(Q(i,:,:)); Qaux = Qaux(:).';
%            Xp(i) = Xp(i)+Qaux*XX;
%        end
%    end
%end
    

gsys = Xp;
end