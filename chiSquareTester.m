n1 = 124; N1 = 921;
n2 = 75; N2 = 921;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)