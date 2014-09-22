function iter1
%%% In this iteration of the model, Cellulose is in abundance with Ammonia
%%% produced by the Oligosaccharide degraders. Oligosaccharides and organic
%%% acids may not be completely used. Proteins are produced by all bacteria.
ac
ao
aoa
cc
co
cn
do
[t,y] = ode45(@diffs,[0,100],[10,10,10,0,0,10]);
    function derivatives = diffs(~,N)
        cellulose_d = N(1);
        oligo_d = N(2);
        org_acid_d = N(3);
        
        oligo = N(4);
        org_acid = N(5);
        ammonia = N(6);
        
        dcellulose_d = ac*(1-p_c/k_c)*cellulose_d;
        doligo_d = ao*(1-p_o/k_o)*oligo_d;
        dorg_acid_d = aoa*(1-p_oa/k_oa)*org_acid_d;
        
        doligo = cc*cellulose_d-do*oligo_d;
        dorg_acid = co*oligo_d-doa*org_acid_d;
        dammonia = cn*n_fixers-do*(cellulose_d+oligo_d+org_acid_d);
             
        derivatives=[dcellulose_d;doligo_d;dorg_acid_d;doligo;dorg_acid;dammonia];
    end
end