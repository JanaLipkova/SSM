model_name = 'lacy_lacz';

% add one space before and after the species name
reaction = {...
   ' PLac + RNAP -> PLacRNAP ',...
   ' PLacRNAP -> PLac + RNAP ',...
   ' PLacRNAP -> TrLacZ1 ',...
   ' TrLacZ1 -> RbsLacZ + PLac + TrLacZ2 ',...
   ' TrLacZ2 -> TrLacY2 ',...
   ' TrLacY1 -> RbsLacY + TrLacY2 ',...
   ' TrLacY2 -> RNAP ',...
   ' Ribosome + RbsLacZ -> RbsribosomeLacZ ',...
   ' RbsribosomeLacZ -> Ribosome + RbsLacZ ',...
   ' Ribosome + RbsLacY -> RbsribosomeLacY ',...
   ' RbsribosomeLacY -> Ribosome + RbsLacY '...
   ' RbsribosomeLacZ -> TrRbsLacZ + RbsLacZ ',...
   ' RbsribsomeLacY -> TrRbsLacY + RbsLacY ',...
   ' TrRbsLacZ -> LacZ ',...
   ' TrRbsLacZ -> LacY ',...
   ' LacZ -> dgrLacZ ',...
   ' LacY -> dgrLacY ',...
   ' RbsLacZ -> dgrLacY ',...
   ' RbsLacZ -> dgrRbsLacY ',...
   ' LacZ + lactose -> LacZlactose ',...
   ' LacZlactose -> product + LacZ ',...
   ' LacY -> lactose + LacY '
};
   

rate = [ 0.17;
        0.17;
        1;
        1;
        0.015;
        1;
        0.36;
        0.17;
        0.17;
         0.45;
        0.45;
        0.4;
        0.4;
        0.015;
        0.036;
        6.42e-05;
        6.42e-05;
        0.3;
        0.3;
        9.52e-05;
        431;
        14; 
  ];


species={
   'PLac',...
   'RNAP',...
   'PLacRNAP',...
   'TrLacZ1',...
   'RbsLacZ',...
   'TrLacZ2',...
   'TrLacY2',...
   'TrLacY1',...
   'RbsLacY',...
   'Ribosome',...
   'RbsribosomeLacZ',...
   'RbsribosomeLacY',...
   'TrRbsLacZ',...
   'RbsribsomeLacY',...
   'TrRbsLacY',...
   'LacZ',...
   'LacY',...
   'dgrLacZ',...
   'dgrLacY',...
   'dgrRbsLacY',...
   'lactose',...
   'LacZlactose',...
   'product',...
   };