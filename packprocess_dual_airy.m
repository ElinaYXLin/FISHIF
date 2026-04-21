cziconvert2('Duallist.xls');
flyimagetile4_quick_bhh('Duallist.xls',[],[],2);
% % % flyimagetile4_quick('Duallist.xls',[],[],1);
Dualprocess0;
emmask_manual;
DAPI_seg3D2('Duallist.xls',[],[],[],[],[],1/3);
DAPI_seg3D2('Duallist.xls',[],0.7,[],3,[1000,50000],1/3,[],0,[2,2,0],[],2.1,1.2);
% % DAPI_seg3D2('Duallist.xls',[],0.65,{1:2,{[1,2,4],1}},[],[1000,30000],1/3,[],0); DAPI_seg3D2('Duallist.xls',[],0.65,{1,3},0,[1000,30000],1/3,[],1.5,15);
% % DAPI_seg3D2('Duallist.xls',[],0.7,{1,4},[],[2000,50000],1/3,[],0); DAPI_seg3D2('Duallist.xls',[],0.7,{2,1},0,[2000,50000],1/3,[],1.5,15);
% % % % DAPI_seg3D2('Duallist.xls',[],0.7,{1,3:6},0,[1500,50000],1/3,[],0,0);
% % % % DAPI_seg3D2('Duallist.xls',[],0.7,{3,1},0,[500,50000],1/3,[],0,[0,5,0]);
% % % % DAPI_seg3D2dense('Duallist.xls',[],0.6,{1,1},0,[500,50000],1/3,[],0,0,0);
% % % % DAPI_seg3D2fractal('Duallist.xls',[],0.7,{2,3},0,[500,50000],1/3,[],0,0,0);
% % % % DAPI_seg3D2('Duallist.xls',[],0.5,{1,5},0,[1000,50000],[],[],0,[3,3,0],[],4);
% % % % DAPI_seg3D2({'\\192.168.1.108\wangjingyaodata\FISHIF_WANG\confocal\S1\12242019_S1_60X\'},[],0.5,{1,1},0,[1000,50000],[],[],0,[5,5,5],[],5);

stack_RNA
stack_protein_new4
stack_protein_enrich
% stack_RNA_check([],1/3)
% foci_alignment_fun
% stack_protein_purify

hist_fit_RNA('Duallist.xls');
hist_fit_protein('Duallist.xls');

% nuclei_manual3D_foci([],1/3)
% nuclei_manual3D_supple('Duallist.xls',[],1/3,500,[],[],[],true);
% nuclei_overlap_detect('Duallist.xls',[],1/3);
% nushape_check;

% Dualprocess7_2RNA;

