%{
str1 = [8 8 7 5 1 9 9 7 9 8 4 4 8 8 7 5 8 8 7 5 1 1 1 1];
str2 = [9 0 9 1 9 6 5 9 0 9 0 1 9 0 9 1 9 0 9 1 4 6 5 2];
str3 = [6 6 0 1 2 5 0 9 6 6 0 6 6 0 1 1 9 6 6 0 1 1 2 3];

plot(str1);
%}

x = 0:0.02:2*pi;
sinx=sin(x);
test_data = sin(x);
test_data(316:510)=0;
test_data(511:800)= sinx(1:290);
test_data(801:904)=2.5;
test_data(905:1200)=sinx(1:296);
test_data(1201:1300)=1;
test_data(1301:1500)=sinx(81:280);
test_data(1501:1700)=-2;
test_data(1701:1850)=sinx(51:200);
test_data(1851:2000)=0;
test_data = test_data';

test_datab = sin(x);
test_datab(316:510)=0;
test_datab(511:800)= sinx(1:290);
test_datab(801:904)=2.5;
test_datab(905:1200)=sinx(1:296);
test_datab(1201:1300)=1;
test_datab(1301:1500)=sinx(81:280);
test_datab(1501:1700)=-2;
test_datab(1701:1850)=sinx(51:200);
test_datab(1851:2000)=0;
test_datab = test_datab';

test_datac = sin(x);
test_datac(316:510)=0;
test_datac(511:800)= sinx(1:290);
test_datac(801:904)=2.5;
test_datac(905:1200)=sinx(1:296);
test_datac(1201:1300)=1;
test_datac(1301:1500)=sinx(81:280);
test_datac(1501:1700)=-2;
test_datac(1701:1850)=sinx(51:200);
test_datac(1851:2000)=0;
test_datac = test_datac';

%plot(test_data);

%a = CK_TAR_BGL(1001:3000,1);
%b = CK_TAR_BOD(1001:3000,1);
%c = CK_TAR_SESKMAC(1001:3000,1);
a_norm = normoliazation(test_data);
a_norm=a_norm';

b_norm = normoliazation(test_datab);
b_norm=b_norm';


c_norm = normoliazation(test_datac);
c_norm=c_norm';


N=2000;
n=100;
win_size=N/n;

str1 = timeseries2symbol(a_norm,N,n,20);
str2 = timeseries2symbol(b_norm,N,n,20);
str3 = timeseries2symbol(c_norm,N,n,20);

%plot(str2);
%hold on;
%str2
[motif1 motif1occ motif_id_pos_end1]= momo_find_motif(str1,win_size);

[motif2 motif2occ motif_id_pos_end2]= momo_find_motif(str2,win_size);

[motif3 motif3occ motif_id_pos_end3]= momo_find_motif(str3,win_size);


[merged_motif1 merged_motif1_occ_length]= momo_mergemotif(motif1,motif1occ,win_size,win_size);
%merged_motif1_occ_length
[merged_motif2 merged_motif2_occ_length]= momo_mergemotif(motif2,motif2occ,win_size,win_size);
%merged_motif2_occ_length
[merged_motif3 merged_motif3_occ_length]= momo_mergemotif(motif3,motif3occ,win_size,win_size);
%merged_motif3_occ_length

%[num_of_motif motif_size]=size(merged_motif1);

%{
for i=1:num_of_motif
    x=merged_motif(i,1):1:(merged_motif(i,1)+motif_size-2);
    y=merged_motif(i,2:motif_size);
    plot(x,y,'LineWidth',2,'color','r');
    hold on;
end

%plot(motif1);
%[motif2 motif2occ] = find_motif(str2,n,win_size);
%motif2
%[motif3 motif3occ] = find_motif(str3,n,win_size);
%motif3
%}


motifs(1).motif = merged_motif1;
motifs(2).motif = merged_motif2;
motifs(3).motif = merged_motif3;

motif_occ(1).motif = merged_motif1_occ_length;
motif_occ(2).motif = merged_motif2_occ_length;
motif_occ(3).motif = merged_motif3_occ_length;

%[motif_collection weighted_incidence_table] = momo_incidence_table(motifs);


[motif_collection incidence_table_weight]=momo_incidence_table(motifs);
[multi_dim_motif] = momo_multi_motif_detection(motif_collection,incidence_table_weight,1);
 
%clustered_motifs = multi_dim_motif_dec(all_motif_length_occ,incidence_table_weight,win_size);
% [num_of_motifs sb] = size(clustered_motifs);
 
%}