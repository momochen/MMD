% Copyright Yu Chen, Tetherless World Constellation,Rensselaer Polytechnic Institute
% Email: cheny18@rpi.edu

% This is the testing programming for the multi dimensional motif discovery
% algorithm

% Generate data for testing
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

% Normalize the data
a_norm = normoliazation(test_data);
a_norm=a_norm';

b_norm = normoliazation(test_datab);
b_norm=b_norm';


c_norm = normoliazation(test_datac);
c_norm=c_norm';

% Configure the length of the time series,number of symbols generated and 
% alphabetic size
N_1=length(a_norm');
n_1=100;
win_size_1=N_1/n_1;

N_2=length(a_norm');
n_2=100;
win_size_2=N_2/n_2;

N_3=length(a_norm');
n_3=100;
win_size_3=N_3/n_3;

alphabetic_size = 20;

% Run SAX to convert time series data to symbolic strings
str1 = timeseries2symbol(a_norm,N_1,n_1,alphabetic_size);
str2 = timeseries2symbol(b_norm,N_2,n_2,alphabetic_size);
str3 = timeseries2symbol(c_norm,N_3,n_3,alphabetic_size);

% Preliminary base motif detection
[motif1 motif1occ motif_id_pos_end1]= momo_find_motif(str1,win_size_1);

[motif2 motif2occ motif_id_pos_end2]= momo_find_motif(str2,win_size_2);

[motif3 motif3occ motif_id_pos_end3]= momo_find_motif(str3,win_size_3);

% Motif merge to obtain the most appropriate motif
[merged_motif1 merged_motif1_occ_length]= momo_mergemotif(motif1,motif1occ,win_size_1,win_size_1);
[merged_motif2 merged_motif2_occ_length]= momo_mergemotif(motif2,motif2occ,win_size_1,win_size_1);
[merged_motif3 merged_motif3_occ_length]= momo_mergemotif(motif3,motif3occ,win_size_1,win_size_1);

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

% Stack the motifs together
motifs(1).motif = merged_motif1;
motifs(2).motif = merged_motif2;
motifs(3).motif = merged_motif3;

motif_occ(1).motif = merged_motif1_occ_length;
motif_occ(2).motif = merged_motif2_occ_length;
motif_occ(3).motif = merged_motif3_occ_length;

% Multi-dimensional motif detection
[motif_collection incidence_table_weight]=momo_incidence_table(motifs);
[multi_dim_motif] = momo_multi_motif_detection(motif_collection,incidence_table_weight,1);
 
% Finally multi-dimensional motif will be output