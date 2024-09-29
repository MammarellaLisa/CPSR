%           C O N D I T I O N A L     P R O B A B I L I T Y       
%                O F     S U R F A C E     R U P T U R E

% This MATLAB script is associated with the following paper:
% Mammarella, L., et al. (2024). "Conditional probability of surface rupture: 
% a numerical approach for principal faulting" Earthquake Spectra 

% Copyright (c) 2024, MammarellaLisa
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors
%    may be used to endorse or promote products derived from this software 
%    without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clear all

% W = rupture width             (down-dip)
% W_z = rupture width           (vertical projaction)
% Z_s = seismogenic depth       
% Z_hypo = hypocentral depth
% W_top = vertical distance from the hypocenter to the top of the Rupture
% HDR = Hypocentral Depth Ratio
% HDD = Hypocentral Depth Distribution

%% INPUT txt

% MAGNITUDE RANGE
Mw = 5:0.1:8;

fid = fopen('INPUT.txt', 'r');
values = textscan(fid, '%d %d %s %f %f %f %f %f %f', 'Delimiter', ',', 'CommentStyle', '%');
fclose(fid);
MSR = values{1};
SoF = values{2};
HDD = values{3}{1};
dip_mu = values{4};
dip_sigma = values{5};
t_d = values{6};
Zs_mu = values{7};
Zs_sigma = values{8};
t_z = values{9};

%% TAB FOR MSR
tab1 = [
    0, 3, 3.63, 2.5, 0.15;
    0, 4, 3.63, 2.5, 0.15;
    0, 5, 3.88, 2.5, 0.15;
    1, 3, 4.14, 2.5, 0.15;
    1, 4, 4.14, 2.5, 0.15;
    1, 5, 4.22, 2.5, 0.15;
    2, 3, -0.829, 0.323, 0.128;
    2, 4, -1.669, 0.435, 0.087;
    2, 5, -0.543, 0.261, 0.105;
];

idx1 = find(tab1(:,1) == MSR & tab1(:,2) == SoF);
if ~isempty(idx1)

    a = tab1(idx1, 3);
    b = tab1(idx1, 4);
    W_sigma = tab1(idx1, 5);
else
    disp('Incorrect choice for MSR or SoF');
end

% truncation MSR distribution (OBS: ± 1 = 68%; ± 2 = 95%; ± 3 = 99.7%)
t_W = 1;           
eps = t_W * W_sigma;
% truncation DIP distribution
eps_dip = t_d * dip_sigma;
% truncation Zs distribution
eps_Z = t_z * Zs_sigma;

%% TAB FOR HDD
tab2 = {
    'ITA_N', 0.67, 0.21;
    'GB_N', 0.64, 0.25;
    'AGG_N', 0.65, 0.23;
    'ITA_R', 0.65, 0.24;
    'TAI_R', 0.60, 0.25;
    'JAP_R', 0.74, 0.20;
    'AGG_R', 0.66, 0.24;
    'CA_S', 0.67, 0.21;
    'NZ_S', 0.72, 0.21;
    'JAP_S', 0.79, 0.15;
    'AGG_S', 0.70, 0.23;
    };

idx2 = find(strcmp(tab2(:,1), HDD));
if ~isempty(idx2)
    HDD_mu = tab2{idx2, 2};
    HDD_sigma = tab2{idx2, 3};
else
    disp('Incorrect choice for HDD');
end

%% MSR

out = [];

if MSR == 0 || MSR == 1
    calc_LogW_mu = @(Mw) (Mw - a) / b;
elseif MSR == 2
    calc_LogW_mu = @(Mw) a + (b * Mw);
else
    error('Invalid Magnitude Scaling Rrelation!');
end


for i = 1:length(Mw)
    LogW_mu = calc_LogW_mu(Mw(i));  

    truncated = truncate(makedist('Normal', LogW_mu, W_sigma), LogW_mu - eps, LogW_mu + eps);
    LogW_troncati = linspace(LogW_mu - eps, LogW_mu + eps, 25);

    p_LogW = truncated.pdf(LogW_troncati);
    p_LogW = p_LogW ./ sum(p_LogW);  

    W = 10.^LogW_troncati;

    out = [out; repmat(Mw(i), length(W), 1), W', p_LogW'];
end


%% FIGURE W DISTRIBUTION

% Mag_fig1 = 6;
% indice = (out(:,1) == Mag_fig1);
% select = out(indice, :);
% figure(1);
% plot(select(:, 2), select(:, 3), '-o', 'Color', [0, 0, 0]);
% xlabel('W [km]'); 
% ylabel('Probability (%)'); 
% grid on
% titolo = sprintf('RUPTURE WIDTH DISTRIBUTION FOR MAGNITUDE %d', Mag_fig1);
% title(titolo);

%%                            AVERAGE DIP  

dip = (dip_mu - eps_dip):(dip_mu + eps_dip);
dip_dist = makedist('normal','mu',dip_mu,'sigma',dip_sigma);
DIP = dip_dist.pdf(dip)./sum(dip_dist.pdf(dip));

% -----------------------------------------------FIGURE DIP DISTRIBUTION 
% figure(2);
% 
% plot(dip, DIP, '-o', 'Color', [0, 0, 0]);
% xlabel('DIP [deg]'); 
% ylabel('Probability (%)'); 
% grid on
% title('Probability Distribution of Dip Angle for all magnitude');
%------------------------------------------------------------------------

out1 = [];
out1 = repelem(out,length(DIP),1);
DIP2 = repmat(dip',size(out,1),1);
W_z = out1(:,2).* sind(DIP2); 
p_Dip = repmat(DIP',size(out,1),1);

out1(:,end+1) = DIP2;
out1(:,end+1) = p_Dip;
out1(:,end+1) = W_z;    


%%             POBABILITY DISTRIBUTION for uniform HDR:

r = 0.1:0.1:0.9;
HDR_unif = ones(1,length(r))./ (length(r));

out2 = [];
out2 = repelem(out1,length(HDR_unif),1);
R = repmat(r',size(out1,1),1);

HDR = repmat(r',size(out1,1),1);
p_HDR_unif = repmat(HDR_unif',size(out1,1),1);
Wtop = out2(:,6).* R;

out2(:,end+1) = HDR;
out2(:,end+1) = p_HDR_unif; 
out2(:,end+1) = Wtop; 

%%                  S E I M O G E N I C   D E P T H 

z = (Zs_mu - eps_Z):0.5:(Zs_mu + eps_Z);

Zseismopdf = makedist('normal','mu',Zs_mu,'sigma',Zs_sigma);
pesi = Zseismopdf.pdf(z)./sum(Zseismopdf.pdf(z));

% -------------------------------------------- FIGURE Zseismo DISTRIBUTION
% figure(3);
% 
% plot(z, pesi, '-o', 'Color', [0, 0, 0]);
% xlabel('Zseismo [km]'); 
% ylabel('Probability (%)'); 
% grid on
% title('Probability Distribution of Seismogenic Depth for all magnitude');
%-------------------------------------------------------------------------

out3 = [];
out3 = repelem(out2,length(pesi),1);

Zseismo = repmat(z',size(out2,1),1);
p_Zseismo = repmat(pesi',size(out2,1),1);

out3(:,end+1) = Zseismo; 
out3(:,end+1) = p_Zseismo;


%%                  POBABILITY DISTRIBUTION for HDD:

r2 = 0.1:0.1:0.9;

HDD_emp = makedist('normal','mu',HDD_mu,'sigma',HDD_sigma);
HDD_pesi = HDD_emp.pdf(r2)./sum(HDD_emp.pdf(r2));

% %----------------------------------------------- FIGURE Zhypo DISTRIBUTION
% figure(4);
% 
% plot(r2, HDD_pesi, '-o', 'Color', [0, 0, 0]);
% xlabel('input HDD [Zhypo / Zseismo]'); 
% ylabel('Probability (%)'); 
% xlim([0, 1]);
% grid on
% title('Probability Distribution of Hypocentral Depth for all magnitude');
%-------------------------------------------------------------------------

out4 = [];
out4 = repelem(out3,length(HDD_pesi),1);
R2 = repmat(r2',size(out3,1),1);

HDD = repmat(r2',size(out3,1),1);
p_HDD = repmat(HDD_pesi',size(out3,1),1);
Zhypo = out4(:,10) .* R2;

out4(:,end+1) = HDD;
out4(:,end+1) = p_HDD;
out4(:,end+1) = Zhypo;


%%                    R E A L L O C A T I O N

NEW_HDR = out4(:, 7);           out4(:, end+1) = NEW_HDR; 

% idx_1) ∀ Wz >= Zseismo  
idx1 = out4(:, 6) >= out4(:, 10);
% we assume Wz = Zseismo
NEW_Wz = out4(:, 6);
NEW_Wz(idx1, 1) = out4(idx1, 10);       
out4(:, end+1) = NEW_Wz;
% & HDR = Zhypo/Zseismo = HDD
out4(idx1, 15) = out4(idx1, 12);


% idx_2) ∀ rupture violating the LOWER BOUNDARY LAYER 
idx_2 = (out4(:, 6) - out4(:, 9)) >= (out4(:, 10) - out4(:, 14));
idx2 = idx_2 & ~idx1;
% adjusted HDR:
out4(idx2, 15) = 1 - ((out4(idx2, 10) - out4(idx2, 14)) ./ out4(idx2, 6));


% idx_3) ∀ rupture violating the UPPER BOUNDARY LAYER 
idx_3 = out4(:, 9) >= out4(:, 14);
idx3 = idx_3 &  ~idx1;
% adjusted HDR:
out4(idx3, 15) = out4(idx3, 14) ./ out4(idx3, 6);


%% Adjusted Wtop

NEW_Wtop = out4(:,16).* out4(:,15);
out4(:,end+1) = NEW_Wtop; 

% figure(5)
% histogram(out4(:, 15), linspace(-0.05, 1.05, 12), 'FaceColor', 'yellow','Normalization', 'pdf');
% xticks(0:0.1:1);
% 
% ylabel('Probability Density Function');
% xlabel('nucleation position in the rupture - HDR');
% title('Recomputed HDR - normal expl ');
% grid on;


%%                   E X P L O R A T I O N    T R E E

% out4(:,1)  = mag;        out4(:,2)  = W;             out4(:,3)  = P(W)
% out4(:,4)  = dip;        out4(:,5)  = P(dip);        out4(:,6)  = W_z
% out4(:,7)  = HDR;        out4(:,8)  = P(HDR);        out4(:,9)  = W_top
% out4(:,10) = Zseismo;    out4(:,11) = P(Zseismo);    out4(:,12) = HDD
% out4(:,13) = P(HDD);     out4(:,14) = Zhypo;         out4(:,15) = NEW_HDR
% out4(:,16) = NEW_Wz;     out4(:,17) = NEW_Wtop

tab = table(out4(:,1), out4(:,2), out4(:,3), out4(:,4), out4(:,5), out4(:,6), ...
                out4(:,7), out4(:,8), out4(:,9), out4(:,10), out4(:,11), out4(:,12), ...
                out4(:,13), out4(:,14), out4(:,15), out4(:,16), out4(:,17), ...
                'VariableNames', {'mag', 'W', 'P_W', 'dip', 'P_dip', 'W_z', ...
                                   'HDR', 'P_HDR', 'W_top', 'Zseismo', 'P_Zseismo', ...
                                   'HDD', 'P_HDD', 'Zhypo', 'NEW_HDR', 'New_Wz' , 'New_Wtop'});

%% surface rupture condition (SRC) : if Wtop >= Zhypo --> yes rupture
SRC = (out4(:,17) >= out4(:,14));

outcome = zeros(size(out4, 1), 1);
outcome(SRC) = 1;       
outcome(~SRC) = 0;

%% VERIFICA 

ind_0 = find(outcome == 0);
outcome_0 = out4(ind_0, [14, 17]); % [Zhypo, New Wtop]

diff_values = outcome_0(:, 1) - outcome_0(:, 2);
outcome_0 = [outcome_0, diff_values];

if any(diff_values < 0); disp('In outcome(0) ci sono casi di rottura. Max value ='); disp(min(diff_values)); end

ind_1 = find(outcome == 1);
outcome_1 = out4(ind_0, [14, 17]); % [Zhypo, New Wtop]

diff_values = outcome_1(:, 1) - outcome_1(:, 2);
outcome_1 = [outcome_1, diff_values];

if any(diff_values < 0); disp('In outcome(1) ci sono casi di non rottura. Max value ='); disp(max(diff_values)); end


%%        CONDITIONAL PROBABILITY OF SURFACE RUPTURE


P = zeros(1,length(Mw));
for i = 1:length(Mw)  
    g = find(out4(:,1) == Mw(i) & outcome(:,1) == 1);
    P(1, i) = sum(out4(g, 3) .* out4(g, 5) .* out4(g, 8) .* out4(g, 11) .* out4(g,13));

end

%%
figure(6)

hold on
plot(Mw, P, '-k', 'LineWidth', 1);
grid on

xlabel('M_w')
ylabel('CPSR')
ylim([0 1])
xticks(min(Mw):0.5:max(Mw))



