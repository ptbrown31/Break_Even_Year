close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% decide on type of experiment

time_desc = 2015:5:2335;

T = 65;

global abatement_cost_on
abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 

global Burke_damage_on;
Burke_damage_on = 0;  %1 means Burke damage through TFP and capital depreciation pathways, 
                      %0 means default DICE, 
                      %2 means no damage

global optimize_on;
optimize_on = 1;  %1 means optimizing, 0 means its reading in a control rate

%     if optimize_on == 0;
%         
         global aopt
%         aopt = ones(T,1);
%         
%     end

global optimize_for_total_utility;
optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility

global utility_function;
utility_function = 1; %1 means utility function in code and 0 means utility function in text,


%% assign paramater values

range_values_2nd_half = linspace(1,2,6);
range_values_1st_half = 1./(range_values_2nd_half(2:end));

range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);

%accordian the variables

%damage function
damage_dep_rate_coeffs_on_T = range_values*0.013;
damage_TFP_growth_coeffs_on_T = range_values*0.0055;

%Elasticity of marginal utility 
%of consumption
elasmu_vals = range_values*1.45;         % 3
elasmu_vals(elasmu_vals <= 1) = 1.001;   

% carbon intensity of output
Sigg0_vals = range_values*-0.0152;         % 4
Siggd_vals = range_values*-0.001;        % 5

%initial cost decline backstop 
%cost per period
pd_vals = range_values*0.025;            % 6

%abatement cost curve exponent
theta2_vals = range_values*2.6;          % 7

%asymptotic population
LA_vals = range_values*11500;            % 8

%total factor productivity
Ag0_vals = range_values*0.076;           % 9
Agd_vals = range_values*0.005;           % 10

%pure rate of time preference
prstp_vals = range_values*0.015;         % 11

%climate sensitivity
deltaT_vals = range_values*3.1;          % 12

%damage convexity

sai3_vals = range_values*2;


num_vars = 12;

var_names = {'damage cap-dep-rate coeff on T',...
             'damage TFP coeff on T',...
             'elasticity of utility',...
             'C intensity, growth rate',...
             'C intensity, change in growth rate',...
             'backstop cost decline',...
             'abatement cost curve exponent',...
             'asymptotic population',...
             'TFP, growth rate',...
             'TFP, change in growth rate',...
             'pure rate of time preference',...
             'climate sensitivity'};

global damage_dep_rate_coeff_on_T
global damage_TFP_growth_coeff_on_T
global elasmu
global Sigg0
global Siggd
global pd
global theta2
global LA
global Ag0
global Agd
global prstp
global deltaT

global damage_scalar
global mitigation_scalar

global sai3

middle_value_i = 6;

        %set everything to middle
        
                           damage_dep_rate_coeff_on_T = damage_dep_rate_coeffs_on_T(middle_value_i);
                           damage_TFP_growth_coeff_on_T = damage_TFP_growth_coeffs_on_T(middle_value_i);
                           elasmu = elasmu_vals(middle_value_i);
                           Sigg0 = Sigg0_vals(middle_value_i);
                           Siggd = Siggd_vals(middle_value_i);
                           pd = pd_vals(middle_value_i);      
                           theta2 = theta2_vals(middle_value_i);
                           LA = LA_vals(middle_value_i);    
                            Ag0 = Ag0_vals(middle_value_i);
                            Agd = Agd_vals(middle_value_i);    
                           deltaT = deltaT_vals(middle_value_i);
                           
                           prstp = prstp_vals(middle_value_i);
                           % prstp = prstp_vals(end);
                           
                           damage_scalar = range_values(middle_value_i);
                            %damage_scalar = range_values(end);
                            
                           mitigation_scalar = range_values(middle_value_i);
                            %mitigation_scalar = range_values(end);
                           
                           sai3 = sai3_vals(middle_value_i);
                               
                                %run DICE
                                addpath(genpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth'))

                                DICE;

                                NPV_utility = fval;

                                %Calculate the Carbon Price
                                opt_price = pbacktime.*(aopt'./pai(1:length(aopt))).^(theta2-1); % Carbon price, $ per ton of CO2

                                overall_plot_series = cat(1,S(16,:),...
                                                    S(17,:),...
                                                    S(7,:),...
                                                    S(19,:),...
                                                    S(9,:),...
                                                    Sig,...
                                                    L,...
                                                    R,...
                                                    S(20,:),...
                                                    S(12,:),...
                                                    S(4,:),...
                                                    S(2,:),...
                                                    S(14,:),...
                                                    aopt',...
                                                    opt_price,...
                                                    S(10,:));
                                                
                                                overall_plot_series_all = overall_plot_series;
                                                
    % S
    % 1 - regular capital X
    % 2 - GMST X
    % 3 - lower ocean temperature X
    % 4 - atmospheric CO2 X
    % 5 - upper ocean CO2 X
    % 6 - lower ocean CO2 X
    % 7 - gross output X
    % 8 - zero-carbon emissions captial (leave NaN)
    % 9 - total capital (same as 1)
    % 10 - abatement cost
    % 11 - epsilon (leave NaN)
    % 12 - Industrial CO2 emissions
    % 13 - CO2 Radiative Forcing
    % 14 - Climate Damage
    % 15 - 
    % 16 - total utility
    % 17 - per capita utility
    % 18 - consumption
    % 19 - consumption per capita
    % 20 - net output
    % 21 - deltaK
    % 22 - TFP growth
    % 23 - TFP
      
% Net output: 9 ($ trillion 2005 USD)
% Population: 7 (million people)
% Output/capita = $ million/person
   
per_capita_GDP_million_per_person = overall_plot_series_all(9,:,:)./overall_plot_series_all(7,:,:);

%convert to thousand / person

per_capita_GDP_thousand_per_person = (1000).*per_capita_GDP_million_per_person;

overall_plot_series_all = cat(1,overall_plot_series_all,per_capita_GDP_thousand_per_person);

titles = {'total utility',...
          'per-capita utility',...
          'gross output',...
          'per-capita consumption',...
          'total capital',...
          'CO_2 intensity of output',...
          'population',...
          'discount factor',...
          'net output (Y_G)',...
          'industrial CO_2 emissions',...
          'ATM CO_2',...
          'global temperature',...
          'damage from warming fraction',...
          'control rate',...
          'carbon tax',...
          'abatement cost fraction',...
          'per-capita GWP'};
      
ylabels = {'utils',...
          'utils / person',...
          '$ trillion 2010 USD',...
          '$ thousand 2010 USD / person',...
          '$ trillion 2010 USD',...
          'MTC/$1000',...
          'million people',...
          'fraction',...
          '$ trillion 2010 USD',...
          'Gt CO_2 / year',...
          'Gt CO_2',...
          '\Delta T above PI',...
          'frac output w.r.t. no warming',...
          'fraction',...
          '$1000 / Ton CO_2',...        
          'fract output w.r.t. no abatement',...
          '$ thousand 2010 USD / person'};
      
overall_plot_series_all_optimal = overall_plot_series_all;

%% plot 

row_num = 4;
col_num = 5;

end_year_plot = 2150;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',6);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
for var_i = 1:length(titles)
    
    subplot(row_num,col_num,var_i)
    hold on
    grid on
    
            plot(time_desc,squeeze(overall_plot_series_all_optimal(var_i,:)),'LineWidth',0.5)
                      
        %xlim([2010 end_year_plot])    

        xlim([2015 end_year_plot])

        xlabel('time')
        ylabel(ylabels{var_i})
        title(titles{var_i})

end

%% run again with no mitigation

%load DICE-default

abatement_cost_on = 1;
Burke_damage_on = 0;
optimize_on = 0;
optimize_for_total_utility = 1;
utility_function = 1;

aopt = zeros(T,1);

                                %run DICE
                                addpath(genpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth'))

                                DICE;

                                NPV_utility = fval;

                                %Calculate the Carbon Price
                                opt_price = pbacktime.*(aopt'./pai(1:length(aopt))).^(theta2-1); % Carbon price, $ per ton of CO2

                                overall_plot_series = cat(1,S(16,:),...
                                                    S(17,:),...
                                                    S(7,:),...
                                                    S(19,:),...
                                                    S(9,:),...
                                                    Sig,...
                                                    L,...
                                                    R,...
                                                    S(20,:),...
                                                    S(12,:),...
                                                    S(4,:),...
                                                    S(2,:),...
                                                    S(14,:),...
                                                    aopt',...
                                                    opt_price,...
                                                    S(10,:));
                                                
                                                overall_plot_series_all = overall_plot_series;
   
per_capita_GDP_million_per_person = overall_plot_series_all(9,:,:)./overall_plot_series_all(7,:,:);

%convert to thousand / person

per_capita_GDP_thousand_per_person = (1000).*per_capita_GDP_million_per_person;

overall_plot_series_all = cat(1,overall_plot_series_all,per_capita_GDP_thousand_per_person);

      
overall_plot_series_all_no_mitigation = overall_plot_series_all;
      
%% optpion to overwrite the optimal with one of the idealized linear control rates
      
% make schematic

plot_var = 4;
%plot_var = 2;

difference_optimal = overall_plot_series_all_optimal - overall_plot_series_all_no_mitigation;
ratio_optimal = overall_plot_series_all_optimal./overall_plot_series_all_no_mitigation;


difference_optimal_cumsum = cumsum(squeeze(difference_optimal(plot_var,2:end)));

%calculate the break-even year

break_even_year_optimal = [];
payback_period_year_optimal = [];

for yeari = 1:length(time_desc)-1

   if difference_optimal(plot_var,yeari) < 0 && ...
      difference_optimal(plot_var,yeari+1) >= 0

       break_even_year_optimal = time_desc(yeari)+2.5;

   end
   if difference_optimal_cumsum(yeari) < 0 && ...
      difference_optimal_cumsum(yeari+1) >= 0

       payback_period_year_optimal = time_desc(yeari)+2.5;

   end   
end

break_even_year_optimal
payback_period_year_optimal

year_max = 2150;

% ymaxes = [50 3.5 300 1];
% ymins = [9 0 -50 0];

ymaxes = [6 9 400 2];
ymins = [0 -0.3 -75 0];

%  ymaxes = [6 9 3 2];
%  ymins = [0 -0.3 -3 0];

FigHandle = figure('Position', [100, 100, 1100, 800]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',6);
set(0,'defaultAxesFontName', 'helvetica')

% subplot(3,1,1)
% hold on
% 
% %     plot(time_desc,squeeze(overall_plot_series_all_default_DICE(plot_var,:,:)),'-m','LineWidth',1)
% %     plot(time_desc,squeeze(overall_plot_series_all_default_DICE(plot_var,:,end)),'-m','LineWidth',2)
% 
% %     plot(time_desc,overall_plot_series_all_default_DICE_no_mitigation(plot_var,:),'-k','LineWidth',2)
% %     plot(time_desc,overall_plot_series_all(plot_var,:),'-b','LineWidth',2)
% 
%     plot(time_desc,overall_plot_series_all_no_mitigation(12,:),'-k','LineWidth',2)
%     plot(time_desc,overall_plot_series_all_optimal(12,:),'-b','LineWidth',2)
% 
% %     if ~isempty(break_even_year_optimal)
% %         plot([break_even_year_optimal break_even_year_optimal],[ymins(1) ymaxes(1)],'--b')
% %     end     
% %     if ~isempty(break_even_year_target)
% %         plot([break_even_year_target break_even_year_target],[ymins(1) ymaxes(1)],'--m')
% %     end
%     
%     legend('no policy',...
%            'optimal policy',...
%            'location','NorthWest')
%            
%     %title(strcat('a, ',titles{plot_var}))
%     %title(strcat('a, per-capita consumption'))
%     title(strcat('a, global temperature')) 
%     
%     %ylabel(ylabels{plot_var})
%     %ylabel('thousand 2010 US$ / person')
%     ylabel('degrees C')
%     
%     xlim([2015 year_max])
%     ylim([ymins(1) ymaxes(1)])
    
subplot(3,2,1)
hold on

damagage_frac_var = 13;
abatement_frac_var = 16;

    plot(time_desc,100*(1-overall_plot_series_all_no_mitigation(damagage_frac_var,:)),'-r','LineWidth',2)
    plot(time_desc,100*(1-overall_plot_series_all_optimal(damagage_frac_var,:)),'--r','LineWidth',2)

    plot(time_desc,100*(overall_plot_series_all_no_mitigation(abatement_frac_var,:)),'-m','LineWidth',2)    
    plot(time_desc,100*(overall_plot_series_all_optimal(abatement_frac_var,:)),'--m','LineWidth',2)
    
        legend('no policy, damage cost',...
               'optimal policy, damage cost',...
               'no policy, mitigation cost',...
               'optimal policy, mitigation cost',...
               'location','NorthWest',...
               'AutoUpdate','off')    
    
%     if ~isempty(break_even_year_optimal)
%         plot([break_even_year_optimal break_even_year_optimal],[ymins(2) ymaxes(2)],'--b')
%     end    
%     if ~isempty(break_even_year_target)
%         plot([break_even_year_target break_even_year_target],[ymins(2) ymaxes(2)],'--m')
%     end    

    
    title(strcat('a, cost of mitigation and cost of damages with and without policy'))
    ylabel('% output lost')
    
    xlim([2015 year_max])
    ylim([ymins(2) ymaxes(2)])
    
subplot(3,2,3)
hold on
    
    plot(time_desc,100*difference_optimal(damagage_frac_var,:),'-g','LineWidth',2)
    plot(time_desc,100*difference_optimal(abatement_frac_var,:),'-r','LineWidth',2)
    
    legend('benefits of avoided damages from mitigation',...
           'costs of policy',...      
           'location','NorthWest',...
           'AutoUpdate','off')    

%     if ~isempty(break_even_year_optimal)
%         plot([break_even_year_optimal break_even_year_optimal],[ymins(4) ymaxes(4)],'--b')
%     end
    if ~isempty(payback_period_year_optimal)
        plot([payback_period_year_optimal payback_period_year_optimal],[ymins(3) ymaxes(3)],'--b')
    end
    
    xlim([2015 year_max])
    title('b, Influence of policy on costs of mitigation and damages')
    ylabel('% of output')    
    
    ylim([ymins(4) ymaxes(4)])    
        
subplot(3,2,5)
hold on

    plot(time_desc,1000*difference_optimal(plot_var,:),'-k','LineWidth',2)
    
    legend('optimal policy w.r.t. no policy',...
           'location','NorthWest',...
           'AutoUpdate','off')
       
    plot([2015 year_max],[0 0],'-k','LineWidth',1)       
    
%     if ~isempty(break_even_year_optimal)
%         plot([break_even_year_optimal break_even_year_optimal],[ymins(3) ymaxes(3)],'--b')
%     end
    if ~isempty(payback_period_year_optimal)
        plot([payback_period_year_optimal payback_period_year_optimal],[ymins(3) ymaxes(3)],'--b')
    end
    
    xlim([2015 year_max])
    title('c, Influence of policy on per-capita consumption')
    ylabel('2010 US$ / person / year')
    
    ylim([ymins(3) ymaxes(3)])
    
    xlabel('year')
    
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/break_even_year_fig_1','-dpdf','-bestfit')

    


%% make schematic - Ken version

close all

plot_var = 4;
%plot_var = 2;

difference_optimal = overall_plot_series_all_optimal - overall_plot_series_all_no_mitigation;

difference_optimal_cumsum = cumsum(squeeze(difference_optimal(plot_var,2:end)));

%calculate the break-even year

break_even_year_optimal = [];
payback_period_year_optimal = [];

for yeari = 1:length(time_desc)-1

   if difference_optimal(plot_var,yeari) < 0 && ...
      difference_optimal(plot_var,yeari+1) >= 0

       break_even_year_optimal = time_desc(yeari)+2.5;

   end
   if difference_optimal_cumsum(yeari) < 0 && ...
      difference_optimal_cumsum(yeari+1) >= 0

       payback_period_year_optimal = time_desc(yeari)+2.5;

   end   
end

break_even_year_optimal
payback_period_year_optimal

year_max = 2150;

% ymaxes = [50 3.5 300 1];
% ymins = [9 0 -50 0];

ymaxes = [6 9 400 2];
ymins = [0 -0.3 -75 0];

%  ymaxes = [6 9 3 2];
%  ymins = [0 -0.3 -3 0];

FigHandle = figure('Position', [100, 100, 1100, 800]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',7);
set(0,'defaultAxesFontName', 'helvetica')

% subplot(3,1,1)
% hold on
% 
% %     plot(time_desc,squeeze(overall_plot_series_all_default_DICE(plot_var,:,:)),'-m','LineWidth',1)
% %     plot(time_desc,squeeze(overall_plot_series_all_default_DICE(plot_var,:,end)),'-m','LineWidth',2)
% 
% %     plot(time_desc,overall_plot_series_all_default_DICE_no_mitigation(plot_var,:),'-k','LineWidth',2)
% %     plot(time_desc,overall_plot_series_all(plot_var,:),'-b','LineWidth',2)
% 
%     plot(time_desc,overall_plot_series_all_no_mitigation(12,:),'-k','LineWidth',2)
%     plot(time_desc,overall_plot_series_all_optimal(12,:),'-b','LineWidth',2)
% 
% %     if ~isempty(break_even_year_optimal)
% %         plot([break_even_year_optimal break_even_year_optimal],[ymins(1) ymaxes(1)],'--b')
% %     end     
% %     if ~isempty(break_even_year_target)
% %         plot([break_even_year_target break_even_year_target],[ymins(1) ymaxes(1)],'--m')
% %     end
%     
%     legend('no policy',...
%            'optimal policy',...
%            'location','NorthWest')
%            
%     %title(strcat('a, ',titles{plot_var}))
%     %title(strcat('a, per-capita consumption'))
%     title(strcat('a, global temperature')) 
%     
%     %ylabel(ylabels{plot_var})
%     %ylabel('thousand 2010 US$ / person')
%     ylabel('degrees C')
%     
%     xlim([2015 year_max])
%     ylim([ymins(1) ymaxes(1)])
    
subplot(3,2,1)
hold on

damagage_frac_var = 13;
abatement_frac_var = 16;

    plot(time_desc,100*(1-overall_plot_series_all_no_mitigation(damagage_frac_var,:)),'-b','LineWidth',1.5)
    plot(time_desc,100*(1-overall_plot_series_all_optimal(damagage_frac_var,:)),'--b','LineWidth',1.5)

    plot(time_desc,100*(overall_plot_series_all_no_mitigation(abatement_frac_var,:)),'-r','LineWidth',1.5)    
    plot(time_desc,100*(overall_plot_series_all_optimal(abatement_frac_var,:)),'--r','LineWidth',1.5)
    
        legend('no policy, damage cost',...
               'optimal policy, damage cost',...
               'no policy, mitigation cost',...
               'optimal policy, mitigation cost',...
               'location','NorthWest',...
               'AutoUpdate','off')    
    
%     if ~isempty(break_even_year_optimal)
%         plot([break_even_year_optimal break_even_year_optimal],[ymins(2) ymaxes(2)],'--b')
%     end    
%     if ~isempty(break_even_year_target)
%         plot([break_even_year_target break_even_year_target],[ymins(2) ymaxes(2)],'--m')
%     end    

    
    title(strcat('a, cost of mitigation and cost of damages with and without policy'))
    ylabel('% output lost')
    
    xlim([2015 year_max])
    ylim([ymins(2) ymaxes(2)])
    
subplot(3,2,3)
hold on
    
    plot(time_desc,100*difference_optimal(damagage_frac_var,:),'-b','LineWidth',0.1)
    plot(time_desc,100*difference_optimal(abatement_frac_var,:),'-r','LineWidth',0.1)
    
    legend('benefits of avoided damages from mitigation',...
           'costs of mitigation',...      
           'location','NorthWest',...
           'AutoUpdate','off')    

%     if ~isempty(break_even_year_optimal)
%         plot([break_even_year_optimal break_even_year_optimal],[ymins(4) ymaxes(4)],'--b')
%     end
%     if ~isempty(payback_period_year_optimal)
%         plot([payback_period_year_optimal payback_period_year_optimal],[ymins(3) ymaxes(3)],'--b')
%     end
    
    xlim([2015 year_max])
    title('b, Influence of policy on costs of mitigation and damages')
    ylabel('% of output')    
    
    ylim([ymins(4) ymaxes(4)])    
        
subplot(3,2,5)
hold on

    plot(time_desc,1000*difference_optimal(plot_var,:),'-k','LineWidth',2)
    
%     legend('optimal policy w.r.t. no policy',...
%            'location','NorthWest',...
%            'AutoUpdate','off')
       
    plot([2015 year_max],[0 0],'-k','LineWidth',1)       
    
%     if ~isempty(break_even_year_optimal)
%         plot([break_even_year_optimal break_even_year_optimal],[ymins(3) ymaxes(3)],'--b')
%     end
%     if ~isempty(payback_period_year_optimal)
%         plot([payback_period_year_optimal payback_period_year_optimal],[ymins(3) ymaxes(3)],'--b')
%     end
    
    xlim([2015 year_max])
    title('c, Influence of policy on per-capita consumption')
    ylabel('2010 US$ per person per year')
    
    ylim([ymins(3) ymaxes(3)])
    
    xlabel('year')
    
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/break_even_year_fig_1_KC_version','-dpdf','-bestfit')

%% ratio

% FigHandle = figure('Position', [100, 100, 700, 600]); %[left bottom width height]
% set(gcf,'color',[1 1 1]);
% set(0, 'DefaultAxesFontSize',16);
% set(0,'defaultAxesFontName', 'helvetica')
% hold on
% 
%     semilogy(time_desc(1:25),ratio_optimal(plot_var,1:25))
% 
%     xlim([2015 year_max])
%     title('consumption ratio')
%     ylabel('fracto')

%% life sampling calculation

consumption_diff_t_series = 1000*difference_optimal(plot_var,2:end); %get rid of NaN
consumption_ratio_t_series = ratio_optimal(plot_var,2:end); %get rid of NaN


time_desc_2 = time_desc(2:end);

life_time_lengths_years = 5:5:120; %this is the "lifetime" index

life_time_lengths_indicies = life_time_lengths_years/5;

birth_years = 1920:5:2100;

num_negative_years = time_desc_2(1) - birth_years(1);

num_negative_indicies = num_negative_years/5;

extended_consumption_diff_t_series = horzcat(zeros(1,num_negative_indicies),consumption_diff_t_series);

extended_time_desc = birth_years(1):5:time_desc_2(end);

personal_sampling_discounted = NaN(length(birth_years),length(life_time_lengths_years));
 
 for lifetime_window_lengths_i = 1:length(life_time_lengths_years)
     for birth_year_i = 1:length(birth_years)
         
         first_year_sample = birth_years(birth_year_i);
         last_year_sample = first_year_sample + life_time_lengths_years(lifetime_window_lengths_i);
         
         first_year_sample_i = find(extended_time_desc == first_year_sample);
         last_year_sample_i = find(extended_time_desc == last_year_sample);
         
         personal_sampling_now_non_discounted = extended_consumption_diff_t_series(first_year_sample_i:last_year_sample_i);
         
         personal_sampling_now_discounted = NaN(length(personal_sampling_now_non_discounted),1);
         
         discount_rate_per_year = 0.05;
         
         for time_i = 1:length(personal_sampling_now_non_discounted)
             
             personal_sampling_now_non_discounted_per_5_years = personal_sampling_now_non_discounted*5; %this is a consumption difference expressed every 5 years but it is per year to need to multiply by 5 to get total consumption difference
             
             total_years = time_i*5 - 5;
             
             personal_sampling_now_discounted(time_i) = personal_sampling_now_non_discounted_per_5_years(time_i).*(1./((1 + discount_rate_per_year)^total_years));
             
         end
         
         personal_sampling_discounted(birth_year_i,lifetime_window_lengths_i) = sum(personal_sampling_now_discounted);
     end
 end

%% plot contour

FigHandle = figure('Position', [100, 100, 700, 600]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',13);
set(0,'defaultAxesFontName', 'helvetica')
hold on

num_contours = 30;
red_vector = horzcat(linspace(1,1,num_contours/2),linspace(1,0,num_contours/2));
green_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,1,num_contours/2));
blue_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,0,num_contours/2));

rgb2 = horzcat(red_vector',green_vector',blue_vector');
colormap(rgb2)
%colormap(winter)

%colormap(redblue)

frac_contours_2nd_half = linspace(1,1.05,20);
frac_contours_1st_half = 1./(frac_contours_2nd_half(2:end));

%contour_range_values = horzcat(fliplr(frac_contours_1st_half),frac_contours_2nd_half);

%contour_range_values = [0.9984 0.9986 0.9988 0.9990 0.9992 0.9994 0.9996 0.9998 1 1.0002 1.0004 1.0006 1.0008 1.0008 1.001 1.0012 1.0014];
%contour_range_values = [0.9984 0.9986 0.9988 0.9990 0.9992 0.9994 0.9996 0.9998 1 1.0002 1.001 1.005 1.01 1.05 1.1 1.15 1.2];

contour_range_values = [-10000, -3000, -1000, -300, -100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000];
% contour_range_values = [-100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100]

%     max_color_range = 5000;
%     min_color_range = -5000;
%      max_color_range = max(contour_range_values);
%      min_color_range = min(contour_range_values);

      max_color_range = 0.8*max(max(abs(personal_sampling_discounted)));
      min_color_range = -1*max_color_range;

      %contour_range_values = linspace(min_color_range,max_color_range,num_contours)';


    %contourf(birth_years,life_time_lengths_years,personal_sampling',contour_range_values, 'LineStyle','none')
    contourf(birth_years,life_time_lengths_years,personal_sampling_discounted',contour_range_values, 'LineStyle','none')
    
    %contour(birth_years,life_time_lengths_years,personal_sampling_fraction',[1.0001 1.00001],'k','LineStyle','-','LineWidth',4);
    contour(birth_years,life_time_lengths_years,personal_sampling_discounted',[0.000 0.000],'k','LineStyle','-','LineWidth',4);

    caxis([min_color_range max_color_range])

    xlabel('Birth Year')
    ylabel('Life Length (years)')
    title('Mitigation Effect on Discounted Cumulative Lifetime Consumption')
    
    h = colorbar;
    %h.Ruler.Scale = 'log';
    
    %set(gca,'colorscale','log')
    
    ylabel(h,'2010 US$ per person');
    %ylabel(h,'log(fraction of no-policy consumption)');

    
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/lifetime_sample','-dpdf','-bestfit')
    
%% plot contour

plot_matrix = personal_sampling_discounted';

plot_matrix_pos = NaN(size(plot_matrix));
plot_matrix_neg = NaN(size(plot_matrix));

for i = 1:size(plot_matrix,1)
    for j = 1:size(plot_matrix,2)
        
        if plot_matrix(i,j) > 0
            
            plot_matrix_pos(i,j) = plot_matrix(i,j);
            
        end
        if plot_matrix(i,j) <= 0
            
            plot_matrix_neg(i,j) = plot_matrix(i,j);
            
        end        
    end
end

contour_range_values_neg = linspace(-1000,1000,20);
contour_range_values_pos = linspace(-20000,20000,20);

% contour_range_values_pos = [-25000 -10000,-7500,-5000,-2500,-1000,-750,-500,-250,-100,-75,-50,-25,-10,-7.5,-5,-2.5,-1,0,...
%                             1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000,2500,5000,7500,10000,25000];
% contour_range_values_neg = contour_range_values_pos(6:end-6+1);


FigHandle = figure('Position', [100, 100, 700, 600]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',13);
set(0,'defaultAxesFontName', 'helvetica')
hold on

num_contours = 30;
red_vector = horzcat(linspace(1,1,num_contours/2),linspace(1,0,num_contours/2));
blue_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,1,num_contours/2));
green_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,0,num_contours/2));

rgb3 = horzcat(red_vector',green_vector',blue_vector');
colormap(rgb3)

%colormap(redblue)

    contourf(birth_years,life_time_lengths_years,plot_matrix_pos,contour_range_values_pos, 'LineStyle','none')
    %sanePColor(birth_years,life_time_lengths_years,plot_matrix_pos)
    
        caxis([min(contour_range_values_pos) max(contour_range_values_pos)])
                
    xlabel('Birth Year')
    ylabel('Life Length (years)')
%    title('Mitigation impact on Avgerage Lifetime Cumulative Consumption')
    
    h = colorbar;
%    h.Ticks = contour_range_values_pos;
%     h.Limits = [min(contour_range_values_pos) max(contour_range_values_pos)];
    ylabel(h,'lifetime cumulative consumption difference');
        
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/lifetime_sample_pos','-dpdf','-bestfit')    
            
FigHandle = figure('Position', [100, 100, 700, 600]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',13);
set(0,'defaultAxesFontName', 'helvetica')
hold on       

colormap(rgb2)

%colormap(redblue)

    contourf(birth_years,life_time_lengths_years,plot_matrix_neg,contour_range_values_neg, 'LineStyle','none')
    %sanePColor(birth_years,life_time_lengths_years,plot_matrix_neg)    
    
        caxis([min(contour_range_values_neg) max(contour_range_values_neg)])

    contour(birth_years,life_time_lengths_years,personal_sampling_discounted',[0.000 0.000],'k','LineStyle','-','LineWidth',3);


    xlabel('Birth Year')
    ylabel('Life Length (years)')
    %title('Mitigation impact on Avgerage Lifetime Cumulative Consumption')
    
    h = colorbar;
    %h.Ruler.Scale = 'log';
    
    %set(gca,'colorscale','log')
    
    %ylabel(h,'2010 US$ / person / year');
    ylabel(h,'lifetime cumulative consumption difference');
    
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/lifetime_sample_neg','-dpdf','-bestfit')

% two_d_plot = optimim_level_of_mitigation_stab_temp;
% x_axis = time_horzs;
% y_axis = discount_rates;
% 
% FigHandle = figure('Position', [100, 100, 700, 600]); %[left bottom width height]
% set(gcf,'color',[1 1 1]);
% set(0, 'DefaultAxesFontSize',16);
% set(0,'defaultAxesFontName', 'helvetica')
% hold on
% grid on
% %colormap(redblue)
% %colormap(redgreencmap)
% %colormap(pink)
% 
%                 max_color_range = 4.0;
%                 min_color_range = 1.3;
%                 num_contours = 20;
% 
%                 contourf(x_axis,y_axis,two_d_plot,linspace(min_color_range,max_color_range,num_contours), 'LineStyle','none');
%                 contour(x_axis,y_axis,two_d_plot,[1.5 2 3],'k','LineStyle','-','LineWidth',4);
%                 caxis([min_color_range max_color_range])
%                 
%                 plot([2100 2100],[min(discount_rates) max(discount_rates)],'--k')
%                 plot([2300 2300],[min(discount_rates) max(discount_rates)],'--k')
% 
%                 plot([min(time_horzs) max(time_horzs)],[0.0 0.0],'--k')
%                 plot([min(time_horzs) max(time_horzs)],[0.03 0.03],'--k')
%                 %plot([min(time_horzs) max(time_horzs)],[0.05 0.05],'--k')
%                 %plot([min(time_horzs) max(time_horzs)],[0.07 0.07],'--k')
%                 
%                 %sanePColor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');
%                 %pcolor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');
% 
%                 h = colorbar;
% 
%                 xlim([min(time_horzs) max(time_horzs)])
%                 ylim([min(discount_rates) max(discount_rates)])
%                 
% %                 h.Ruler.Scale = 'log';
% %                 h.Ruler.MinorTick = 'on';
% 
%                 ylabel(h,'Optimal level of mitigation effort (global warming, C)');
%                 
%                 ylabel('Discount rate')
%                 xlabel('Time horizon (year)')
%                 
%                 %title('optimal mitigation level f(discount rate,time horizon)')