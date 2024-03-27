clear
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex'); % For all other text elements
% InitialState = 0; % Start the voting from an unfished state (i.e., underfished)
InitialState = 1; % Start the voting from an open access equilibrium (i.e., overfished)

N = 7; % How many countries are there?

T = round(N/2); % This is the number of votes you need to pass by majority
Results_majority = Solve_for_voting(N,T,InitialState);
disp(['Majority = ' num2str(T) ' votes'])

T = ceil(N*0.75); % This is the number of votes you need to pass by supermajority
Results_supermajority = Solve_for_voting(N,T,InitialState);
disp(['Supermajority = ' num2str(T) ' votes'])

T = N; % This is the number of votes you need to pass by consensus
[Results_consensus,W] = Solve_for_voting(N,T,InitialState);
disp(['Consensus = ' num2str(T) ' votes'])

for k = 1:N
    XT{k} = ['$w_' num2str(k) '= ' num2str(W(k),2) '$'];
end
if InitialState == 0 
    XTL = {'Ex','C','S','M'};
    Filename = '../Manuscript/Figures/Figure_results_Ex.tiff';
else
    XTL = {'OA','C','S','M'};
    Filename = '../Manuscript/Figures/Figure_results_OA.tiff';
end

% Construct a utility matrix
U = [Results_consensus{1,1} Results_consensus{2,1} Results_supermajority{2,1} Results_majority{2,1}];
U = U./repmat(max(U,[],2),1,4);

figure(1), clf; FS = 15;
CL = viridis(5);
subplot(3,1,1), hold on
bb = bar(U); set(bb,'edgecolor','none');
for k = 1:4; set(bb(k),'facecolor',CL(k,:)); end
set(gca,'xtick',[1:N],'xticklabel',XT,'fontsize',FS-2)
ylabel({'National utility','($U_i$)'},'fontsize',FS)
text(1,1.15,{'Economic','focus'},'fontsize',FS-2,'horizontalalignment','center')
text(4,1.15,{'Balanced'},'fontsize',FS-2,'horizontalalignment','center')
text(7,1.15,{'Conservation','focus'},'fontsize',FS-2,'horizontalalignment','center')
% ;'Balanced';'Conservation'}
ylim([0 1])
CornerLetterLabel('(A)',[-0.1 1.05],FS)


subplot(3,1,2), hold on
bb = bar([  Results_consensus{1,2} ...
            Results_consensus{2,2} ...
            Results_supermajority{2,2} ...
            Results_majority{2,2}]); set(bb,'edgecolor','none');
for k = 1:4; set(bb(k),'facecolor',CL(k,:)); end
ylabel({'National harvest','effort ($h_i$)'},'fontsize',FS)
set(gca,'xtick',[1:N],'xticklabel',XT,'fontsize',FS-2,'ytick',[0:0.1:0.5])
CornerLetterLabel('(B)',[-0.1 1.05],FS)

subplot(3,3,7), hold on
bb=bar(1,Results_consensus{1,3},0.8); set(bb,'edgecolor','none','facecolor',CL(1,:));
bb=bar(2,Results_consensus{2,3},0.8); set(bb,'edgecolor','none','facecolor',CL(2,:));
bb=bar(3,Results_supermajority{2,3}, 0.8); set(bb,'edgecolor','none','facecolor',CL(3,:));
bb=bar(4,Results_majority{2,3}, 0.8); set(bb,'edgecolor','none','facecolor',CL(4,:));
ylabel({'Abundance','($N^*$)'},'fontsize',FS)
set(gca,'fontsize',FS-2,'xtick',[1:4],'xticklabel',XTL)
CornerLetterLabel('(C)',[-0.35 1.05],FS)

subplot(3,3,8), hold on
bb=bar(1,sum(Results_consensus{1,2}),0.8); set(bb,'edgecolor','none','facecolor',CL(1,:));
bb=bar(2,sum(Results_consensus{2,2}),0.8); set(bb,'edgecolor','none','facecolor',CL(2,:));
bb=bar(3,sum(Results_supermajority{2,2} ),0.8); set(bb,'edgecolor','none','facecolor',CL(3,:));
bb=bar(4,sum(Results_majority{2,2} ),0.8); set(bb,'edgecolor','none','facecolor',CL(4,:));
ylabel({'Total effort','($H$)'},'fontsize',FS)
set(gca,'fontsize',FS-2,'xtick',[1:4],'xticklabel',XTL)
CornerLetterLabel('(D)',[0.9 0.95],FS)

subplot(3,3,9), hold on
bb = patch([0 1 1 0],[0 0 1 1],'k'); set(bb,'edgecolor','none','facecolor',CL(1,:));
bb = patch([0 1 1 0],[0 0 1 1],'k'); set(bb,'edgecolor','none','facecolor',CL(2,:));
bb = patch([0 1 1 0],[0 0 1 1],'k'); set(bb,'edgecolor','none','facecolor',CL(3,:));
bb = patch([0 1 1 0],[0 0 1 1],'k'); set(bb,'edgecolor','none','facecolor',CL(4,:));
if InitialState == 1; L = legend('Open access','Consensus','Supermajority','Majority');
else InitialState == 0; L = legend('Exploratory','Consensus','Supermajority','Majority');
end
set(L,'box','off','fontsize',FS-2,'location','west')
set(gca,'visible','off')
xlim([6 10])

%Make_TIFF(Filename,[0 0 30 20])

function [Results,w] = Solve_for_voting(N,VotingThreshold,InitialState)

% Choose system parameters
r = 1;
cost = 0;

% Define their utility functions
w = linspace(0.05,0.95,N)'; % Low values of w don't care about conservation as much
R = 1; % This is a downweighting for conservation (which tends to be too high)

% Search algorithm tuning parameters
init_delta = 0.05; % Start with large changes
delta_threshold = 1e-6; % End with small changes
DecreaseCount = 500*N; % How many searches before reducing delta
shrink = 0.9; % Proportional reduction in delta
LowestHarvest = 5e-3; % Go no lower than this harvest rate

% Choose an initial set of harvests
h = ones(N,1).*LowestHarvest;
[y,n] = EquilPop(r,h); % Solve for the resulting harvest
c = ConservationFunction(n,R);
p = ProfitFunction(y,h,cost);
b = UtilityFunction(p,c,w);

if InitialState == 1 % If we're modelling a stock which starts at the OA state

    % Iterate unilateral harvest rates until we reach the Nash open access
    count = 0; delta = init_delta;
    while delta > delta_threshold

        % Choose a country at random
        this_country = randi(N);

        % Choose a direction to change at random
        this_direction = (randi(2)-1.5)*2;
        h_i = h; h_i(this_country) = max(LowestHarvest, h_i(this_country) + this_direction*delta);

        % Forecast the effect on equilibrium outcomes
        [y_i,n_i] = EquilPop(r,h_i);
        c_i = ConservationFunction(n_i,R);
        p_i = ProfitFunction(y_i,h_i,cost);
        b_i = UtilityFunction(p_i,c_i,w);

        % Accept the change if it improves the benefit for this particular country
        if b_i(this_country) > b(this_country)
            y = y_i; h = h_i; n = n_i; b = b_i; p = p_i; c = c_i;
            count = 0; % Reset the counter
        else count = count + 1;
        end

        if count == DecreaseCount
            delta = delta*shrink;
            count = 0;
        end
    end
end

Initial_Benefit = b;
Initial_h = h;
Initial_H = sum(h);
Initial_H_prop = h./Initial_H;

Results{1,1} = Initial_Benefit;
Results{1,2} = Initial_h;
Results{1,3} = n;
Results{1,4} = sum(y);
Results{1,5} = sum(y)./sum(h);
Results{1,6} = Initial_H;

% Randomly implement quota changes and accept solutions that pass voting
count = 0; delta = init_delta; TotalH = Initial_H;
while delta > delta_threshold

    % Change the overall quota, keeping the proportions 
    this_direction = (randi(2)-1.5)*2;
    TotalHi = TotalH;
    TotalHi = max(N*LowestHarvest, TotalHi + this_direction*delta);
    h_i = Initial_H_prop.*TotalHi;
    
    % Forecast the effect on equilibrium outcomes
    [y_i,n_i] = EquilPop(r,h_i);
    c_i = ConservationFunction(n_i,R);
    p_i = ProfitFunction(y_i,h_i,cost);
    b_i = UtilityFunction(p_i,c_i,w);

    % If enough countries vote for the measure, it passes
    if sum(b_i > b) >= VotingThreshold
        y = y_i; h = h_i; n = n_i; b = b_i; p = p_i; c = c_i; TotalH = TotalHi;
    end
    count = count + 1;

    if count == DecreaseCount % If we've tried a lot of things and nothing has improved
        delta = delta*shrink;
        count = 0;
    end
end
Voting_Benefit = b;
Voting_h = h;

Results{2,1} = Voting_Benefit;
Results{2,2} = Voting_h;
Results{2,3} = n;
Results{2,4} = sum(y);
Results{2,5} = sum(y)./sum(h);
Results{2,6} = sum(Voting_h);

end

% ================================
% === FUNCTIONS ==================
% ================================

function B = UtilityFunction(P,C,w)

% Cobb-Douglas utility function
B = C.^w .* P.^(1-w);

% % Linear additive utility function
% B = C*w + P.*(1-w);

end

function C = ConservationFunction(n,R)

% Conservation benefits accrue nonlinearly
% C = R.*(exp(1.*n)-1);

% Conservation benefits accrue linearly
C = R.*n;

end

function p_i = ProfitFunction(y_i,h_i,cost)

% Linear cost function
p_i = max(0,y_i - cost.*h_i);

end

function [y,n_star] = EquilPop(r,h)

% What's the total harvest?
H = sum(h);

% What are the proportional harvest rates?
h_prop = h./H;

% What's the equilibrium population
n_star = max(0,1 - H./r);

% What's the resulting SY?
Y = H.*n_star;

% What do each country get?
y = Y.*h_prop;

end