function [F ,S_LocMatrix] = ObjFun_calcu_Myfq(veh,t,R)
[t_loc nt_loc] = ObjFun_calcu_V(veh,t);
veh1.member = t_loc;
S_loc = ObjFun_GetCdrty(veh1);
location = [S_loc.member{1} S_loc.member{2} S_loc.member{3} ...
            S_loc.member{4} S_loc.member{5}];
index = find(imag(location) ~= 0); 
location1 = location(index);
S_LocMatrix = ObjFun_CreateLocMatrix(location1,R);
tt_loc = S_LocMatrix.loc;
veh1.member = nt_loc;
nS_loc = ObjFun_GetCdrty(veh1);
nlocation = [nS_loc.member{1} nS_loc.member{2} nS_loc.member{3} ...
            nS_loc.member{4} nS_loc.member{5}];
index1 = find(imag(nlocation) ~= 0); 
 nlocation1 = nlocation(index1);    
 nS_LocMatrix = ObjFun_CreateLocMatrix(nlocation1,R);
 Ddist = nS_LocMatrix.dist - S_LocMatrix.dist ;
 F = ObjFun_CF(Ddist,S_LocMatrix.dist,R,S_LocMatrix);
end

function F = ObjFun_CF(Ddist,dist,R,S_LocMatrix)
% F = ObjFun_OneCF(Ddist,dist,R);

% F = ObjFun_testCF();
% F = ObjFun_thrCF(Ddist,dist,R,S_LocMatrix);
% F = ObjFun_fourCF(Ddist,dist,R,S_LocMatrix);
%   F = ObjFun_fiveCF(Ddist,dist,R,S_LocMatrix);
%  F = ObjFun_CF1(Ddist,dist,R,S_LocMatrix);%不行
F = ObjFun_SevenCF(Ddist,dist,R,S_LocMatrix);

end

function F = ObjFun_OneCF(Ddist,dist,R)
 k = Ddist;
 QiQj = k;
 index1 = find(Ddist>0);
%  k(index1) = exp(-Ddist(index1));
 k(index1) = 1./(1+(Ddist(index1)));
 QiQj(index1) = (2/3)*(2/3);
%  k(index1 ) = -Ddist(index1);
 index2 = find(Ddist<0);
%  k(index2) = exp(Ddist(index2));
  k(index2) = 1./(1-(Ddist(index2)));
 QiQj(index2) = (3/2)*(3/2);
 index3 = find(Ddist==0);
%  k(index3) = exp(Ddist(index3));
 k(index3) = 1./(1+(Ddist(index3)));
 QiQj(index3) = 2*2;
%  k(index2 ) = -1./Ddist(index2);
 index = find(dist == 0);
 t_F = k.*QiQj./(dist.*dist);
 t_F(index) =0;
% %超出通信距离
 index = find(dist>R);
 t_F(index) =0;
 F = ObjFun_OneV(t_F);
end
function F = ObjFun_TwoCF(Ddist,dist)
 k = Ddist;
 QiQj = k;
 index1 = find(Ddist>0);
 k(index1) = 1./(1+(Ddist(index1)));
 QiQj(index1) = (2/3)*(2/3);
 index2 = find(Ddist<0);
 k(index2) = 1./(1-(Ddist(index2)));
 QiQj(index2) = (3/2)*(3/2);
 index3 = find(Ddist==0);
 k(index3) = 1./(1+(Ddist(index3)));
 QiQj(index3) = 2*2;
 index = find(dist == 0);
 t_F = 1000*k.*QiQj./(dist.*dist);
 t_F(index) =0;
 F = ObjFun_OneV(t_F);
end

function F = ObjFun_OneV(FQ)
 t_sum = sum(FQ');
 len = length(t_sum);
 for n = 1:len
      FQ(n,n) = 0;
 end
 index = find(FQ>0);
 m_FO = min(FQ(index));
 F = FQ./m_FO;
 max_F = max(F(index));
 %对角元素给0
 for n = 1:len
      F(n,n) =1;
 end
end
function F = ObjFun_testCF()
% F = [1 10 10 10 8 8;10 8 0 0 0 0;10 0 8 0 0 0;10 0 0 8 0 0;8 0 0 0 8 0;8 0 0 0 0 8];
% F = [0 10 10 10 8 8;10 0 0 0 0 0;10 0 0 0 0 0;10 0 0 0 0 0;8 0 0 0 0 0;8 0 0 0 0 0];
% F = [1 10 10 10 0 0;10 1 0 0 0 0;10 0 1 0 0 0;10 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
% F = [1 1/10 1/10 1/10 0 0;1/10 1 0 0 0 0;1/10 0 1 0 0 0;1/10 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
%  F = [1 1/1001 1/1001 0 0 0 0 0 0 0;
%      1/1001 1 0 0 0 0 0 0 0 0;
%      1/1001 0 1 0 0 0 0 0 0 0; 
%      0 0 0 1 1/1000 1/1000 0 0 0 0;
%      0 0 0 1/1000 1 0 0 0 0 0;
%      0 0 0 1/1000 0 1 0 0 0 0;
%      0 0 0 0 0 0 1 0 0 0;
%      0 0 0 0 0 0 0 1 0 0;
%      0 0 0 0 0 0 0 0 1 0;
%      0 0 0 0 0 0 0 0 0 1;];
 t_F = [6 6 6 0 0 0 0 0 0;
        0 6 6 0 0 0 0 0 0;
        0 0 6 0 0 0 0 0 0;
        0 0 0 6 6 6 0 0 0;
        0 0 0 0 6 6 0 0 0;
        0 0 0 0 0 6 0 0 0;
        0 0 0 0 0 0 6 6 6;
        0 0 0 0 0 0 0 6 6;  
        0 0 0 0 0 0 0 0 6;
        ]
    F = t_F+t_F';
end

function [t_loc nt_loc] = ObjFun_calcu_V(veh,t)
    for n =1:5
        speed = veh.speed{n};
        vehicle = veh.member{n};
        vehicle = vehicle(1:t);
        index = find(vehicle>0);
        v_loc = vehicle;
        v_loc(index) = 1;
        t_loc{n} = v_loc.*((fliplr(1:t))')*(speed/4);
        nt_loc{n} = v_loc.*((fliplr(2:(t+1)))')*(speed/4) ;
    end
end

function F = ObjFun_thrCF(Ddist,dist,R,S_LocMatrix)
 ab_dist = abs(S_LocMatrix.loc) ;
 len = length(ab_dist);
 AB_dist = zeros(len);
 for n = 1:len;
     AB_dist(n,:) = ab_dist;
 end
 AB_dist = (AB_dist+ AB_dist')./2;
 k = Ddist;
 QiQj = k;
 index1 = find(Ddist>0);
 k(index1) = 1./(1+(Ddist(index1)));
 QiQj(index1) = (2/3)*(2/3);
 index2 = find(Ddist<0);
  k(index2) = 1./(1-(Ddist(index2)));
 QiQj(index2) = (3/2)*(3/2);
 index3 = find(Ddist==0);
 k(index3) = 1./(1+(Ddist(index3)));
 QiQj(index3) = 2*2;
 index = find(dist == 0);
 t_F = k.*QiQj./(1+dist);
 t_F(index) =0;
% %超出通信距离
 index = find(dist>R);
 t_F(index) =0;
%  offset = (AB_dist').*(1/R);
%  offset(index) = 0;
   F = t_F +eye(length(t_F)); 
end
function F = ObjFun_fourCF(Ddist,dist,R,S_LocMatrix)
 k = Ddist;
 QiQj = k;
 index1 = find(Ddist>0);
 k(index1) = 1./(1+(Ddist(index1)));
 QiQj(index1) = (2/3)*(2/3);
 index2 = find(Ddist<0);
 k(index2) = 1./(1-(Ddist(index2)));
 QiQj(index2) = (3/2)*(3/2);
 index3 = find(Ddist==0);
 k(index3) = 1./(1+(Ddist(index3)));
 QiQj(index3) = 2*2;
 index = find(dist == 0);
%  t_F = 100.*k.*QiQj./((1+dist));
%  t_F = 100./(dist.*dist);
 t_F(index) = 0;
% %超出通信距离
%  index = find(dist>R);
%  t_F(index) =0;
 %% 孤立点处理
 F = t_F +eye(length(t_F)); 
end
function F = ObjFun_sixCF(Ddist,dist,R,S_LocMatrix)

 k = Ddist;
 QiQj = k;
 index1 = find(Ddist>0);
 k(index1) = Ddist(index1);
 QiQj(index1) = (2/3)*(2/3);
 index2 = find(Ddist<0);
 k(index2) = -Ddist(index2);
 QiQj(index2) = (3/2)*(3/2);
 index3 = find(Ddist==0);
 k(index3) = Ddist(index3);
 QiQj(index3) = (3/2)*(3/2);
 index = find(dist == 0);
 t_F = zeros(length(dist));
 index1 = find(dist>R);
 t_F(index1) = 0.*QiQj(index1)./((1+k(index1)).*(1+dist(index1)));
 index2 = find(dist<R);
%  t_F(index2) = 1.*QiQj(index2)./((1+k(index2)/R).*(1+dist(index2)/R));
 t_F(index2) = 100.*QiQj(index2)./((1+k(index2)).*(1+(dist(index2).^2)));
  index3 = find(dist==R);
 t_F(index3) = 100.*QiQj(index3)./((1+k(index3)).*(1+(dist(index3)).^2));
 t_F(index) = 0;
 %% 孤立点处理
 F = t_F+eye(length(t_F));
end
function F = ObjFun_fiveCF(Ddist,dist,R,S_LocMatrix)

 k = Ddist;
 QiQj = k;
 index1 = find(Ddist>0);
 k(index1) = Ddist(index1);
 QiQj(index1) = (2/3)*(2/3);
 index2 = find(Ddist<0);
 k(index2) = -Ddist(index2);
 QiQj(index2) = (3/2)*(3/2);
 index3 = find(Ddist==0);
 k(index3) = Ddist(index3);
 QiQj(index3) = (3/2)*(3/2);
 index = find(dist == 0);
 t_F = zeros(length(dist));
 index1 = find(dist>R);
 t_F(index1) = 0.*QiQj(index1)./((1+k(index1)).*(1+dist(index1)));
 index2 = find(dist<R);
%  t_F(index2) = 1.*QiQj(index2)./((1+k(index2)/R).*(1+dist(index2)/R));
 t_F(index2) = 1000.*QiQj(index2)./((1+k(index2)).*((1+dist(index2).^2)));
  index3 = find(dist==R);
 t_F(index3) = 1000.*QiQj(index3)./((1+k(index3)).*((1+dist(index3)).^2));
 t_F(index) = 0;
 %% 孤立点处理
 F = t_F+eye(length(t_F));
end
function F = ObjFun_CF1(Ddist,dist,R,S_LocMatrix)
% 1/(1+r+d)^2  不行
RD = Ddist./2 + dist;
t_F = 1000./(1+RD).^2
% more than R
index = find(dist>R);
t_F(index) = 0;
%self == 1
index = find(dist == 0);
t_F(index) = 1;
F =  t_F;
end

function F = ObjFun_SevenCF(Ddist,dist,R,S_LocMatrix)

%%% 1/(1+exp(-x))

t_F = 1./(1+exp(-Ddist));
index = find(dist>R);
t_F(index) = 0;
F = t_F;
end
