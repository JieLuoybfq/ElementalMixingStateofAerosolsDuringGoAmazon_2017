function systemout = ODsystem(f,OD,mass,eleflaglist)

%this function is meant to represent a single particle

%this function is setup to be solved by fsolve, which takes only a single
%variable input.  f(1:14) stands for mass fractions of the elements
%C,N,O,Na,Mg,P,S,Cl,K,Ca,Mn,Fe,Ni,Zn in that order.
%
%
%f(15) currently stands for rho*t because of the inability to solve the
%equation with those two variables separate
%f(15) is a stand in for rho, the density of the particle, and f(16) stands
%for t, the thickness of the particle.
%
%OD is a constant
%
%Even though this will almost always produce a system of equations with
%more equations than unknown (aka an overdetermined system) it is good.
%While an overdetermined system, mathematically speaking, has conflicting,
%duplicate, or linear combinations of existing equations making up the
%extra equations and thus will often have no analytical solution, solving a
%least squares problem will have a minimum which will not be ruined by
%additional equations.
%
% The addt'l eqn's will probably help the problem converge faster.

ulist = STXMEDXulist;

systemout(1) = -OD(1) + f(15).*(...
    f(1) .*ulist{1,1} + ...
    f(2) .*ulist{2,1} + ...
    f(3) .*ulist{3,1} + ...
    f(4) .*ulist{4,1} + ...
    f(5) .*ulist{5,1} + ...
    f(6) .*ulist{6,1} + ...
    f(7) .*ulist{7,1} + ...
    f(8) .*ulist{8,1} + ...
    f(9) .*ulist{9,1} + ...
    f(10).*ulist{10,1} + ...
    f(11).*ulist{11,1} + ...
    f(12).*ulist{12,1} + ...
    f(13).*ulist{13,1} + ...
    f(14).*ulist{14,1});


systemout(2) = -OD(2) + f(15).*(...
    f(1) .*ulist{1,2} + ...
    f(2) .*ulist{2,2} + ...
    f(3) .*ulist{3,2} + ...
    f(4) .*ulist{4,2} + ...
    f(5) .*ulist{5,2} + ...
    f(6) .*ulist{6,2} + ...
    f(7) .*ulist{7,2} + ...
    f(8) .*ulist{8,2} + ...
    f(9) .*ulist{9,2} + ...
    f(10).*ulist{10,2} + ...
    f(11).*ulist{11,2} + ...
    f(12).*ulist{12,2} + ...
    f(13).*ulist{13,2} + ...
    f(14).*ulist{14,2});

systemout(3) = -OD(3) + f(15).*(...
    f(1) .*ulist{1,3} + ...
    f(2) .*ulist{2,3} + ...
    f(3) .*ulist{3,3} + ...
    f(4) .*ulist{4,3} + ...
    f(5) .*ulist{5,3} + ...
    f(6) .*ulist{6,3} + ...
    f(7) .*ulist{7,3} + ...
    f(8) .*ulist{8,3} + ...
    f(9) .*ulist{9,3} + ...
    f(10).*ulist{10,3} + ...
    f(11).*ulist{11,3} + ...
    f(12).*ulist{12,3} + ...
    f(13).*ulist{13,3} + ...
    f(14).*ulist{14,3});

systemout(4) = -OD(4) + f(15).*(...
    f(1) .*ulist{1,4} + ...
    f(2) .*ulist{2,4} + ...
    f(3) .*ulist{3,4} + ...
    f(4) .*ulist{4,4} + ...
    f(5) .*ulist{5,4} + ...
    f(6) .*ulist{6,4} + ...
    f(7) .*ulist{7,4} + ...
    f(8) .*ulist{8,4} + ...
    f(9) .*ulist{9,4} + ...
    f(10).*ulist{10,4} + ...
    f(11).*ulist{11,4} + ...
    f(12).*ulist{12,4} + ...
    f(13).*ulist{13,4} + ...
    f(14).*ulist{14,4});

systemout(5) = -OD(5) + f(15).*(...
    f(1) .*ulist{1,5} + ...
    f(2) .*ulist{2,5} + ...
    f(3) .*ulist{3,5} + ...
    f(4) .*ulist{4,5} + ...
    f(5) .*ulist{5,5} + ...
    f(6) .*ulist{6,5} + ...
    f(7) .*ulist{7,5} + ...
    f(8) .*ulist{8,5} + ...
    f(9) .*ulist{9,5} + ...
    f(10).*ulist{10,5} + ...
    f(11).*ulist{11,5} + ...
    f(12).*ulist{12,5} + ...
    f(13).*ulist{13,5} + ...
    f(14).*ulist{14,5});

systemout(6) = -OD(6) + f(15).*(...
    f(1) .*ulist{1,6} + ...
    f(2) .*ulist{2,6} + ...
    f(3) .*ulist{3,6} + ...
    f(4) .*ulist{4,6} + ...
    f(5) .*ulist{5,6} + ...
    f(6) .*ulist{6,6} + ...
    f(7) .*ulist{7,6} + ...
    f(8) .*ulist{8,6} + ...
    f(9) .*ulist{9,6} + ...
    f(10).*ulist{10,6} + ...
    f(11).*ulist{11,6} + ...
    f(12).*ulist{12,6} + ...
    f(13).*ulist{13,6} + ...
    f(14).*ulist{14,6});


% 
% systemout(7) = -f(1) + f(2).*(mass(1)./mass(2));
% 
% systemout(8) = -f(2) + f(3).*(mass(2)./mass(3));

eleidx_zeros = find(eleflaglist==0);
eleidx = find(eleflaglist==1); %retrieves linear index of all non-zero elements

%This, much simpler thing, only works if C, N, and O STXM data is present
%(which it usually is for this study)
% systemout(7) = -f(1) + f(2) .* mass(1) ./ mass(2) ;
% systemout(8) = -f(1) + f(3) .* mass(1) ./ mass(3) ;

% %this mess is adding 2 equations to our system of unkn's and determines how
% %to get those two equations using C,N, and O data.  Either 2 useful
% %equations come out, or 2 trivial equations come out depending on if C,N,
% %and O data exists.
if ~isempty(eleidx(eleidx == 1)); %if eleidx(lindex of nonzero elements) is 1, C is present
    if ~isempty(eleidx(eleidx ==2)); %if eleidx(lindex of nonzero elements) is 2, N is present
        
        systemout(7) = -f(1) + f(2).*(mass(1)./mass(2));
        
        if ~isempty(eleidx(eleidx == 3)); %if eleidx(lindex of nonzero elements) is 3, O is present
            systemout(8) = -f(2) + f(3) .*(mass(2)./mass(3));
        else %No O
            systemout(8) = f(3);
        end
        
    else %No N
        
        systemout(7) = f(2);
        
        if ~isempty(eleidx(eleidx == 3)); %O present
            systemout(8) = -f(1) + f(3) .*(mass(1)./mass(3));
        else %No O
            systemout(8) = f(3);
        end
    end
    
else %No C
    systemout(7) = f(1);
    
    if ~isempty(eleidx(eleidx == 2)) && ~isempty(eleidx(eleidx==3)); %N AND O present
        systemout(8) = -f(2) + f(3) .*(mass(2)./mass(3));
        
    else %Either N or O or neither is present.  This doesn't check the case where only N or only O is present, but then a system with 2 mass parameters cannot be made.
        systemout(8) = f(2);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%



semeleidx = eleidx(eleidx>3); %this gets the indexes of all non-zero elements that were analyzed with SEM
for i = 1:(length(semeleidx)-1);
    systemout(8+i) = -f(semeleidx(i)) + f(semeleidx(i+1)).*(mass(semeleidx(i))./mass(semeleidx(i+1))); 
end

numeqns_sem = 8+length(semeleidx)-1;
for j = 1:length(eleidx_zeros) %this explicitly makes f of any zero element, equal to zero
    systemout(numeqns_sem+j) = f(eleidx_zeros(j));
end
numeqns_sem_zeros = numeqns_sem+length(eleidx_zeros);

systemout(numeqns_sem_zeros+1) = 1- sum(f(1:14));

% numeqns_mfracs = numeqns_sem + length(eleidx_zeros);
% 
% systemout(numeqns_mfracs+1) = f(15).*f(16).*A - totM;




end
    