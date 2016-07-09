function P = BSofHexagonNetwork(n,r)
%
%INPUT:
%n:         which tier to be constructed, n=1,2,3
%r:         radius of hexagonal cell
%
%OUTPUT:
%P:         BS positions
%           P(1,:): x positions of cells
%           P(2,:): y postions of cells
%
%
%This function is used in the paper:
%
%Xueru Li, Emil Björnson, Shidong Zhou, Jing Wang, “Massive MIMO with
%Multi-Antenna Users: When are Additional User Antennas Beneficial?,”
%Proceedings of International Conference on Telecommunications (ICT),
%Thessaloniki, Greece, May 2016.
%
%Download article: http://arxiv.org/pdf/1603.09052
%
%This is version 1.0 (Last edited: 2016-07-09)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original paper listed above.



a=[3/2; sqrt(3)/2]*r;
b=[0;sqrt(3)]*r;

if(n==1)
    m=6;
    alfa1=[1 0 -1 -1 0 1];
    alfa2=[0 1 1 0 -1 -1];
    PP=repmat(a,1,m).*repmat(alfa1,2,1) + repmat(b,1,m).*repmat(alfa2,2,1); 
elseif(n==2)
    m=12;
    alfa1=[2 1 0 -1 -2 -2 -2 -1 0 1 2 2];
    alfa2=[0 1 2 2 2 1 0 -1 -2 -2 -2 -1];
    PP=repmat(a,1,m).*repmat(alfa1,2,1) + repmat(b,1,m).*repmat(alfa2,2,1); 
elseif(n==3)
    m=18;
    alfa1=[3 2 1 0 -1 -2 -3 -3 -3 -3 -2 -1 0 1 2 3 3 3];
    alfa2=[0 1 2 3 3 3 3 2 1 0 -1 -2 -3 -3 -3 -3 -2 -1];
    PP=repmat(a,1,m).*repmat(alfa1,2,1) + repmat(b,1,m).*repmat(alfa2,2,1); 
end

P=PP(1,:) + 1i* PP(2,:);
 