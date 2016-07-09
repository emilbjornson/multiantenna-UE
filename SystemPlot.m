function Distance = SystemPlot(nbrBSs, K, radius)
%
%INPUT:
%nbrBSs:    number of BSs
%K:         number of UEs per cell
%radius:    radius of a cell
%
%OUTPUT:
%Distance:  squared distance from the users to BSs
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



%Generate the BS locations in the area, Positions of BSs of the central

%Compute alternative BS locations by using wrap around.

if (nbrBSs ==1)
    
    % 1 cell structure
    BSpositions = [0 0 ];
    wrapHorizontal=0;
    wrapVertical = 0;
    
elseif (nbrBSs == 7)
    
    %cell and the one tier around
    BSpositions = [0 BSofHexagonNetwork(1,radius) ];
    
    % 7 cells structure
     wrapHorizontal = radius * [0 4.5 1.5 -3 -4.5 -1.5 3];
    wrapVertical = radius*sqrt(3)/2*[0 1 5 4 -1 -5 -4];
end

wrapLocations = wrapHorizontal+ 1i*wrapVertical;  % Move the original nbrBSs BSs to the left, right, up, down, leftup corner, leftdown corner, rightup corner, rightdown corner.
BSpositionsWrapped = repmat(BSpositions.',[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);  % Each row of BSpositionsWrapped represents the wrapped positions of a BS, a nbrBss*7 matrix

%Prepare to put out UEs in the cells
UEpositions = zeros(K,nbrBSs);

distancesSquaredBSj = zeros(K*nbrBSs, nbrBSs);

%Go through all the cells
for l = 1:nbrBSs
     
    %Put out K users in each cell
    vertices=exp(1j*pi/3*(1:6));
    UEpositions(:,l) = DropUEinHexagonCell(K,vertices,radius) + BSpositions(l);
 
    for j = 1:nbrBSs
        %Compute the distance from the users to BS j
        distancesSquaredBSj(1+(l-1)*K:l*K, j) = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);

    end

end

Distance=distancesSquaredBSj;



