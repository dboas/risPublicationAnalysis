function pos = graphForces( pos, AuthorAdj, lstA, nIterMax );

nAuthorsA = length(lstA);

nIter = 0;
while nIter < nIterMax
    
    nIter = nIter + 1;
    
    % get distances between all nodes
    rhoM = zeros(nAuthorsA,nAuthorsA);
    for ii = 1:nAuthorsA
%         for jj = (ii+1):nAuthorsA
%             rhoM(ii,jj) = max(norm(pos(lstA(ii),:)-pos(lstA(jj),:)),0.1);
%         end
        rhoM(ii,(ii+1):nAuthorsA) = max( sum( (ones(nAuthorsA-ii,1)*pos(lstA(ii),:)-pos(lstA((ii+1):nAuthorsA),:)).^2,2).^0.5, 0.1 );
    end
    
    % get repulsion forces on each node
    Fx = zeros(nAuthorsA,1);
    Fy = zeros(nAuthorsA,1);
    for ii = 1:(nAuthorsA-1)
        for jj = (ii+1):nAuthorsA
            R = 2 / rhoM(ii,jj)^2;
            Fx(ii) = Fx(ii) + R * (pos(lstA(ii),1)-pos(lstA(jj),1))/rhoM(ii,jj);
            Fy(ii) = Fy(ii) + R * (pos(lstA(ii),2)-pos(lstA(jj),2))/rhoM(ii,jj);

            Fx(jj) = Fx(jj) - R * (pos(lstA(ii),1)-pos(lstA(jj),1))/rhoM(ii,jj);
            Fy(jj) = Fy(jj) - R * (pos(lstA(ii),2)-pos(lstA(jj),2))/rhoM(ii,jj);
        end
        
%         R = (1 ./ rhoM(ii,(ii+1):nAuthorsA).^2)';
%         Fx(ii) = Fx(ii) + sum( R .* (ones(nAuthorsA-ii,1)*pos(lstA(ii),1)-pos(lstA((ii+1):nAuthorsA),1))./rhoM(ii,(ii+1):nAuthorsA)' );
%         Fy(ii) = Fy(ii) + sum( R .* (ones(nAuthorsA-ii,1)*pos(lstA(ii),2)-pos(lstA((ii+1):nAuthorsA),2))./rhoM(ii,(ii+1):nAuthorsA)' );
% 
%         Fx((ii+1):nAuthorsA) = Fx((ii+1):nAuthorsA) - sum( R .* (ones(nAuthorsA-ii,1)*pos(lstA(ii),1)-pos(lstA((ii+1):nAuthorsA),1))./rhoM(ii,(ii+1):nAuthorsA)' );
%         Fy((ii+1):nAuthorsA) = Fy((ii+1):nAuthorsA) - sum( R .* (ones(nAuthorsA-ii,1)*pos(lstA(ii),2)-pos(lstA((ii+1):nAuthorsA),2))./rhoM(ii,(ii+1):nAuthorsA)' );
    end
    
    % now add the attractive forces
    for ii = 1:nAuthorsA
        for jj = (ii+1):nAuthorsA
            A = 0.5 * ( ( 1+log2(1+AuthorAdj(lstA(ii),lstA(jj))) ) ) * rhoM(ii,jj);
            Fx(ii) = Fx(ii) - A * (pos(lstA(ii),1)-pos(lstA(jj),1))/rhoM(ii,jj);
            Fy(ii) = Fy(ii) - A * (pos(lstA(ii),2)-pos(lstA(jj),2))/rhoM(ii,jj);

            Fx(jj) = Fx(jj) + A * (pos(lstA(ii),1)-pos(lstA(jj),1))/rhoM(ii,jj);
            Fy(jj) = Fy(jj) + A * (pos(lstA(ii),2)-pos(lstA(jj),2))/rhoM(ii,jj);
        end
    end
    
    % and add the attraction to the origin
    for ii = 1:nAuthorsA
%        Fx(ii) = Fx(ii) - 1e-2 * pos(lstA(ii),1) / norm(pos(lstA(ii),:))^2;
%        Fy(ii) = Fy(ii) - 1e-2 * pos(lstA(ii),2) / norm(pos(lstA(ii),:))^2;
        Fx(ii) = Fx(ii) - 2e-1 * pos(lstA(ii),1) / norm(pos(lstA(ii),:));
        Fy(ii) = Fy(ii) - 2e-1 * pos(lstA(ii),2) / norm(pos(lstA(ii),:));
    end
    
    % update positions
    scl = 0.1 / max(max(Fx),max(Fy));
    
    pos(lstA,1) = pos(lstA,1) + Fx * scl;
    pos(lstA,2) = pos(lstA,2) + Fy * scl;
    
end