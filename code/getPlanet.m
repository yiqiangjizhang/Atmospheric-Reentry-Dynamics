function [Planet] = getPlanet()

% Function that returns the planet for Atmospheric Re-Entry

a = true;
    while(a)
        disp('Insert a planet to analyse the Atmospheric Re-Entry');
        disp('You can choose between Mars and Earth');
        disp('FYI: Mars will only include atmospheric thermodynamical properties');
        prompt = '[Answer (M(Mars)/E(Earth)]:';
        atmosphere = input(prompt,'s');
        
        % Get planet from user input
        if ((atmosphere == "E") || (atmosphere == "Earth") || (atmosphere == "e"))
            Planet = '1';
            a = false;
        elseif ((atmosphere == "M") || (atmosphere == "Mars") || (atmosphere == "m"))
            Planet = '2';
            a = false;
        else
            fprintf('Invalid answer. Please re-enter the answer.')
        end
    end

end

