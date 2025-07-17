function[a] = draw_structure(G)


f = figure('visible','on')
plot(G(:,1), G(:,2),'*','MarkerSize',2);
hold on

%draw atomic bonds
for  site = 1:length(G)
 text( G(site,1),G(site,2),num2str(site),'FontSize',18,'FontWeight','bold')
   for  site2 = site+1:length(G)


      if (  norm(G(site,:)-G(site2,:)) < 1.8  )  %neighbors


         a = G(site,1);
         b = G(site2,1);
         c = G(site,2);
         d = G(site2,2);

         line([a b],[c d],"linestyle", "-", "color", "k")

      end % end if neighbors


   end % site2

end %site

axis("equal")
axis off





end %end function draw_orbital
