function[a] = draw_orbital_save(G,C,Nel,title_fig,varargin)

Ccil = C;  nstate = Nel;
f = figure('visible','off')
plot(G(:,1), G(:,2),'*','MarkerSize',2);
hold on
for  site = 1:length(G)
   if ( Ccil(site,nstate) < 0)
      color_circle = 'b';
   else
      color_circle = 'r';
   end
   %plot(G(site,1), G(site,2),color_circle, 'MarkerSize',5);
   %plot(G(site,1), G(site,2), '.r', 'MarkerSize',max(abs(Ccil(site,nstate)),0.000001))
   %circles([G(site,1), G(site,2)],max(abs(Ccil(site,nstate)),0.000001),'Color',color_circle);
  filledCircle([G(site,1), G(site,2)],max(2.5*abs(Ccil(site,nstate)),0.000001),1000,color_circle);
end % site

%draw atomic bonds
for  site = 1:length(G)

   for  site2 = site+1:length(G)


      if (  norm(G(site,:)-G(site2,:)) < 1.7  )  %neighbors


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

narginv = length(varargin);

if (narginv >= 1)
 title(varargin{1})
end %if narginv >= 1

%title(title_fig)
saveas (f, title_fig,'jpg')

end %end function draw_orbital
