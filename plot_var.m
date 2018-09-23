function plot_var(x,l,u,c)

    x = [ x x(end:-1:1)];
    y = [ l u(end:-1:1)];

    id=~isnan(x)&~isnan(y);
    
    patch(x(id),y(id),c,'FaceAlpha',0.5);
   