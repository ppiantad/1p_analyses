function [hl,ht] = xlin(x,txt,ylo,yhi,ytxt) 
% Documentation: 
%   x       = x-Position
%   txt     = Text String
%   ylo     = Low y-Value (Start)
%   yhi     = High y-Value (End)
%   ytxt    = Text Starting Position
hold on
hl = plot([x x],[ylo yhi],'DisplayName',txt, 'LineWidth',1);
ht = text(x,ytxt, txt, 'Horiz','left', 'Vert','top', 'Rotation',90);
hold off
end