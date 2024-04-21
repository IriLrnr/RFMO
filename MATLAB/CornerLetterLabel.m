function tx = CornerLetterLabel(TXT,Prop,FS)

if length(Prop) == 1
   Prop = [Prop Prop];
end

if nargin == 3
   LT = 0;
end

XL = xlim; YL = ylim; 
XR = XL(2)-XL(1);
YR = YL(2)-YL(1);

tx = text(XL(1)+Prop(1)*XR,YL(1)+Prop(2)*YR,TXT);

set(tx,'fontsize',FS);
