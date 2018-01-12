h = [35; 36; 37; 38];
g = [60;61;62; 63];
v0 = [5; 6; 7; 8];
for i = 1:4
    v_body(:,i) = [v0(i)*cosd(h(i)) + v0(i)*cosd(g(i));...
                   v0(i)*sind(h(i)); v0(i)*sind(g(i))];
R1=[cosd(h(i)) -sind(h(i)) 0; sind(h(i)) cosd(h(i)) 0;0 0 1];
v(:,i) = R1*v_body(:,i); 
end

