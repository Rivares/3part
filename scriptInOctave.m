SCH = zeros(100001, 5);

[a, b] = textread ("Tv.txt", "%f %s");

z = 1;
for i = 1:100001
  for j = 1:5

  if z+1 > 600000
        return;
    endif;
    
    if rem(z,6) != 0
      SCH(i, j) = a(z);
    else
      ++z;
      SCH(i, j) = a(z);
    endif;
    ++z;
    
  endfor;
endfor;

plot(SCH(1:100000,2))
hold on
plot(SCH(1:100000,3))