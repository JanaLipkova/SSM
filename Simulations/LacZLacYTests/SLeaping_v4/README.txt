S Leaping + if L = 0, recompute until L>0 and set dt = dt + #repeats*dt

dt = tau
L = pois(a0*dt)
count = 0;

while (L == 0)
{
   L = pois(a0*dt);
   count++
}

dt = dt + dt * count

