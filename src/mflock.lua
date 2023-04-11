-- Lua 5.3.5

function getConfig(iter, newx)
  -- iter is current iteration
  -- newx is 1: if points were placed randomly or,
  --         0: if points were loaded from a file

  maxiter = 8000

  ---- proportion of steps taken
  q = iter/maxiter

  ---- Compression force, kCom
  -- Compacts chromosomes towards their centre of mass
  -- q <= 0.2         OFF, so that the beads can relax a little
  -- 0.2 <= q <= 0.8  ON
  -- q >= 0.8         OFF again, to avoid overlapping beads etc
  kCom = 0
  if newx == 1 then
    if q > 0.2 and q < 0.8 then
      kCom = 1
    end
  end

  ---- Domain force, kDom
  -- Force that keeps bead in domain
  kDom = 1.0

  ---- Volumetric interaction, kVol
  -- Restrain beads from overlapping
  kVol = 1

  -- Radial force, kRad,
  -- Only use if GPSeq is included
  -- In the GPSeq paper we used 0.005 for 10k Haploid structures at 1 MB
  kRad = 0

  ---- Brownian force, fBrown
  -- Default: Decreases linearly from 0.7 to 0 along the iterations
  fBrown = 0.7*(1-iter/maxiter)

  ---- Interaction force, kInt
  -- attracts beads that should be in contact (according to W)
  kInt = 0.5*(1.0 + math.sin( (q - .5)*math.pi )) -- from 0 to 1

  ---- Interaction distance, dInteraction
  dInteraction = (1.1 + 0.9*q)*2 -- from 1.8 to 4 diameters radi

  ---- Exit condition
  quit = 0
  if iter >= maxiter then
    quit = 1
  end

  -- To slow down the simulations:
  -- usleep(1000000/60)


end
