function psi = za_geom2other_point( r1, r2, dalpha )
  if( dalpha == 0 )
    psi = 0;
  else
    global DEG2RAD RAD2DEG
    r1s   = r1 * r1;
    r2s   = r2 * r2;
    l     = sqrt( r1s + r2s - 2*r1*r2*cos(DEG2RAD*abs(dalpha)) );
    psi   = RAD2DEG*acos( (r2s-r1s-l*l) / (2*r1*l)  ); 
    if( dalpha < 0 )
      psi = -psi;
    end
  end
return
