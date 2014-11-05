; stolen from IDL's online HELP on "Three-Dimensional Graphics"
FUNCTION CVT_TO_2D, x, y, z
  P = [!X.S[0] + !X.S[1] * X,	$
       !Y.S[0] + !Y.S[1] * Y,	$
       !Z.S[0] + !Z.S[1] * Z,	$
       1]

  P = P # !P.T

  RETURN, [P[0] / P[3], P[1] / P[3]]
END


PRO plot_axes, frame, pos, et, _Extra=_e
  FOR axis=0,2 DO BEGIN
    xyz = FltArr(3)
    xyz[axis] = 0.5

    ; XYZ is in our desired frame... convert to 'HAE', Heliocentric Aries 
    ; Ecliptic, which is pretty close to what we're plotting in.
    vec = cxform(xyz * 1.5D8, frame, 'HAE', et) / 1.5D8

    PlotS, [pos[0],vec[0]], [pos[1],vec[1]], [pos[2],vec[2]], /T3D, $
      _Extra=_e, Thick=1 + (axis EQ 0)*2
  ENDFOR
END

PRO coord_anim, frame1, frame2, step=step, _Extra=_e
  IF N_Elements(frame1) EQ 0 THEN BEGIN
    MESSAGE, 'Usage: ' + myname() + ', "FRAME1 [,FRAME2] [,step=]"', /Continue
    RETURN
  ENDIF
    
  spice_init_generic

  ; Set up the display.  If we don't already have a window up, make one...
  Set_Plot, 'X'
  IF !D.Window EQ -1 THEN ww, 'Coord_Anim'
  ; ...and make a Z buffer of the same dimensions
  xysize = [ !D.X_Size, !D.Y_Size ]
  Set_Plot, 'Z'
  Device, Set_Resolution=xysize
  Set_Plot, 'X'

  ; Define colors
  c_earth = 3	& TvLct,   0,   0, 255, c_earth		; blue, for Earth
  c_sun   = 4	& TvLct, 255, 255,   0, c_sun		; yellow, for Sun
  c_frame1= 10	& TvLct, 128, 128, 128, c_frame1	; gray, for frame 1
  c_frame2= 11	& TvLct, 200,   0,   0, c_frame2	; Frame 2

  ; Define the user symbol ass a circle, for drawing Earth & Sun
  UserSym_Circle, radius=2

  ; Set up the starting date (Vernal Equinox), and go for one year
  cspice_str2et, '1998 mar 21 0:00:00', et0
  et1 = et0 + (365D * 86400D)

  smooth = 1

  ; Assume one-hours steps, but that can be overriden from the command line
  stepsize = 3600D
  IF N_Elements(step) NE 0 THEN stepsize = stepsize * step

  ; Set up plot range
  Surface, dist(10), /save, /NoData,				$
    XRange=[-1.5, 1.5], XStyle=5, XMargin=[0,0],		$
    YRange=[-1.5, 1.5], YStyle=5, YMargin=[0,0],		$
    ZRange=[-1.5, 1.5], ZStyle=5, ZMargin=[0,0], _Extra=_e

  FOR et=et0, et1, stepsize DO BEGIN
    IF smooth THEN Set_Plot, 'Z'
    erase

    ; Draw Mister Sun
    xy = Cvt_to_2D(0,0,0)
    PlotS, xy[0], xy[1], PSym=8, Color=c_sun, SymSize=2, /Normal

    ; Find out where Earth is, and draw it
    cspice_spkez, 399, et, 'ECLIPJ2000', 'NONE', 10, st, cl
    st = st / 1.5D8

    xy = Cvt_To_2D(st[0],st[1],st[2])
    PlotS, xy[0], xy[1], PSym=8, Color=c_earth, /Normal

    ; Draw X, Y, and Z axes at Earth, showing the pointing of this frame
    Plot_Axes, frame1, st[0:2], et, Color=c_frame1
    IF N_Elements(frame2) NE 0 THEN $
      Plot_Axes, frame2, st[0:2], et, Color=c_frame2

    cspice_et2utc, et, 'C', 0, tout
    XYOutS, 0.01, 0.95, tout, /Normal

    XYOutS, 0.01, 0.91, frame1, /Normal, Color=c_frame1, CharThick=2
    IF N_Elements(frame2) NE 0 THEN $
      XYOutS, 0.01, 0.87, frame2, /Normal, Color=c_frame2, CharThick=2

    IF smooth THEN BEGIN
      pic = TvRd()
      Set_Plot, 'X'
      Tv, pic
    ENDIF
  ENDFOR
END
