; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
;--------------------------
FUNCTION aii_plotsymbols, i
;--------------------------
;
;==========================================================================
;+
;NAME:
;       aii_plotsymbols
; PURPOSE:
;       defines user defined plotting symbol USERSYM (PSYM=8) 
; EXPLANATION:
;       Position PSYM=8 is the free position in IDL to define 
;       a user defined plotting symbol called USERSYM. Because
;       IDL has only seven build in symbols, it might be necessary
;       to have more than these for ARTS. This function provides
;       13 different symbols, including those of the standart IDL ones.
;
; CALLING EXAMPLES:
;       psym=aii_plotsymbols(5) e.g. in the keyword definition of 
;       the command 'plot' or 'oplot' to define a triangle for 
;       USERSYM (PSYM=7)
;
; INPUTS:
;       i   (integer)  selection of the symbol.
;                      possible value: 1 to 13.
;
; OUTPUTS:
;       no own output, but temporarily definition of the
;       plotting symbol USERSYM (PSYM=8).
;       
;
;
; MODIFICATION HISTORY:
;       03/04/01  TKS  alpha version created 
;-
; ==========================================================================
;
;
CASE i OF
    1: begin
;       plus sign:
        return, i
       end
;
    2: begin
;       asterisk:
        return, i
       end
;
    3: begin
;       open circle
        px = FINDGEN(48) * (!PI*2.0/48.0) ;Make a vector of 16 points, A[i] = 2pi/16.
        USERSYM, COS(px), SIN(px), /THICK
        return, 8
       end
;
    4: begin
;       open diamond:
        return, i
       end
;
    5: begin
;       triangle
        return, i
       end
;
    6: begin
;       open square:
        return, i
       end
;
    7: begin
;       X:
        return, i
       end
;
    8: begin
;       filled circle
        px = FINDGEN(48) * (!PI*2.0/48.0) ;Make a vector of 16 points, A[i] = 2pi/16.
        USERSYM, COS(px), SIN(px), /FILL
        return, 8
       end
;
    9: begin
;       open triangle  on the edge
        px = 1.0*[-1.0,  0.0, 1.0, -1.0]
        py = 2.0*[ 0.5, -0.5, 0.5,  0.5]
        USERSYM, px, py, /THICK
        return, 8
       end
    10: begin
;        filled triangle  on the edge
         px = 1.0*[-1.0,  0.0, 1.0, -1.0]
         py = 2.0*[ 0.5, -0.5, 0.5,  0.5]
         USERSYM, px, py, /FILL
         return, 8
        end
    11: begin
;        filled triangle on the basis:
         px = 1.0*[-1.0, 0.0, 1.0, -1.0]
         py = 2.0*[ -0.5, 0.5, -0.5,  -0.5]
         USERSYM, px, py, /FILL
         return, 8
        end
    12: begin
;        filled diamond:
         px = 1.0*[-1.0,  0.0, 1.0, 0.0, -1.0]
         py = 2.0*[ 0.0, -0.5, 0.0, 0.5,  0.0]
         USERSYM, px, py, /FILL
         return, 8
        end
    13: begin
;       filled square:
        px = 1.0*[-1.0,  1.0,  1.0, -1.0]
        py = 1.0*[-1.0, -1.0,  1.0,  1.0]
        USERSYM, px, py, /FILL
        return, 8
       end
ELSE:  begin
;       filled square:
        px = 1.0*[-1.0, 1.0, 1.0, -1.0]
        py = 2.0*[ 0.0, 0.0, 1.0,  1.0]
        USERSYM, px, py, /FILL
        return, 8
       end
ENDCASE
;
return, 8
;
END
; ###############################################################################
