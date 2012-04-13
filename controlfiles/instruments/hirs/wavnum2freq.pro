FUNCTION wavnum2freq, wavnum
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;
;; Convert wavenumber to frequency. 
;;
;; Input  :  Wavenumber ( per cm)
;; Output :  Frequency  ( Hz )
;;
;; 16-01-2005 VOJ<vojohn@uni-bremen.de>
;;
;----------------------------------------------------------------

wavnum  =  DOUBLE( wavnum )
c       =  3.0*1e10

freq  =  c * wavnum

RETURN, freq
END  
